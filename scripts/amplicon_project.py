#!/usr/bin/env python

import os, sys
import argparse
import shutil
import pandas as pd
import glob
import yaml
import re
import json
import urllib.request
import subprocess
from collections import defaultdict
from pathlib import Path

# ==================================
#           UPDATE CHECKER
# ==================================

def check_for_updates(repo_owner: str, repo_name: str):
    """
    Checks for a new version of the software on GitHub.

    This function compares the local version (from /NAAM_VERSION) with the latest
    release on GitHub. If a newer version is available, it notifies the user
    and asks if they want to abort the script.

    Args:
        repo_owner (str): The owner of the GitHub repository (e.g., 'your_username').
        repo_name (str): The name of the GitHub repository (e.g., 'naam-workflow').
    """
    # --- Define some colors for the output message for better visibility ---
    class Colors:
        YELLOW = '\033[93m'
        GREEN = '\033[92m'
        RED = '\033[91m'
        ENDC = '\033[0m'

    # Define the two possible locations for the version file.
    container_path = Path("/NAAM_VERSION")

    # When running locally, the script is in '.../scripts/amplicon_project.py'.
    # The version file is in the project root, so we go up two directories.
    script_location = Path(__file__).resolve() # Absolute path to this script
    project_root = script_location.parent.parent # .../scripts/ -> .../
    local_path = project_root / "NAAM_VERSION"

    # Now, figure out which one to use.
    if container_path.exists():
        local_version_file = container_path
    elif local_path.exists():
        local_version_file = local_path
    else:
        # If neither file exists, we can't check, so we just skip.
        print("Warning: NAAM_VERSION file not found in standard locations. Skipping update check.", file=sys.stderr)
        return

    local_version_str = local_version_file.read_text().strip()

    # 2. Fetch the latest release from the GitHub API
    api_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/releases/latest"
    try:
        with urllib.request.urlopen(api_url, timeout=5) as response:
            if response.status != 200:
                # Handle cases like private repos or API rate limiting
                return
            data = json.loads(response.read().decode())
            latest_version_str = data['tag_name']

    except Exception:
        # If there's any issue (no internet, API down), just fail silently.
        # The user's workflow should not be blocked by the version check.
        return

    # 3. Compare the versions
    # A simple helper to parse versions like 'v1.2.3' into comparable tuples (1, 2, 3)
    def parse_version(v_str):
        return tuple(map(int, v_str.lstrip('v').split('.')))

    local_v = parse_version(local_version_str)
    latest_v = parse_version(latest_version_str)

    if latest_v > local_v:
        # A newer version is available. Prompt the user.
        # --- CHANGE 1: Add ENDC to the separator lines ---
        print(f"\n{Colors.YELLOW}{'='*70}{Colors.ENDC}", file=sys.stderr)
        print(f"{Colors.YELLOW}  WARNING: A new version of this workflow is available!{Colors.ENDC}", file=sys.stderr)
        print(f"{Colors.YELLOW}{'-'*70}{Colors.ENDC}", file=sys.stderr)
        print(f"  You are using version: {Colors.RED}{local_version_str}{Colors.ENDC}", file=sys.stderr)
        print(f"  Latest version is:     {Colors.GREEN}{latest_version_str}{Colors.ENDC}", file=sys.stderr)
        print(f"  Please update at: https://github.com/{repo_owner}/{repo_name}/releases", file=sys.stderr)
        # --- Also add ENDC to the final separator line ---
        print(f"{Colors.YELLOW}{'='*70}{Colors.ENDC}", file=sys.stderr)
        
        try:
            # --- CHANGE 2: Color the input prompt string itself, ending with ENDC ---
            prompt = f"  Do you want to abort the current setup to update? (y/n): {Colors.ENDC}"
            answer = input(prompt).lower().strip()

            if answer == 'y':
                print(f"\n{Colors.RED}Aborting script. Please pull the latest container version.{Colors.ENDC}", file=sys.stderr)
                sys.exit(1)
            else:
                # --- CHANGE 3: Add ENDC to the final message for a clean exit ---
                print(f"Continuing with the current version...\n{Colors.ENDC}", file=sys.stderr)
        except (EOFError, KeyboardInterrupt):
            # Handle Ctrl+D or Ctrl+C during input as an abort
            print(f"\n\n{Colors.RED}Aborting script.{Colors.ENDC}", file=sys.stderr)
            sys.exit(1)

# ==================================
#        HELPER FUNCTIONS
# ==================================

def load_sample_map(sample_map_path):
    """
    Load a sample map file with automatic delimiter detection.
    Supports CSV/TSV/semicolon/pipe-separated files and Excel files.
    """
    # Detect Excel formats
    excel_ext = (".xlsx", ".xls", ".xlsm", ".xlsb", ".odf", ".ods", ".odt")
    try:
        # Excel formats
        if sample_map_path.lower().endswith(excel_ext):
            try:
                df = pd.read_excel(sample_map_path)
                print(f"Loaded manual sample map (Excel): {sample_map_path}")
            except Exception as e:
                raise ValueError(f"Failed to read Excel file '{sample_map_path}': {e}")
        # Or delimited text formats (CSV, TSV, etc.)
        else:
            try:
                df = pd.read_csv(
                    sample_map_path,
                    sep=None,               # autodetect delimiter
                    engine='python',        # needed for SEP autodetect
                    skipinitialspace=True   # trim spaces after delimiters
                )
                print(f"Loaded manual sample map (text file): {sample_map_path}")
            except Exception as e:
                raise ValueError(f"Failed to detect delimiter or read text file '{sample_map_path}': {e}")
        # Validate required columns
        required = {'barcode_dir', 'virus_id'}
        missing = required - set(df.columns)

        if missing:
            raise ValueError(
                f"Missing required column(s): {', '.join(missing)}.\n"
                f"Your sample map must contain: {', '.join(required)}"
            )
        return df
    except Exception as e:
        print(f"\nERROR while loading sample map '{sample_map_path}':\n{e}\n", file=sys.stderr)
        sys.exit(1)

def extract_fasta_header(fasta_path):
    """Gets the first sequence ID from a FASTA file."""
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith(">"):
                return line[1:].split()[0]
    return None

def auto_generate_sample_map(reads_dir, config, reads_to_test, min_reads, output_csv_path, project_dir):
    """Subsamples reads and uses minimap2 to auto-detect the virus in each barcode."""
    print("\nBuilding reference database for auto-detection...")
    ref_to_virus_id = {}
    master_fasta_path = os.path.join(project_dir, "temp_master_db.fasta")
    
    with open(master_fasta_path, 'w') as master_fasta:
        for virus_id, params in config.items():
            ref_path = params.get('reference_genome')
            if not ref_path or not os.path.isfile(ref_path):
                continue
            
            header = extract_fasta_header(ref_path)
            if header:
                ref_to_virus_id[header] = virus_id
            
            with open(ref_path, 'r') as ref_file:
                master_fasta.write(ref_file.read())
                master_fasta.write("\n")

    results = []
    search_pattern = os.path.join(os.path.abspath(reads_dir), 'barcode*')
    
    print("Scanning barcodes to auto-assign viruses...")
    for barcode_dir in sorted(glob.glob(search_pattern)):
        if not os.path.isdir(barcode_dir):
            continue
            
        barcode_name = os.path.basename(barcode_dir)
        
        fastq_files = glob.glob(os.path.join(barcode_dir, "*.fastq.gz"))
        is_gzipped = True
        if not fastq_files:
            fastq_files = glob.glob(os.path.join(barcode_dir, "*.fastq"))
            is_gzipped = False
            if not fastq_files:
                continue
                
        lines_to_extract = reads_to_test * 4
        file_list_str = " ".join(fastq_files)
        
        if is_gzipped:
            extract_cmd = f"zcat {file_list_str} 2>/dev/null | head -n {lines_to_extract}"
        else:
            extract_cmd = f"cat {file_list_str} 2>/dev/null | head -n {lines_to_extract}"

        minimap_cmd = f"minimap2 -x map-ont --secondary=no {master_fasta_path} - 2>/dev/null"
        full_cmd = f"{extract_cmd} | {minimap_cmd}"
        
        process = subprocess.Popen(full_cmd, shell=True, stdout=subprocess.PIPE, text=True)
        paf_output, _ = process.communicate()
        
        votes = defaultdict(int)
        seen_reads = set()
        
        for line in paf_output.strip().split('\n'):
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) >= 6:
                read_id = cols[0]
                target_header = cols[5]
                
                if read_id not in seen_reads:
                    seen_reads.add(read_id)
                    votes[target_header] += 1

        if votes:
            winning_header = max(votes, key=votes.get)
            winning_votes = votes[winning_header]
            winning_virus_id = ref_to_virus_id.get(winning_header, "Unknown")
            total_votes = sum(votes.values())
            confidence = (winning_votes / total_votes) * 100
            
            if winning_votes >= min_reads:
                print(f"  {barcode_name}: Assigned to -> {winning_virus_id} ({confidence:.1f}% confidence based on {total_votes} unique mapped reads)")
            else:
                print(f"  {barcode_name}: FAILED (Top match '{winning_virus_id}' had only {winning_votes} reads. Threshold is {min_reads})")
                winning_virus_id = "unassigned"
                
            results.append({"barcode_dir": barcode_name, "virus_id": winning_virus_id})
        else:
            print(f"  {barcode_name}: FAILED (No viral reads detected in subsample.)")
            results.append({"barcode_dir": barcode_name, "virus_id": "unassigned"})

    if os.path.exists(master_fasta_path):
        os.remove(master_fasta_path)

    df = pd.DataFrame(results)
    df.to_csv(output_csv_path, index=False)
    print(f"\nAuto-generated sample map saved to: {output_csv_path}")
    return df

# ==================================
#           MAIN LOGIC
# ==================================

# Argument parsing
parser = argparse.ArgumentParser(
    description="Interactive tool for setting up a multi-virus amplicon analysis project.",
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("-p", "--project_dir", help="Project directory path (default: current directory)", default=".") # current directory
parser.add_argument("-n", "--study_name", required=True, help="Name of the study")
parser.add_argument("-d", "--raw_fastq_dir", required=True, help="Directory containing raw FASTQ barcode subdirectories")
parser.add_argument("--virus-config", required=True, help="Path to the virus configuration YAML file.")
parser.add_argument("--sample-map", required=False, help="Path to a manual sample map CSV file. If omitted, it will be auto-generated.")
parser.add_argument("--reads-to-test", type=int, default=2000, help="Number of reads to sample per barcode for auto-detection. (default: 2000)")
parser.add_argument("--min-reads", type=int, default=50, help="Minimum mapped reads required to assign a virus during auto-detection. (default: 50)")

args = parser.parse_args()

check_for_updates(repo_owner="EMC-Viroscience", repo_name="nanopore-amplicon-analysis-manual")


### Project directory path
# Use abspath to ensure it's absolute before using it to join other paths
project_dir = os.path.abspath(args.project_dir)
if not os.path.exists(project_dir):
    try:
        os.makedirs(project_dir)
        print(f"Created project directory: {project_dir}")
    except OSError as e:
        print(f"Error creating project directory '{project_dir}': {e}")
        sys.exit(1)
else:
     # Check if the specified path is actually a directory
    if not os.path.isdir(project_dir):
        print(f"Error: Specified project path '{project_dir}' exists but is not a directory.")
        sys.exit(1)
print(f"Using project directory: {project_dir}")


### Copy Snakemake file and update STUDY_NAME
# --- Robustly find the source Snakefile ---
container_snakefile_path = Path("/snakefile_naam.smk")
# When local, it's in '<project_root>/workflow/snakefile_naam.smk'
script_location = Path(__file__).resolve()
project_root = script_location.parent.parent
local_snakefile_path = project_root / "workflow" / "snakefile_naam.smk"

if container_snakefile_path.exists():
    src_snakemake = container_snakefile_path
elif local_snakefile_path.exists():
    src_snakemake = local_snakefile_path
else:
    print(f"Error: Source Snakemake file not found at '{container_snakefile_path}' or '{local_snakefile_path}'.")
    sys.exit(1)

dest_snakemake = os.path.join(project_dir, "Snakefile")

try:
    shutil.copy2(src_snakemake, dest_snakemake)
    with open(dest_snakemake, 'r') as file:
        filedata = file.read()

    # Replace the STUDY_NAME variable
    study_name = args.study_name
    # Make the replacement more robust (e.g., handle spaces around =)
    filedata = re.sub(r'STUDY_NAME\s*=\s*""', f'STUDY_NAME = "{study_name}"', filedata)
    # filedata = filedata.replace('STUDY_NAME = ""', f'STUDY_NAME = "{study_name}"') # Original less robust way

    with open(dest_snakemake, 'w') as file:
        file.write(filedata)

    print(f"Copied and modified Snakemake file to: {dest_snakemake}")

except Exception as e:
    print(f"Error processing Snakemake file: {e}")
    sys.exit(1)


# --- Load Virus Config YAML ---
try:
    with open(args.virus_config, 'r') as f:
        virus_config = yaml.safe_load(f)
    if not isinstance(virus_config, dict):
        raise TypeError("YAML content is not a dictionary.")
    print(f"Successfully loaded virus config from: {args.virus_config}")
except Exception as e:
    print(f"Error: Failed to load or parse virus config file '{args.virus_config}': {e}", file=sys.stderr)
    sys.exit(1)

# --- Resolve Sample Map (Manual or Auto) ---
if args.sample_map:
    sample_map_df = load_sample_map(args.sample_map)
else:
    print("\nNo manual sample map provided. Generating auto sample map...")
    auto_map_path = os.path.join(project_dir, "sample_map.csv")
    sample_map_df = auto_generate_sample_map(
        reads_dir=args.raw_fastq_dir,
        config=virus_config,
        reads_to_test=args.reads_to_test,
        min_reads=args.min_reads,
        output_csv_path=auto_map_path,
        project_dir=project_dir
    )

print("\nGenerating final sample sheet for Snakemake...")
sample_data = []
unassigned_samples = []
raw_fastq_dir_abs = os.path.abspath(args.raw_fastq_dir)

# Find all potential barcode directories
search_pattern = os.path.join(raw_fastq_dir_abs, 'barcode*')
potential_barcode_paths = glob.glob(search_pattern)

for item_path in potential_barcode_paths:
    if os.path.isdir(item_path):
        barcode_dir_name = os.path.basename(item_path)

        # Extract the base barcode
        number_part = barcode_dir_name[len("barcode"):]
        if not number_part.isdigit():
            print(f"  - Warning: Directory name '{barcode_dir_name}' is not in 'barcode<number>' format. Skipping.")
            continue
        base_barcode = f"BC{int(number_part):02d}"

        # Look up matching rows for this barcode from the sample map
        matching_rows = sample_map_df[sample_map_df['barcode_dir'] == barcode_dir_name]
        
        if matching_rows.empty:
            print(f"  - Warning: Directory '{barcode_dir_name}' found on disk but not in sample map. Skipping.")
            continue

        for _, row in matching_rows.iterrows():
            virus_id = row['virus_id']

            if virus_id not in virus_config:
                if virus_id == "unassigned":
                    print(f"  - Note: '{barcode_dir_name}' did not meet read thresholds. Bypassing pipeline, but will report in final evaluation.")
                else:
                    print(f"  - Warning: virus_id '{virus_id}' not found in virus config. Skipping.")
                
                unassigned_samples.append({"barcode": base_barcode, "virus_id": "Failed / Unassigned"})
                continue
                
            params = virus_config[virus_id]
            unique_id = f"{base_barcode}_{virus_id}"
            sequence_name = f"{args.study_name}_{unique_id}"

            try:
                sample_row = {
                    "unique_id": unique_id,
                    "sequence_name": sequence_name,
                    "fastq_path": item_path,
                    "virus_id": virus_id,
                    "reference_genome": os.path.abspath(params.get('reference_genome')),
                    "primer": os.path.abspath(params.get('primer')),
                    "primer_reference": os.path.abspath(params.get('primer_reference')),
                    "min_length": params.get('min_length'),
                    "coverage": params.get('coverage'),
                    "run_nextclade": params.get('run_nextclade', False), 
                    "nextclade_dataset": params.get('nextclade_dataset'), 
                    "primer_allowed_mismatch": params.get('primer_allowed_mismatch')
                }
                for key in ['reference_genome', 'primer', 'primer_reference']:
                    if not os.path.isfile(sample_row[key]):
                        raise FileNotFoundError(f"File for '{key}' not found at '{sample_row[key]}'")

                sample_data.append(sample_row)
            except (KeyError, TypeError) as e:
                print(f"  - Error: Missing a required key for virus '{virus_id}'. Skipping '{barcode_dir_name}'. Details: {e}", file=sys.stderr)
            except FileNotFoundError as e:
                print(f"  - Error: {e}. Skipping '{barcode_dir_name}'. Please check paths in virus config.", file=sys.stderr)

if not sample_data:
    print("\nWarning: No valid samples were processed. The generated sample.tsv will be empty.", file=sys.stderr)

expected_columns = [
    'unique_id', 'sequence_name', 'fastq_path', 'virus_id', 'reference_genome', 'primer',
    'primer_reference', 'min_length', 'coverage', 'run_nextclade', 'nextclade_dataset',
    'primer_allowed_mismatch'
]
samples_df = pd.DataFrame(sample_data, columns=expected_columns)

if not samples_df.empty:
    samples_df = samples_df.sort_values(by="unique_id").reset_index(drop=True)

samples_tsv_path = os.path.join(project_dir, "sample.tsv")
try:
    samples_df.to_csv(samples_tsv_path, sep='\t', index=False)
    print(f"\nGenerated sample sheet with {len(samples_df)} valid samples: {samples_tsv_path}")

    if unassigned_samples:
        unassigned_df = pd.DataFrame(unassigned_samples)
        unassigned_path = os.path.join(project_dir, "unassigned_samples.tsv")
        unassigned_df.to_csv(unassigned_path, sep='\t', index=False)
        print(f"Saved {len(unassigned_samples)} bypassed/failed samples to: {unassigned_path}")
        
except Exception as e:
    print(f"Error writing sample sheet to {samples_tsv_path}: {e}", file=sys.stderr)
    sys.exit(1)

# Victory lap
print("\nProject setup complete.")
print(f"Project Directory: {project_dir}")
print(f"Snakemake File: {dest_snakemake}")
print(f"Sample Sheet: {samples_tsv_path}")
print("\nYou can now navigate to the project directory and run Snakemake.")
