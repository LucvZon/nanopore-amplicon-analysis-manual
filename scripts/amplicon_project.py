#!/usr/bin/env python

import os, sys
import argparse
import shutil
import pandas as pd
import glob
import yaml
import re
import json # version checker
import urllib.request # version checker
from pathlib import Path # version checker

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
                print(f"Loaded sample map (Excel): {sample_map_path}")
            except Exception as e:
                raise ValueError(
                    f"Failed to read Excel file '{sample_map_path}': {e}"
                )
        # Or delimited text formats (CSV, TSV, etc.)
        else:
            try:
                df = pd.read_csv(
                    sample_map_path,
                    sep=None,               # autodetect delimiter
                    engine='python',        # needed for SEP autodetect
                    skipinitialspace=True   # trim spaces after delimiters
                )
                print(f"Loaded sample map (text file): {sample_map_path}")
            except Exception as e:
                raise ValueError(
                    f"Failed to detect delimiter or read text file '{sample_map_path}': {e}"
                )
        # Validate required columns
        required = {'barcode_dir', 'virus_id'}
        missing = required - set(df.columns)

        if missing:
            raise ValueError(
                f"Missing required column(s): {', '.join(missing)}.\n"
                f"Your sample map must contain: {', '.join(required)}"
            )
        # Use barcode_dir as index
        df.set_index('barcode_dir', inplace=True)
        return df
    except Exception as e:
        print(f"\nERROR while loading sample map '{sample_map_path}':\n{e}\n", file=sys.stderr)
        sys.exit(1)

# Argument parsing
parser = argparse.ArgumentParser(
    description="Interactive tool for setting up a multi-virus amplicon analysis project.",
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("-p", "--project_dir", help="Project directory path (default: current directory)", default=".") # current directory
parser.add_argument("-n", "--study_name", required=True, help="Name of the study")
parser.add_argument("-d", "--raw_fastq_dir", required=True, help="Directory containing raw FASTQ barcode subdirectories")
parser.add_argument("--virus-config", required=True, help="Path to the virus configuration YAML file.")
parser.add_argument("--sample-map", required=True, help="Path to the sample map CSV file.\n(CSV must have 'barcode_dir' and 'virus_id' columns)")

args = parser.parse_args() # reads command from user

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

# --- Load Sample Map CSV ---
sample_map_df = load_sample_map(args.sample_map)

print("\nScanning for barcode directories and generating sample sheet...")
sample_data = []
raw_fastq_dir_abs = os.path.abspath(args.raw_fastq_dir)


# Find all potential barcode directories
search_pattern = os.path.join(raw_fastq_dir_abs, 'barcode*')
potential_barcode_paths = glob.glob(search_pattern)

for item_path in potential_barcode_paths:
    if os.path.isdir(item_path):
        barcode_dir_name = os.path.basename(item_path)

        # 1. Look up the virus_id from the sample map
        if barcode_dir_name not in sample_map_df.index:
            print(f"  - Warning: Directory '{barcode_dir_name}' found on disk but not in sample map. Skipping.")
            continue
        virus_id = sample_map_df.loc[barcode_dir_name, 'virus_id']

        # 2. Get the parameters for this virus from the virus config
        if virus_id not in virus_config:
            print(f"  - Warning: virus_id '{virus_id}' (from sample map for '{barcode_dir_name}') not found in virus config. Skipping.")
            continue
        params = virus_config[virus_id]

        # 3. Create the unique ID and sequence name
        number_part = barcode_dir_name[len("barcode"):]
        if not number_part.isdigit():
            print(f"  - Warning: Directory name '{barcode_dir_name}' is not in 'barcode<number>' format. Skipping.")
            continue
        
        unique_id = f"BC{int(number_part):02d}"
        sequence_name = f"{args.study_name}_{unique_id}"

        # 4. Build the data dictionary for this sample, making paths absolute
        try:
            sample_row = {
                "unique_id": unique_id,
                "sequence_name": sequence_name,
                "fastq_path": item_path,
                "virus_id": virus_id,
                # Use .get() for safety and make all paths absolute for robustness
                "reference_genome": os.path.abspath(params.get('reference_genome')),
                "primer": os.path.abspath(params.get('primer')),
                "primer_reference": os.path.abspath(params.get('primer_reference')),
                "min_length": params.get('min_length'),
                "coverage": params.get('coverage'),
                "run_nextclade": params.get('run_nextclade', False), # Default to False if not specified
                "nextclade_dataset": params.get('nextclade_dataset'), # Can be None
                "primer_allowed_mismatch": params.get('primer_allowed_mismatch')
            }
            # Final check that essential file paths exist
            for key in ['reference_genome', 'primer', 'primer_reference']:
                if not os.path.isfile(sample_row[key]):
                    raise FileNotFoundError(f"File for '{key}' not found at '{sample_row[key]}'")

            sample_data.append(sample_row)
        except (KeyError, TypeError) as e:
            print(f"  - Error: Missing a required key (e.g., 'reference_genome') for virus '{virus_id}' in virus config. Skipping sample '{barcode_dir_name}'. Details: {e}", file=sys.stderr)
        except FileNotFoundError as e:
            print(f"  - Error: {e}. Skipping sample '{barcode_dir_name}'. Please check paths in virus config.", file=sys.stderr)

# --- Create and Save the DataFrame ---
if not sample_data:
    print("\nWarning: No valid samples were processed. The generated sample.tsv will be empty.", file=sys.stderr)

# Define column order for consistency
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
    print(f"\nGenerated sample sheet with {len(samples_df)} samples: {samples_tsv_path}")
except Exception as e:
    print(f"Error writing sample sheet to {samples_tsv_path}: {e}", file=sys.stderr)
    sys.exit(1)


# Victory lap
print("\nProject setup complete.")
print(f"Project Directory: {project_dir}")
print(f"Snakemake File: {dest_snakemake}")
print(f"Sample Sheet: {samples_tsv_path}")
print("\nYou can now navigate to the project directory and run Snakemake.")
