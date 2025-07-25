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

    # 1. Get the local version from the file inside the container
    local_version_file = Path("/NAAM_VERSION")
    if not local_version_file.exists():
        # If the file doesn't exist, we can't check, so we just skip.
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

parser = argparse.ArgumentParser(description="Interactive tool for setting up a Snakemake project.")
parser.add_argument("-p", "--project_dir", help="Project directory path (default: current directory)", default=".") # current directory
parser.add_argument("-n", "--study_name", required=True, help="Name of the study")
parser.add_argument("-d", "--raw_fastq_dir", required=True, help="Directory containing raw FASTQ files")
parser.add_argument("-P", "--primer", required=True, help="Fasta file containing primer sequences")
parser.add_argument("-r", "--primer_reference", required=True, help="Fasta file containing primer reference sequence")
parser.add_argument("-R", "--reference_genome", required=True, help="Fasta file containing reference genome")
parser.add_argument("-m", "--min_length", required=False, help="Minimum read length (default: 1000)", default=1000, type=int)
parser.add_argument("-c", "--coverage", required=False, help="Minimum coverage required for consensus (default: 30)", default=30, type=int)
parser.add_argument("-t", "--threads", required=False, help="Maximum number of threads for the Snakefile (default: 8)", default=8, type=int)
parser.add_argument("--use_sars_cov_2_workflow", action='store_true', help="Add this parameter if you want to analyze SARS-CoV-2 data")
parser.add_argument("--nextclade_dataset", required=False, help="Path to a custom Nextclade dataset directory, OR an official Nextclade dataset name (e.g., 'nextstrain/sars-cov-2/wuhan-hu-1/orfs'). Check official nextclade datasets with `nextclade dataset list`.")

args = parser.parse_args() # reads command from user

check_for_updates(repo_owner="LucvZon", repo_name="nanopore-amplicon-analysis-manual")

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
src_snakemake = "/snakefile_naam.smk" # This path is set for the singularity image
if not os.path.exists(src_snakemake):
    print(f"Error: Source Snakemake file '{src_snakemake}' not found.")
    # Potentially exit here, or provide a way to specify it
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

### Path validations
# Ensure optional paths are stored as absolute paths
ref_genome_path = os.path.abspath(args.reference_genome)
primer_path = os.path.abspath(args.primer)
primer_ref_path = os.path.abspath(args.primer_reference)

# Check existence of optional files if paths were provided
for path, name in [(ref_genome_path, "Reference genome"), (primer_path, "Primer FASTA"), (primer_ref_path, "Primer reference FASTA")]:
    if path and not os.path.isfile(path):
         print(f"Error: Specified {name} file '{path}' does not exist or is not a file.")
         sys.exit(1)

# Handle nextclade dataset
nextclade_dataset_identifier = args.nextclade_dataset
if nextclade_dataset_identifier:
    print(f"Nextclade dataset identifier provided: {nextclade_dataset_identifier}")
    print("  Note: If this is a path to a custom dataset, ensure it's an absolute path.")
    print("  If it's an official dataset name, the workflow will attempt to download it.")


### Set up raw FASTQ files
raw_fastq_dir_arg = args.raw_fastq_dir
# Check existence
if not os.path.exists(raw_fastq_dir_arg):
    print(f"Error: Raw FASTQ directory '{raw_fastq_dir_arg}' does not exist.")
    sys.exit(1)
# Check if it's a directory
if not os.path.isdir(raw_fastq_dir_arg):
     print(f"Error: Provided raw FASTQ path '{raw_fastq_dir_arg}' is not a directory.")
     sys.exit(1)

# Ensure raw_fastq_dir itself is absolute for robust path joining and globbing
raw_fastq_dir_abs = os.path.abspath(raw_fastq_dir_arg)

print(f"Scanning for barcode directories in: {raw_fastq_dir_abs}")

sample_data = []
# Use glob to find items starting with 'barcode' inside the raw fastq directory
search_pattern = os.path.join(raw_fastq_dir_abs, 'barcode*')
potential_barcode_paths = glob.glob(search_pattern)

for item_path in potential_barcode_paths:
    # Check if the found item is actually a directory
    if os.path.isdir(item_path):
        dir_name = os.path.basename(item_path) # e.g., "barcode49"
        # Extract the number part after "barcode"
        if dir_name.startswith("barcode"):
            number_part = dir_name[len("barcode"):] # Get the part after "barcode"
            if number_part.isdigit(): # Check if it's purely numeric
                barcode_num = int(number_part)
                # Format the unique_id with leading zero if needed (e.g., BC01, BC49)
                unique_id = f"BC{barcode_num:02d}"
                # The fastq_path should be the absolute path to the barcode directory
                fastq_path = item_path # glob with absolute base returns absolute paths
                # Format sequence name {study_name}_{unique_id}  
                sequence_name = f"{args.study_name}_{unique_id}"
                sample_data.append({
                    "unique_id": unique_id,
                    "sequence_name": sequence_name,
                    "fastq_path": fastq_path,
                    "reference": ref_genome_path,
                    "primers": primer_path,
                    "primer_reference": primer_ref_path,
                    "coverage": args.coverage,
                    "min_length": args.min_length,
                    "nextclade_db": nextclade_dataset_identifier
                })
            else:
                 print(f"Warning: Directory '{dir_name}' in '{raw_fastq_dir_abs}' starts with 'barcode' but is not followed by only digits. Skipping.")
        # No else needed as glob already filtered for 'barcode*'

# Check if any valid barcode directories were found
if not sample_data:
    print(f"Warning: No valid barcode directories (e.g., 'barcode1', 'barcode49') found in '{raw_fastq_dir_abs}'.")
    print("The generated sample.tsv file will be empty or contain only headers.")
    # Depending on Snakemake requirements, you might want to exit here:
    # print("Exiting because no samples were found.")
    # sys.exit(1)

# Create a pandas DataFrame
samples_df = pd.DataFrame(sample_data)

# Sort by UniqueID for consistency (optional but nice)
if not samples_df.empty:
    samples_df = samples_df.sort_values(by="unique_id").reset_index(drop=True)

# Define the output path for sample.tsv within the project directory
samples_tsv_path = os.path.join(project_dir, "sample.tsv")

# Write the DataFrame to a TSV file
try:
    samples_df.to_csv(samples_tsv_path, sep='\t', index=False)
    print(f"Generated sample sheet: {samples_tsv_path}")
    if samples_df.empty:
        print("Note: The sample sheet is empty as no valid barcode directories were found.")
except Exception as e:
    print(f"Error writing sample sheet to {samples_tsv_path}: {e}")
    sys.exit(1)

### Set numeric values
min_length = args.min_length
coverage = args.coverage
threads = args.threads

### Need to expand this with "use_nextclade"
config_data = {
    "threads": threads,
	"use_sars_cov_2_workflow": args.use_sars_cov_2_workflow # Either 'True' or 'False'. 
}

# Add nextclade to config.yaml
if nextclade_dataset_identifier:  # If --nextclade_dataset argument was provided
    config_data["run_nextclade"] = True
else: # If --nextclade_dataset argument was NOT provided
    config_data["run_nextclade"] = False

# Define config file path within project directory
yaml_file = os.path.join(project_dir, "config.yaml")

# Write the config data to YAML
try:
    with open(yaml_file, "w") as outfile:
        # Use sort_keys=False to maintain the order defined in the dictionary
        yaml.dump(config_data, outfile, default_flow_style=False, sort_keys=False)
    print(f"Created configuration file: {yaml_file}")
except Exception as e:
    print(f"Error writing configuration file {yaml_file}: {e}")
    sys.exit(1)

# Victory lap
print("\nProject setup complete.")
print(f"Project Directory: {project_dir}")
print(f"Snakemake File: {dest_snakemake}")
print(f"Sample Sheet: {samples_tsv_path}")
print(f"Config File: {yaml_file}")
print("\nYou can now navigate to the project directory and run Snakemake.")
