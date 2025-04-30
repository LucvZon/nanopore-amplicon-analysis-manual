#!/usr/bin/env python

import os, sys
import argparse
import shutil
import pandas as pd
import glob
import yaml
import re

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

args = parser.parse_args() # reads command from user

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
src_snakemake = "/home/r050834/scripts/workflows/snakefile_naam.smk" # Temporary hard coded path
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

### Maybe introduce checks here?
# Ensure optional paths are stored as absolute paths
ref_genome_path = os.path.abspath(args.reference_genome)
primer_path = os.path.abspath(args.primer)
primer_ref_path = os.path.abspath(args.primer_reference)

# Check existence of optional files if paths were provided
for path, name in [(ref_genome_path, "Reference genome"), (primer_path, "Primer FASTA"), (primer_ref_path, "Primer reference FASTA")]:
    if path and not os.path.isfile(path):
         print(f"Error: Specified {name} file '{path}' does not exist or is not a file.")
         sys.exit(1)


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
                sample_data.append({"unique_id": unique_id, "sequence_name": sequence_name,"fastq_path": fastq_path, "reference": ref_genome_path, "primers": primer_path,
                                    "primer_reference": primer_ref_path, "coverage": args.coverage, "min_length": args.min_length})
                # print(f"Found: {unique_id} -> {fastq_path}") # Optional: for debugging
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

config_data = {
    "threads": threads,
	"use_sars_cov_2_workflow": args.use_sars_cov_2_workflow # Either 'True' or 'False'. 
}

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