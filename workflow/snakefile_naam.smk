# Generic workflow for amplicon based sequencing


STUDY_NAME = ""


import sys, string, shutil, glob
from pandas import read_table
from collections import defaultdict
import pandas as pd
import os

# Check for an environment variable that is ONLY set inside the Singularity container.
# If it exists, we are in the container. Otherwise, we are running locally.
if "SINGULARITY_NAME" in os.environ:
    # We are in the container, scripts are in the root.
    SCRIPT_PATH = "/" 
else:
    # We are running locally, scripts are in the 'scripts' sub-directory.
    SCRIPT_PATH = "scripts/"

# Load in sample data
sample_data = pd.read_csv("sample.tsv", sep="\t")
SAMPLES = sample_data["unique_id"].tolist()

# Get a list of all unique viruses in this run
VIRUSES = sample_data["virus_id"].unique().tolist()

sample_data = sample_data.set_index('unique_id')
sample_data = sample_data.to_dict("index")


# Helper functions:
def get_virus_param(wildcards, param_name):
    """
    Finds the first sample belonging to a virus group and returns a specific parameter for it.
    """
    virus = wildcards.virus
    # Use next() with a generator for efficiency. It finds the first match.
    sample_for_virus = next(s for s, data in sample_data.items() if data['virus_id'] == virus)
    # Return the requested parameter from that sample's data.
    return sample_data[sample_for_virus][param_name]

def get_required_nextclade_datasets():
    """
    Scans the sample sheet and returns a unique list of all nextclade_dataset
    identifiers that are needed for this run and are not local paths.
    """
    datasets = set()
    for sample, data in sample_data.items():
        # Only consider samples where nextclade is set to run
        if data.get('run_nextclade', False):
            dataset_id = data.get('nextclade_dataset')
            # Add to the set if it's not None and doesn't look like an absolute/existing path
            if dataset_id and not os.path.isdir(dataset_id):
                datasets.add(dataset_id)
    return list(datasets)

def get_nextclade_dataset_path(wildcards):
    """
    Determines the correct path for the nextclade dataset for a given virus.
    If it's a local path in the config, it returns that.
    If it's an official name, it returns the path where it WILL BE downloaded.
    """
    dataset_id = get_virus_param(wildcards, 'nextclade_dataset')
    if os.path.isdir(dataset_id):
        return dataset_id # It's a pre-existing local path
    else:
        # It's an official name, so return the path where we will download it.
        # We replace '/' with '_' to create a valid directory name.
        safe_name = dataset_id.replace('/', '_')
        return f"result/nextclade_downloads/{safe_name}"


# The 'rule all' needs to request the final outputs based on the chosen workflow
def get_final_outputs():
    outputs = []
    required_dataset_paths = []
    for virus in VIRUSES:
        # Aggragted alignment for this virus
        outputs.append(f"result/{virus}/alignment/all_consensus_aligned.fasta")
        outputs.append(f"result/{virus}/variants/all_variants.tsv")
    
        if any(d.get('run_nextclade', False) for s, d in sample_data.items() if d['virus_id'] == virus):
             outputs.append(f"result/{virus}/nextclade/success.txt")
             # Determine the required dataset path for this virus and add it to our list
             dataset_id = next(d['nextclade_dataset'] for s, d in sample_data.items() if d['virus_id'] == virus and d.get('run_nextclade'))
             if not os.path.isdir(dataset_id):
                 safe_name = dataset_id.replace('/', '_')
                 required_dataset_paths.append(f"result/nextclade_downloads/{safe_name}")

    # Non-grouped outputs
    outputs.append("result/readstats/raw.tsv")
    outputs.append("result/readstats/QC.tsv")
    outputs.append("result/readstats/trimmed.tsv")
    outputs.append("result/readstats/mapped.tsv")
    
    outputs.extend(list(set(required_dataset_paths)))
    # Ensure return is correctly indented at the function level
    return outputs

rule all:
    input:
        get_final_outputs()

rule merge_barcodes:
    input:
        lambda wildcards: sample_data[wildcards.sample_id]["fastq_path"]
    output:
        "result/raw/{sample_id}.fastq"
    threads: 1
    shell:
        """
        zcat {input}/*.fastq.gz > {output}
        """

rule QC_Nanopore_reads:
    input:
        "result/raw/{sample_id}.fastq"
    output:
        fastq="result/QC/{sample_id}.fastq",
        report="result/QC/{sample_id}.html"
    log: "logs/{sample_id}_fastp.log"
    threads: 2
    params:
        min_length=lambda wildcards: sample_data[wildcards.sample_id]["min_length"]
    shell:
        """
        fastp -i {input} -o {output.fastq} -j /dev/null -h {output.report} \
        --disable_trim_poly_g \
        --disable_adapter_trimming \
        --qualified_quality_phred 10 \
        --unqualified_percent_limit 50 \
        --length_required {params.min_length} \
        -w {threads} > {log} 2>&1
        """

rule map_to_primer_reference:
	input:
		fastq="result/QC/{sample_id}.fastq"
	output:
		mapped="result/trimmed/{sample_id}_mapped.bam"
	params:
		reference=lambda wildcards: sample_data[wildcards.sample_id]["primer_reference"]
	threads: 4
	shell:
		"""
		minimap2 -Y -t {threads} -x map-ont -a {params.reference} {input.fastq} 2> /dev/null | samtools view -bF 4 - | samtools sort -@ {threads} - > {output.mapped}
		"""

rule trim_primers:
    input:
        mapped="result/trimmed/{sample_id}_mapped.bam"
    output:
        clipped="result/trimmed/{sample_id}_clipped.bam",
        trimmed="result/trimmed/{sample_id}_trimmed.fastq"
    log:
        "logs/{sample_id}_ampliclip.log"
    params:
        reference=lambda wildcards: sample_data[wildcards.sample_id]["primer_reference"],
        primers=lambda wildcards: sample_data[wildcards.sample_id]["primer"],
        min_length=lambda wildcards: sample_data[wildcards.sample_id]["min_length"],
        mismatch=lambda wildcards: sample_data[wildcards.sample_id]["primer_allowed_mismatch"]
    threads: 1
    shell:
        """
        samtools index {input.mapped}

        ampliclip \
        --infile {input.mapped} \
        --outfile {output.clipped}_ \
        --outfastq {output.trimmed} \
        --primerfile {params.primers} \
        --referencefile {params.reference}\
        -fwd LEFT -rev RIGHT \
        --padding 20 --mismatch {params.mismatch} --minlength {params.min_length} > {log} 2>&1

        samtools sort {output.clipped}_ > {output.clipped}
        rm {output.clipped}_
        """

rule map_to_reference:
    input:
        fastq="result/trimmed/{sample_id}_trimmed.fastq"
    output:
        bam="result/mapped/{sample_id}_mapped.bam",
        filtered="result/filtered/{sample_id}_filtered.fastq"
    threads: 8
    params:
        reference=lambda wildcards: sample_data[wildcards.sample_id]["reference_genome"]
    shell:
        """
        minimap2 -Y -t {threads} -x map-ont -a {params.reference} {input.fastq} 2> /dev/null | samtools view -bF 4 - | samtools sort -@ {threads} - > {output.bam}
        samtools index -@ {threads} {output.bam}
        samtools fastq {output.bam} 2> /dev/null > {output.filtered}
        """

### CONSENSUS CALLING ###
rule create_consensus:
    input:
        "result/mapped/{sample_id}_mapped.bam"
    output:
        fa="result/consensus/{sample_id}_consensus.fasta",
        vf_raw=temp("result/variants/{sample_id}_variants_raw.tsv")
    threads: 8
    params:
        reference=lambda wildcards: sample_data[wildcards.sample_id]["reference_genome"],
        coverage=lambda wildcards: sample_data[wildcards.sample_id]["coverage"],
        sequence_name=lambda wildcards: sample_data[wildcards.sample_id]["sequence_name"]
    shell:
        """
        samtools index {input}

        virconsens \
        -b {input} \
        -o {output.fa} \
        -vf {output.vf_raw} \
        -n {params.sequence_name} \
        -r {params.reference} \
        -d {params.coverage} \
        -af 0.1 \
        -c {threads}
        """

rule add_samplename_to_variants:
    input:
        vf_raw="result/variants/{sample_id}_variants_raw.tsv"
    output:
        vf="result/variants/{sample_id}_variants.tsv"
    params:
        sequence_name=lambda wildcards: sample_data[wildcards.sample_id]["sequence_name"]
    shell:
        """
        awk -v sample_name="{params.sequence_name}" \
            'BEGIN {{ OFS="\\t" }} \
            NR==1 {{ print "sample_name", $0; next }} \
            $3 != $4 {{ print sample_name, $0 }}' \
            {input.vf_raw} > {output.vf}
        """

# Helper function to get all VARIANTS for a specific virus
def get_variants_for_virus(wildcards):
    samples_for_this_virus = [s for s, data in sample_data.items() if data['virus_id'] == wildcards.virus]
    return expand("result/variants/{sample_id}_variants.tsv", sample_id=samples_for_this_virus)

rule aggregate_variants:
    input:
        #expand("result/variants/{sample_id}_variants.tsv", sample_id=SAMPLES)
        get_variants_for_virus
    output:
        all_vf="result/{virus}/variants/all_variants.tsv"
    threads: 1
    shell:
        """
        # Check if there are any input files to prevent errors with head/tail on empty list
        if [ -n "{input}" ]; then
            # Print header from the first file
            head -n 1 $(echo {input} | cut -d' ' -f1) > {output.all_vf}
            # Append content (without header) from all files
            for f in {input}; do
                tail -n +2 "$f" >> {output.all_vf}
            done
        else
            # Create an empty file or a file with just a header if no inputs
            echo "seqName\tref_pos\tnum_aln\tREF\tALT\tALT_count\tALT_AF" > {output.all_vf}
        fi
        """

# Helper function to get all CONSENSUS files for a specific virus
def get_consensus_for_virus(wildcards):
    samples_for_this_virus = [s for s, data in sample_data.items() if data['virus_id'] == wildcards.virus]
    return expand("result/consensus/{sample_id}_consensus.fasta", sample_id=samples_for_this_virus)

rule aggregate_consensus:
    input:
        #expand("result/consensus/{sample_id}_consensus.fasta",sample_id=SAMPLES)
        get_consensus_for_virus
    output:
        "result/{virus}/alignment/all_consensus.fasta"
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule prepare_nextclade_dataset:
    output:
        # The output is a directory that will contain the dataset.
        # The wildcard 'dataset_path' will match the safe name, e.g., 'nextstrain_sars-cov-2_wuhan-hu-1_orfs'
        directory("result/nextclade_downloads/{dataset_path}")
    params:
        # The input to this rule is the original, unsafe dataset name with slashes
        # We reverse the transformation to get the original name.
        original_name=lambda wildcards: wildcards.dataset_path.replace('_', '/')
    log:
        "logs/nextclade_download_{dataset_path}.log"
    shell:
        """
        echo "Preparing Nextclade dataset: {params.original_name}" > {log}
        # This rule only runs for official datasets. The `rule all` logic prevents it from
        # running for local paths. Therefore, we can just run the download command.
        nextclade dataset get --name '{params.original_name}' --output-dir '{output}' >> {log} 2>&1
        """


rule run_nextclade:
    input:
        consensus="result/{virus}/alignment/all_consensus.fasta",
        dataset_dir=lambda wildcards: get_nextclade_dataset_path(wildcards)
    output:
        outdir=directory("result/{virus}/nextclade/"),
        json="result/{virus}/nextclade/nextclade.json"
    shell:
       """
       echo "Running Nextclade for {wildcards.virus} using dataset at: {input.dataset_dir}"
       nextclade run \
       --input-dataset {input.dataset_dir} \
       --output-all={output.outdir} \
       --output-json={output.json} \
       {input.consensus}
       """

rule visualize_mutation_table:
    input:
        nextclade_json="result/{virus}/nextclade/nextclade.json",
        variants_tsv="result/{virus}/variants/all_variants.tsv"
    output:
        success="result/{virus}/nextclade/success.txt"
    log:
        "logs/{virus}_visualize_mutation_table.log"
    params:
        r_script=os.path.join(SCRIPT_PATH, "viz_nextclade_cli.R"),
        nextclade_input_dir="result/{virus}/nextclade",
        plotly_output_dir="result/{virus}/nextclade/plotly",
        ggplotly_output_dir="result/{virus}/nextclade/ggplotly"
    shell:
        """
        set -euo pipefail

        Rscript {params.r_script} \
        --nextclade-input-dir {params.nextclade_input_dir} \
        --json-file {input.nextclade_json} \
        --virconsens-variants {input.variants_tsv} \
        --plotly-output-dir {params.plotly_output_dir} \
        --ggplotly-output-dir {params.ggplotly_output_dir} 2>&1 | tee {log}
        """


rule align_consensus:
    input:
        "result/{virus}/alignment/all_consensus.fasta"
    output:
        "result/{virus}/alignment/all_consensus_aligned.fasta"
    threads: 8
    ##
    log:
        "logs/{virus}_align_consensus.log"
    params:
        reference=lambda wildcards: get_virus_param(wildcards, 'reference_genome')
    shell:
        """
        minimap2 -t {threads} -a -x asm20 --sam-hit-only --secondary=no --score-N=0 {params.reference} {input} -o result/{wildcards.virus}/alignment/all_consensus_aligned.sam > {log} 2>&1
        gofasta sam toMultiAlign -s result/{wildcards.virus}/alignment/all_consensus_aligned.sam -o {output}
        """

### STATS GENERATION ###
rule generate_readstats_raw:
    input:
        expand("result/raw/{sample_id}.fastq", sample_id=SAMPLES)
    output:
        "result/readstats/raw.tsv"
    threads: 1
    shell:
        """
        seqkit stats -T {input} > {output}
        """

rule generate_readstats_QC:
    input:
        expand("result/QC/{sample_id}.fastq", sample_id=SAMPLES)
    output:
        "result/readstats/QC.tsv"
    threads: 1
    shell:
        """
        seqkit stats -T {input} > {output}
        """

rule generate_readstats_trimmed:
    input:
        expand("result/trimmed/{sample_id}_trimmed.fastq", sample_id=SAMPLES)
    output:
        "result/readstats/trimmed.tsv"
    threads: 1
    shell:
        """
        seqkit stats -T {input} > {output}
        """

rule generate_readstats_mapped:
    input:
        expand("result/mapped/{sample_id}_mapped.bam", sample_id=SAMPLES)
    output:
        "result/readstats/mapped.tsv"
    threads: 1
    shell:
        """
        for file in {input}; do
            samtools fastq $file | seqkit stats -T --stdin-label $file
        done > {output}
        """

### DESTROY .snakemake/ AFTER THE WORKFLOW HAS SUCCESFULLY FINISHED
onsuccess:
    """
    This code runs only after the entire workflow completes successfully.
    It launches a detached, background process to remove the .snakemake
    directory a few seconds after the main workflow process exits.
    This avoids the "file handle in use" race condition.
    """
    import subprocess
    import shlex

    print("Workflow finished successfully.")
    print("Spawning a detached cleanup process to remove .snakemake directory in a few seconds...")

    cleanup_script = "import time, shutil; time.sleep(2); shutil.rmtree('.snakemake', ignore_errors=True)"
    command = shlex.split(f'{sys.executable} -c "{cleanup_script}"')
    
    subprocess.Popen(
        command,
        start_new_session=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    
onerror:
    """
    Optional: This code runs if the workflow fails at any point.
    Useful for explicitly stating that cleanup WON'T happen.
    """
    print("Workflow failed. The .snakemake directory will NOT be removed for debugging.")
