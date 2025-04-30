# Generic workflow for amplicon based sequencing


STUDY_NAME = ""


import sys, string, shutil, glob
from pandas import read_table
from collections import defaultdict
import pandas as pd

configfile: "config.yaml"

# Load in sample data
sample_data = pd.read_csv("sample.tsv", sep="\t")
SAMPLES = sample_data["unique_id"].tolist()

ALIGNMENT_REFERENCE = sample_data.iloc[0]["reference"]

sample_data = sample_data.set_index('unique_id')
sample_data = sample_data.to_dict("index")

# The 'rule all' needs to request the final outputs based on the chosen workflow
def get_final_outputs(wildcards):
    outputs = []
    # Standard outputs - just append the file path strings
    outputs.append("result/alignment/all_consensus_aligned.fasta") # Removed keyword= and comma
    outputs.append("result/readstats/raw.tsv")
    outputs.append("result/readstats/QC.tsv")
    outputs.append("result/readstats/trimmed.tsv")
    outputs.append("result/readstats/mapped.tsv")

    if config.get("use_sars_cov_2_workflow", False):
        # Use consistent path 'result/' and correct append syntax
        outputs.append("result/pangolin.csv")

    # Ensure return is correctly indented at the function level
    return outputs

rule all:
    input:
        get_final_outputs

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
        primers=lambda wildcards: sample_data[wildcards.sample_id]["primers"],
        min_length=lambda wildcards: sample_data[wildcards.sample_id]["min_length"]
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
        --padding 20 --mismatch 2 --minlength {params.min_length} > {log} 2>&1

        samtools sort {output.clipped}_ > {output.clipped}
        rm {output.clipped}_
        """

rule map_to_reference:
    input:
        fastq="result/trimmed/{sample_id}_trimmed.fastq"
    output:
        "result/mapped/{sample_id}_mapped.bam"
    threads: 8
    params:
        reference=lambda wildcards: sample_data[wildcards.sample_id]["reference"]
    shell:
        """
        minimap2 -Y -t {threads} -x map-ont -a {params.reference} {input.fastq} 2> /dev/null | samtools view -bF 4 - | samtools sort -@ {threads} - > {output}
        """

### CONSENSUS CALLING ###
rule create_consensus:
    input:
        input="result/mapped/{sample_id}_mapped.bam"
    output:
        "result/consensus/{sample_id}_consensus.fasta"
    threads: 8
    params:
        reference=lambda wildcards: sample_data[wildcards.sample_id]["reference"],
        coverage=lambda wildcards: sample_data[wildcards.sample_id]["coverage"],
        sequence_name=lambda wildcards: sample_data[wildcards.sample_id]["sequence_name"]
    shell:
        """
        samtools index {input}

        virconsens \
        -b {input} \
        -o {output} \
        -n {params.sequence_name} \
        -r {params.reference} \
        -d {params.coverage} \
        -af 0.1 \
        -c {threads}
        """

rule aggregate_consensus:
    input:
        expand("result/consensus/{sample_id}_consensus.fasta",sample_id=SAMPLES)
    output:
        "result/alignment/all_consensus.fasta"
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

if config.get("use_sars_cov_2_workflow", False):
    rule run_pangolin:
        input:
            "result/alignment/all_consensus.fasta"
        output:
            "result/pangolin.csv"
        threads: 1
        shell:
            """
            pangolin {input} --outfile {output}
            """
		
rule align_consensus:
    input:
        "result/alignment/all_consensus.fasta"
    output:
        "result/alignment/all_consensus_aligned.fasta"
    threads: 8
    ##
    log:
        "logs/align_consensus.log"
    params:
        reference=ALIGNMENT_REFERENCE
    shell:
        """
        minimap2 -t {threads} -a -x asm20 --sam-hit-only --secondary=no --score-N=0 {params.reference} {input} -o result/alignment/all_consensus_aligned.sam > {log} 2>&1
        gofasta sam toMultiAlign -s result/alignment/all_consensus_aligned.sam -o {output}
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
            samtools fastq $file | seqkit stats -T --stdin-label $file | tail -1
        done > {output}
        """

### DESTROY .snakemake/ AFTER THE WORKFLOW HAS SUCCESFULLY FINISHED
onsuccess:
    """
    This code runs only after the entire workflow completes successfully.
    It removes the .snakemake directory.
    """
    print("Workflow finished successfully.")
    snakemake_dir = ".snakemake" # The directory Snakemake creates

    if os.path.exists(snakemake_dir):
        try:
            print(f"Attempting to remove {snakemake_dir} directory...")
            shutil.rmtree(snakemake_dir)
            print(f"Successfully removed {snakemake_dir}.")
        except OSError as e:
            print(f"Error removing {snakemake_dir}: {e}")
    else:
        print(f"{snakemake_dir} directory not found. Skipping removal.")

onerror:
    """
    Optional: This code runs if the workflow fails at any point.
    Useful for explicitly stating that cleanup WON'T happen.
    """
    print("Workflow failed. The .snakemake directory will NOT be removed for debugging.")
