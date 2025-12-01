## Nanopore amplicon analysis manual (NAAM).

[![Deploy Quarto Book](https://github.com/EMC-Viroscience/nanopore-amplicon-analysis-manual/actions/workflows/deploy.yml/badge.svg)](https://github.com/EMC-Viroscience/nanopore-amplicon-analysis-manual/actions/workflows/deploy.yml)

To view the fully fledged manual, please visit the [EMC-Viroscience website](https://EMC-Viroscience.github.io/workflows.html)!

Interesting files can be found here:

-   **workflow/snakefile_naam.smk**: This is the heart of the workflow.
-   **scripts/amplicon_project.py**: Script to setup your project directory.
-   **scripts/viz_nextclade_cli.R**: Script to visualize Nextclade mutations.
-   **envs/environment.yml**: The conda environment used for the workflow.
-   **envs/imam_workflow.def**: The definition file used for container creation.

## Installation and quick start

1.  **Download pre-built image**
```         
wget https://github.com/EMC-Viroscience/nanopore-amplicon-analysis-manual/releases/latest/download/naam_workflow.sif
```

2.  **Prepare project directory**
```         
singularity exec \
  --bind /mnt/viro0002-data:/mnt/viro0002-data \ # Use whatever --bind is required for your context
  --bind $HOME:$HOME \
  --bind $PWD:$PWD \
  naam_workflow.sif \
  python /amplicon_project.py \
    -p {project.folder} \
    -n {name} \
    -d {reads} \
    --virus-config {virus_config.yaml} \
    --sample-map {sample_map.csv}
```

Example of virus_config.yaml and sample_map.csv is available in config/.

3.  **Run Snakemake workflow**
```         
singularity exec \
  --bind /mnt/viro0002-data:/mnt/viro0002-data \
  --bind $HOME:$HOME \
  --bind $PWD:$PWD \
  naam_workflow.sif \
  snakemake --snakefile Snakefile \
  --cores {threads}
```

**Note**: If you want to run the Snakemake workflow _without_ the singularity container, then you have to install the required Conda environment manually. The **envs/environment.yml** file in this repo will create a Conda environment called 'naam', which has most of the tools required. In addition, you need to install Virconsens and Ampliclip manually into your new conda environment. 

Here are all the steps outlined:
```         
# Create the environment
conda env create -f envs/environment.yml

# Installation of Virconsens and Ampliclip need to be performed with the 'naam' conda environment activated
conda activate naam

git clone https://github.com/dnieuw/Virconsens.git
cd Virconsens
pip install .
cd ..

git clone https://github.com/dnieuw/Ampliclip.git
cd Ampliclip
pip install .

cd ..
```

## Brief workflow summary:

1.  **Quality control**
    - Merging and decompressing with zcat
    - Quality trimming with [fastp](https://github.com/OpenGene/fastp)
    - Primer trimming with [Ampliclip](https://github.com/dnieuw/Ampliclip)

2.  Generating consensus sequences
    - Mapping trimmed reads to a reference with [minimap2](https://github.com/lh3/minimap2)
    - Generate consensus sequences with [Virconsens](https://github.com/dnieuw/Virconsens)

3.  Comparing sequences to references
    - Create multiple sequence alignment between consensus sequences and a reference with minimap2 and [gofasta](https://github.com/virus-evolution/gofasta)
    - Create some read statistics for raw, QC, trimmed and mapped reads. 

4.  Optional: [Nextclade CLI](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/index.html).
