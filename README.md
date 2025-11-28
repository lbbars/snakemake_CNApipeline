# Copy Number Alteration Analysis Pipeline
A Snakemake workflow for ONT circulating tumor DNA (ctDNA) analysis

This pipeline processes low-pass Whole Genome Sequencing with Oxford Nanopore to estimate copy-number alterations (CNAs) and tumor fraction from cell-free DNA (cfDNA). It performs:

1. Basecalling using Dorado (GPU)

2. Alignment & QC (cramino, pycoQC, samtools)

3. Genome binning (QDNAseq)

4. Segmentation & CN inference (ACE)

5. Tumor fraction & CN confidence estimation (ichorCNA)

## Pipeline overview
```scss
Raw POD5 files
        ↓ (GPU basecalling: Dorado)
Basecalled FASTQ
        ↓ (alignment + QC)
Aligned BAM + QC reports
        ↓ (binning: QDNAseq)
Binned copy number tracks
        ↓ (segmentation: ACE)
Segment-level CN estimates
        ↓ (ichorCNA modeling)
Tumor fraction & CNA profiles
```

## Repository Layout
```bash
snakemake_CNApipeline/
├── Snakefile
├── config.yaml     # Pipeline configuration (edit this)
├── scripts/        # R analysis scripts
├── extdata/        # Reference genomes, models & panels
├── data/           # Input POD5 files go here
└── results/        # Outputs generated here

```


# 1. Download/clone pipeline repo:
```bash
git clone https://gitlab.com/hazarsandakly/snakemake_cnapipeline
cd snakemake_CNApipeline
```
# Place Required External Data:
dorado model, reference genome and panel of normals
- Add them to this directory:
```bash
mkdir path/to/snakemake_CNApipeline/extdata
```
Place the following inside extdata:

|Required	|Example Location|
|-----------|----------------|
|Dorado model	  |extdata/dorado_model/…|
|Reference genome |(FASTA)	extdata/reference/…|
|Panel of normals | (RDS)	extdata/|

The expected paths are already defined in config.yaml, so no edits are needed if you follow the folder layout.

## Environment setup
### Build Tool Containers (Apptainer/Singularity)
Activate apptainer on your environment: 
```bash
# For genotoul cluster
module load containers/Apptainer
```
Each tool is built into a reproducible .sif image:
```bash
cd path/to/snakemake_CNApipeline

make sif TOOL=toolname VER=XXX
```
The required tools are the following:
- Dorado for basecalling: 
```bash
make sif TOOL=dorado VER=0.9.1
```
- Quality control checks with cramino and pycoQC:
```bash
make sif TOOL=cramino VER=1.1.0

make sif TOOL=pycoqc VER=2.5.2
```
- 3 tools for CNA analysis:
```bash
make sif TOOL=ichorCNA VER=0.2.0

make sif TOOL=QDNAseq VER=1.9.2

make sif TOOL=ACE VER=1.27.1
```
This only needs to be done once.

## Load snakemake on cluster:
```bash
module load bioinfo/Snakemake
```
## Create a SLURM profile for snakemake (Required):
```bash
mkdir -p ~/.config/snakemake/profiles/default
cd ~/.config/snakemake/profiles/default
# Create file called config.yaml
nano config.yaml
```
Inside the config.yaml, we will specify the parameters used by snakemake to execute slurm jobs. This is a working example of this file: 

```bash
executor: slurm

executor-settings:
  default:
    # Add the slurm partition you are working on
    slurm_partition: workq
    slurm_extra: "--mem=64G"

# These ressources are specific to the basecalling rule which requires gpu access
set-resources: 
  basecalling:
    # Change the slurm partition you are working on for gpu access: 
    slurm_partition: gpuq
    # Change the gpu nodes you are going to use for basecalling: 
    slurm_extra: '"--gres=gpu:nvidia_a100" "--mem=128G"' 

use-apptainer: true
# Bind the directories you're working in, so they can be identifiable by the containers:
apptainer-args: "--nv --bind /home/username:/home/username --bind /work:/work --writable-tmpfs"

```
# Data placement:
For easier execution, all data files should be in the same repo as the pipeline. For example, the control pod5 files should be added to 
```bash
path/to/snakemake_CNApipeline/data/control
```
These control samples will be used to generate a panel of normals (PoN). It is used to normalize noise and accurately adjust to low pass sequencing.

As for the tumor samples, they should be added to: 
```bash
path/to/snakemake_CNApipeline/data/samples
```
- These paths should be specified in the "config.yaml" present in the pipeline's directory:
```yaml
# Input samples
samples_dir: "path/to/snakemake_CNApipeline/data/samples"

# Directory for Panel of Normals (PoN) creation
# If this is empty or commented out, PoN generation will be skipped
controls_dir: "path/to/snakemake_CNApipeline/data/control"
```

# To run the pipeline: 
```bash
cd path/to/snakemake_CNApipeline

# Execute the bash script containing the snakemake run command: 
./run_pipeline.sh

# OR run for specific samples
./run_pipeline.sh sample1 sample2 
```

- If at some point you had to kill the job mid-execution, and you get a snakemake locked error. Go to the **run_pipeline.sh** script and uncomment this line: 
```bash
# Uncomment this line to unlock directory:
snakemake --unlock
```
Sometimes, snakemake requires you to clean the directories that were created during the cancelled run. While trying to run, it will generate an error with all uncomplete files specified. You can manually delete them, or uncomment this line in the **run_pipeline.sh** script.
```bash
# Uncomment this line to clean incomplete files:
snakemake --rerun-incomplete
``` 
NB: Make sure to comment the snakemake execution command when you uncomment one of these lines and run the pipeline. 

This script automatically:

- Creates a timestamped log directory
- Runs Snakemake with your default cluster profile
- Tracks execution time
- Generates a summary report

 ## Pipeline Output Directory Structure (results/)
```bash
results/
├── basecalling/  # Bams and fastq files 
├── QC/        # Per-sample QC summaries (pycoQC & cramino)
├── CNAs/      # Copy-number alteration results
│   ├── ACE/
|   ├── QDNAseq/ 
│   └── IchorCNA/
├── logs/   # Execution logs per rule & per sample
```

## Most Important Final Result Files
| Step                           | File to report              | Folder                                           |
| ------------------------------ | --------------------------- | ------------------------------------------------ |
| **QC**                         | `pycoQC_<sample>.html`      | `results/QC/<sample>/`                           |
| **Binned CN**                  | `<sample>_QDNAseq_plot.png` | `results/CNAs/QDNAseq/plots/`                    |
| **Ploidy estimation**          | `summary_<sample>.png`      | `results/CNAs/ACE/<sample>/1000kbp/2N/`          |
| **Final tumor fraction & CNA**        | `<sample>_genomeWide.pdf`   | `results/CNAs/IchorCNA/<sample>_ichor/<sample>/` |

## Personalizing ichorCNA Parameters via config.yaml
The behavior of the ichorCNA algorithm in this pipeline is fully customizable through the "ichorcna_params" section of the config file. 
These parameters adjust how tumor fraction, ploidy, and copy number states are inferred for each sample.
| Parameter              | Description                                                                                | Typical Use                                                         |
| ---------------------- | ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------- |
| `genomeBuild`          | Reference genome build used in copy number calling. Must match the reference `.wig` files. | Set this to `hg19` or `hg38` depending on your reference.           |
| `normalizeMaleX`       | Whether to apply male X chromosome normalization.                                          | Set `True` for male samples, `False` for female.                    |
| `chrs`                 | Chromosomes included in segmentation.                                                      | Typically autosomes (`1–22`).                                       |
| `chrTrain`             | Chromosomes used to train the emission models.                                             | Usually same as `chrs`.                                             |
| `chrNormalize`         | Chromosomes used for GC/bias normalization.                                                | Usually autosomal chromosomes.                                      |
| `ploidy`               | Vector of ploidy states to evaluate (e.g., diploid vs. triploid).                          | Expand if tumors are known to be highly aneuploid.                  |
| `normal`               | Expected **normal cell fraction** range (tumor purity search space).                       | Adjust range if tumor fraction is very high or low.                 |
| `maxCN`                | Maximum copy number state allowed in the model.                                            | Increase for highly amplified tumors.                               |
| `includeHOMD`          | Whether to allow homozygous deletions.                                                     | Enable (`True`) if deep deletions are expected.                     |
| `estimateNormal`       | Whether ichorCNA should estimate normal fraction.                                          | Turn off (`False`) if purity is already known.                      |
| `estimatePloidy`       | Whether ichorCNA should infer tumor ploidy.                                                | Disable in uniform ploidy datasets.                                 |
| `estimateScPrevalence` | Estimate subclonal prevalence.                                                             | Set `False` if the tumor is assumed clonal.                         |
| `scStates`             | Subclonal CNA states tested.                                                               | Leave empty (`c()`) unless analyzing strongly heterogeneous tumors. |

### How the Pipeline Uses These Settings

The parameters in config.yaml are passed directly to the runIchorCNA.R script during the IchorCNA rule execution.
This means **you do not need** to modify any workflow code. Adjusting parameters in the config is enough to reconfigure the entire IchorCNA model behavior.
