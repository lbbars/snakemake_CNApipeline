import os 
import glob

configfile:"config.yaml"

# ─────────────────────────────
# Resolve paths dynamically
# ─────────────────────────────
root = config["paths"]["root"]

def resolve_path(path):
    """Replace {root} placeholders and make absolute paths."""
    if isinstance(path, str):
        # if it uses the {root} placeholder
        path = path.replace("{root}", root)
        # if it's relative, join it with root
        if not path.startswith("/"):
            path = os.path.join(root, path)
        return os.path.normpath(path)
    elif isinstance(path, dict):
        return {k: resolve_path(v) for k, v in path.items()}
    else:
        return path

# Expand all config sections that may contain {root}
paths        = resolve_path(config["paths"])
scripts      = resolve_path(config["scripts"])
ichor_refs   = resolve_path(config["ichorcna_refs"])
# samples      = config["samples"]
ichor_params = config["ichorcna_params"]

results  = paths["results"]
logs_dir = paths["logs"]

# (optional) debug check
print("=== DEBUG PATHS ===")
for k,v in paths.items(): print(f"{k}: {v}")
for k,v in scripts.items(): print(f"script {k}: {v}")
for k,v in ichor_refs.items(): print(f"ref {k}: {v}")
print("====================")


# 1. Find CASE samples (the main analysis targets)
SAMPLES_DIR = config.get("samples_dir")
if SAMPLES_DIR is None:
    raise ValueError("samples_dir must be defined in config.yaml")

case_pod5_files = sorted(glob.glob(os.path.join(SAMPLES_DIR, "*.pod5")))
CASE_SAMPLES_MAP = {
    os.path.splitext(os.path.basename(f))[0]: f for f in case_pod5_files
}
print(f"Detected {len(CASE_SAMPLES_MAP)} CASE samples: {', '.join(CASE_SAMPLES_MAP.keys())}")

# 2. Find CONTROL samples (for the Panel of Normals)
CONTROLS_DIR = config.get("controls_dir")
CONTROL_SAMPLES_MAP = {} # Initialize as empty

if CONTROLS_DIR and os.path.exists(CONTROLS_DIR):
    control_pod5_files = sorted(glob.glob(os.path.join(CONTROLS_DIR, "*.pod5")))
    CONTROL_SAMPLES_MAP = {
        os.path.splitext(os.path.basename(f))[0]: f for f in control_pod5_files
    }
print(f"Detected {len(CONTROL_SAMPLES_MAP)} CONTROL samples: {', '.join(CONTROL_SAMPLES_MAP.keys())}")

if not CONTROL_SAMPLES_MAP:
    print("WARNING: No control samples found. Skipping Panel of Normals generation.")

# 3. Create ONE map of ALL samples for basecalling/filtering rules
ALL_SAMPLES_MAP = {**CASE_SAMPLES_MAP, **CONTROL_SAMPLES_MAP}

# 4. Handle --config sample=... override
SAMPLES_TO_RUN_MAP = {} 

TARGET_SAMPLES_STR = config.get("target_samples", None)

if TARGET_SAMPLES_STR:
    requested_list = TARGET_SAMPLES_STR.split(",")
    
    print(f"Running for {len(requested_list)} requested samples: {', '.join(requested_list)}")
    
    for sample in requested_list:
        if sample not in ALL_SAMPLES_MAP:
            raise ValueError(f"Requested sample '{sample}' not found in either samples_dir or controls_dir")
        SAMPLES_TO_RUN_MAP[sample] = ALL_SAMPLES_MAP[sample]
        
else:
    SAMPLES_TO_RUN_MAP = ALL_SAMPLES_MAP
    print(f"Running for ALL {len(SAMPLES_TO_RUN_MAP)} samples found in data directories.")

# 5. Define final path for the Panel of Normals
PANEL_OF_NORMALS_RDS = f"{results}/CNAs/IchorCNA/generated_panel_of_normals_median.rds"


# --- *** HELPER FUNCTIONS *** ---

def get_ichorcna_inputs(wildcards):
    """
    Dynamically determines the input files for the ichorcna rule.
    """
    inputs = {
        "cleanwig": f"{results}/CNAs/IchorCNA/wig/{wildcards.sample}_filtered.wig"
    }
    # Check if the CONTROLS map is not empty
    if CONTROL_SAMPLES_MAP: 
        inputs["panel"] = PANEL_OF_NORMALS_RDS
    return inputs

def get_ichorcna_panel_arg(wildcards, input, **kwargs):
    """
    Dynamically creates the --normalPanel command line argument.
    """
    if "panel" in input:
        return f"--normalPanel {input.panel}"
    else:
        return "" # Return an empty string if no panel


# ────────────────────────────────
# Global expected outputs
# ────────────────────────────────
rule all:
    input:
        # Utilisez SAMPLES_TO_RUN_MAP ici
        expand(f"{results}/basecalling/{{sample}}/{{sample}}.bam", sample=SAMPLES_TO_RUN_MAP.keys()),
        expand(f"{results}/basecalling/{{sample}}/{{sample}}.fasta", sample=SAMPLES_TO_RUN_MAP.keys()),
        expand(f"{results}/QC/{{sample}}/cramino_{{sample}}.txt",  sample=SAMPLES_TO_RUN_MAP.keys()),
        expand(f"{results}/QC/{{sample}}/pycoQC_{{sample}}.html", sample=SAMPLES_TO_RUN_MAP.keys()),
        expand(f"{results}/CNAs/QDNAseq/plots/{{sample}}_QDNAseq_plot.png", sample=SAMPLES_TO_RUN_MAP.keys()),
        expand(f"{results}/CNAs/ACE/{{sample}}/1000kbp/2N/{{sample}}/summary_{{sample}}.png", sample=SAMPLES_TO_RUN_MAP.keys()),
        expand(f"{results}/CNAs/IchorCNA/wig/{{sample}}_filtered.wig", sample=SAMPLES_TO_RUN_MAP.keys()),
        expand(f"{results}/CNAs/IchorCNA/{{sample}}_ichor/{{sample}}/{{sample}}_genomeWide.pdf", sample=SAMPLES_TO_RUN_MAP.keys())
# ────────────────────────────────
# Basecalling
# ────────────────────────────────
rule basecalling:
    input:
        pod5 = lambda wc: ALL_SAMPLES_MAP[wc.sample],
        model    = paths["model"],
        reference= paths["reference"]
    output:
        bam     = f"{results}/basecalling/{{sample}}/{{sample}}.bam",
        fastq   = f"{results}/basecalling/{{sample}}/{{sample}}.fasta",
        summary = f"{results}/basecalling/{{sample}}/summary_{{sample}}.txt"
    container: 
        "containerfiles/dorado/0.9.1/dorado_0.9.1.sif"
    log:
        f"{logs_dir}/basecalling_{{sample}}.log"

    shell:
        """
        mkdir -p $(dirname {output.bam})

        dorado basecaller {input.model} {input.pod5} \
            --device cuda:auto --min-qscore 8 \
            --modified-bases 5mCG_5hmCG \
            --reference {input.reference} > {output.bam}

        dorado summary --verbose {output.bam} > {output.summary}
        samtools sort -o {output.bam} {output.bam}
        samtools index {output.bam}
        samtools fastq {output.bam} > {output.fastq}
        """
# ────────────────────────────────
# Quality control
# ────────────────────────────────
rule cramino:
    input:
        bam = f"{results}/basecalling/{{sample}}/{{sample}}.bam",
        summary = f"{results}/basecalling/{{sample}}/summary_{{sample}}.txt"
    output:
        cramino = f"{results}/QC/{{sample}}/cramino_{{sample}}.txt"
    container: 
        "containerfiles/cramino/1.1.0/cramino_1.1.0.sif"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.cramino})
        
        cramino {input.bam} > {output.cramino}
        """

rule pycoQC:
    input: 
        summary = f"{results}/basecalling/{{sample}}/summary_{{sample}}.txt"
    output: 
        html = f"{results}/QC/{{sample}}/pycoQC_{{sample}}.html"
    container: 
        "containerfiles/pycoqc/2.5.2/pycoqc_2.5.2.sif"
    resources:
        mem_mb = 16000
    shell: 
        """      
        mkdir -p $(dirname {output.html})
        
        pycoQC -f {input.summary} -o {output.html}
        """    

# ────────────────────────────────
# QDNAseq
# ────────────────────────────────
rule qdnaseq:
    input:
        bam = f"{results}/basecalling/{{sample}}/{{sample}}.bam"
    output:
        png = f"{results}/CNAs/QDNAseq/plots/{{sample}}_QDNAseq_plot.png"
    params:
        outdir = lambda wc: f"{results}/CNAs/QDNAseq/{wc.sample}"
    container:
        "containerfiles/QDNAseq/1.9.2/QDNAseq_1.9.2.sif"
    log:
        f"{logs_dir}/qdnaseq_{{sample}}.log"
    shell:
        """
        Rscript {scripts[qdnaseq]} {input.bam} {output.png}
        """
# ────────────────────────────────
# ACE
# ────────────────────────────────
rule ACE:
    input:
        bam = f"{results}/basecalling/{{sample}}/{{sample}}.bam"
    output:
        png = f"{results}/CNAs/ACE/{{sample}}/1000kbp/2N/{{sample}}/summary_{{sample}}.png"
    params:
        outdir = lambda wc: f"{results}/CNAs/ACE/{wc.sample}",
        bamdir = lambda wc: os.path.dirname(f"{results}/basecalling/{wc.sample}/{wc.sample}.bam")
    container: 
        "containerfiles/ACE/1.27.1/ACE_1.27.1.sif"
    log:
        f"{logs_dir}/ace_{{sample}}.log"
    shell:
        """
        Rscript {scripts[ace]} {params.bamdir} {params.outdir}
        """

rule filterwigs:
    input: 
        bam = f"{results}/basecalling/{{sample}}/{{sample}}.bam"
    output: 
        cleanwig = f"{results}/CNAs/IchorCNA/wig/{{sample}}_filtered.wig"
    params:
        rawwig = f"{results}/CNAs/IchorCNA/wig/{{sample}}.wig"
    container:
        "containerfiles/ichorCNA/0.2.0/ichorCNA_0.2.0.sif"
    shell:
        r"""
        mkdir -p $(dirname {params.rawwig})

        # Step 1: Create raw wig file
        readCounter --w 1000000 --q 8 {input.bam} > {params.rawwig}

        # Step 2: Filter wig to keep only autosomes + chrX/Y
        mawk '
        BEGIN {{
            # keep only autosomes + chrX/Y
            for (i=1;i<=22;i++) CHR["chr"i]=1
            CHR["chrX"]=1
            CHR["chrY"]=1
            keep=0
        }}

        /^fixedStep/ {{
            # Extract chrom from header
            chrom=""
            n=split($0,f," ")
            for (i=1;i<=n;i++) {{
                if (index(f[i],"chrom=")==1) {{
                    split(f[i],g,"=")
                    chrom=g[2]
                    if (substr(chrom,1,3)!="chr") chrom="chr"chrom
                    break
                }}
            }}

            # Decide whether to keep this chromosome
            keep = (chrom in CHR)

            # Print the full header as-is if we keep it
            if (keep) print $0
            next
        }}

        # Print numeric values if we are keeping the chromosome
        /^[0-9eE.+-]+$/ {{
            if (keep) print
        }}' {params.rawwig} > {output.cleanwig}
        """

# ────────────────────────────────
# Panel of Normals
# ────────────────────────────────
# ────────────────────────────────
# Panel of Normals
# ────────────────────────────────
rule create_ichorcna_pon:
    input:
        wigs = expand(f"{results}/CNAs/IchorCNA/wig/{{sample}}_filtered.wig", sample=CONTROL_SAMPLES_MAP.keys()),
    output:
        rds = PANEL_OF_NORMALS_RDS 
    params:
        gcwig = ichor_refs["gcwig"],
        mapwig = ichor_refs["mapwig"],
        centromere = ichor_refs["centromere"],

        outfile_base = f"{results}/CNAs/IchorCNA/generated_panel_of_normals"
    container:
        "containerfiles/ichorCNA/0.2.0/ichorCNA_0.2.0.sif"
    log:
        f"{logs_dir}/create_ichorcna_pon.log"
    shell:
        r"""
        # 1. Créer le fichier temporaire (inchangé)
        printf "%s\n" {input.wigs} > {params.outfile_base}.filelist.txt
        
        # 2. Appeler le script R en lui passant le NOM DE BASE
        Rscript /opt/ichorCNA/scripts/createPanelOfNormals.R \
            --filelist {params.outfile_base}.filelist.txt \
            --gcWig {params.gcwig} \
            --mapWig {params.mapwig} \
            --centromere {params.centromere} \
            --outfile {params.outfile_base} # Le script ajoutera _median.rds à cela
        
        # 3. Nettoyer (inchangé)
        rm {params.outfile_base}.filelist.txt
        """
        
rule ichorcna:
    input:
        # Use the helper function to get dynamic inputs
        unpack(get_ichorcna_inputs)
    output:
        pdf = f"{results}/CNAs/IchorCNA/{{sample}}_ichor/{{sample}}/{{sample}}_genomeWide.pdf"
    params:
        outdir = lambda wc: f"{results}/CNAs/IchorCNA/{wc.sample}_ichor",
        centromere = ichor_refs["centromere"],
        gcwig = ichor_refs["gcwig"],
        mapwig = ichor_refs["mapwig"],

        # Use the helper function to get the dynamic argument string
        panel_arg = get_ichorcna_panel_arg,
        
        # (All other params are unchanged)
        genomeBuild = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("genomeBuild", config["ichorcna_params"]["genomeBuild"]),
        normalizeMaleX = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("normalizeMaleX", config["ichorcna_params"]["normalizeMaleX"]),
        chrs = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("chrs", config["ichorcna_params"]["chrs"]),
        chrTrain = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("chrTrain", config["ichorcna_params"]["chrTrain"]),
        chrNormalize = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("chrNormalize", config["ichorcna_params"]["chrNormalize"]),
        ploidy = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("ploidy", config["ichorcna_params"]["ploidy"]),
        normal = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("normal", config["ichorcna_params"]["normal"]),
        maxCN = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("maxCN", config["ichorcna_params"]["maxCN"]),
        includeHOMD = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("includeHOMD", config["ichorcna_params"]["includeHOMD"]),
        estimateNormal = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("estimateNormal", config["ichorcna_params"]["estimateNormal"]),
        estimatePloidy = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("estimatePloidy", config["ichorcna_params"]["estimatePloidy"]),
        estimateScPrev = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("estimateScPrevalence", config["ichorcna_params"]["estimateScPrevalence"]),
        scStates = lambda wc: config["ichorcna_params"].get(wc.sample, {}).get("scStates", config["ichorcna_params"]["scStates"])
    container:
        "containerfiles/ichorCNA/0.2.0/ichorCNA_0.2.0.sif"
    threads: 16
    log: 
        f"{logs_dir}/ichor_{{sample}}.log"
    shell:
        r"""
        mkdir -p {params.outdir}

        Rscript /usr/local/bin/runIchorCNA \
            --id {wildcards.sample} \
            --WIG {input.cleanwig} \
            --gcWig {params.gcwig} \
            --mapWig {params.mapwig} \
            --outDir {params.outdir} \
            --centromere {params.centromere} \
            {params.panel_arg} \
            --genomeBuild '{params.genomeBuild}' \
            --normalizeMaleX '{params.normalizeMaleX}' \
            --chrs '{params.chrs}' \
            --chrTrain '{params.chrTrain}' \
            --chrNormalize '{params.chrNormalize}' \
            --ploidy '{params.ploidy}' \
            --normal '{params.normal}' \
            --maxCN '{params.maxCN}' \
            --includeHOMD '{params.includeHOMD}' \
            --estimateNormal '{params.estimateNormal}' \
            --estimatePloidy '{params.estimatePloidy}' \
            --estimateScPrevalence '{params.estimateScPrev}' \
            --scStates '{params.scStates}'
        """
