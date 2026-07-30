# =============================================================================
# longpipe — long-read WGS + fiberseq/FIRE orchestrator
#
# Replaces the step-by-step runbook. One sample sheet + one config = the whole
# DAG. `snakemake -n` previews it; `--rerun-incomplete --keep-going` resumes it.
#
# TOPOLOGY (important):
#   Orchestrator steps (miniwdl driver, FIRE pixi) run as *local* rules on the
#   submit node (rhino) and dispatch their own cluster jobs. Leaf tools
#   (add-nucleosomes, sawfish, svtopo, hifihla, bedtools, kent) are submitted to
#   SLURM by the executor. => Run the top-level snakemake ON rhino, not inside an
#   salloc. Set --local-cores to how many WGS drivers you want running at once.
#
# KICKOFF: fill config/runsheet.csv with columns
#   core_name,lab_name,run,cell,barcode,sex,hia5
#   - label      = lab_name if given else core_name (downstream sample id)
#   - hifi bam   = derived: {rawdir}/{run}/{cell}/hifi_reads/*hifi_reads.{barcode}.bam
#   - hia5=no    => WGS/SV/HLA only; skips add-nucleosomes + FIRE branch
#                  (toggle with hia5_gates_fiberseq in config)
# =============================================================================

import os, glob, json
import pandas as pd

configfile: "config/config.yaml"

NAME     = config["name"]
OUTDIR   = config["outdir"].format(name=NAME)
RAWDIR   = config["rawdir"]                     # dir CONTAINING the r84265_* run folders
RUNSHEET = config["runsheet"]
HIA5_GATES = config.get("hia5_gates_fiberseq", True)

FASTA    = config["ref_fasta"]
HLA_REF  = config["hla_ref_fasta"]
CHROMSZ  = config["chrom_sizes"]
PEAKS    = config["ccre_peaks"]                 # rules use ${PEAKS}.gz (tabix'd)
REFMAP   = config["ref_map_file"]
TERTMAP  = config["tertiary_map_file"]
WDL      = config["wdl_singleton"]
FIRE_PIXI= config["fire_pixi"]
FIRE_REF = config["fire_ref"]
CNV_EXCL = config["cnv_excluded_regions"]
EXP_CN   = {"MALE": config["expected_cn_xy"], "FEMALE": config["expected_cn_xx"]}

# CONFIRM: gizmo Lmod init so `ml` works in Snakemake's non-login shells.
MODINIT  = config.get("module_init", "/etc/profile.d/modules.sh")
shell.prefix(f'set -euo pipefail; source {MODINIT} >/dev/null 2>&1 || true; ')

def act(env):
    return f'set +u; eval "$(conda shell.bash hook)"; conda activate {env}; set -u; '

workdir: OUTDIR

# ── sample sheet: single source of truth ────────────────────────────────────
samples = pd.read_csv(RUNSHEET, dtype=str).fillna("")

def _label(r):
    lab, core = str(r.get("lab_name", "")).strip(), str(r.get("core_name", "")).strip()
    if not (lab or core):
        raise ValueError("each runsheet row needs core_name (lab_name optional)")
    return lab or core

samples["label"] = samples.apply(_label, axis=1)
samples = samples.set_index("label", drop=False)
if samples.index.duplicated().any():
    raise ValueError(f"duplicate labels: {list(samples.index[samples.index.duplicated()])}")
SAMPLES = list(samples.index)

def _truthy(x):
    return str(x).strip().lower() in ("y", "yes", "true", "t", "1")

def is_hia5(label):
    return _truthy(samples.loc[label, "hia5"])

def fiberseq(label):
    # a sample gets the fiberseq path (nucleosomes + FIRE) if Hia5-treated
    # (or unconditionally, if the gate is turned off)
    return (not HIA5_GATES) or is_hia5(label)

FIRE_SAMPLES = [s for s in SAMPLES if fiberseq(s)]

def raw_hifi(label):
    r  = samples.loc[label]
    bc = r["barcode"].replace("bc-", "bc")
    # os.path.join drops RAWDIR automatically if `run` is an absolute path
    pat = os.path.join(RAWDIR, r["run"], r["cell"], "hifi_reads", f"*hifi_reads.{bc}.bam")
    hits = glob.glob(pat)
    if len(hits) != 1:
        raise ValueError(f"{label}: expected exactly 1 hifi bam matching {pat}, found {hits}")
    return hits[0]

def sex_of(label):
    s = str(samples.loc[label, "sex"]).strip().upper()
    if s in ("M", "MALE"):   return "MALE"
    if s in ("F", "FEMALE"): return "FEMALE"
    return config.get("default_sex", "FEMALE")

def wgs_hifi(w):
    # fiberseq samples feed the nucleosome bam to WGS; WGS-only samples feed raw hifi
    return f"fibertools/{w.s}.nucleosomes.bam" if fiberseq(w.s) else raw_hifi(w.s)

HAPBAM = "pacbiowdlR/{s}/{s}.GRCh38.haplotagged.bam"   # the fan-out point

localrules: all, wdl_input, hifi_wgs, fire_config, fire

# ── final targets ───────────────────────────────────────────────────────────
rule all:
    input:
        # FIRE branch — Hia5 samples only
        expand("fibertools/FIRE/results/{s}/{s}-fire-v0.1.1-pileup.bed.gz", s=FIRE_SAMPLES),
        expand("fibertools/FIRE/results/bigwigs/{s}-fire-v0.1.1-pileup_shrunk.bw", s=FIRE_SAMPLES),
        expand("fibertools/FIRE/scATAC_peak/{s}_scATAC_max_scores.bed", s=FIRE_SAMPLES),
        # SV branch — all samples
        expand("sawfish/{s}_joint_call_dir/genotyped.sv.vcf.gz", s=SAMPLES),
        expand("sawfish/{s}.svtopo.done", s=SAMPLES),
        # HLA branch — all samples
        expand("hla/{s}/{s}_hifihla.json", s=SAMPLES),

# ── 1. add nucleosomes (GPU) — fiberseq samples only ────────────────────────
rule add_nucleosomes:
    input:  lambda w: raw_hifi(w.s)
    output: "fibertools/{s}.nucleosomes.bam"
    threads: 30
    resources:
        mem_mb_per_cpu=21000, runtime=180,
        slurm_partition="campus-new",
        slurm_extra="'--gres=gpu:1'",
    shell:
        "ml Python/3.11.5-GCCcore-13.2.0 CUDA/12.4.0; " + act("fibertools") +
        "ft add-nucleosomes -t {threads} {input} > {output}"

# ── 2. HiFi-human-WGS-WDL (miniwdl driver = local; dispatches to cluster) ────
rule wdl_input:
    input:  wgs_hifi
    output: "pbWGS/{s}/input.json"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        d = {
            "humanwgs_singleton.sample_id": wildcards.s,
            "humanwgs_singleton.sex": sex_of(wildcards.s),
            "humanwgs_singleton.hifi_reads": [os.path.abspath(input[0])],
            "humanwgs_singleton.ref_map_file": REFMAP,
            "humanwgs_singleton.tertiary_map_file": TERTMAP,
            "humanwgs_singleton.backend": "HPC",
            "humanwgs_singleton.preemptible": True,
        }
        with open(output[0], "w") as fh:
            json.dump(d, fh, indent=2)

rule hifi_wgs:
    input:  bam=wgs_hifi, cfg="pbWGS/{s}/input.json"
    # CONFIRM outputs.json location: some miniwdl versions nest a run subdir and
    # the real file is pbWGS/{s}/_LAST/outputs.json — adjust here + copy_pbwdl if so.
    output: "pbWGS/{s}/outputs.json"
    shell:
        "ml Python/3.9.6-GCCcore-11.2.0; "
        "source ~/.virtualenvs/hifi/bin/activate; "
        "miniwdl run {WDL} -i {input.cfg} --dir pbWGS/{wildcards.s} --verbose"

# ── 3. extract WDL outputs -> haplotagged bam ───────────────────────────────
rule copy_pbwdl:
    input:  "pbWGS/{s}/outputs.json"
    output: HAPBAM
    threads: 20
    resources: mem_mb_per_cpu=21000, runtime=60, slurm_partition="campus-new"
    shell:
        "~/develop/pacbiowdlR/copy_pbWDL -i {input} -o pacbiowdlR/{wildcards.s} -t {threads} -m"

# ── 4. FIRE (Hia5 samples at once; pixi/snakemake driver = local) ───────────
rule fire_config:
    input:  expand(HAPBAM, s=FIRE_SAMPLES)
    output: tbl="fibertools/FIRE/config.tbl", yaml="fibertools/FIRE/config.yaml"
    run:
        os.makedirs("fibertools/FIRE", exist_ok=True)
        with open(output.tbl, "w") as t:
            t.write("sample\tbam\n")
            for s in FIRE_SAMPLES:
                t.write(f"{s}\t{os.path.abspath(f'pacbiowdlR/{s}/{s}.GRCh38.haplotagged.bam')}\n")
        with open(output.yaml, "w") as y:
            y.write(f"ref: {FIRE_REF}\nref_name: hg38\nmanifest: {os.path.abspath(output.tbl)}\n")

rule fire:
    input:  tbl="fibertools/FIRE/config.tbl", yaml="fibertools/FIRE/config.yaml"
    output: expand("fibertools/FIRE/results/{s}/{s}-fire-v0.1.1-pileup.bed.gz", s=FIRE_SAMPLES)
    shell:
        "cd fibertools/FIRE && "
        "pixi run --manifest-path {FIRE_PIXI} slurm --configfile config.yaml -k"

# ── 5a. bigwigs from FIRE pileups ───────────────────────────────────────────
rule bigwig:
    input:  "fibertools/FIRE/results/{s}/{s}-fire-v0.1.1-pileup.bed.gz"
    output: "fibertools/FIRE/results/bigwigs/{s}-fire-v0.1.1-pileup_shrunk.bw"
    threads: 10
    resources: mem_mb_per_cpu=8000, runtime=30, slurm_partition="campus-new"
    shell:
        "ml Kent_tools/468-GCC-12.3.0; "
        "mkdir -p fibertools/FIRE/results/bigwigs; "
        "BG=fibertools/FIRE/results/bigwigs/{wildcards.s}-fire-v0.1.1-pileup_shrunk.bg; "
        "~/develop/pacbiowdlR/shrink_fire_bed -c {CHROMSZ} {input} $BG -t {threads}; "
        "bedGraphToBigWig $BG {CHROMSZ} {output}; rm $BG"

# ── 5b. scATAC peak mapping from FIRE pileups ───────────────────────────────
rule scatac:
    input:  "fibertools/FIRE/results/{s}/{s}-fire-v0.1.1-pileup.bed.gz"
    output:
        mx="fibertools/FIRE/scATAC_peak/{s}_scATAC_max_scores.bed",
        md="fibertools/FIRE/scATAC_peak/{s}_scATAC_median_scores.bed",
        sm="fibertools/FIRE/scATAC_peak/{s}_scATAC_sum_scores.bed",
    threads: 8
    resources: mem_mb_per_cpu=21000, runtime=30, slurm_partition="campus-new"
    shell:
        "ml BEDTools/2.31.0-GCC-12.3.0 HTSlib/1.19.1-GCC-13.2.0; "
        "mkdir -p fibertools/FIRE/scATAC_peak; "
        "bedtools map -a {PEAKS}.gz -b {input} -c 6 -o max    > {output.mx}; "
        "bedtools map -a {PEAKS}.gz -b {input} -c 6 -o median > {output.md}; "
        "bedtools map -a {PEAKS}.gz -b {input} -c 6 -o sum    > {output.sm}"

# ── 6. sawfish discover -> joint-call -> svtopo (all samples) ───────────────
rule sawfish_discover:
    input:  HAPBAM
    output: directory("sawfish/{s}_discover")
    params: cn=lambda w: EXP_CN[sex_of(w.s)]
    threads: 30
    resources: mem_mb_per_cpu=5000, runtime=120, slurm_partition="campus-new"
    shell:
        "sawfish discover --threads {threads} --ref {FASTA} --bam {input} "
        "--expected-cn {params.cn} --cnv-excluded-regions {CNV_EXCL} "
        "--output-dir {output}"

rule sawfish_jointcall:
    input:  "sawfish/{s}_discover"
    output: "sawfish/{s}_joint_call_dir/genotyped.sv.vcf.gz"
    threads: 30
    resources: mem_mb_per_cpu=5000, runtime=120, slurm_partition="campus-new"
    shell:
        "sawfish joint-call --threads {threads} --report-supporting-reads "
        "--sample {input} --output-dir sawfish/{wildcards.s}_joint_call_dir"

rule svtopo:
    input:  bam=HAPBAM, vcf="sawfish/{s}_joint_call_dir/genotyped.sv.vcf.gz"
    output: touch("sawfish/{s}.svtopo.done")   # CONFIRM svtopo's real output names
    resources: mem_mb_per_cpu=5000, runtime=30, slurm_partition="campus-new"
    shell:
        "svtopo --bam {input.bam} --svtopo-dir sawfish/{wildcards.s}_joint_call_dir "
        "--prefix {wildcards.s} "
        "--variant-readnames sawfish/{wildcards.s}_joint_call_dir/supporting_reads.json.gz "
        "--vcf {input.vcf} --exclude-regions {CNV_EXCL}"

# ── 7. HLA (per-sample DAG — native deps replace hand-rolled afterok) ────────
rule hla_pbindex:
    input:  HAPBAM
    output: "pacbiowdlR/{s}/{s}.GRCh38.haplotagged.bam.pbi"
    threads: 30
    resources: mem_mb_per_cpu=21000, runtime=60, slurm_partition="campus-new"
    shell:  act("hifi-hla") + "pbindex -j {threads} {input}"

rule hla_fastq:
    input:  bam=HAPBAM, pbi="pacbiowdlR/{s}/{s}.GRCh38.haplotagged.bam.pbi"
    output: "fastq/{s}/{s}.fastq.gz"
    threads: 30
    resources: mem_mb_per_cpu=21000, runtime=120, slurm_partition="campus-new"
    shell:  act("hifi-hla") + "bam2fastq -j {threads} {input.bam} -o fastq/{wildcards.s}/{wildcards.s}"

rule hla_hifiasm:
    input:  "fastq/{s}/{s}.fastq.gz"
    output: h1="hifiasm/{s}/{s}.asm.bp.hap1.p_ctg.gfa",
            h2="hifiasm/{s}/{s}.asm.bp.hap2.p_ctg.gfa"
    threads: 30
    resources: mem_mb_per_cpu=21000, runtime=480, slurm_partition="campus-new"
    shell:  act("hifi-hla") + "hifiasm -o hifiasm/{wildcards.s}/{wildcards.s}.asm -t {threads} {input}"

rule hla_gfa2fa:
    input:  h1="hifiasm/{s}/{s}.asm.bp.hap1.p_ctg.gfa",
            h2="hifiasm/{s}/{s}.asm.bp.hap2.p_ctg.gfa"
    output: f1="hifiasm/{s}/{s}.asm.bp.hap1.p_ctg.fa",
            f2="hifiasm/{s}/{s}.asm.bp.hap2.p_ctg.fa"
    resources: mem_mb=4000, runtime=20, slurm_partition="campus-new"
    shell:
        r"""awk '/^S/{{print ">"$2;print $3}}' {input.h1} > {output.f1}
awk '/^S/{{print ">"$2;print $3}}' {input.h2} > {output.f2}"""

rule hla_align:
    input:  f1="hifiasm/{s}/{s}.asm.bp.hap1.p_ctg.fa",
            f2="hifiasm/{s}/{s}.asm.bp.hap2.p_ctg.fa"
    output: "hifiasm/{s}/{s}.GRCh38.asm_sorted.bam"
    threads: 30
    resources: mem_mb_per_cpu=21000, runtime=180, slurm_partition="campus-new"
    shell:
        "ml SAMtools/1.11-GCC-10.2.0 minimap2/2.26-GCCcore-12.3.0; "
        r'minimap2 -t {threads} -L --secondary=no --eqx -ax asm5 '
        r'-R "@RG\tID:{wildcards.s}\tSM:{wildcards.s}" {HLA_REF} {input.f1} {input.f2} | '
        "samtools sort -@ {threads} -o {output}"

rule hla_index:
    input:  "hifiasm/{s}/{s}.GRCh38.asm_sorted.bam"
    output: "hifiasm/{s}/{s}.GRCh38.asm_sorted.bam.bai"
    resources: mem_mb_per_cpu=21000, runtime=30, slurm_partition="campus-new"
    shell:  "ml SAMtools/1.11-GCC-10.2.0; samtools index {input}"

rule hla_call:
    input:  bam="hifiasm/{s}/{s}.GRCh38.asm_sorted.bam",
            bai="hifiasm/{s}/{s}.GRCh38.asm_sorted.bam.bai",
            f1="hifiasm/{s}/{s}.asm.bp.hap1.p_ctg.fa",
            f2="hifiasm/{s}/{s}.asm.bp.hap2.p_ctg.fa"
    output: "hla/{s}/{s}_hifihla.json"   # CONFIRM hifihla's exact output suffix
    threads: 30
    resources: mem_mb_per_cpu=21000, runtime=120, slurm_partition="campus-new"
    shell:
        act("hifi-hla") +
        "hifihla call-contigs --abam {input.bam} --hap1 {input.f1} --hap2 {input.f2} "
        "--threads {threads} --out_prefix hla/{wildcards.s}/{wildcards.s}_hifihla"
