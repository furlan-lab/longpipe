FIRE_CORE_ROOT = "intermediate/fire/core"
FIRE_FILTERED_CRAM = (
    f"{FIRE_CORE_ROOT}/{{sample}}/"
    f"{{sample}}-fire-{FIRE_VERSION}-filtered.cram"
)
FIRE_FILTERED_CRAM_INDEX = f"{FIRE_FILTERED_CRAM}.crai"
FIRE_BEDGRAPH = (
    f"{FIRE_CORE_ROOT}/{{sample}}/coverage/"
    f"{{sample}}-{FIRE_VERSION}.bed.gz"
)
FIRE_BEDGRAPH_INDEX = f"{FIRE_BEDGRAPH}.tbi"
FIRE_MEDIAN_COVERAGE = (
    f"{FIRE_CORE_ROOT}/{{sample}}/coverage/"
    f"{{sample}}-{FIRE_VERSION}-median-coverage.txt"
)
FIRE_MINIMUM_COVERAGE = (
    f"{FIRE_CORE_ROOT}/{{sample}}/coverage/"
    f"{{sample}}-{FIRE_VERSION}-minimum-coverage.txt"
)
FIRE_MAXIMUM_COVERAGE = (
    f"{FIRE_CORE_ROOT}/{{sample}}/coverage/"
    f"{{sample}}-{FIRE_VERSION}-maximum-coverage.txt"
)
FIRE_CHROMOSOME_FIBERS = (
    f"{FIRE_CORE_ROOT}/{{sample}}/chrom/"
    f"{FIRE_VERSION}-{{fire_chrom}}.fiber-locations.bed.gz"
)
FIRE_FIBERS = (
    f"{FIRE_CORE_ROOT}/{{sample}}/coverage/"
    f"{FIRE_VERSION}-fiber-locations.bed.gz"
)
FIRE_FIBERS_INDEX = f"{FIRE_FIBERS}.tbi"
FIRE_FILTERED_FIBERS = (
    f"{FIRE_CORE_ROOT}/{{sample}}/coverage/filtered-for-coverage/"
    f"{FIRE_VERSION}-fiber-locations.bed.gz"
)
FIRE_FILTERED_FIBERS_INDEX = f"{FIRE_FILTERED_FIBERS}.tbi"
FIRE_SHUFFLE_EXCLUDES = (
    f"{FIRE_CORE_ROOT}/{{sample}}/coverage/"
    "exclude-from-shuffles.bed.gz"
)
FIRE_SHUFFLED_FIBERS = (
    f"{FIRE_CORE_ROOT}/{{sample}}/shuffle/"
    f"{FIRE_VERSION}-{{fire_chrom}}.fiber-locations-shuffled.bed.gz"
)
FIRE_CHROMOSOME_SHUFFLED_PILEUP = (
    f"{FIRE_CORE_ROOT}/{{sample}}/shuffle/"
    f"{FIRE_VERSION}-{{fire_chrom}}.pileup.bed.gz"
)
FIRE_SHUFFLED_PILEUP = (
    f"{FIRE_CORE_ROOT}/{{sample}}/shuffle/"
    f"{FIRE_VERSION}-shuffled-pileup.bed.gz"
)
FIRE_SHUFFLED_PILEUP_INDEX = f"{FIRE_SHUFFLED_PILEUP}.tbi"
FIRE_FDR_TABLE = (
    f"{FIRE_CORE_ROOT}/{{sample}}/fire-peaks/"
    f"{{sample}}-{FIRE_VERSION}-fire-score-to-fdr.tbl"
)
FIRE_CHROMOSOME_PILEUP = (
    f"{FIRE_CORE_ROOT}/{{sample}}/pileup/"
    f"{FIRE_VERSION}-{{fire_chrom}}.pileup.bed.gz"
)
FIRE_CHROMOSOME_FDR_TRACK = (
    f"{FIRE_CORE_ROOT}/{{sample}}/fire-peaks/"
    f"{FIRE_VERSION}-{{fire_chrom}}-FDR.track.bed"
)


rule fire_filtered_cram:
    input:
        bam=eligible_fire_bam,
        bai=HAPLOTAGGED_BAM_INDEX,
        reference=FIRE_REF,
        reference_index=f"{FIRE_REF}.fai"
    output:
        cram=temp(FIRE_FILTERED_CRAM),
        crai=temp(FIRE_FILTERED_CRAM_INDEX)
    threads: 11
    resources:
        mem_mb=32768,
        runtime=600,
        io_heavy=1
    log:
        "logs/fire/filtered_cram/{sample}.log"
    benchmark:
        "benchmarks/fire/filtered_cram/{sample}.tsv"
    params:
        flag=FIRE_FILTER_FLAG,
        min_msp=FIRE_MIN_MSP,
        min_average_msp_size=FIRE_MIN_AVERAGE_MSP_SIZE
    conda:
        "../../envs/fire-core.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.cram:q})" "$(dirname {log:q})"
        (
            samtools view -@ {threads} -u -F {params.flag} {input.bam:q} \
                | ft fire -F {params.flag} -t {threads} \
                    --min-msp {params.min_msp} \
                    --min-ave-msp-size {params.min_average_msp_size} \
                    --skip-no-m6a - - \
                | samtools view -C -@ {threads} -T {input.reference:q} \
                    --output-fmt-option embed_ref=1 \
                | samtools view -C -@ {threads} -T {input.reference:q} \
                    --output-fmt-option embed_ref=1 \
                    --input-fmt-option required_fields=0x1bff \
                    --write-index -o {output.cram:q}
            reads_in_cram="$(
                samtools view {output.cram:q} | head | wc -l || true
            )"
            if [[ "$reads_in_cram" -eq 0 ]]; then
                printf '%s\n' 'No reads passed the FIRE filters.' >&2
                exit 1
            fi
        ) > {log:q} 2>&1
        test -s {output.cram:q}
        test -s {output.crai:q}
        """


rule fire_genome_bedgraph:
    input:
        reference=FIRE_REF,
        reference_index=f"{FIRE_REF}.fai",
        cram=rules.fire_filtered_cram.output.cram,
        crai=rules.fire_filtered_cram.output.crai
    output:
        bedgraph=temp(FIRE_BEDGRAPH),
        index=temp(FIRE_BEDGRAPH_INDEX)
    threads: 11
    resources:
        mem_mb=32768,
        runtime=240,
        io_heavy=1
    log:
        "logs/fire/genome_bedgraph/{sample}.log"
    benchmark:
        "benchmarks/fire/genome_bedgraph/{sample}.tsv"
    conda:
        "../../envs/fire-core.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.bedgraph:q})" "$(dirname {log:q})"
        temporary_directory="$(
            mktemp -d "${{TMPDIR:-/tmp}}/{wildcards.sample}.mosdepth.XXXXXX"
        )"
        trap 'rm -rf "$temporary_directory"' EXIT
        (
            mosdepth -F 4 -f {input.reference:q} -t {threads} \
                "$temporary_directory/coverage" {input.cram:q}
            bgzip -cd "$temporary_directory/coverage.per-base.bed.gz" \
                | LC_ALL=C sort --parallel={threads} \
                    -k1,1 -k2,2n -k3,3n -k4,4 \
                | bgzip -@ {threads} \
                > {output.bedgraph:q}
            tabix -f -p bed {output.bedgraph:q}
        ) > {log:q} 2>&1
        test -s {output.bedgraph:q}
        test -s {output.index:q}
        """


rule fire_coverage:
    input:
        bedgraph=rules.fire_genome_bedgraph.output.bedgraph
    output:
        median=temp(FIRE_MEDIAN_COVERAGE),
        minimum=temp(FIRE_MINIMUM_COVERAGE),
        maximum=temp(FIRE_MAXIMUM_COVERAGE)
    threads: 1
    resources:
        mem_mb=65536,
        runtime=120
    log:
        "logs/fire/coverage/{sample}.log"
    benchmark:
        "benchmarks/fire/coverage/{sample}.tsv"
    params:
        script=workflow.source_path("../scripts/fire_coverage.py"),
        chromosomes=FIRE_CHROMS,
        coverage_within_n_sd=FIRE_COVERAGE_WITHIN_N_SD,
        min_coverage=FIRE_MIN_COVERAGE
    conda:
        "../../envs/fire-python.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.median:q})" "$(dirname {log:q})"
        chromosome_arguments=()
        for chromosome in {params.chromosomes:q}; do
            chromosome_arguments+=(--chromosome "$chromosome")
        done
        POLARS_MAX_THREADS={threads} python {params.script:q} \
            --coverage-within-n-sd {params.coverage_within_n_sd} \
            --min-coverage {params.min_coverage} \
            "${{chromosome_arguments[@]}}" \
            {input.bedgraph:q} \
            {output.median:q} \
            {output.minimum:q} \
            {output.maximum:q} \
            > {log:q} 2>&1
        """


rule fire_fiber_locations_chromosome:
    input:
        cram=rules.fire_filtered_cram.output.cram,
        crai=rules.fire_filtered_cram.output.crai
    output:
        bed=temp(FIRE_CHROMOSOME_FIBERS)
    threads: 4
    resources:
        mem_mb=16384,
        runtime=120
    log:
        "logs/fire/fiber_locations/{sample}.{fire_chrom}.log"
    benchmark:
        "benchmarks/fire/fiber_locations/{sample}.{fire_chrom}.tsv"
    conda:
        "../../envs/fire-core.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.bed:q})" "$(dirname {log:q})"
        (
            samtools view -@ {threads} -u {input.cram:q} \
                {wildcards.fire_chrom:q} \
                | ft extract -t {threads} -s --all - \
                | hck -F '#ct' -F st -F en -F fiber -F strand -F HP \
                | (grep -v '^#' || true) \
                | bgzip -@ {threads} \
                > {output.bed:q}
        ) > {log:q} 2>&1
        test -s {output.bed:q}
        """


rule fire_fiber_locations:
    input:
        fibers=lambda wildcards: expand(
            FIRE_CHROMOSOME_FIBERS,
            sample=wildcards.sample,
            fire_chrom=FIRE_CHROMS,
        ),
        bedgraph=rules.fire_genome_bedgraph.output.bedgraph,
        minimum=rules.fire_coverage.output.minimum,
        maximum=rules.fire_coverage.output.maximum
    output:
        bed=temp(FIRE_FIBERS),
        index=temp(FIRE_FIBERS_INDEX),
        filtered=temp(FIRE_FILTERED_FIBERS),
        filtered_index=temp(FIRE_FILTERED_FIBERS_INDEX)
    threads: 4
    resources:
        mem_mb=32768,
        runtime=120,
        io_heavy=1
    log:
        "logs/fire/fiber_locations/{sample}.merge.log"
    benchmark:
        "benchmarks/fire/fiber_locations/{sample}.merge.tsv"
    params:
        max_fraction_overlap=0.2
    conda:
        "../../envs/fire-core.yaml"
    shell:
        """
        mkdir -p \
            "$(dirname {output.bed:q})" \
            "$(dirname {output.filtered:q})" \
            "$(dirname {log:q})"
        (
            cat {input.fibers:q} > {output.bed:q}
            tabix -f -p bed {output.bed:q}
            minimum_coverage="$(cat {input.minimum:q})"
            maximum_coverage="$(cat {input.maximum:q})"
            bedtools intersect -header -sorted -v \
                -f {params.max_fraction_overlap} \
                -a {output.bed:q} \
                -b <(
                    bgzip -cd {input.bedgraph:q} \
                        | awk \
                            -v MAX="$maximum_coverage" \
                            -v MIN="$minimum_coverage" \
                            '$4 <= MIN || $4 >= MAX'
                ) \
                | bgzip -@ {threads} \
                > {output.filtered:q}
            tabix -f -p bed {output.filtered:q}
        ) > {log:q} 2>&1
        test -s {output.bed:q}
        test -s {output.index:q}
        test -s {output.filtered:q}
        test -s {output.filtered_index:q}
        """


rule fire_exclude_from_shuffle:
    input:
        filtered=rules.fire_fiber_locations.output.filtered,
        reference_index=f"{FIRE_REF}.fai",
        excludes=FIRE_EXCLUDE_REGIONS
    output:
        bed=temp(FIRE_SHUFFLE_EXCLUDES)
    threads: 4
    resources:
        mem_mb=16384,
        runtime=120
    log:
        "logs/fire/exclude_from_shuffle/{sample}.log"
    benchmark:
        "benchmarks/fire/exclude_from_shuffle/{sample}.tsv"
    conda:
        "../../envs/fire-core.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.bed:q})" "$(dirname {log:q})"
        exclude_files=({input.excludes:q})
        (
            (
                bedtools genomecov -bga \
                    -i {input.filtered:q} \
                    -g {input.reference_index:q} \
                    | awk '$4 == 0'
                for exclude_file in "${{exclude_files[@]}}"; do
                    if [[ "$exclude_file" == *.gz ]]; then
                        bgzip -cd "$exclude_file"
                    else
                        cat "$exclude_file"
                    fi
                done
            ) \
                | cut -f 1-3 \
                | bedtools sort \
                | bedtools merge \
                | bgzip -@ {threads} \
                > {output.bed:q}
        ) > {log:q} 2>&1
        test -s {output.bed:q}
        """


rule fire_shuffled_fiber_locations_chromosome:
    input:
        filtered=rules.fire_fiber_locations.output.filtered,
        filtered_index=rules.fire_fiber_locations.output.filtered_index,
        excludes=rules.fire_exclude_from_shuffle.output.bed,
        reference_index=f"{FIRE_REF}.fai"
    output:
        bed=temp(FIRE_SHUFFLED_FIBERS)
    threads: 4
    resources:
        mem_mb=16384,
        runtime=120
    log:
        "logs/fire/shuffle/{sample}.{fire_chrom}.log"
    benchmark:
        "benchmarks/fire/shuffle/{sample}.{fire_chrom}.tsv"
    conda:
        "../../envs/fire-core.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.bed:q})" "$(dirname {log:q})"
        (
            tabix {input.filtered:q} {wildcards.fire_chrom:q} \
                | bioawk -t '{{print $1,$2,$3,$4,$2}}' \
                | bedtools shuffle -chrom -seed 42 \
                    -excl {input.excludes:q} \
                    -i - \
                    -g {input.reference_index:q} \
                | sort -k1,1 -k2,2n -k3,3n -k4,4 \
                | bgzip -@ {threads} \
                > {output.bed:q}
        ) > {log:q} 2>&1
        test -s {output.bed:q}
        """


rule fire_shuffled_pileup_chromosome:
    input:
        cram=rules.fire_filtered_cram.output.cram,
        shuffled=rules.fire_shuffled_fiber_locations_chromosome.output.bed
    output:
        bed=temp(FIRE_CHROMOSOME_SHUFFLED_PILEUP)
    threads: 4
    resources:
        mem_mb=32768,
        runtime=120
    log:
        "logs/fire/shuffled_pileup/{sample}.{fire_chrom}.log"
    benchmark:
        "benchmarks/fire/shuffled_pileup/{sample}.{fire_chrom}.tsv"
    conda:
        "../../envs/fire-core.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.bed:q})" "$(dirname {log:q})"
        (
            ft pileup {input.cram:q} {wildcards.fire_chrom:q} \
                -t {threads} \
                --fiber-coverage \
                --shuffle {input.shuffled:q} \
                --no-msp \
                --no-nuc \
                | bgzip -@ {threads} \
                > {output.bed:q}
        ) > {log:q} 2>&1
        test -s {output.bed:q}
        """


rule fire_shuffled_pileup:
    input:
        beds=lambda wildcards: expand(
            FIRE_CHROMOSOME_SHUFFLED_PILEUP,
            sample=wildcards.sample,
            fire_chrom=FIRE_CHROMS,
        )
    output:
        bed=temp(FIRE_SHUFFLED_PILEUP),
        index=temp(FIRE_SHUFFLED_PILEUP_INDEX)
    threads: 4
    resources:
        mem_mb=16384,
        runtime=60,
        io_heavy=1
    log:
        "logs/fire/shuffled_pileup/{sample}.merge.log"
    benchmark:
        "benchmarks/fire/shuffled_pileup/{sample}.merge.tsv"
    conda:
        "../../envs/fire-core.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.bed:q})" "$(dirname {log:q})"
        (
            cat {input.beds:q} > {output.bed:q}
            tabix -f -p bed {output.bed:q}
        ) > {log:q} 2>&1
        test -s {output.bed:q}
        test -s {output.index:q}
        """


rule fire_fdr_table:
    input:
        shuffled=rules.fire_shuffled_pileup.output.bed,
        minimum=rules.fire_coverage.output.minimum,
        maximum=rules.fire_coverage.output.maximum
    output:
        table=temp(FIRE_FDR_TABLE)
    threads: 11
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 92 * 1024,
        runtime=240
    log:
        "logs/fire/fdr_table/{sample}.log"
    benchmark:
        "benchmarks/fire/fdr_table/{sample}.tsv"
    params:
        script=workflow.source_path("../scripts/fire_fdr_table.py")
    conda:
        "../../envs/fire-python.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.table:q})" "$(dirname {log:q})"
        minimum_coverage="$(cat {input.minimum:q})"
        maximum_coverage="$(cat {input.maximum:q})"
        POLARS_MAX_THREADS={threads} python {params.script:q} -v 1 \
            {input.shuffled:q} \
            {output.table:q} \
            --max-cov "$maximum_coverage" \
            --min-cov "$minimum_coverage" \
            > {log:q} 2>&1
        test -s {output.table:q}
        """


rule fire_pileup_chromosome:
    input:
        cram=rules.fire_filtered_cram.output.cram
    output:
        bed=temp(FIRE_CHROMOSOME_PILEUP)
    threads: 4
    resources:
        mem_mb=32768,
        runtime=120
    log:
        "logs/fire/pileup/{sample}.{fire_chrom}.log"
    benchmark:
        "benchmarks/fire/pileup/{sample}.{fire_chrom}.tsv"
    conda:
        "../../envs/fire-core.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.bed:q})" "$(dirname {log:q})"
        (
            ft pileup -t {threads} \
                --haps \
                --fiber-coverage \
                {input.cram:q} \
                {wildcards.fire_chrom:q} \
                | bgzip -@ {threads} \
                > {output.bed:q}
        ) > {log:q} 2>&1
        test -s {output.bed:q}
        """


rule fire_fdr_track_chromosome:
    input:
        pileup=rules.fire_pileup_chromosome.output.bed,
        fdr_table=rules.fire_fdr_table.output.table
    output:
        bed=temp(FIRE_CHROMOSOME_FDR_TRACK)
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 92 * 1024,
        runtime=240
    log:
        "logs/fire/fdr_track/{sample}.{fire_chrom}.log"
    benchmark:
        "benchmarks/fire/fdr_track/{sample}.{fire_chrom}.tsv"
    params:
        script=workflow.source_path("../scripts/fire_fdr_table.py")
    conda:
        "../../envs/fire-python.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.bed:q})" "$(dirname {log:q})"
        POLARS_MAX_THREADS={threads} python {params.script:q} -v 1 \
            --fdr-table {input.fdr_table:q} \
            {input.pileup:q} \
            {output.bed:q} \
            > {log:q} 2>&1
        test -e {output.bed:q}
        """


rule fire_pileup:
    input:
        beds=lambda wildcards: expand(
            FIRE_CHROMOSOME_FDR_TRACK,
            sample=wildcards.sample,
            fire_chrom=FIRE_CHROMS,
        )
    output:
        pileup=FIRE_PILEUP,
        index=FIRE_PILEUP_INDEX
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
        io_heavy=1
    log:
        "logs/fire/pileup/{sample}.merge.log"
    benchmark:
        "benchmarks/fire/pileup/{sample}.merge.tsv"
    conda:
        "../../envs/fire-core.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.pileup:q})" "$(dirname {log:q})"
        track_files=({input.beds:q})
        (
            (
                cat "${{track_files[@]}}" \
                    | grep '^#' \
                    | head -n 1
            ) || true
            (
                cat "${{track_files[@]}}" \
                    | grep -v '^#'
            ) || true
        ) \
            | bgzip -@ {threads} \
            > {output.pileup:q} 2> {log:q}
        tabix -f -p bed {output.pileup:q} >> {log:q} 2>&1
        test -s {output.pileup:q}
        test -s {output.index:q}
        """


rule fire_haplotype_tracks:
    input:
        pileup=FIRE_PILEUP,
        index=FIRE_PILEUP_INDEX,
        chrom_sizes=CHROM_SIZES
    output:
        h1=FIRE_H1_BIGWIG,
        h2=FIRE_H2_BIGWIG,
        difference=FIRE_HAPLOTYPE_DIFFERENCE_BIGWIG
    threads: 4
    resources:
        mem_mb=16000,
        runtime=60,
        io_heavy=1
    log:
        "logs/fire/haplotypes/{sample}.tracks.log"
    benchmark:
        "benchmarks/fire/haplotypes/{sample}.tracks.tsv"
    params:
        script=workflow.source_path("../scripts/fire_haplotype_tracks.py")
    conda:
        "../../envs/intervals.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.h1:q})" "$(dirname {log:q})"
        temporary_directory="$(
            mktemp -d "${{TMPDIR:-/tmp}}/{wildcards.sample}.haptracks.XXXXXX"
        )"
        trap 'rm -rf "$temporary_directory"' EXIT
        (
            python {params.script:q} \
                {input.pileup:q} \
                "$temporary_directory/H1.bedgraph" \
                "$temporary_directory/H2.bedgraph" \
                "$temporary_directory/difference.bedgraph"
            first_chromosome="$(awk 'NF >= 2 {{print $1; exit}}' {input.chrom_sizes:q})"
            if [[ -z "$first_chromosome" ]]; then
                printf '%s\n' 'Chromosome sizes file is empty.' >&2
                exit 1
            fi
            raw_tracks=(
                "$temporary_directory/H1.bedgraph"
                "$temporary_directory/H2.bedgraph"
                "$temporary_directory/difference.bedgraph"
            )
            sorted_tracks=(
                "$temporary_directory/H1.sorted.bedgraph"
                "$temporary_directory/H2.sorted.bedgraph"
                "$temporary_directory/difference.sorted.bedgraph"
            )
            bigwigs=({output.h1:q} {output.h2:q} {output.difference:q})
            for track_index in "${{!raw_tracks[@]}}"; do
                if [[ ! -s "${{raw_tracks[$track_index]}}" ]]; then
                    printf '%s\t0\t1\t0\n' "$first_chromosome" \
                        > "${{raw_tracks[$track_index]}}"
                fi
                bedtools sort \
                    -g {input.chrom_sizes:q} \
                    -i "${{raw_tracks[$track_index]}}" \
                    > "${{sorted_tracks[$track_index]}}"
                bedGraphToBigWig \
                    "${{sorted_tracks[$track_index]}}" \
                    {input.chrom_sizes:q} \
                    "${{bigwigs[$track_index]}}"
            done
        ) > {log:q} 2>&1
        test -s {output.h1:q}
        test -s {output.h2:q}
        test -s {output.difference:q}
        """


rule fire_haplotype_differences:
    input:
        pileup=FIRE_PILEUP,
        index=FIRE_PILEUP_INDEX
    output:
        bed=FIRE_HAPLOTYPE_DIFFERENCES,
        index=FIRE_HAPLOTYPE_DIFFERENCES_INDEX,
        summary=FIRE_HAPLOTYPE_SUMMARY
    threads: 4
    resources:
        mem_mb=16000,
        runtime=60
    log:
        "logs/fire/haplotypes/{sample}.differences.log"
    benchmark:
        "benchmarks/fire/haplotypes/{sample}.differences.tsv"
    params:
        script=workflow.source_path("../scripts/fire_haplotype_differences.py"),
        minimum_coverage=FIRE_HAPLOTYPE_MIN_COVERAGE,
        maximum_fire_fdr=FIRE_HAPLOTYPE_MAX_FDR,
        qvalue_threshold=FIRE_HAPLOTYPE_QVALUE
    conda:
        "../../envs/intervals.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.bed:q})" "$(dirname {log:q})"
        temporary_bed="$(
            mktemp "${{TMPDIR:-/tmp}}/{wildcards.sample}.hapdiff.XXXXXX.bed"
        )"
        trap 'rm -f "$temporary_bed"' EXIT
        (
            python {params.script:q} \
                --minimum-coverage {params.minimum_coverage} \
                --maximum-fire-fdr {params.maximum_fire_fdr} \
                --qvalue-threshold {params.qvalue_threshold} \
                {input.pileup:q} \
                "$temporary_bed" \
                {output.summary:q}
            bgzip -@ {threads} -c "$temporary_bed" > {output.bed:q}
            tabix -f -p bed {output.bed:q}
        ) > {log:q} 2>&1
        test -s {output.bed:q}
        test -s {output.index:q}
        test -s {output.summary:q}
        """
