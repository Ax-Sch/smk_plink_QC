configfile: '././config/config.yaml'

rule all_pca38:
    input:
        expand("results/1000G/1000G_chr{contig}.bcf",contig=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                                     11, 12, 13, 14, 15, 16, 17, 18,
                                                     19, 20, 21, 22]),
        "resources/fasta/GRCh38_full_analysis_set_plus_decoy_hla.dict",
        "resources/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai",
        "resources/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann"

rule download_chromosomes:
    output:
        vcf1000G= config["location_1000G"]+ "ALL.chr{contig}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz",
        vcf1000G_tbi=config["location_1000G"]+ "ALL.chr{contig}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi",
    resources: cpus=1, mem_mb=3000, time_job=720
    params:
        partition='batch'
    shell:
        """
        wget -P resources/1000G/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr{wildcards.contig}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
        wget -P resources/1000G/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr{wildcards.contig}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
        """

rule download_fasta_files:
    output:
        "resources/fasta/GRCh38_full_analysis_set_plus_decoy_hla.dict",
        "resources/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai",
        "resources/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann",
        "resources/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa",
    resources: cpus=1, mem_mb=3000, time_job=720
    params:
        partition='batch',
    shell:
        """
        wget -P resources/fasta/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.dict
        wget -P resources/fasta/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
        wget -P resources/fasta/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann
        wget -P resources/fasta/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
        """

rule prepare_1000G_for_ancestry_PCA_step1:
    input:
        vcf1000G=config["location_1000G"]+ "ALL.chr{contig}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz",
        fasta="resources/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    output:
        bcf="results/1000G/1000G_chr{contig}.bcf",
    resources: cpus=1, mem_mb=18000, time_job=720
    conda: "../envs/bcftools.yaml"
    params:
        partition='batch',
        maf1=config["pca_1000G"]["bcf_maf"],
    shell:
        """
        if bcftools view -q {params.maf1}:minor "{input.vcf1000G}" | \
        bcftools norm -m-any --check-ref w -f "{input.fasta}" | \
        bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        sed 's/^chr//; s/\\tchr/\\t/' | \
        sed 's/##contig=<ID=chr/##contig=<ID=/' | \
        bcftools norm -Ob --rm-dup both \
        > {output.bcf} ; then
        echo "no error"
        fi

        bcftools index {output.bcf}
        """
