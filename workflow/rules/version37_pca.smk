rule all:
    input:
        expand("results/1000G/1000G_chr{contig}.bcf",contig=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                                     11, 12, 13, 14, 15, 16, 17, 18,
                                                     19, 20, 21, 22]),
        "resources/fasta/README.human_g1k_v37.fasta.txt",
        "resources/fasta/human_g1k_v37.fasta.fai",

rule download_chromosomes:
    output:
        bcf1000G= config["location_1000G"]+"ALL.chr{contig}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf",
        bcf1000G_csi="resources/1000G/ALL.chr{contig}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf.csi",
    resources: cpus=1, mem_mb=3000, time_job=720
    params:
        partition='batch'
    shell:
        """
        wget -P resources/1000G/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files/ALL.chr{wildcards.contig}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf
        wget -P resources/1000G/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files/ALL.chr{wildcards.contig}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf.csi
        """

rule download_fasta_files:
    output:
        "resources/fasta/README.human_g1k_v37.fasta.txt",
        "resources/fasta/human_g1k_v37.fasta.fai",
        fasta_gz="resources/fasta/human_g1k_v37.fasta.gz"
        fasta_file= "resources/fasta/human_g1k_v37.fasta"
    resources: cpus=1, mem_mb=3000, time_job=720
    params:
        partition='batch',
    shell:
        """
        wget -P resources/fasta/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/README.human_g1k_v37.fasta.txt
        wget -P resources/fasta/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
        wget -P resources/fasta/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz

        gunzip -c {output.fasta_gz} > {output.fasta_file}
        """

rule prepare_1000G_for_ancestry_PCA_step1:
    input:
        bcf1000G=config["location_1000G"]+"ALL.chr{contig}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf",
        fasta=config["fasta"]
    output:
        bcf="results/1000G/1000G_chr{contig}.bcf",
    resources: cpus=1, mem_mb=18000, time_job=720
    conda: "envs/bcftools.yaml"
    conda: "envs/bcftools.yaml"
    params:
        partition='batch',
        maf1=config["pca_1000G"]["bcf_maf"],
    shell:
        """
        if bcftools view -q {params.maf1}:minor "{input.bcf1000G}" | \
        bcftools norm -m-any --check-ref w -f "{input.fasta}" | \
        bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        bcftools norm -Ob --rm-dup both \
        > {output.bcf} ; then
        echo "no error"
        fi

        bcftools index {output.bcf}
        """