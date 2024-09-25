configfile: '././config/config.yaml'

rule all_pca37:
    input:
        expand("results/H_1000G_PCA/H_step1/1000G_chr{contig}.bcf",contig=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                                     11, 12, 13, 14, 15, 16, 17, 18,
                                                     19, 20, 21, 22]),
        "resources/H_1000G/H_Fasta/README.human_g1k_v37.fasta.txt",
        "resources/H_1000G/H_Fasta/human_g1k_v37.fasta.fai",

rule H_Download_1000G_chromosomes:
    output:
        bcf1000G= config["location_1000G"]+"ALL.chr{contig}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf",
        bcf1000G_csi="resources/H_1000G/H_Genotypes/ALL.chr{contig}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf.csi",
    resources: cpus=1, mem_mb=3000, time_job=720
    params:
        partition='batch'
    shell:
        """
        wget -P resources/H_1000G/H_Genotypes/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files/ALL.chr{wildcards.contig}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf
        wget -P resources/H_1000G/H_Genotypes/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files/ALL.chr{wildcards.contig}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf.csi
        """

rule H_Download_fasta_files:
    output:
        "resources/H_1000G/H_Fasta/README.human_g1k_v37.fasta.txt",
        "resources/H_1000G/H_Fasta/human_g1k_v37.fasta.fai",
        fasta_gz="resources/H_1000G/H_Fasta/human_g1k_v37.fasta.gz",
    resources: cpus=1, mem_mb=3000, time_job=720
    params:
        partition='batch',
    shell:
        """
        wget -P resources/H_1000G/H_Fasta ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/README.human_g1k_v37.fasta.txt
        wget -P resources/H_1000G/H_Fasta ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
        wget -P resources/H_1000G/H_Fasta ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
        """

rule H_Unzip_fasta:
    input:
        fasta_gz="resources/H_1000G/H_Fasta/human_g1k_v37.fasta.gz"
    output:
        fasta_file="resources/H_1000G/H_Fasta/human_g1k_v37.fasta"
    shell:
        """
        set +e
        gunzip -c {input.fasta_gz} > {output.fasta_file}
        exit 0
        """


rule H_Prepare_1000G_for_ancestry_PCA_step1:
    input:
        bcf1000G=config["location_1000G"]+"ALL.chr{contig}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf",
        fasta="resources/H_1000G/H_Fasta/human_g1k_v37.fasta"
    output:
        bcf="results/H_1000G_PCA/H_step1/1000G_chr{contig}.bcf",
    resources: cpus=1, mem_mb=18000, time_job=720
    conda: "../envs/bcftools.yaml"
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
