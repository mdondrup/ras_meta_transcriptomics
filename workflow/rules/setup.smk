# Rule: Setup (create directories)
rule setup:
    input: "rRNA_databases_v4/"
    output: ".setup_done"
    shell:
        r"""
        mkdir -p fastp sortmerna aligned misc multiqc_report \
        featurecounts results benchmarks logs \
        && touch .setup_done
        """


rule dowload_sortmeRNA_reference:
    output:
        directory("rRNA_databases_v4")
    shell:
        r"""
        wget -c https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
        mkdir -p rRNA_databases_v4
        tar -xf database.tar.gz -C rRNA_databases_v4
        
        """
