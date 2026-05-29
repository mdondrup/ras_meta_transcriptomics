REFS = config['metaquast_refs']

# Rule: Setup (create directories)



rule setup:
    input: dir="rRNA_databases_v4/smr_v4.3_default_db.fasta",
           data=expand("ncbi_dataset/data/{ref}.fna", ref=REFS)
    output: ".setup_done"
    shell:
        r"""
        echo {input.data}
        mkdir -p fastp sortmerna aligned misc multiqc_report \
        featurecounts results benchmarks logs 
        """


rule dowload_sortmeRNA_reference:
    output:
        "rRNA_databases_v4/smr_v4.3_default_db.fasta"
    shell:
        r"""
        wget -c https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
        mkdir -p rRNA_databases_v4
        tar -xf database.tar.gz -C rRNA_databases_v4
        
        """


        
        
rule download_metaquast_references:
    conda: "../envs/setup.yaml"
    output:
        expand("ncbi_dataset/data/{ref}.fna", ref=REFS)
    params:
        refs=REFS
    message: "downloading NCBI reference datasets ({params.refs}) for metaQUAST"     
    shell:
        r"""
        #mkdir -p ncbi_dataset/
        datasets download genome accession {params.refs}
        unzip -o ncbi_dataset.zip
        # We cannot know exact filenames beforehand 
        for ACC in {params.refs}
        do
          rm -f ncbi_dataset/data/$ACC.fna
          mv -nv ncbi_dataset/data/$ACC/$ACC*_genomic.fna ncbi_dataset/data/$ACC.fna
        done
        rm ncbi_dataset.zip
        """
        
