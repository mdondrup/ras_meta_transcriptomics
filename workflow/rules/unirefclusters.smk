
rule map_uniref:
    input:
        "results/uniref90_mapping.tsv",
        "results/reference_annotated.gff"

# 1. Download UniRef90 FASTA from UniProt
rule download_uniref90:
    output:
        "db/uniref90.fasta.gz"
    shell:
        """
        wget -O {output} ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
        """

# 2. Unzip UniRef90 FASTA
rule unzip_uniref90:
    input:
        "db/uniref90.fasta.gz"
    output:
        "db/uniref90.fasta"
    shell:
        """
        gunzip -c {input} > {output}
        """

# 3. Create MMseqs2 database from UniRef90 FASTA
rule create_uniref_db:
    conda: "../envs/mmseqs2.yaml"
    input:
        "db/uniref90.fasta"
    output:
        db="db/uniref90_mmseqs_db"
    shell:
        """
        mmseqs createdb {input} {output.db}
        """

# 4. Create MMseqs2 database from your query FASTA
rule create_query_db:
    conda: "../envs/mmseqs2.yaml"
    input:
        fasta="input/proteins.faa"
    output:
        db="results/query_mmseqs_db"
    shell:
        """
        mmseqs createdb {input.fasta} {output.db}
        """

# 5. Search query proteins against UniRef90
rule search_mmseqs:
    conda: "../envs/mmseqs2.yaml"
    input:
        query_db="results/query_mmseqs_db",
        target_db="db/uniref90_mmseqs_db"
    output:
        search_db_dbtype="results/search_db.dbtype",
        aln="results/mmseqs_result.m8"
    params:
        tmpdir="results/tmp",
        minident=0.9,
        coverage=0.9
    threads: 80
    log: "logs/mmseqs.log"
    message: "Mapping to UniRef90 with MMseqs2. See the log: {log}"     
    shell:
        """
        
        mmseqs search {input.query_db} {input.target_db} results/search_db {params.tmpdir} \
            --min-seq-id {params.minident} \
            --cov-mode 1 -c {params.coverage} \
            --max-seqs 1 \
            --threads {threads} > {log} 2>&1
        mmseqs convertalis {input.query_db} {input.target_db} results/search_db {output.aln} \
        --format-output "query,target,fident,qalnlen,qlen" > {log} 2>&1
        """

# 6. Format the results as a simple mapping table
rule format_mmseqs_output:
    input:
        "results/mmseqs_result.m8"
    output:
        "results/uniref90_mapping.tsv"
    shell:
        """
        awk '{{print $1 "\\t" $2}}' {input} > {output}
        """

#7. Add uniref cluster to the GFF file        
rule annotate_gff_with_uniref:
    input:
        gff="input/reference/reference.gff",
        mapping="results/uniref90_mapping.tsv"
    output:
        gff="results/reference_annotated.gff"
    message: "Adding uniref attribute to gff file"
    run:
        # Load UniRef90 mapping into a dictionary
        mapping = {}
        with open(input.mapping) as map_file:
            for line in map_file:
                prot_id, uniref = line.strip().split("\t")
                mapping[prot_id] = uniref

        with open(input.gff) as fin, open(output.gff, "w") as fout:
            for line in fin:
                if line.startswith("#"):
                    fout.write(line)
                    continue
                cols = line.strip().split("\t")
                if len(cols) != 9 or cols[2] != "CDS":
                    fout.write(line)
                    continue

                # Extract ID attribute
                attributes = cols[8]
                attrs = {k: v for k, v in [field.split("=", 1) for field in attributes.split(";") if "=" in field]}
                prot_id = attrs.get("ID") or attrs.get("protein_id")

                # Annotate with UniRef90 if available
                if prot_id and prot_id in mapping:
                    uniref_tag = f"uniref90={mapping[prot_id]}"
                    if "uniref90=" not in attributes:
                        attributes += ";" + uniref_tag
                cols[8] = attributes
                fout.write("\t".join(cols) + "\n")
