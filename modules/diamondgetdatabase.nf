process diamond_download_db {
        if (params.cloudProcess) { publishDir "${params.cloudDatabase}/diamond", mode: 'copy', pattern: "database_uniprot.dmnd" }
        else { storeDir 'nextflow-autodownload-databases/diamond' }  
        label 'diamond' 
      output:
        file("database_uniprot.dmnd")
      script:
        """
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        diamond makedb --in uniprot_sprot.fasta.gz -d database_uniprot
        rm -f uniprot_sprot.fasta.gz
        """
    }

