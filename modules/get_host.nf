/*Comment section: */

process get_host {
  label 'basics'
  if (params.cloudProcess) { 
      publishDir "${params.cloudDatabase}/hosts/${params.species}", mode: 'copy', pattern: "*.fa.gz" 
  }
  else { 
      storeDir "nextflow-autodownload-databases/hosts/${params.species}" 
  }  

  output:
      file("*.fa.gz")

  script:
    """
    if [ ${params.species} == 'hsa' ]; then
      wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
      mv *.fa.gz ${params.species}.fa.gz
    fi
    if [ ${params.species} == 'mmu' ]; then
      wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
      mv *.fa.gz ${params.species}.fa.gz
    fi
    if [ ${params.species} == 'eco' ]; then
      wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//fasta/bacteria_90_collection/escherichia_coli_k_12/dna/Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz
      mv *.fa.gz ${params.species}.fa.gz
    fi


    """
}

