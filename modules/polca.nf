process polca {
    label 'masurca'  
    publishDir "${params.output}/${name}/", mode: 'copy', pattern: "${name}_polca.fasta"

    /*errorStrategy { 'retry' }
    cpus { 24 }
    memory { 36 * task.attempt }
    clusterOptions { '-P bigmem' }
    maxRetries 3*/

  input:
    tuple val(name), file(assembly), file(shortRead)

  output:
    tuple val(name), file("${name}_polca.fasta") 

  script:
    """
    polca.sh -a ${assembly} -r '${shortRead}' -t ${task.cpus} -m ${task.memory}G
    mv *.PolcaCorrected.fa ${name}_polca.fasta
    """
  }