process medaka {
  label 'medaka'
  publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "${name}_polished.fasta.gz"
      input:
        tuple val(name), file(read), file(consensus) 
      output:
  	    tuple val(name), file("${name}_polished.fasta.gz") 
      script:
        """
      	medaka_consensus -i ${read} -d ${consensus} -o polished -t ${task.cpus} -m ${params.model}
        mv polished/consensus.fasta ${name}_polished.fasta
        gzip -f ${name}_polished.fasta
      	"""
}

