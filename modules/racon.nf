process racon {
      label 'racon'

    errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
    maxRetries 2

   input:
      tuple val(name), file(read), file(assembly), file(mapping) 
   output:
   	tuple val(name), file(read), file("${name}_consensus.fasta") 
   shell:
      """
      racon -t ${task.cpus} ${read} ${mapping} ${assembly} > ${name}_consensus.fasta
      """
  }