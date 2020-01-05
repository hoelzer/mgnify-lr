process racon {
      label 'racon'
   input:
      tuple val(name), file(read), file(assembly), file(mapping) 
   output:
   	tuple val(name), file(read), file("${name}_consensus.fasta") 
   shell:
      """
      racon -t ${task.cpus} ${read} ${mapping} ${assembly} > ${name}_consensus.fasta
      """
  }