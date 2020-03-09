process racon {
      label 'racon'

    errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
    cpus { 24 }
    memory { 60.GB * task.attempt }
    clusterOptions { '-P bigmem' }
    maxRetries 2

   input:
      tuple val(name), file(read), file(assembly), file(mapping) 
   output:
   	tuple val(name), file(read), file("${name}_racon.fasta") 
   shell:
      """
      racon -t ${task.cpus} ${read} ${mapping} ${assembly} > ${name}_racon.fasta
      """
  }