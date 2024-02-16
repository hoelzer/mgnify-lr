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
      racon -m 8 -x -6 -g -8 -w 500 -q -1 -t ${task.cpus} ${read} ${mapping} ${assembly} > ${name}_consensus.fasta
      """
  }

// From: https://ncbi.nlm.nih.gov/pmc/articles/PMC9861289/
// For Racon, we used the ONT suggested parameters: score for matching bases (-m 8), 
// score for mismatching bases (-x -6), gap penalty (-g -8), window size (-w 500), and mean quality threshold for each window (-q -1).