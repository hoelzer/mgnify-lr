process minimap2_to_polish {
  label 'minimap2'
    input:
  	  tuple val(name), file(read), file(assembly) 
    output:
      tuple val(name), file(read), file(assembly), file("${name}.paf") 
    script:
      """
      minimap2 -x map-ont -t ${task.cpus} ${assembly} ${read} > ${name}.paf
      """
}
