  
process diamond {
      label 'diamond'
    input:
      tuple val(name), val(ASSEMBLY_STATUS), file(proteins)
      file(database)
    output:
      tuple val(name), val(ASSEMBLY_STATUS), file("${name}.data")
    script:
      """
      diamond blastp --threads ${task.cpus} --db ${database} \
      --query ${proteins} --outfmt 6 qlen slen --out ${name}.data
      """
}