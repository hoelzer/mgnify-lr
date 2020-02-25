  
process diamond {
      label 'diamond'
    input:
      tuple val(name), file(proteins)
      file(database)
    output:
      tuple val(name), file("${name}.data")
    script:
      """
      diamond blastp --threads ${task.cpus} --max-target-seqs 1 --db ${database} \
      --query ${proteins} --outfmt 6 qlen slen --out ${name}.data
      """
}