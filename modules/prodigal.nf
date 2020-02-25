process prodigal {
    label 'prodigal'

    input:
      tuple val(name), file(fasta)
    
    output:
      tuple val(name), file("${name}.faa")
    
    script:
      """
      prodigal -p meta -a ${name}.faa -q -i ${fasta} > ${name}.log
      """
} 