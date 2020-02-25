process prodigal {
    label 'prodigal'

    input:
      tuple val(name), file(fasta)
    
    output:
      tuple val(name), env(ASSEMBLY_STATUS), file("${name}.faa")
    
    script:
      """
      #ERR3662306_1_small_raw_assembly.fasta
      #ERR3662306_1_small_polished.fasta

      BN=\$(basename ${fasta} .fasta)
      ASSEMBLY_STATUS=\$(echo \$BN | sed 's/${name}_//g') # [raw_assembly, polished]

      prodigal -p meta -a ${name}.faa -q -i ${fasta} > ${name}.log
      """
} 