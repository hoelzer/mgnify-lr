process prodigal {
    label 'prodigal'

    input:
      tuple val(name), file(fasta)
    
    output:
      tuple val(name), env(ASSEMBLY_STATUS), file("${name}.faa")
    
    script:
      """
      #possible input channels
      #ERR3662306_1_small_raw_assembly.fasta
      #ERR3662306_1_small_polished.fasta
      #ERR3662306_1_small_pilon_polished.fasta
      #ERR3662306_1_small_raw_assembly_clean.fasta
      #ERR3662306_1_small_polished_clean.fasta
      #ERR3662306_1_small_pilon_polished_clean.fasta

      BN=\$(basename ${fasta} .fasta)
      ASSEMBLY_STATUS=\$(echo \$BN | sed 's/${name}_//g') # [raw_assembly, polished, pilon_polished, raw_assembly_clean, ...]

      prodigal -p meta -a ${name}.faa -q -i ${fasta} > ${name}.log
      """
} 