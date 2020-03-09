process polca {
    label 'masurca'  
    publishDir "${params.output}/${name}/", mode: 'copy', pattern: "${name}_polca.fasta"

    /*errorStrategy { 'retry' }
    cpus { 24 }
    memory { 36.GB * task.attempt }
    clusterOptions { '-P bigmem' }
    maxRetries 3*/

  input:
    tuple val(name), file(assembly), file(shortRead)

  output:
    tuple val(name), file("${name}_polca.fasta") 

  script:
    """
    MEM=\$(echo ${task.memory} | awk '{print \$2}' | sed 's/ GB//g')

    # replace the samtools part in the shell script
    #SCRIPT=\$(which polca.sh)
    #cp \$SCRIPT polca_refine.sh
    #OLD="\$SAMTOOLS sort -m \$MEM -@ \$NUM_THREADS <(samtools view -uhS \$BASM.unSorted.sam) \$BASM.alignSorted 2>>samtools.err && \\"
    #NEW="\$SAMTOOLS view -bS \$BASM.unSorted.sam -o \$BASM.unSorted.bam 2>>samtools.err && \$SAMTOOLS sort -m \$MEM -@ \$NUM_THREADS \$BASM.unSorted.bam \$BASM.alignSorted 2>>samtools.err && \\"

    polca_refine.sh -a ${assembly} -r '${shortRead}' -t ${task.cpus} -m \${MEM}G
    mv *.PolcaCorrected.fa ${name}_polca.fasta
    """
  }



  /*
  REPLACE
$SAMTOOLS sort -m $MEM -@ $NUM_THREADS <(samtools view -uhS $BASM.unSorted.sam) $BASM.alignSorted 2>>samtools.err && \
  WITH
$SAMTOOLS view -bS $BASM.unSorted.sam -o $BASM.unSorted.bam 2>>samtools.err && $SAMTOOLS sort -m $MEM -@ $NUM_THREADS $BASM.unSorted.bam $BASM.alignSorted 2>>samtools.err && \
  */