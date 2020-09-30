/*
Sort BAM and index
*/
process samtools {
        label 'samtools'  
    input:
        tuple val(name), path(bam)
    output:
        tuple val(name), path("${name}.sorted.bam"), path("${name}.sorted.bam.bai")
    script:
        """
        samtools sort -@ ${task.cpus} ${bam} > ${name}.sorted.bam
        samtools index -@ ${task.cpus} ${name}.sorted.bam
        """
  }