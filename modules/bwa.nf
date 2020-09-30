/* general bwa command for paired-end reads, generating a BAM */
process bwa {
    label "bwa"
    input:
        tuple val(name), path(reads), path(fasta)
    output:
        tuple val(name), path("${name}.bam")
    script:
        """
        bwa index ${fasta}
        bwa mem -t ${task.cpus} ${fasta} ${reads[0]} ${reads[1]} | samtools view -bS - > ${name}.bam
        """
}