process spades {
        label 'spades'
        publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "${name}_assembly.fasta.gz"
    input:
        tuple val(name), file(ont), file(illumina)
    output:
        tuple val(name), file("${name}_assembly.fasta.gz")
    script:
        """
        spades.py --only-assembler -1 ${illumina[0]} -2 ${illumina[1]} --meta --nanopore ${ont} -o spades_output -t ${task.cpus}
        mv spades_output/contigs.fasta  ${name}_assembly.fasta
        gzip -f ${name}_assembly.fasta
        """
}