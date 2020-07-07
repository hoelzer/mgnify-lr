process spades {
    label 'spades'
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "${name}_raw_assembly.fasta"

    input:
        tuple val(name), file(ont), file(illumina)
    output:
        tuple val(name), file("${name}_raw_assembly.fasta")
    script:
        """
        MEM=\$(echo ${task.memory} | sed 's/ GB//g')
        spades.py --only-assembler -1 ${illumina[0]} -2 ${illumina[1]} --nanopore ${ont} --meta -o spades_output -t ${task.cpus} -m \${MEM}
        mv spades_output/contigs.fasta  ${name}_raw_assembly.fasta
        """
}