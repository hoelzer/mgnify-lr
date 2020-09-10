process spades {
    label 'spades'
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "${name}_raw_assembly.fasta"

    input:
        tuple val(name), file(ont), file(illumina)
    output:
        tuple val(name), file("${name}_raw_assembly.fasta")
    script:
        """
	SIZE=\$(echo ${task.memory} | awk 'BEGIN{FS=" "};{print \$2}')
	if [ \$SIZE == "TB" ]; then
          MEM=\$(echo ${task.memory} | sed 's/ TB//g') 
	  MEM=\$(echo \$MEM*1000 | bc | awk 'BEGIN{FS="."};{print \$1}')
	else
          MEM=\$(echo ${task.memory} | sed 's/ GB//g') 
	fi
        spades.py --only-assembler -1 ${illumina[0]} -2 ${illumina[1]} --nanopore ${ont} --meta -o spades_output -t ${task.cpus} -m \${MEM}
        mv spades_output/contigs.fasta  ${name}_raw_assembly.fasta
        """
}