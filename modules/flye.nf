process flye {
    label 'flye'
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "${name}_raw_assembly.fasta"
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "flye.log"
    input:
    tuple val(name), file(ont), file(genome_size)
    output:
    tuple val(name), file(ont), file("${name}_raw_assembly.fasta")
    tuple val(name), file("flye.log")
    shell:
    """
    size=\$(cat !{genome_size})
    if [[ "\$(echo \$size | awk 'BEGIN{FS="."}{print \$1}')" == "0" ]]; then
        size=1m
    fi

    flye --nano-corr !{ont} -o flye_output -t !{task.cpus} --plasmids --meta --genome-size \$size
    mv flye_output/assembly.fasta ${name}_raw_assembly.fasta
    mv flye_output/flye.log flye.log
    """

}