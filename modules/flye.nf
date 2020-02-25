process flye {
    label 'flye'
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "${name}_raw_assembly.fasta"
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "flye.log"

    //errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
    //cpus { task.cpus }
    //memory { task.memory * task.attempt }
    //maxRetries 3
    
    input:
    tuple val(name), file(ont), file(genome_size)
    
    output:
    tuple val(name), file(ont), file("${name}_raw_assembly.fasta")
    tuple val(name), file("flye.log")
    
    shell:
    """
    size=\$(cat !{genome_size})
    flye --nano-corr !{ont} -o flye_output -t !{task.cpus} --plasmids --meta --genome-size \$size
    mv flye_output/assembly.fasta ${name}_raw_assembly.fasta
    mv flye_output/flye.log flye.log
    """

}