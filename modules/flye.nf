process flye {
    label 'flye'
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "${name}_raw_flye.fasta"
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "flye.log"
    
    input:
    tuple val(name), file(ont)
    
    output:
    tuple val(name), file(ont), file("${name}_raw_flye.fasta"), emit: assembly
    tuple val(name), file("flye.log"), emit: log
    
    script:
    """
    flye --nano-raw ${ont} -o flye_output -t ${task.cpus} --plasmids --meta
    mv flye_output/assembly.fasta ${name}_raw_flye.fasta
    mv flye_output/flye.log flye.log
    """

}