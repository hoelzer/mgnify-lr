process flye {
    label 'flye'
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "${name}_raw_assembly.fasta"
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "flye.log"
    
    input:
    tuple val(name), file(ont)
    
    output:
    tuple val(name), file(ont), file("${name}_raw_assembly.fasta")
    tuple val(name), file("flye.log")
    
    shell:
    """
    # --nano-raw should be used for R9 data, --nano-hq for R10!
    
    flye --nano-raw ${ont} -o flye_output -t ${task.cpus} --meta
    mv flye_output/assembly.fasta ${name}_raw_assembly.fasta
    mv flye_output/flye.log flye.log
    """
}

// Genome size parameter is no longer required since the version 2.8.
