// not used at the moment, low abundace kmers are also trimmed via the estimate_gsize process
process trim_low_abund {
    label 'khmer'

    input:
        tuple val(name), file(reads)

    output:
        tuple val(name), file("*.abundtrim")
    
    shell:
    """
    trim-low-abund.py -k 25 -C 3 -M ${params.maxmem} ${reads}
    """
}

process estimate_gsize {
    label 'kmc'
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "genome_size.txt"

    // gess.py seems to throw an error that is not recognized here and so the workflow continues with wrong gsize estimation
    // in that case, simply execute a second time
    errorStrategy 'retry'
    maxRetries 1

    input:
        tuple val(name), file(reads)

    output:
        tuple val(name), file('genome_size.txt')
        

    shell:
    """
    if [[ "${params.gsize}" != "" ]]; then
        echo "${params.gsize}m" > genome_size.txt
    else
        # maybe change this to /scratch on HPC
        TMP=/tmp
        gess.py --threads ${task.cpus} --cutoff 3 ${reads} -t \$TMP | awk 'BEGIN{FS=" "};{print \$5"m"}' > genome_size.txt

        cp genome_size.txt estimated_genome_size.txt
        size=\$(cat genome_size.txt)
        if [ \$(echo \$size | awk 'BEGIN{FS="."}{print \$1}') -lt 1 ]; then
            echo '1m' > genome_size.txt
        fi
    fi
    """
}