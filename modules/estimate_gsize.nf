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

process gess_gsize {
    label 'kmc'
    publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "genome_size.txt"

    input:
        tuple val(name), file(reads)

    output:
        tuple val(name), file('genome_size.txt')
        

    shell:
    """
    TMP=/tmp
    if [[ "${workflow.profile}" == "ebi" ]]; then 
        TMP=/scratch
    fi
    gess.py --threads ${task.cpus} --cutoff 3 ${reads} -t \$TMP | awk 'BEGIN{FS=" "};{print \$5"m"}' > genome_size.txt

    size=\$(cat genome_size.txt)
    if [[ \$(echo \$size | awk 'BEGIN{FS="."}{print \$1}') < 101 ]]; then
        echo '100m' > genome_size.txt
    fi
    """
}