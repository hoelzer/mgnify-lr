process sourmash_metagenome_size {
    label 'sourmash'
    publishDir "${params.output}/${name}_assembly/", mode: 'copy', pattern: "genome_size.txt"
    input:
    tuple val(name), file(ont)
    file(database)
    output:
    tuple val(name), file('genome_size.txt')
    shell:
    """
    #sourmash compute -p !{task.cpus} --scaled 10000 -k 31 !{ont} -o !{name}.sig 
    #sourmash lca gather  !{name}.sig !{database} --ignore-abundance -o metagenomic-composition.txt
    #sum_ont=\$(cat metagenomic-composition.txt | cut -d ',' -f 1 | paste -sd+ | bc)
    #total_m_ont=\$(bc -l <<< "scale=2 ; \$sum_ont /10^6")
    #echo \$total_m_ont"M" > genome_size.txt
    echo "100m" > genome_size.txt
    """
}