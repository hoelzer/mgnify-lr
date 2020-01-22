process minimap2_to_polish {
  label 'minimap2'
    input:
  	  tuple val(name), file(read), file(assembly) 
    output:
      tuple val(name), file(read), file(assembly), file("${name}.paf") 
    script:
      """
      minimap2 -x map-ont -t ${task.cpus} ${assembly} ${read} > ${name}.paf
      """
}

process minimap2_to_decontaminate {
  label 'minimap2'
  publishDir "${params.output}/${name}/decontamination/", mode: 'copy', pattern: "${name}.*.fastq.gz"  

  input: 
    tuple val(name), file(fastq)
    file(db)

  output:
    tuple val(name), file("*.clean.fastq.gz")

  script:
    """

    # remove spaces in read IDs to keep them in the later cleaned output
    if [[ ${fastq} =~ \\.gz\$ ]]; then
      zcat ${fastq} | sed 's/ /DECONTAMINATE/g' > ${name}.id.fastq
    else
      sed 's/ /DECONTAMINATE/g' ${fastq} > ${name}.id.fastq
    fi

    minimap2 -ax map-ont -t ${task.cpus} -o ${name}.sam ${db} ${name}.id.fastq
    samtools fasta -f 4 -0 ${name}.clean.id.fastq ${name}.sam
    samtools fasta -F 4 -0 ${name}.contamination.id.fastq ${name}.sam

    sed 's/DECONTAMINATE/ /g' ${name}.clean.id.fastq | gzip > ${name}.clean.fastq.gz
    sed 's/DECONTAMINATE/ /g' ${name}.contamination.id.fastq | gzip > ${name}.contamination.fastq.gz
     
    rm ${name}.sam ${name}.clean.id.fastq ${name}.contamination.id.fastq ${name}.id.fastq
    """
}

