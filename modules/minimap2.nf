process minimap2_index_ont {
  label 'minimap2'
  if (params.cloudProcess) { 
      publishDir "${params.cloudDatabase}/minimap2/", mode: 'copy', pattern: "*.mmi" 
  }
  else { 
      storeDir "nextflow-autodownload-databases/minimap2/" 
  }  
    input:
  	  tuple val(name), file(fasta) 
    output:
      file("${name}.ont.mmi") 
    script:
      """
      minimap2 -x map-ont -t ${task.cpus} -d ${name}.ont.mmi ${fasta}
      """
}

process minimap2_index_ill {
  label 'minimap2'
  if (params.cloudProcess) { 
      publishDir "${params.cloudDatabase}/minimap2/", mode: 'copy', pattern: "*.mmi" 
  }
  else { 
      storeDir "nextflow-autodownload-databases/minimap2/" 
  }  
    input:
  	  tuple val(name), file(fasta) 
    output:
      file("${name}.ill.mmi") 
    script:
      """
      minimap2 -x sr -t ${task.cpus} -d ${name}.ill.mmi ${fasta}
      """
}

process minimap2_index_fna {
  label 'minimap2'
  if (params.cloudProcess) { 
      publishDir "${params.cloudDatabase}/minimap2/", mode: 'copy', pattern: "*.mmi" 
  }
  else { 
      storeDir "nextflow-autodownload-databases/minimap2/" 
  }  
    input:
  	  tuple val(name), file(fasta) 
    output:
      file("${name}.fna.mmi") 
    script:
      """
      minimap2 -x asm5 -t ${task.cpus} -d ${name}.fna.mmi ${fasta}
      """
}


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

process minimap2_to_decontaminate_fastq {
  label 'minimap2'
  publishDir "${params.output}/${name}/decontamination/", mode: 'copy', pattern: "${name}.*.fastq.gz"  

  input: 
    tuple val(name), file(fastq)
    file(db)

  output:
    tuple val(name), file("*.clean.fastq.gz")
    tuple val(name), file("*.contamination.fastq.gz")

  script:
    """

    # remove spaces in read IDs to keep them in the later cleaned output
    if [[ ${fastq} =~ \\.gz\$ ]]; then
      zcat ${fastq} | sed 's/ /DECONTAMINATE/g' > ${name}.id.fastq
    else
      sed 's/ /DECONTAMINATE/g' ${fastq} > ${name}.id.fastq
    fi

    minimap2 -a -t ${task.cpus} -o ${name}.sam ${db} ${name}.id.fastq
    samtools fastq -f 4 -0 ${name}.clean.id.fastq ${name}.sam
    samtools fastq -F 4 -0 ${name}.contamination.id.fastq ${name}.sam

    sed 's/DECONTAMINATE/ /g' ${name}.clean.id.fastq | gzip > ${name}.clean.fastq.gz
    sed 's/DECONTAMINATE/ /g' ${name}.contamination.id.fastq | gzip > ${name}.contamination.fastq.gz
     
    rm ${name}.sam ${name}.clean.id.fastq ${name}.contamination.id.fastq ${name}.id.fastq
    """
}


process minimap2_to_decontaminate_fasta {
  label 'minimap2'
  publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "*.fasta.gz"  

  input: 
    tuple val(name), file(fasta)
    file(db)

  output:
    tuple val(name), file("*clean.fasta.gz")
    tuple val(name), file("*contamination.fasta.gz")

  script:
    """

    # remove spaces in fasta IDs to keep them in the later cleaned output
    if [[ ${fasta} =~ \\.gz\$ ]]; then
      zcat ${fasta} | sed 's/ /DECONTAMINATE/g' > ${name}.id.fasta
    else
      sed 's/ /DECONTAMINATE/g' ${fasta} > ${name}.id.fasta
    fi

    minimap2 -a -t ${task.cpus} -o ${name}.sam ${db} ${name}.id.fasta
    samtools fasta -f 4 -0 ${name}.clean.id.fasta ${name}.sam
    samtools fasta -F 4 -0 ${name}.contamination.id.fasta ${name}.sam


    sed 's/DECONTAMINATE/ /g' ${name}.clean.id.fasta | gzip > ${name}_polished_clean.fasta.gz
    sed 's/DECONTAMINATE/ /g' ${name}.contamination.id.fasta | gzip > ${name}_polished_contamination.fasta.gz
     
    rm ${name}.sam ${name}.clean.id.fasta ${name}.contamination.id.fasta ${name}.id.fasta
    """
}


process minimap2_to_decontaminate_ill {
  label 'minimap2'
  publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "*.fastq.gz"  

  input: 
    tuple val(name), file(reads)
    file(db)

  output:
    tuple val(name), file("*clean*.fastq.gz")
    tuple val(name), file("*contamination*.fastq.gz")

  script:
    """
    # replace the space in the header to retain the full read IDs after mapping (the mapper would split the ID otherwise after the first space)
    # this is working for ENA reads that have at the end of a read id '/1' or '/2'
    EXAMPLE_ID=\$(zcat ${reads[0]} | head -1)
    if [[ \$EXAMPLE_ID == */1 ]]; then 
      if [[ ${reads[0]} =~ \\.gz\$ ]]; then
        zcat ${reads[0]} | sed 's/ /DECONTAMINATE/g' > ${name}.R1.id.fastq
      else
       sed 's/ /DECONTAMINATE/g' ${reads[0]} > ${name}.R1.id.fastq
     fi
     if [[ ${reads[1]} =~ \\.gz\$ ]]; then
       zcat ${reads[1]} | sed 's/ /DECONTAMINATE/g' > ${name}.R2.id.fastq
     else
       sed 's/ /DECONTAMINATE/g' ${reads[1]} > ${name}.R2.id.fastq
     fi
    else
      # this is for paried-end SRA reads that don't follow the ENA pattern
      if [[ ${reads[0]} =~ \\.gz\$ ]]; then
        zcat ${reads[0]} > ${name}.R1.id.fastq
        zcat ${reads[1]} > ${name}.R2.id.fastq
      else
        cp ${reads[0]} ${name}.R1.id.fastq
        cp ${reads[1]} ${name}.R2.id.fastq
      fi
    fi

    # Use samtools -F 2 to discard only reads mapped in proper pair:
    minimap2 -ax sr -t ${task.cpus} -o ${name}.sam ${db} ${name}.R1.id.fastq ${name}.R2.id.fastq
    samtools fastq -F 2 -1 ${name}.clean.R1.id.fastq -2 ${name}.clean.R2.id.fastq ${name}.sam
    samtools fastq -f 2 -1 ${name}.contamination.R1.id.fastq -2 ${name}.contamination.R2.id.fastq ${name}.sam
    
    # restore the original read IDs
    sed 's/DECONTAMINATE/ /g' ${name}.clean.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | gzip > ${name}.clean.R1.fastq.gz 
    sed 's/DECONTAMINATE/ /g' ${name}.clean.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | gzip > ${name}.clean.R2.fastq.gz
    sed 's/DECONTAMINATE/ /g' ${name}.contamination.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | gzip > ${name}.contamination.R1.fastq.gz 
    sed 's/DECONTAMINATE/ /g' ${name}.contamination.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | gzip > ${name}.contamination.R2.fastq.gz

    # remove intermediate files
    rm ${name}.R1.id.fastq ${name}.R2.id.fastq ${name}.clean.R1.id.fastq ${name}.clean.R2.id.fastq ${name}.contamination.R1.id.fastq ${name}.contamination.R2.id.fastq ${name}.sam

    """
}
