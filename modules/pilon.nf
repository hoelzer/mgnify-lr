process pilon {
  label 'pilon'
  publishDir "${params.output}/${name}/assembly/", mode: 'copy', pattern: "${name}_polished_${round}.fasta"
      input:
        tuple val(name), path(assembly), path(bam), path(index) 
        val(round)
      output:
  	    tuple val(name), file("${name}_polished_${round}.fasta") 
      script:
        """
        VALUE=\$(echo ${task.memory} | awk '{print \$2}')
        if [[ \$VALUE == "GB" ]]; then
          MEM=\$(echo ${task.memory} | sed 's/ GB//g')
          VALUE=g
        else
          MEM=\$(echo ${task.memory} | sed 's/ TB//g')
          VALUE=t
        fi
        #only allow pilon to use 80% of the available memory
        ADJUSTEDMEM=\$(echo 0.8*\$MEM | bc | awk 'BEGIN{FS="."};{print \$1}')
        pilon -Xmx\${ADJUSTEDMEM}\${VALUE} --threads ${task.cpus} --genome ${assembly} --frags ${bam} --output ${name}_polished_${round}
      	"""
}

/*
MEM=$(echo 1.1 TB | sed 's/ TB//g')
MEM=$(echo 720 GB | sed 's/ GB//g')
*/