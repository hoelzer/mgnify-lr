process removeSmallReads {
    label 'basics'
  input:
    tuple val(name), file(reads) 
  output:
	  tuple val(name), file("${name}_filtered.fastq.gz") 
  shell:
    """
    case "!{reads}" in
      *.fastq.gz ) 
        zcat !{reads} | paste - - - - | awk -F"\\t" 'length(\$2)  >= 500' | sed 's/\\t/\\n/g' | gzip > "!{name}_filtered.fastq.gz"
        ;;
      *.fastq)
        cat !{reads} | paste - - - - | awk -F"\\t" 'length(\$2)  >= 500' | sed 's/\\t/\\n/g' | gzip > "!{name}_filtered.fastq.gz"
        ;;
    esac   
    """
}

/* Comments:
This is a super fast process to remove short reads.

it can take .fastq or .fastq.gz
*/