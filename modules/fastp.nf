process fastp {
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "${name}.fastp.R*.fastq"
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "${name}.qc.html"
      label 'fastp'
    input:
      tuple val(name), file(sread)
    output:
      tuple val(name), file("${name}.fastp.R*.fastq")
      tuple val(name), file("${name}.qc.html") 
    script:
      """
      fastp -w ${task.cpus} --in1 ${sread[0]} --in2 ${sread[1]} --out1 "${name}.fastp.R1.fastq" \
        --out2 "${name}.fastp.R2.fastq" --json "${name}.qc.json" --html "${name}.qc.html"
      """
}