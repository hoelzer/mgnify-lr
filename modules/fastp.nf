process fastp {
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "${name}.fastp.R*.fastq.gz"
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "${name}.qc.html"
      label 'fastp'
    input:
      tuple val(name), file(sread)
    output:
      tuple val(name), file("${name}.fastp.R*.fastq.gz")
      tuple val(name), file("${name}.qc.html") 
    script:
      """
      fastp -w ${task.cpus} --in1 ${sread[0]} --in2 ${sread[1]} --out1 "${name}.fastp.R1.fastq.gz" \
        --out2 "${name}.fastp.R2.fastq.gz" --json "${name}.qc.json" --html "${name}.qc.html" ${params.fastp_additional_params}
      """
}