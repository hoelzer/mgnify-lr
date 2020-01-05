process nanoplot {
    label 'nanoplot'
      publishDir "${params.output}/${name}/", mode: 'copy', pattern: "*.html"
      publishDir "${params.output}/${name}/readQCdir/", mode: 'copy', pattern: "*_read_quality.txt"
      publishDir "${params.output}/${name}/readQCdir/", mode: 'copy', pattern: "*.png"
      publishDir "${params.output}/${name}/readQCdir/", mode: 'copy', pattern: "*.pdf"
    input:
      tuple val(name), file(reads)
    output:
      tuple val(name), file("*.html"), file("*.png"), file("*.pdf"), file("${name}_read_quality.txt") 
    script:
      """
      NanoPlot -t ${task.cpus} --fastq ${reads} --title ${name} --color darkslategrey --N50 --plots hex --loglength -f png --store
      NanoPlot -t ${task.cpus} --pickle NanoPlot-data.pickle --title ${name} --color darkslategrey --N50 --plots hex --loglength -f pdf
      mv *.html ${name}_read_quality_report.html
      mv NanoStats.txt ${name}_read_quality.txt
      """
}

/* Comments:
We run nanoplot 2 times to get png and pdf files.
The second time its done via the pickle file of the previous run, to save computing time
*/