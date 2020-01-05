#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
Nextflow -- Metagenomics Pipeline
Author: christian@nanozoo.com
*/

/************************** 
* META & HELP MESSAGES 
**************************/

/* 
Comment section: First part is a terminal print for additional user information,
followed by some help statements (e.g. missing input) Second part is file
channel input. This allows via --list to alter the input of --nano & --illumina
to add csv instead. name,path   or name,pathR1,pathR2 in case of illumina 
*/

// terminal prints
if (params.help) { exit 0, helpMSG() }

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "
if (workflow.profile == 'standard') {
println "\033[2mCPUs to use: $params.cores"
println "Output dir name: $params.output\u001B[0m"
println " "}

if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if (((params.nano == '' &&  params.illumina == '') || params.illumina == '' ) && params.bins == '' ) 
    { exit 1, "input missing, use [--nano] and/or [--illumina]; or provide a bin directory [--bins]"}

/************************** 
* INPUT CHANNELS 
**************************/

// nanopore reads input & --list support
if (params.nano && params.list) { nano_input_ch = Channel
  .fromPath( params.nano, checkIfExists: true )
  .splitCsv()
  .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)] }
  .view() }
else if (params.nano) { nano_input_ch = Channel
  .fromPath( params.nano, checkIfExists: true)
  .map { file -> tuple(file.simpleName, file) }
  .view() }

// illumina reads input & --list support
if (params.illumina && params.list) { illumina_input_ch = Channel
  .fromPath( params.illumina, checkIfExists: true )
  .splitCsv()
  .map { row -> ["${row[0]}", [file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true)]] }
  .view() }
else if (params.illumina) { illumina_input_ch = Channel
  .fromFilePairs( params.illumina , checkIfExists: true )
  .view() }

// bin input channel (this skips over the whole metagenomic assembly part) -- list input not implemented
  /*
  if (params.bins && params.list) { bins_input_ch = Channel
    .fromPath( params.bins, checkIfExists: true )
    .splitCsv()
    .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)] }
    .view() }
  */
if (params.bins) { bins_input_ch = Channel
  .fromPath( "${params.bins}/*", checkIfExists: true)
  .map { file -> tuple(file.getParent().getName(), file) }
  .groupTuple() }


/************************** 
* NON-SAMPLE-READ CHANNELS 
**************************/

/* 
This special input section contains a csv input for other metagenomic reads
These reads are only to improve binning quality heavily
ont input is:  val, read
illumina is: val, R1, R2
IMPORTANT: val refers here to the actual sample to map against and NOT the name of these reads
*/

// extra ont reads
if (params.extra_ont) { extra_ont_ch = Channel
  .fromPath( params.extra_ont, checkIfExists: true )
  .splitCsv()
  .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)] }
  }
else { extra_ont_ch = file("noinput.fastq") }

// extra illumina
if (params.extra_ill) { extra_ill_ch = Channel
  .fromPath( params.extra_ill, checkIfExists: true )
  .splitCsv()
  .map { row -> ["${row[0]}", [file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true)]] }
  }
else { extra_ill_ch = file("noinput.fastq") }

/************************** 
* MODULES
**************************/

    include bwa_to_bam as bwa_bin from './modules/bwa'  
    include bwa_to_bam as bwa_to_bam_extra from './modules/bwa'
    include bwa_to_bam from './modules/bwa'
    include cat_fasta from './modules/cat_fasta'
    include checkm from './modules/checkm' params(output : params.output)
    include concoct from './modules/concoct' params(output : params.output)
    include concoct_extra from './modules/concoct' params(output : params.output)
    include concoct_extra_ill_only from './modules/concoct' params(output : params.output)
    include concoct_ill_only from './modules/concoct' params(output : params.output)
    include contig_IDs_in_bins from './modules/contig_IDs_in_bins' params( output : params.output)
    include fastp from './modules/fastp' params(output: params.output)
    include flye from './modules/flye' params(output: params.output, gsize: params.gsize)
    include maxbin2 from './modules/maxbin2' params(output : params.output)
    include maxbin2_ill_only from './modules/maxbin2' params(output : params.output)
    include medaka from './modules/medaka' params(output: params.output, model: params.model)
    include megahit from './modules/megahit' params( output : params.output)
    include metabat2 from './modules/metabat2' params(output : params.output)
    include metabat2_extra from './modules/metabat2' params(output : params.output)
    include metabat2_extra_ill_only from './modules/metabat2' params(output : params.output)
    include metabat2_ill_only from './modules/metabat2' params(output : params.output)
    include metawrap from './modules/metawrap' params(output : params.output)
    include minimap2_to_bam as minimap2_bin from './modules/minimap2'
    include minimap2_to_bam as minimap2_to_bam_extra from './modules/minimap2'
    include minimap2_to_bam from './modules/minimap2'
    include minimap2_to_polish from './modules/minimap2'
    include nanoplot from './modules/nanoplot' params(output: params.output)
    include pilon from './modules/pilon' params(output: params.output)
    include racon from './modules/racon'
    include reads_retrieval from './modules/seqtk_retrieve_reads' params(output : params.output)
    include reads_retrieval_unmapped from './modules/seqtk_retrieve_reads'params(output : params.output)
    include reads_retrieval_unmapped_ill_only from './modules/seqtk_retrieve_reads'params(output : params.output)
    include removeSmallReads from './modules/removeSmallReads' params(output: params.output)
    include sourmash_checkm_parser from './modules/parser/checkm_sourmash_parser'params(output: params.output)
    include sourmash_download_db from './modules/sourmashgetdatabase' params(cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
    include sourmash_metagenome_size from './modules/sourmash_metagenome_size' params(output: params.output)
    include sourmash_tax_classification from './modules/sourmash_tax_classification' params(output : params.output)
    include spades from './modules/spades' params( output : params.output)
    include spades_ill_only from './modules/spades' params( output : params.output)
    include unicycler from './modules/unicycler' params(output : params.output)

    include eggnog_download_db from './modules/eggnoggetdatabase' params(cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
    include eggnog_bin from './modules/eggnog'params(output : params.output)
    include kegg_parser from './modules/parser/KEGG_parser_no_rna'params(output: params.output)

/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow download_sourmash {
    main:
        if (params.sour_db) { database_sourmash = file(params.sour_db) }
        else if (!params.cloudProcess) { sourmash_download_db() ; database_sourmash = sourmash_download_db.out }
        else if (params.cloudProcess) { 
            sour_db_preload = file("${params.cloudDatabase}/sourmash/gtdb.lca.json")
            if (sour_db_preload.exists()) { database_sourmash = sour_db_preload }    
            else  { sourmash_download_db() ; database_sourmash = sourmash_download_db.out }
        }

    emit: database_sourmash
} 

workflow download_eggnog {
    main:
        if (params.egg_db) { database_eggnog = file(params.egg_db) }
        else if (!params.cloudProcess) { eggnog_download_db() ; database_eggnog = eggnog_download_db.out }
        else if (params.cloudProcess) { 
            eggnog_db_preload = file("${params.cloudDatabase}/eggnog/eggnog-db")
            if (eggnog_db_preload.exists()) { database_eggnog = eggnog_db_preload }
            else  { eggnog_download_db(); database_eggnog = eggnog_download_db.out }
        }

    emit: database_eggnog
} 

/************************** 
* SUB WORKFLOWS
**************************/

/* Hybrid Assembly Workflow */

workflow hybrid_assembly_wf {
  get:   nano_input_ch
          illumina_input_ch
          extra_ont_ch
          extra_ill_ch
          database_sourmash

  main:
      // trimming and QC of reads
        removeSmallReads(nano_input_ch)
        fastp(illumina_input_ch)
        nanoplot(nano_input_ch)
      // size estimation for flye // not working well - bc not installed n sourmash nanozoo container
        if (params.assemblerHybrid == 'flye') { sourmash_metagenome_size(nano_input_ch,database_sourmash) }
      // assembly with assembler choice via --assemblerHybrid; assemblerOutput should be the emiting channel
        if (params.assemblerHybrid == 'flye') { flye(removeSmallReads.out.join(sourmash_metagenome_size.out)) ; assemblerUnpolished = flye.out[0] ; graphOutput = flye.out[1]}
        if (params.assemblerHybrid == 'flye') { medaka(racon(minimap2_to_polish(assemblerUnpolished))) }
        if (params.assemblerHybrid == 'flye') { pilon(medaka.out.join(fastp.out[0])) }
        if (params.assemblerHybrid == 'flye') { assemblerOutput = pilon.out }
        if (params.assemblerHybrid == 'spades') { spades(removeSmallReads.out.join(fastp.out[0])) ; assemblerOutput = spades.out[0] ; graphOutput = spades.out[1]}
      // read mapping for the binner  
        minimap2_to_bam(assemblerOutput.join(nano_input_ch)) 
        bwa_to_bam(assemblerOutput.join(illumina_input_ch)) 
      // extra "not-sample-reads" mapping to improve binning     
        if (params.extra_ont) { minimap2_to_bam_extra(assemblerOutput.combine(extra_ont_ch, by: 0)) }
        if (params.extra_ill) { bwa_to_bam_extra(assemblerOutput.combine(extra_ill_ch, by: 0)) } 
        // merge "not-sample-reads-bam-files" into extra_bam channel
            if (params.extra_ont && params.extra_ill) { extra_bam = bwa_to_bam_extra.out.concat(minimap2_to_bam_extra.out).groupTuple() } 
            if (params.extra_ont && !params.extra_ill) { extra_bam = minimap2_to_bam_extra.out }
            if (!params.extra_ont && params.extra_ill) { extra_bam = bwa_to_bam_extra.out } 
      // Multiple Binning apporaches
        // metabat2
        if (params.extra_ont || params.extra_ill) { metabat2_extra(assemblerOutput.join(minimap2_to_bam.out).join(bwa_to_bam.out).join(extra_bam)) ; metabat2_output = metabat2_extra.out }
        if (!params.extra_ont && !params.extra_ill) { metabat2(assemblerOutput.join(minimap2_to_bam.out).join(bwa_to_bam.out)) ; metabat2_output = metabat2.out }
        // maxbin2
        maxbin2(assemblerOutput.join(nano_input_ch).join(illumina_input_ch))
        // concoct
        if (params.extra_ont || params.extra_ill) { concoct_extra(assemblerOutput.join(minimap2_to_bam.out).join(bwa_to_bam.out).join(extra_bam)) ; concoct_output = concoct_extra.out }
        if (!params.extra_ont && !params.extra_ill)  { concoct(assemblerOutput.join(minimap2_to_bam.out).join(bwa_to_bam.out)) ; concoct_output = concoct.out }
      // Bin refinment of all bins
        metawrap(metabat2_output.join(maxbin2.out).join(concoct_output))
      // Bin re-assembly module
        // get contig ID for each bin and collect all bins into one file (cat_fasta) for mapping
          contig_IDs_in_bins(metawrap.out[0])
          cat_fasta(metawrap.out[0])
        // map reads against all bins 
          bwa_bin(cat_fasta.out.join(illumina_input_ch))    
          minimap2_bin(cat_fasta.out.join(nano_input_ch))
        //retrieve reads per bin and emit unmapped read
          reads_retrieval(contig_IDs_in_bins.out.transpose().combine((bwa_bin.out.join(minimap2_bin.out).join(illumina_input_ch).join(nano_input_ch)), by: 0))
        // Unicyler re-assembly for each bin
          unicycler(reads_retrieval.out)
      // retrive unmapped reads
        reads_retrieval_unmapped(bwa_bin.out.join(minimap2_bin.out).join(illumina_input_ch).join(nano_input_ch))

  emit:   
        clean_bins = unicycler.out[0].groupTuple()

}

/* 
Illumina only Assembly Workflow 
UNTESTED MODULES NEED TO BE ADJUSTED, because the in and output are different in some cases!!!
*/

workflow illumina_assembly_wf {
  get:    illumina_input_ch
          extra_ont_ch
          extra_ill_ch

  main:
      // trimming and QC of reads
        fastp(illumina_input_ch)
      // assembly
        if (params.assemblerShort == 'spades') { spades_ill_only(fastp.out[0]) ; assemblerOutput = spades_ill_only.out[0] ; graphOutput = spades_ill_only.out[1]}
        if (params.assemblerShort == 'megahit') { megahit(fastp.out[0]) ; assemblerOutput = megahit.out[0] ; graphOutput = megahit.out[1]}
      // read mapping for the binner  
        bwa_to_bam(assemblerOutput.join(illumina_input_ch)) 
      // extra "not-sample-reads" mapping to improve binning     
        if (params.extra_ont) { minimap2_to_bam_extra(assemblerOutput.combine(extra_ont_ch, by: 0)) }
        if (params.extra_ill) { bwa_to_bam_extra(assemblerOutput.combine(extra_ill_ch, by: 0)) } 
        // merge "not-sample-reads-bam-files" into extra_bam channel
            if (params.extra_ont && params.extra_ill) { extra_bam = bwa_to_bam_extra.out.concat(minimap2_to_bam_extra.out).groupTuple() } 
            if (params.extra_ont && !params.extra_ill) { extra_bam = minimap2_to_bam_extra.out }
            if (!params.extra_ont && params.extra_ill) { extra_bam = bwa_to_bam_extra.out } 
      // Multiple Binning apporaches
        // metabat2
        if (params.extra_ont || params.extra_ill) { metabat2_extra_ill_only(assemblerOutput.join(bwa_to_bam.out).join(extra_bam)) ; metabat2_output = metabat2_extra_ill_only.out }
        if (!params.extra_ont && !params.extra_ill) { metabat2_ill_only(assemblerOutput.join(bwa_to_bam.out)) ; metabat2_output = metabat2_ill_only.out }
        // maxbin2
        maxbin2_ill_only(assemblerOutput.join(illumina_input_ch))
        // concoct
        if (params.extra_ont || params.extra_ill) { concoct_extra_ill_only(assemblerOutput.join(bwa_to_bam.out).join(extra_bam)) ; concoct_output = concoct_extra_ill_only.out }
        if (!params.extra_ont && !params.extra_ill)  { concoct_ill_only(assemblerOutput.join(bwa_to_bam.out)) ; concoct_output = concoct_ill_only.out }
      // Bin refinment of all bins
        metawrap(metabat2_output.join(maxbin2_ill_only.out).join(concoct_output))
      // retrive unmapped reads - UNTESTED
        //reads_retrieval_unmapped_ill_only(bwa_bin.out.join(illumina_input_ch))

  emit:   
      clean_bins = metawrap.out[2]
}

/* General analysis workflow for a dir with lots of fasta files - bins
Its just separated so one could also start a workflow starting the whole assembly/binning procedure
 */


workflow bin_analysis_wf {
    get:  clean_bins
          database_sourmash

    main:   
        checkm(clean_bins)
        sourmash_tax_classification(clean_bins.transpose(), database_sourmash) 
        sourmash_checkm_parser(checkm.out[0].join(sourmash_tax_classification.out.groupTuple())) 

    emit:   
        overview_bins = sourmash_checkm_parser.out
}


workflow kegg_pathways_no_rna_wf {
    get:  clean_bins
          database_eggnogg

    main:
        eggnog_bin(clean_bins.transpose().map { it, file -> tuple( it, file.baseName, file) }, database_eggnogg)     
   
        //bin_annotated_ch=eggnog_bin.out[0].groupTuple(by:0).view()
        //parser(rna_annot_ch,bin_annotated_ch)  
    emit: clean_bins
                        
}




/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {
      // assembly workflows
      if (params.nano && params.illumina ) { hybrid_assembly_wf(nano_input_ch, illumina_input_ch, extra_ont_ch, extra_ill_ch, download_sourmash()) ; bins_input_ch = hybrid_assembly_wf.out }
      if (!params.nano && params.illumina ) { illumina_assembly_wf(illumina_input_ch, extra_ont_ch, extra_ill_ch) ; bins_input_ch = illumina_assembly_wf.out }

      // analysis workflows
      bin_analysis_wf(bins_input_ch, download_sourmash())
      //kegg_pathways_no_rna_wf(bins_input_ch, download_eggnog())

}


/**************************  
* --help
**************************/
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    Workflow: Template
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run main.nf --nano '*/*.fastq' --illumina '*.R{1,2}.fastq.gz' --extra_ill illumina_reads.csv

    ${c_yellow}Input:${c_reset}
    ${c_green} --nano ${c_reset}            '*.fasta' or '*.fastq.gz'   -> one sample per file
    ${c_green} --illumina ${c_reset}        '*.R{1,2}.fastq.gz'         -> file pairs
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset} 

    ${c_green} --bins ${c_reset}            'sample1/' or 'sample*/'   -> one bin folder per sample, containing .fa files     

    ${c_yellow}Input differential Binning:${c_reset}
    ${c_green} --extra_ont ${c_reset}        'nanopore_reads.csv'  -> read list sparated with ,
    ${c_green} --extra_ill ${c_reset}        'illumina_reads.csv'  -> read list sparated with ,
    ${c_dim} See differential binning explanation below${c_reset}


    ${c_yellow}Options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --output            name of the result folder [default: $params.output]
    --assemblerHybrid   hybrid assembly tool used [spades | flye, default: $params.assemblerHybrid]
    --assemblerShort    illumina assembly tool used [spades | megahit, default: $params.assemblerShort]


    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}LSF computing:${c_reset}
    For execution of the workflow on a HPC with LSF adjust the following parameters:
    --databases         defines the path where databases are stored [default: $params.cloudDatabase]
    --workdir           defines the path where nextflow writes tmp files [default: $params.workdir]
    --cachedir          defines the path where images (singularity) are cached [default: $params.cachedir] 

    Profile:
    -profile                 standard (local, pure docker) [default]
                             conda (mixes conda and docker)
                             lsf (HPC w/ LSF, singularity/docker)
                             nanozoo (googlegenomics and docker)  
                             gcloudAdrian (googlegenomics and docker)
                             gcloudChris (googlegenomics and docker)
                             gcloudMartin (googlegenomics and docker)
                             ebi (EBI cluster specific, singularity and docker)
                             ${c_reset}

    ${c_dim}Differential Binning: 
    a Create a csv file with additional reads to guide binning.
    b These additional Reads are only used for "binning guidance" - nothing more
    c Structure of the csv looks like this for illumina:

      sampleID-1,path-to-read-1,path-to-read-2
      sampleID-2,path-to-read-1,path-to-read-2
      ...

    IMPORTANT: "sampleID" refers to the actual sample that gets assembled, NOT the name of
    the reads in this csv. This number is used to know which of the additional reads should be used
    for each sample. The sampleID refers to the actual samples from --nano / --illumina 
      ${c_reset}
    """.stripIndent()
}


