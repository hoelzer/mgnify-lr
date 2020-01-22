#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
Nextflow -- MGnify long-read support
Author: hoelzer.martin@gmail.com
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
if (params.nano == '' || (params.nano == '' && params.illumina == '')) 
    { exit 1, "input missing, use [--nano] or [--nano] and [--illumina]"}

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

/************************** 
* MODULES
**************************/

    include sourmash_download_db from './modules/sourmashgetdatabase' params(cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
    include get_host from './modules/get_host' params(species: params.species, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)

    include removeSmallReads from './modules/removeSmallReads' params(output: params.output)
    include nanoplot from './modules/nanoplot' params(output: params.output)
    include sourmash_metagenome_size from './modules/sourmash_metagenome_size' params(output: params.output, gsize: params.gsize)
    include flye from './modules/flye' params(output: params.output)
    include minimap2_to_polish from './modules/minimap2'
    include minimap2_to_decontaminate from './modules/minimap2' params(output:params.output)
    include racon from './modules/racon'
    include medaka from './modules/medaka' params(output: params.output, model: params.model)
    include ena_manifest from './modules/ena_manifest' params(output: params.output, model: params.model, assemblerLong: params.assemblerLong, study: params.study, sample: params.sample, run: params.run)
    include ena_project_xml from './modules/ena_project_xml' params(output: params.output, model: params.model, assemblerLong: params.assemblerLong, study: params.study, sample: params.sample, run: params.run)

/*
    include bwa_to_bam as bwa_bin from './modules/bwa'  
    include bwa_to_bam as bwa_to_bam_extra from './modules/bwa'
    include bwa_to_bam from './modules/bwa'
    include cat_fasta from './modules/cat_fasta'
    include fastp from './modules/fastp' params(output: params.output)
    include megahit from './modules/megahit' params( output : params.output)
    include minimap2_to_bam as minimap2_bin from './modules/minimap2'
    include minimap2_to_bam as minimap2_to_bam_extra from './modules/minimap2'
    include minimap2_to_bam from './modules/minimap2'
    include pilon from './modules/pilon' params(output: params.output)
    include sourmash_checkm_parser from './modules/parser/checkm_sourmash_parser'params(output: params.output)
    include sourmash_tax_classification from './modules/sourmash_tax_classification' params(output : params.output)
    include spades from './modules/spades' params( output : params.output)
    include spades_ill_only from './modules/spades' params( output : params.output)
    include unicycler from './modules/unicycler' params(output : params.output)
*/


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

workflow download_host_genome {
  main:
    // local storage via storeDir
    if (!params.cloudProcess) { get_host(); db = get_host.out }
    // cloud storage via db_preload.exists()
    if (params.cloudProcess) {
      if (params.phix) {
        db_preload = file("${params.cloudDatabase}/hosts/${params.species}_phix/${params.species}_phix.fa.gz")
      } else {
        db_preload = file("${params.cloudDatabase}/hosts/${params.species}/${params.species}.fa.gz")
      }
      if (db_preload.exists()) { db = db_preload }
      else  { get_host(); db = get_host.out } 
    }
  emit: db
}


/************************** 
* SUB WORKFLOWS
**************************/

/* Hybrid Assembly Workflow */

workflow hybrid_assembly_wf {
  get:  nano_input_ch
        illumina_input_ch
        database_sourmash

  main:
      // trimming and QC of reads
        removeSmallReads(nano_input_ch)
        fastp(illumina_input_ch)
        nanoplot(nano_input_ch)

      // size estimation for flye // not working well - bc not installed n sourmash nanozoo container
        if (params.assemblerHybrid == 'flye') { sourmash_metagenome_size(nano_input_ch, database_sourmash) }

      // assembly with assembler choice via --assemblerHybrid; assemblerOutput should be the emiting channel
        if (params.assemblerHybrid == 'flye') { flye(removeSmallReads.out.join(sourmash_metagenome_size.out)) ; assemblerUnpolished = flye.out[0] ; graphOutput = flye.out[1]}
        if (params.assemblerHybrid == 'flye') { medaka(racon(minimap2_to_polish(assemblerUnpolished))) }
        if (params.assemblerHybrid == 'flye') { pilon(medaka.out.join(fastp.out[0])) }
        if (params.assemblerHybrid == 'flye') { assemblerOutput = pilon.out }
        if (params.assemblerHybrid == 'spades') { spades(removeSmallReads.out.join(fastp.out[0])) ; assemblerOutput = spades.out[0] ; graphOutput = spades.out[1]}

  emit:   
        assembly = assemblerOutput.out
}


/* Nanopore-only Assembly Workflow */

workflow nanopore_assembly_wf {
  get:  nano_input_ch
        database_sourmash
        host_genome

  main:
      // decontaminate reads if a host genome is provided
      if (host_genome) {
        minimap2_to_decontaminate(nano_input_ch, host_genome)
        nano_input_ch = minimap2_to_decontaminate.out
      }

      // trimming and QC of reads
        removeSmallReads(nano_input_ch)
        nanoplot(nano_input_ch)

      // size estimation for flye // not working well - bc not installed n sourmash nanozoo container
        if (params.assemblerLong == 'flye') { sourmash_metagenome_size(nano_input_ch, database_sourmash) }

      // assembly with assembler choice via --assemblerLong; assemblerOutput should be the emiting channel
        if (params.assemblerLong == 'flye') { flye(removeSmallReads.out.join(sourmash_metagenome_size.out)) ; assemblerUnpolished = flye.out[0]}
        if (params.assemblerLong == 'flye') { medaka(racon(minimap2_to_polish(assemblerUnpolished))) }
        if (params.assemblerLong == 'flye') { assemblerOutput = medaka.out }

  emit:   
        assemblerOutput
        flye.out[1] // the flye.log
        sourmash_metagenome_size.out
}


/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {

      // nanopore read decontamination
      genome = false
      if (params.species) { 
        download_host_genome()
        genome = download_host_genome.out
      }
      
      // assembly workflows
      if (params.nano && !params.illumina ) { 
        nanopore_assembly_wf(nano_input_ch, download_sourmash(), genome)
        if (params.study || params.sample || params.run) {
          ena_manifest(nanopore_assembly_wf.out[0], nanopore_assembly_wf.out[1], nanopore_assembly_wf.out[2])
          ena_project_xml(nanopore_assembly_wf.out[0], nanopore_assembly_wf.out[1], nanopore_assembly_wf.out[2])
        }
      }
      if (params.nano && params.illumina ) { 
        hybrid_assembly_wf(nano_input_ch, illumina_input_ch, extra_ont_ch, extra_ill_ch, download_sourmash())
      }

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
    
    MGnify-LR: build long-read and hybrid assemblies for input of the MGnify pipeline. 
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run main.nf --nano '*/*.fastq' --illumina '*.R{1,2}.fastq.gz'

    ${c_yellow}Input:${c_reset}
    ${c_green} --nano ${c_reset}            '*.fasta' or '*.fastq.gz'   -> one sample per file
    ${c_green} --illumina ${c_reset}        '*.R{1,2}.fastq.gz'         -> file pairs
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset} 

    ${c_yellow}Options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --output            name of the result folder [default: $params.output]
    --gsize            	estimated genome size for flye assembly [default: $params.gsize]
    --assemblerHybrid   hybrid assembly tool used [spades | flye, default: $params.assemblerHybrid]
    --assemblerLong     nanopore assembly tool used [flye, default: $params.assemblerLong]

    ${c_yellow}Decontamination:${c_reset}
    --species       reference genome for decontamination is selected based on this parameter [default: $params.species]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel]${c_reset}


    ${c_yellow}ENA parameters:${c_reset}
    --study             ENA study ID [default: $params.study]
    --sample            ENA sample ID [default: $params.sample]
    --run               ENA run ID [default: $params.run]
    ${c_dim}Use this when you assemble an ENA run to automatically produce the correct
    manifest.txt and project.xml files for ENA upload. Provide either of the three 
    and the files will be generated, missing values you not provided.${c_reset}
    
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
                             ebi (EBI cluster specific, singularity and docker)
                             gcloudMartin (googlegenomics and docker)
                             ${c_reset}

    """.stripIndent()
}


