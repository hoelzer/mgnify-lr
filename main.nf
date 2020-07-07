#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
Nextflow -- MGnify long-read support
Author: hoelzer.martin@gmail.com
*/

/************************** 
* META & HELP MESSAGES 
**************************/

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

if( !nextflow.version.matches('20.01+') ) {
    println "This workflow requires Nextflow version 20.01 or greater -- You are running version $nextflow.version"
    exit 1
}

if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if (params.fasta == '') {
  if (params.nano == '' || (params.nano == '' && params.illumina == '')) { 
    if (params.sra == '') {
      exit 1, "input missing, use [--nano] or [--sra] or [--nano] and [--illumina] or [--fasta]"
    }
  }
}

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

// fasta assembly input & --list support
if (params.fasta && params.list) { fasta_input_ch = Channel
  .fromPath( params.fasta, checkIfExists: true )
  .splitCsv()
  .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)] }
  .view() }
else if (params.fasta) { fasta_input_ch = Channel
  .fromPath( params.fasta, checkIfExists: true)
  .map { file -> tuple(file.simpleName, file) }
  .view() }


// SRA reads input 
if (params.sra) {
nano_input_ch = Channel
  .fromSRA(params.sra, apiKey: params.key)
  .view()
}

// user defined host genome fasta
if (params.host) {
  host_input_ch = Channel
    .fromPath( params.host, checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
    .view()
}

/************************** 
* MODULES
**************************/

// databases
include get_host from './modules/get_host'
include diamond_download_db from './modules/diamondgetdatabase'

// read preprocessing and qc
include removeSmallReads from './modules/removeSmallReads'
include fastp from './modules/fastp' 
include nanoplot from './modules/nanoplot'
   
// estimate genome size
include estimate_gsize from './modules/estimate_gsize'

// assembly & polishing
include flye from './modules/flye'
include spades from './modules/spades'
include minimap2_to_polish from './modules/minimap2'
include racon from './modules/racon'
include medaka from './modules/medaka' 
include pilon from './modules/pilon' 

// decontamination
include bbduk from './modules/bbduk' 
include minimap2_index_ont from './modules/minimap2' 
//include minimap2_index_ill from './modules/minimap2' 
include minimap2_index_assembly from './modules/minimap2' 
include minimap2_clean_ont as clean_ont from './modules/minimap2' 
include minimap2_clean_assembly as clean_assembly from './modules/minimap2' 
//include minimap2_clean_ill as clean_ill from './modules/minimap2' 

// analysis
include prodigal from './modules/prodigal'
include diamond from './modules/diamond'
include ideel from './modules/ideel'

// ENA submission
include ena_manifest from './modules/ena_manifest' 
include ena_manifest_hybrid from './modules/ena_manifest'
include ena_project_xml from './modules/ena_project_xml'
include ena_project_xml_hybrid from './modules/ena_project_xml' 


/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow download_host_genome {
  main:
    // local storage via storeDir
    if (!params.cloudProcess) { get_host(); db = get_host.out }
    // cloud storage via db_preload.exists()
    else {
        db_preload = file("${params.databases}/hosts/${params.species}/${params.species}.fa.gz")
      }
    if (db_preload.exists()) { db = db_preload }
    else  { get_host(); db = get_host.out } 

  emit: db
}

workflow download_diamond {
    main:
        if (!params.cloudProcess) { diamond_download_db() ; database_diamond = diamond_download_db.out}
        else { 
            dia_db_preload = file("${params.databases}/diamond/database_uniprot.dmnd")
            if (dia_db_preload.exists()) { database_diamond = dia_db_preload }    
            else  { diamond_download_db() ; database_diamond = diamond_download_db.out }
        }
    emit: database_diamond
}  


/************************** 
* SUB WORKFLOWS
**************************/

/**********************************************************************/
/* Nanopore Preprocessing Workflow 
/**********************************************************************/
workflow nanopore_preprocess_wf {
  take:  nano_input_ch
         clean_ont_ch

  main:
      // trimming and QC of reads
        nanoplot(nano_input_ch)
        removeSmallReads(nano_input_ch)
        nano_output_ch = removeSmallReads.out

      // decontaminate reads if a host genome is provided
      if (clean_ont_ch) {
        clean_ont(nano_output_ch, clean_ont_ch)
        nano_output_ch = clean_ont.out[0]
      }

  emit:   
      nano_output_ch
}

/**********************************************************************/
/* Illumina Preprocessing Workflow 
/**********************************************************************/
workflow illumina_preprocess_wf {
  take:  illumina_input_ch
         clean_ill_ch

  main:
      // trimming and QC of reads
        fastp(illumina_input_ch)
        illumina_output_ch = fastp.out[0]

        bbduk(illumina_output_ch, clean_ill_ch)
        illumina_output_ch = bbduk.out[0]

  emit:   
        illumina_output_ch
}

/**********************************************************************/
/* Nanopore-only Assembly Workflow 
/**********************************************************************/
workflow nanopore_assembly_wf {
  take:  nano_input_ch

  main:
        if (params.assemblerLong == 'flye') { 
          // size estimation for flye 
          estimate_gsize(nano_input_ch) 
          // assembly raw
          flye(nano_input_ch.join(estimate_gsize.out))
          // polish
          medaka(racon(minimap2_to_polish(flye.out[0])))
          // collect assembly outputs
          assemblyUnpolished = flye.out[0].map {name, reads, assembly -> [name, assembly]}
          assemblyRacon = racon.out.map {name, reads, assembly -> [name, assembly]}
          assemblyMedaka = medaka.out
        }

      //assemblies = assemblyUnpolished.concat(assemblyRacon).concat(assemblyMedaka)

  emit:   
        assemblyUnpolished
        assemblyRacon
        assemblyMedaka
        flye.out[1] // the flye.log
        estimate_gsize.out
}

/**********************************************************************/
/* Clean Assembly Workflow 
/**********************************************************************/
workflow clean_assembly_wf {
  take:  assembly_input_ch
         clean_assembly_ch

  main:
        clean_assembly(assembly_input_ch, clean_assembly_ch)

  emit:   
        clean_assembly.out[0]
}

/**********************************************************************/
/* SR polishing Workflow 
/**********************************************************************/
workflow illumina_polishing_wf {
  take:  assembly_input_ch
         illumina_input_ch
  main:
        pilon(assembly_input_ch, illumina_input_ch)
  emit:   
        pilon.out
}


/**********************************************************************/
/* SPAdes Hybrid Assembly Workflow 
/**********************************************************************/
workflow hybrid_assembly_wf {
  take:  nano_input_ch
         illumina_input_ch

  main:

      if (params.assemblerHybrid == 'spades') {
        spades(nano_input_ch.join(illumina_input_ch))
        assemblyUnpolished= spades.out[0]
        graphOutput = spades.out[1]
      }

  emit:   
        assemblyUnpolished
}


/**********************************************************************/
/* Analysis Workflow 
/**********************************************************************/
workflow analysis_wf {
  take: assembly
        db_diamond

  main:
        ideel(diamond(prodigal(assembly),db_diamond))
}


/**********************************************************************/
/* Indices Workflow 
/**********************************************************************/
workflow index_wf {
  take: host

  main:
    minimap2_index_ont(host)
    minimap2_index_assembly(host)
    //minimap2_index_ill(host)

  emit:
    minimap2_index_ont.out
    minimap2_index_assembly.out
    //minimap2_index_ill.out
}


/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {

      if (params.fasta) {
        assemblies = fasta_input_ch
      }
      else {
        clean_ont_ch = false
        clean_ill_ch = false // uses bbduk per default
        clean_assembly_ch = false
        if (params.illumina && params.nano) {
          clean_assembly_ch = file("${baseDir}/clean/assembly/NC_001422_DCS.mmi", checkIfExists: true)
        } else {
          if (params.illumina) {
            clean_assembly_ch = file("${baseDir}/clean/assembly/NC_001422_FNA.mmi", checkIfExists: true)          
          }
          if (params.nano) {
            clean_assembly_ch = file("${baseDir}/clean/assembly/DCS_FNA.mmi", checkIfExists: true)          
          }
        }

        // 1) check for user defined minimap2 indices 
        if (params.clean_ont) { clean_ont_ch = file(params.clean_ont, checkIfExists: true) }
        if (params.clean_ill) { clean_ill_ch = file(params.clean_ill, checkIfExists: true) }
        if (params.clean_assembly) { 

          clean_assembly_ch = file(params.clean_assembly, checkIfExists: true) 
        }
 
        // 2) build indices if just a fasta is provided
        // WIP
        if (params.host) {
          index_wf(host_input_ch)
          clean_ont_ch = index_wf.out[0]
          clean_assembly_ch = index_wf.out[1]
        }

        // 3) download genome and build indices
        // WIP
        if (params.species) { 
          download_host_genome()
          genome = download_host_genome.out
          host_input_ch = Channel.of( [params.species, genome] )
          index_wf(host_input_ch)
          clean_ont_ch = index_wf.out[0]
          clean_assembly_ch = index_wf.out[1]
        }
      
        // ONT preprocess
        if (params.nano) {
          nanopore_preprocess_wf(nano_input_ch, clean_ont_ch)
        }

        // ILL preprocess
        if (params.illumina ) { 
          illumina_preprocess_wf(illumina_input_ch, clean_ill_ch)
        }

        // Flye-based assembly followed by optional short-read polishing
        if (!params.illumina || params.assemblerHybrid != 'spades') { 
          // assembly w/ flye
          nanopore_assembly_wf(nanopore_preprocess_wf.out)
          assemblyRaw = nanopore_assembly_wf.out[0]
          assemblyRacon = nanopore_assembly_wf.out[1]
          assemblyReady = nanopore_assembly_wf.out[2]

          // combine for analysis step 
          assemblies = assemblyRaw.concat(assemblyRacon).concat(assemblyReady)

          // polish with short reads
          if (params.illumina) {
            illumina_polishing_wf(assemblyReady, illumina_preprocess_wf.out)
            assemblyReady = illumina_polishing_wf.out
            assemblies = assemblies.concat(assemblyReady)
          }
        }

        // Hybrid SPAdes
        if (params.nano && params.illumina && params.assemblerHybrid == 'spades') { 
          // assembly w/ spades
          hybrid_assembly_wf(nanopore_preprocess_wf.out, illumina_preprocess_wf.out)
          assemblyReady = hybrid_assembly_wf.out

          // collect for analysis step 
          assemblies = assemblyReady
        }

        // clean assembly
        if (clean_assembly_ch) { 
          clean_assembly_wf(assemblyReady, clean_assembly_ch)
          assemblyReady = clean_assembly_wf.out
          assemblies = assemblies.concat(assemblyReady)
        }
      }

      // analysis workflow
      if (params.dia_db) { database_diamond = file(params.dia_db) } 
      else { download_diamond(); database_diamond = download_diamond.out }
      analysis_wf(assemblies, database_diamond)

      // ENA submission
      if (!params.illumina || params.assemblerHybrid != 'spades') { 
        if (params.study || params.sample || params.run) {
          ena_manifest(assemblyReady, nanopore_assembly_wf.out[3], nanopore_assembly_wf.out[4])
          ena_project_xml(assemblyReady, nanopore_assembly_wf.out[3], nanopore_assembly_wf.out[4])
        }
      } else {
        if (params.study || params.sample || params.run) {
          ena_manifest_hybrid(assemblyReady)
          ena_project_xml_hybrid(assemblyReady)
        }
      }
}


/**************************  
* --help
*
*     --maxmem            max memory for kmer operations [default: $params.maxmem]
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
    ${c_green} --fasta ${c_reset}           '*.fasta' or '*.fasta.gz'   -> one sample per file
    ${c_green} --sra ${c_reset}             ERR3407986                  -> Run acc, currently only for ONT data supported
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset} 

    ${c_yellow}Options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --memory            max memory for local use [default: $params.cores]
    --gsize            	will be estimated if not provided, genome size for flye assembly [default: $params.gsize]
    --length            cutoff for ONT read length filtering [default: $params.length]
    --assemblerHybrid   hybrid assembly tool used [spades, flye default: $params.assemblerHybrid]
    --assemblerLong     nanopore assembly tool used [flye, default: $params.assemblerLong]
    --output            name of the result folder [default: $params.output]

    ${c_yellow}Custom databases:${c_reset}
     --dia_db             input for diamond database e.g.: 'databases/database_uniprot.dmnd

    ${c_yellow}Decontamination:${c_reset}
    You have three options to provide references for decontamination. 
    Per default phiX/DCS controls will be used to clean your Illumina/ONT data.
    Deactivate cleaning via '--clean_ont false', etc..

    1) Provide prepared minimap2 indices...
    --clean_ont         minimap2 index prepared with the ``-x map-ont`` flag; clean ONT [default: $params.clean_ont]
    --clean_ill         FASTA file for BBDUK Illumina read decontamination; clean ILLUMINA [default: $params.clean_ill]
    --clean_assembly    minimap2 index prepared with the ``-x asm5`` flag; clean FASTA [default: $params.clean_assembly]
                        during runtime either DCS (for ONT), phiX (for Illumina), or DCS+phiX (hybrid) will be automatically 
                        selected based on your input

    2) Or use your own FASTA...
    --host          use your own FASTA sequence for decontamination, e.g., host.fasta.gz. minimap2 indices will be calculated for you. [default: $params.host]

    3) Or let me download a defined host. 
    Per default phiX and the ONT DCS positive controls are added to the index.
    You can remove controls from the index by additionally specifying the parameters below.   
    --species       reference genome for decontamination is selected based on this parameter [default: $params.species]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel]${c_reset}
    --phix       do NOT use phix in decontamination [Illumina: enterobacteria_phage_phix174_sensu_lato_uid14015, NC_001422]
    --dcs        do NOT use DCS in decontamination [ONT DNA-Seq: 3.6 kb standard amplicon mapping the 3' end of the Lambda genome]

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

    ${c_yellow}Computing:${c_reset}
    For execution of the workflow on a HPC with LSF adjust the following parameters:
    --databases             defines the path where databases are stored [default: $params.dbs]
    --workdir               defines the path where nextflow writes tmp files [default: $params.workdir]
    --condaCacheDir         defines the path where environments (conda) are cached [default: $params.condaCacheDir]
    --singularityCacheDir   defines the path where images (singularity) are cached [default: $params.singularityCacheDir] 


    ${c_yellow}Profile:${c_reset}
    You can merge different profiles for different setups, e.g.

        -profile local,docker
        -profile lsf,singularity
        -profile slurm,singularity

    -profile                 standard (local,docker) [default]

                             local
                             lsf
                             slurm

                             docker
                             singularity
                             conda

                             ebi (lsf,singularity; preconfigured for the EBI cluster)
                             yoda (lsf,singularity; preconfigured for the EBI YODA cluster)
                             ara (slurm,conda; preconfigured for the ARA cluster)
                             nih (slurm,singularity; preconfigured for the NIH cluster)
                             gcloud (use this as template for your own GCP setup)
                             ${c_reset}

    """.stripIndent()
}


