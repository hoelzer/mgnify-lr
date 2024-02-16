process ena_manifest {
    label 'basics'
    publishDir "${params.output}/${name}/", mode: 'copy', pattern: "manifest.txt"
    
    input:
    tuple val(name), file(assembly)
    tuple val(name), file(flye_log)
    //tuple val(name), file(genome_size)
    
    output:
    file("manifest.txt")
    
    shell:
    """
    MD5=\$(md5sum ${assembly} | awk '{print \$1}')
    #SIZE=\$(cat !{genome_size})
    SIZE=\$(grep -v ">" ${assembly} | tr -d '\\n' | wc -m)
    COVERAGE=\$(grep 'Mean coverage' !{flye_log} | awk '{print \$3}')

    FLYE_VERSION=2.5
    RACON_VERSION=1.4.10
    MEDAKA_VERSION=0.10.0

    STUDY=${params.study}
    SAMPLE=${params.sample}
    RUN=${params.run}

    touch manifest.txt
    cat <<EOF >> manifest.txt
STUDY   \${STUDY}_${workflow.scriptId}
SAMPLE  \${SAMPLE}
RUN_REF \${RUN} 
ASSEMBLYNAME    \${RUN}_\${MD5}
ASSEMBLY_TYPE   primary metagenome 
COVERAGE        \${COVERAGE}
PROGRAM         ${params.assemblerLong} v\${FLYE_VERSION} 
PLATFORM        Oxford Nanopore Technologies MinION 
FASTA           ${params.output}/${name}/assembly/${assembly}
DESCRIPTION     Reads < 500 nt removed prior assembly. Draft flye assembly polished with racon v\${RACON_VERSION} and medaka v\${MEDAKA_VERSION} (model: ${params.model}).
EOF
    """
}


process ena_manifest_hybrid {
    label 'basics'
    publishDir "${params.output}/${name}/", mode: 'copy', pattern: "manifest.txt"
    
    input:
    tuple val(name), file(assembly)
    
    output:
    file("manifest.txt")
    
    shell:
    """
    MD5=\$(md5sum ${assembly} | awk '{print \$1}')
    COVERAGE=XXX

    SPADES_VERSION=3.13.1

    STUDY=${params.study}
    SAMPLE=${params.sample}
    RUN=${params.run}

    touch manifest.txt
    cat <<EOF >> manifest.txt
STUDY   \${STUDY}_${workflow.scriptId}
SAMPLE  \${SAMPLE}
RUN_REF \${RUN} 
ASSEMBLYNAME    \$( echo \${RUN} | sed 's/,/_/g')_\${MD5}
ASSEMBLY_TYPE   primary metagenome 
COVERAGE        \${COVERAGE}
PROGRAM         ${params.assemblerHybrid} v\${SPADES_VERSION} 
PLATFORM        Oxford Nanopore Technologies MinION, Illumina MiSeq
FASTA           ${params.output}/${name}/assembly/${assembly}
DESCRIPTION     Reads were quality controlled with fastp prior assembly. Assembly done with SPAdes v\${SPADES_VERSION}.
EOF
    """
}


/*
TUDY: Study accession or unique name (alias)
SAMPLE: Environmental sample accession or unique name (alias)
ASSEMBLYNAME: Unique assembly name
ASSEMBLY_TYPE: ‘primary metagenome’
COVERAGE: The estimated depth of sequencing coverage
PROGRAM: The assembly program
PLATFORM: The sequencing platform, or comma-separated list of platforms
MINGAPLENGTH: Minimum length of consecutive Ns to be considered a gap (optional)
MOLECULETYPE: ‘genomic DNA’, ‘genomic RNA’ or ‘viral cRNA’ (optional)
DESCRIPTION: Free text description of the genome assembly (optional)
RUN_REF: Comma separated list of run accession(s) (optional)
*/
