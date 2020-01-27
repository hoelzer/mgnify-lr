process ena_project_xml {
    label 'basics'
    publishDir "${params.output}/${name}/", mode: 'copy', pattern: "*.xml"
    
    input:
    tuple val(name), file(assembly)
    tuple val(name), file(flye_log)
    tuple val(name), file(genome_size)
    
    output:
    file("*.xml")
    
    shell:
    date = Date().format( 'yyyy-MM-dd' )
    """
    MD5=\$(md5sum ${assembly} | awk '{print \$1}')
    SIZE=\$(cat !{genome_size})
    COVERAGE=\$(grep 'Mean coverage' !{flye_log} | awk '{print \$3}')

    FLYE_VERSION=2.5
    RACON_VERSION=1.4.10
    MEDAKA_VERSION=0.10.0

    STUDY=${params.study}
    SAMPLE=${params.sample}
    RUN=${params.run}

    touch project.xml
    cat <<EOF >> project.xml
<PROJECT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
    <PROJECT alias="\${STUDY}_${workflow.scriptId}">
        <TITLE>Metagenomics assembly of study: \${STUDY}, sample: \${SAMPLE}, run: \${RUN}.</TITLE>
        <DESCRIPTION>This assembly was derived from the Oxford Nanopore Technologies MinION data set \${RUN} and was assembled with ${params.assemblerLong} (v\${FLYE_VERSION}) using the --meta and --plasmids options and an estimated genome size of \${SIZE}. Before assembly, reads smaller 500 nt were filtered and the draft flye assembly was polished by one round of racon (v\$RACON_VERSION) and medaka (v\$MEDAKA_VERSION) using the ${params.model} model. </DESCRIPTION>
        <SUBMISSION_PROJECT>
            <SEQUENCING_PROJECT/>
         </SUBMISSION_PROJECT>
         <PROJECT_LINKS>
            <PROJECT_LINK>
                <XREF_LINK>
                    <DB></DB>
                    <ID></ID>
                </XREF_LINK>
            </PROJECT_LINK>
        </PROJECT_LINKS>
        <PROJECT_ATTRIBUTES>
            <PROJECT_ATTRIBUTE>
                <TAG>new_study_type</TAG>
                <VALUE>Metagenomic assembly</VALUE>
            </PROJECT_ATTRIBUTE>
        </PROJECT_ATTRIBUTES>
    </PROJECT>
</PROJECT_SET>
EOF

    touch submission.xml
    cat <<EOF >> submission.xml
<SUBMISSION>
   <ACTIONS>
      <ACTION>
         <ADD/>
      </ACTION>
      <ACTION>
          <HOLD HoldUntilDate="${date}"/>
       </ACTION>
   </ACTIONS>
</SUBMISSION>
EOF
    """
}

