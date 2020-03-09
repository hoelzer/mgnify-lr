#!/bin/bash

MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
export NUM_THREADS=4
export MEM=1G
export FIX=1
export BATCH_SIZE=0
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

log () {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}


function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}

#parsing arguments
if [[ $# -eq 0 ]];then
echo "Usage:  polca.sh -a <assembly contigs or scaffolds> -r <'Illumina_reads_fastq1 Illumina_reads_fastq'> -t <number of threads> [-n] <optional:do not fix errors that are found> [-m] <optional: memory per thread to use in samtools sort>"
exit 1
fi

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -b|--batch)
            export BATCH_SIZE="$2"
            shift
            ;;
        -t|--threads)
            export NUM_THREADS="$2"
            shift
            ;;
        -n|--nofix)
            export FIX=0
            ;;
        -a|--assembly)
            export ASM="$2"
            shift
            ;;
        -r|--reads)
            READS="$2";
            shift
            ;;
        -m|--memory)
            export MEM="$2";
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "Usage:  polca.sh -a <assembly contigs or scaffolds> -r <'Illumina_reads_fastq1 Illumina_reads_fastq'> -t <number of threads> [-n] <optional:do not fix errors that are found> [-m] <optional: memory per thread to use in samtools sort>"
            echo "Must have bwa, samtools and freebayes available on the PATH"
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

if [[ ! -e $ASM ]];then
echo "Input file $ASM not found or not specified!"
echo "Usage:  polca.sh -a <assembly contigs or scaffolds> -r <'Illumina_reads_fastq1 Illumina_reads_fastq'> -t <number of threads> [-n] <optional:do not fix errors that are found> [-m] <optional: memory per thread to use in samtools sort>"
exit 1
fi

which bwa || error_exit "bwa not found on the PATH, please install bwa aligner"
which freebayes || error_exit "freebayes not found in MaSuRCA bin, please check your MaSuRCA install"
which samtools || error_exit "samtools not found in MaSuRCA bin, please check your MaSuRCA install"

export BASM=`basename $ASM`
export BWA=`which bwa` 
export FREEBAYES=`which freebayes`
export SAMTOOLS=`which samtools`
rm -f bwa.err samtools.err

if [ ! -e $BASM.index.success ];then 
log "Creating BWA index for $ASM"
rm -f $BASM.map.success
$BWA index $ASM -p $BASM.bwa 2>>bwa.err && touch $BASM.index.success || error_exit "Creating BWA index for $ASM failed"
fi

if [ ! -e $BASM.map.success ];then
log "Aligning reads to $ASM"
rm -f $BASM.sort.success
zcat -f $(echo $READS) | $BWA mem -SP -t $NUM_THREADS $BASM.bwa /dev/stdin 1>$BASM.unSorted.sam 2>>bwa.err && \
touch $BASM.map.success
if [ ! -e $BASM.map.success ];then
  error_exit "Aligning reads to $ASM failed"
fi
fi

if [ ! -e $BASM.sort.success ];then
log "Sorting and indexing alignment file"
rm -f $BASM.vc.success
$SAMTOOLS view -bS $BASM.unSorted.sam -o $BASM.unSorted.bam 2>>samtools.err && $SAMTOOLS sort -m $MEM -@ $NUM_THREADS $BASM.unSorted.bam $BASM.alignSorted 2>>samtools.err && \
$SAMTOOLS index $BASM.alignSorted.bam 2>>samtools.err && \
$SAMTOOLS faidx $ASM  2>>samtools.err && \
touch  $BASM.sort.success
if [ ! -e $BASM.sort.success ];then
  error_exit "Sorting and indexing alignment file failed"
fi
fi

#here we are doing variant calling in parallel, per input contig/scaffold
if [ ! -e $BASM.vc.success ];then
rm -f  $BASM.report.success $BASM.fix.success
log "Calling variants in $BASM"
#I do this to mix scaffolds up to equalize batch sizes
ufasta sizes -H $ASM |sort -S 10% -k2 |awk '{print $1}' > $BASM.names
mkdir -p $BASM.work
#begin subshell execution
(
  cd $BASM.work
  rm -f $BASM.vc.success $BASM.fix.success
  CONTIGS=`wc -l ../$BASM.names | awk '{print $1}'`;
  if [ $BATCH_SIZE -lt 1 ];then
    BATCH_SIZE=$(($CONTIGS / $NUM_THREADS+1));
  fi
  if [ $BATCH_SIZE -gt 1000 ];then
    BATCH_SIZE=1000
  fi
  BATCH=1;
  echo "Processing $BATCH_SIZE scaffold(s) per batch"
  LIST="";
  INDEX=1;
  for f in $(cat ../$BASM.names);do
    LIST="$LIST $f"
    let INDEX=$INDEX+1;
      if [ $INDEX -gt $BATCH_SIZE ];then
          echo $LIST | tr " " "\n" > $BATCH.names
          echo $LIST > $BATCH.listnames
          LIST=""
          INDEX=1
          let BATCH=$BATCH+1
      fi
  done
  if [ $INDEX -gt 1 ];then
      echo $LIST | tr " " "\n" > $BATCH.names
      echo $LIST > $BATCH.listnames
  else
    let BATCH=$BATCH-1
  fi

  echo "#!/bin/bash" > commands.sh
  echo "if [ ! -e \$1.vc.success ];then" >> commands.sh
  echo "  $FREEBAYES -C 3 -R 0 -p 1 -F 0.2 -E 0 -b <($SAMTOOLS view -h ../$BASM.alignSorted.bam \`head -n 1 \$1.listnames\` 2>>\$1.samtools.err |$SAMTOOLS view -S -b /dev/stdin 2>>\$1.samtools.err)  -v \$1.vcf -f ../$ASM && touch \$1.vc.success" >> commands.sh 
  echo 'fi' >> commands.sh
  echo "if [ $FIX -gt 0 ];then" >> commands.sh
  echo "  if [ ! -e \$1.fix.success ] && [ -e \$1.vc.success ];then" >> commands.sh
  echo "    $MYPATH/fix_consensus_from_vcf.pl <($MYPATH/ufasta extract -f \$1.names ../$ASM) < \$1.vcf > \$1.fixed.tmp && mv \$1.fixed.tmp \$1.fixed && touch \$1.fix.success"  >> commands.sh
  echo '  fi' >> commands.sh
  echo 'fi' >> commands.sh
  chmod 0755 commands.sh && \

  seq 1 $BATCH |xargs -P $NUM_THREADS -I % ./commands.sh % 

#checking if jobs finished properly

  for f in $(seq 1 $BATCH);do
    if [ ! -e $f.vc.success ];then
      error_exit "Variant calling failed on batch $f in $BASM.work"
    fi
  done
  touch $BASM.vc.success
  if [ $FIX -gt 0 ];then
  for f in $(seq 1 $BATCH);do
    if [ ! -e $f.fix.success ];then
      error_exit "Fixing errors failed on batch $f in $BASM.work"
    fi
  done
  touch $BASM.fix.success
  fi
);
#end subshell execution
if [ -e ./$BASM.work/$BASM.vc.success ];then
  cat ./$BASM.work/*.vcf > $BASM.vcf.tmp && mv $BASM.vcf.tmp $BASM.vcf && touch $BASM.vc.success
else
  error_exit "Variant calling failed in ./$BASM.work"
fi
if [ $FIX -gt 0 ];then
if [ -e ./$BASM.work/$BASM.fix.success ];then
  cat ./$BASM.work/*.fixed  | ufasta format > $BASM.masurca.tmp && mv $BASM.masurca.tmp $BASM.PolcaCorrected.fa && touch $BASM.fix.success
else
  error_exit "Fixing consensus failed in ./$BASM.work"
fi
fi
rm -rf $BASM.work;
fi

if [ ! -e $BASM.report.success ];then
log "Creating report file"
NUMSUB=`grep --text -v '^#'  $BASM.vcf  |perl -ane '{if(length($F[3])==1 && length($F[4])==1){ print "$F[9]:1\n";}}' | awk -F ':' 'BEGIN{nerr=0}{if($6>=3 && $4<=1) nerr+=$NF}END{print nerr}'` 
NUMIND=`grep --text -v '^#' $BASM.vcf  |perl -ane '{if(length($F[3])>1 || length($F[4])>1){$nerr=abs(length($F[3])-length($F[4]));print "$F[9]:$nerr\n";}}' | awk -F ':' 'BEGIN{nerr=0}{if($6>=3 && $4<=1) nerr+=$NF}END{print nerr}'` 
ASMSIZE=`ufasta n50 -S $ASM | awk '{print $2}'` 
NUMERR=$(($NUMSUB+$NUMIND)) 
QUAL=`echo $NUMERR $ASMSIZE | awk '{print 100-$1/$2*100}'`
echo "Substitution Errors: $NUMSUB" > $BASM.report 
echo "Insertion/Deletion Errors: $NUMIND" >> $BASM.report 
echo "Assembly Size: $ASMSIZE" >> $BASM.report 
echo "Consensus Quality: $QUAL" >> $BASM.report 
touch $BASM.report.success
fi
cat $BASM.report
if [ $FIX -gt 0 ];then
log "Success! Final report is in $BASM.report; polished assembly is in $BASM.PolcaCorrected.fa"
else
log "Success! Final report is in $BASM.report"
fi
