#!/bin/bash 
#$ -l local_free=50G
#$ -l h_vmem=4G
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 8
#$ -N {experiment}
#$ -o /cluster/majf_lab/mtinti/polisome/{experiment}
#$ -e /cluster/majf_lab/mtinti/polisome/{experiment}


set -e
set -u 
set -o pipefail

experiment='{experiment}'
genome='{g_version}'
base_fastq='{base_fastq}'


function timer()
{
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local  stime=$1
        etime=$(date '+%s')

        if [[ -z "$stime" ]]; then stime=$etime; fi

        dt=$((etime - stime))
        ds=$((dt % 60))
        dm=$(((dt / 60) % 60))
        dh=$((dt / 3600))
        printf '%d:%02d:%02d' $dh $dm $ds
    fi
}


echo 'copy data in $TMPDIR'
##variables

tmr=$(timer)

#mkdir -p $TMPDIR'/genomes/'$genome'/' && cp -avr  'genomes/'$genome'/' $TMPDIR'/genomes/' 


#mkdir -p $TMPDIR'/'$experiment'/data/' && cp -avr  $experiment'/data/' $TMPDIR'/'$experiment'/'
#mkdir -p $TMPDIR'/FastQ_Screen_Genomes/' && cp -avr  'FastQ_Screen_Genomes/' $TMPDIR
#path_genome_index=$TMPDIR'/genomes/'$genome'/'$genome
#path_trascriptome_index=$TMPDIR'/trascriptome/'$genome'/'$genome
#printf 'Elapsed time: %s\n' $(timer $tmr) 

echo 'folder $TMPDIR'
ls -l $TMPDIR

#path_out=$TMPDIR'/'$experiment'/'
#mkdir $path_out'fastqc'
#mkdir $path_out'fastp'
#mkdir $path_out'fastqs'
#mkdir $path_out'qual_out'
#mkdir $path_out'qual_out/'$experiment
#ls $path_out -l
##qc of fastq 
echo 'run 0.1' 
#mkdir -p $TMPDIR'/trinity_out_'$experiment

cp /cluster/majf_lab/mtinti/polisome/$experiment'/res2/'$experiment/$experiment'_sorted.bam' $TMPDIR
cp /cluster/majf_lab/mtinti/polisome/$experiment'/res2/'$experiment/$experiment'_sorted.bam.bai' $TMPDIR


scallop --library_type unstranded --min_transcript_coverage 10 \
--min_splice_bundary_hits 1200 -i $TMPDIR'/'$experiment'_sorted.bam' \
-o $TMPDIR/all_scallop2.gtf --min_bundle_gap 5

mkdir -p $experiment'/scallop2'
cp $TMPDIR/all_scallop2.gtf /cluster/majf_lab/mtinti/polisome/$experiment'/scallop2'

