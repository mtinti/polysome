#!/bin/bash 
#$ -l local_free=200G
#$ -l h_vmem=24G
#$ -cwd
#$ -V
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

mkdir -p $TMPDIR'/genomes/'$genome'/' && cp -avr  'genomes/'$genome'/' $TMPDIR'/genomes/' 
mkdir -p $TMPDIR'/genomes/'$genome'gmap/' && cp -avr  'genomes/'$genome'gmap/' $TMPDIR'/genomes/' 


#mkdir -p $TMPDIR'/'$experiment'/data/' && cp -avr  $experiment'/data/' $TMPDIR'/'$experiment'/'
#mkdir -p $TMPDIR'/FastQ_Screen_Genomes/' && cp -avr  'FastQ_Screen_Genomes/' $TMPDIR
path_genome_index=$TMPDIR'/genomes/'$genome'/'$genome
path_trascriptome_index=$TMPDIR'/trascriptome/'$genome'/'$genome
printf 'Elapsed time: %s\n' $(timer $tmr) 

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
mkdir $TMPDIR'/trinity_out_'$experiment

cp $experiment'/res2/'$experiment/$experiment'_sorted.bam' $TMPDIR
cp $experiment'/res2/'$experiment/$experiment'_sorted.bam.bai' $TMPDIR

Trinity --genome_guided_bam $TMPDIR'/'$experiment'_sorted.bam' --full_cleanup \
--genome_guided_max_intron 10 --genome_guided_min_coverage 10 \
--max_memory 23G --CPU 8 --jaccard_clip --output $TMPDIR'/trinity_out_'$experiment

gmap -n 0 -D $TMPDIR'/genomes/'$genome'gmap/' \
-d $genome $TMPDIR'/trinity_out_'$experiment'/Trinity-GG.fasta' -f gff3_gene > $TMPDIR'/'$experiment'_trinity_gmap.gff'

gffread -E  $TMPDIR'/'$experiment'_trinity_gmap.gff' -T -o $TMPDIR'/'$experiment'_trinity_gmap.gtf'

mkdir $experiment'/trinity2'
cp $TMPDIR'/trinity_out_'$experiment'/Trinity-GG.fasta' $experiment'/trinity2'
cp $TMPDIR'/'$experiment'_trinity_gmap.gtf' $experiment'/trinity2'
gzip $experiment'/trinity2/Trinity-GG.fasta'











