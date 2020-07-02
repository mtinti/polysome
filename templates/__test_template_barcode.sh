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

mkdir -p $TMPDIR'/genomes/'$genome'/' && cp -avr  'genomes/'$genome'/' $TMPDIR'/genomes/' 
mkdir -p $TMPDIR'/'$experiment'/data/' && cp -avr  $experiment'/data/' $TMPDIR'/'$experiment'/'
#mkdir -p $TMPDIR'/FastQ_Screen_Genomes/' && cp -avr  'FastQ_Screen_Genomes/' $TMPDIR
path_genome_index=$TMPDIR'/genomes/'$genome'/'$genome
path_trascriptome_index=$TMPDIR'/trascriptome/'$genome'/'$genome
printf 'Elapsed time: %s\n' $(timer $tmr) 

echo 'folder $TMPDIR'
ls -l $TMPDIR

path_out=$TMPDIR'/'$experiment'/'
mkdir $path_out'fastqc'
mkdir $path_out'fastp'
#mkdir $path_out'fastqs'
mkdir $path_out'qual_out'
mkdir $path_out'qual_out/'$experiment
ls $path_out -l
##qc of fastq 
echo 'run 0.1' 


fastq_1=$path_out'data/'$base_fastq'1.fastq.gz'
fastq_2=$path_out'data/'$base_fastq'2.fastq.gz'


echo 'find size of fastq file'
inFileLength=$(echo $(zcat $fastq_1 | wc -l)/4|bc)
echo 'inFileLength all: '$inFileLength


#echo 'run fastqc' 
#fastqc  $fastq_1 -o $path_out'fastqc'
#fastqc  $fastq_2 -o $path_out'fastqc'

#remove problematic reads
#this step do not remove polyA unless specified with the --polyX option
echo 'run fastp' 
out_1=$path_out'data/'$base_fastq'f1.fastq.gz'
out_2=$path_out'data/'$base_fastq'f2.fastq.gz'
fastp -i $fastq_1 -I $fastq_2 -o $out_1 -O $out_2 \
-h $path_out'fastp/fastp_'$base_fastq'.html' -j $path_out'fastp/'$base_fastq'_fastp.json'

echo 'run bowtie2'
(bowtie2 --very-sensitive-local -p 8 -x $path_genome_index -1 $out_1 -2 $out_2) \
2>$path_out'qual_out/'$experiment'/'$experiment'.log' | samtools view -bSu | \
samtools sort -@ 8 -o $path_out$base_fastq'sorted.bam'

echo 'run index'
samtools index $path_out$base_fastq'sorted.bam'

#rm $out_1
#rm $out_2


#unset DISPLAY
#echo 'run qualimap bamqc' 
#qualimap bamqc --java-mem-size=4G -bam $path_out$base_fastq'sorted.bam' \
#-outdir $path_out'qual_out/'$experiment \
#-outfile $experiment'.bamqc.html' -outformat 'HTML'

#echo 'run qualimap rnaseq' 
#qualimap rnaseq --java-mem-size=4G -bam $path_out$base_fastq'sorted.bam'  \
#-gtf $path_genome_index'.gtf' \
#-outdir $path_out'qual_out/'$experiment \
#-pe -outfile $experiment'.rnaseq.html' -outformat 'HTML'


## sorting because picard MarkDuplicates
## needs name sorted to filter out secondary aligment
##maybe not
#echo 'run samsort by name' 
#samtools sort -o $path_out$base_fastq'sorted.bam' -n -@ 8 $path_out$base_fastq'sorted.bam'
#rm $path_out$base_fastq'sorted.bam'
#rm $path_out$base_fastq'sorted.bam.bai'

#bedtools genomecov -ibam $path_out$base_fastq'sorted.bam' \
#-bg -pc > $path_out$base_fastq'sorted_pc_bg.bed'
#gzip $path_out$base_fastq'sorted_pc_bg.bed'

#https://academic.oup.com/nar/article/38/15/4946/2409341
#splice leader sequence last 14
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4829303/
#polyA 250 bases
echo 'run filter'
samtools view  -f 2 -F 768 -b \
-o $path_out$base_fastq'sorted_.bam' $path_out$base_fastq'sorted.bam'

echo 'run rename'
mv $path_out$base_fastq'sorted_.bam' $path_out$base_fastq'sorted.bam'

echo 'run extract_barcode' 
python mylib/extract_barcodes_def2.py \
$path_out$base_fastq'sorted.bam' TCTGTACTATATTG CAATATAGTACAGA 1 True

mv $path_out$base_fastq'sorted_F_plus_R.bam' $path_out$base_fastq'sorted_F_plus_R_SL.bam'

bedtools genomecov -ibam $path_out$base_fastq'sorted_F_plus_R_SL.bam' \
-bg > $path_out$base_fastq'sorted_F_plus_R_SL.bed'
gzip $path_out$base_fastq'sorted_F_plus_R_SL.bed'

cp $path_out$base_fastq'sorted_F_plus_R_SL.bam' $experiment'/res2/'$experiment
cp $path_out$base_fastq'sorted_F_plus_R_SL.bed.gz' $experiment'/res2/'$experiment

echo 'run extract_polyA' 
python mylib/extract_barcodes_def2.py \
    $path_out$base_fastq'sorted.bam' AAAAAAAAAA TTTTTTTTTT 1 True

mv $path_out$base_fastq'sorted_F_plus_R.bam' $path_out$base_fastq'sorted_F_plus_R_PoliA.bam'

bedtools genomecov -ibam $path_out$base_fastq'sorted_F_plus_R_PoliA.bam' \
-bg > $path_out$base_fastq'sorted_F_plus_R_PoliA.bed'
gzip $path_out$base_fastq'sorted_F_plus_R_PoliA.bed'

cp $path_out$base_fastq'sorted_F_plus_R_PoliA.bam' $experiment'/res2/'$experiment
cp $path_out$base_fastq'sorted_F_plus_R_PoliA.bed.gz' $experiment'/res2/'$experiment

