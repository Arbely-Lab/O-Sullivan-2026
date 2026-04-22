module load cutadapt/2.10
module load gcc/8.2.0
module load bowtie2/2.4.1

cd /directory/path

mkdir -p /directory/path/[sample]

for f in *.gz; do
  STEM=$(basename "${f}" .gz)
  gunzip -c "${f}" > /directory/path/"${STEM}"
done

cd ./[sample]

cutadapt \
	-m 20 \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o sample.trimmed.R1.fastq -p sample.trimmed.R2.fastq \
    raw_sample_data_1.fq raw_sample_data_2.fq

bowtie2 \
	--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700 --dovetail -p 4 \
	-x /path/to/reference \
	-1 sample.trimmed.R1.fastq \
	-2 sample.trimmed.R2.fastq \
	-S sample.sam

bowtie2 \
	--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700 --dovetail -p 4 \
	-x /path/to/Ecoli/reference \
	-1 sample.trimmed.R1.fastq \
	-2 sample.trimmed.R2.fastq \
	-S sample_ecoli.sam

module load samtools/1.9

samtools stats sample.sam > sample_stats.txt
samtools stats sample_ecoli.sam > sample_ecoli_stats.txt

samtools view -b -h -F 3852 sample.sam > sample.bam
samtools sort -o sample_sorted.bam sample.bam
samtools index sample_sorted.bam

module load gcc/8.2.0
module load bedtools/2.29.0
module load seacr/1.3
module load deeptools/3.3.0

bedtools bamtobed -bedpe -i sample_sorted.bam > sample_seacr.bed

cut -f 1,2,6,8 sample_seacr.clean.bed | sort -k1,1 -k2,2n -k3,3n > sample_seacr.clean.bed
awk -v OFS="\t" '$1=$1' sample_seacr.clean.bed > sample_seacr.final.bed

bedtools genomecov -bg -scale [scale factor] -i sample_seacr.final.bed -g /path/to/ref/genome > sample_seacr_sc.bedgraph

SEACR_1.3.sh sample_seacr_sc.bedgraph 0.00001 non stringent sample_sc_str_0.00001

module purge
module load macs/2.2.7.1
module load gcc/8.2.0
module load samtools/1.9
samtools view -h -b -s [scale_factor] sample_sorted.bam > sample_macs.bam

macs2 callpeak -t sample_sorted.bam -n sample_name -f BAMPE -g 3.1e9 -q 0.00001

module load deeptools/3.3.0

bamCoverage -b sample_sorted.bam -o sample.bw --scaleFactor [scaling factor] -p max/2





