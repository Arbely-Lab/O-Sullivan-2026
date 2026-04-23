#Basecall with guppy or Dorado
#For Dorado, add option to output fastq rather than automatically map with minimap
module load guppy/6.1.2_gpu

#Guppy for R9.x flow cells
guppy_basecaller \
-i /path/to/fast5_pass \
-s /path/to/output \
-c /path/to/config/res_dna_r941_min_modbases_5mC_CpG_v001.cfg \
-m /path/to/json/res_dna_r941_min_modbases_5mC_CpG_v001.jsn \
--device auto \
--bam_out --bam_methylation_threshold 0.05 \
--num_callers 4 \
--gpu_runners_per_device 8 \
--chunks_per_runner 1024 \
--chunk_size 1000 \
--disable_pings


###OR Dorado for R10.x flow cells
module load dorado/0.5.3
dorado basecaller hac,5mC_5hmC,6mA \
./pod5/ \
--chunksize 1000 \
--verbose \
--emit-fastq > sample_basecalls.fastq

##Map with winnowmap
meryl count k=15 output merylDB /path/to/reference.fna
meryl print greater-than distinct=0.9998 merylDB > reference_repetitive_k15.txt

winnowmap -W /path/to/reference_repetitive_k15.txt -ax map-ont --cs --eqx -Y -L -p 0.1 -I8g \
/path/to/reference.fna \
/path/to/fastq/pass/*.fastq > sample.sam

module purge 

module load gcc/8.2.0
module load samtools/1.9
module load bedtools/2.30.0

samtools view -b -h -o sample.bam sample.sam
samtools sort sample.bam > sample.sorted.bam
samtools index sample.sorted.bam
bedtools bamtobed -i sample.sorted.bam > sample.sorted.bed

#Combine guppy mod info in unaligned bam outputs with winnowmap alignment info
#Merging commands adapted from Nick Altemose lab
module purge

cd /path/to/sample/pass/
for file in bam*bam;do
	echo $file >> bam_list_orig.txt
done

module load gcc/8.2.0
module load samtools/1.14

#merge all bams in list (change @ to # cores in instance)
ulimit -n 4096
#check hard user n-limit 
split --verbose -l2000 bam_list_orig.txt bamlist.

#samtools merge -b bamlist.aa -@ 64 sample_merge_aa.bam
#samtools merge -b bamlist.ab -@ 64 sample_merge_ab.bam

echo -e "sample_merge_aa.bam\nsample_merge_ab.bam" > bam_list.txt

#samtools merge -b bam_list.txt -@ 64 sample_complete.bam

#Extract guppy mod_basecalls and read ids (readid, length, Mm, Ml)
samtools view sample_complete.bam | awk -F '\t' 'BEGIN {OFS=FS}{print $1, $9, $12, $13}' > sample_complete_extract.txt

#Get winnowmap bam header
samtools view -H sample.sorted.bam > sample.header.txt

##Combine by read id
#Remove unmapped (4), secondary (256), and supplemental (2048) alignments, 4+256+2048=2308
samtools view -b -@ 96 -F 2308 sample.sorted.bam > sample.sorted.clean.bam 

#Join winnowmap alignment with guppy/dorado Mm and Ml, **here also give TLEN from guppy (2.2 in field 9)***
#Double check field placements between guppy/dorado
#Sort by read id before joining
join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.2,1.10,1.11,2.3,2.4 <(samtools view sample.sorted.clean.bam | sort -k1 -T /path/to/temp/dir) <(sort -k1 -T /path/to/temp/dir /path/to/sample_complete_extract.txt) > sample_merge.txt

#Remove alignments with empty Mm or Ml
#Check syntac for dorado alignments; MM ML
grep -v 'Mm:Z:A+a;C+m;' sample_merge.txt > sample_merge_clean.txt

#Build combined sam by adding header to hybrid, clean read info
cat sample.header.txt sample_merge_clean.txt > sample_final.sam

#Convert to bam and index
samtools view -b sample_final.sam > sample_final.bam
samtools sort sample_final.bam > sample_final_sorted.bam
samtools index sample_final_sorted.bam

#Can be visualized with dimelo package

#For cleaning or splitting before visualization:

#To isolate a single chromosome
samtools view -b -h sample_final_sorted.bam \
   -o sample_chr7.bam NC_060931.1

#To isolate single type of modification calls
module load python/3.10

#To remove modA or modC calls, set to 256; example below
#ModBamClean_ChrParallel_v0.21.py script provided by Dan Xu, Miga Lab
python ModBamClean_ChrParallel_v0.21.py \
 --bam sample.bam --modqA 205 --modqC 256 --out sample_NoC

