https://pmc.ncbi.nlm.nih.gov/articles/PMC8196797/#sec4-ijms-22-05426

# 1. WGBS data processing
Step 1 - TrimGalore (trim 18 bp from 3' and 5' to eliminate adaptase tail, reads <20bp were discarded)  
Step 2 - Mapping to Arabidopsis genome (TAIR10) with BSMAP  
Step 3 - Sambamba (sort, index divided by PCR duplicates)  
Step 4 - Deduplicated read were passed for methylation call    
Step 5 - Qualimap to evaluate alignment  
Step 6 - BSMAP methratio.py to extract methylation call (cutoff of 1 count)  
Step 7 - Annotate cytosine location with table browser function  
Step 8 - Methylation analysis with Bismark  
Step 9 - Deeptools to analyze and plot methylation patterns across gene bodies and 2kb upstream/ downstream (reference bed files from UCSC genome browser and ReMap Regulatory Atlas)

#extract fastq from SRA
```'/home/user/Downloads/sratoolkit.current-ubuntu64/sratoolkit.2.9.4-ubuntu64/bin/fastq-dump' --defline-seq '@$ac[_$sn]/$ri' --split-files '/mnt/D/Liaw/WGS-Bi/SRR8100584.sra' -O '/mnt/D/Liaw/WGS-Bi'```

Step 1 - Trimming
Trim fastq files (remove adapter sequences and also trim 10 bp from 5'-end of reads)
```'/mnt/D/TrimGalore-0.6.0/trim_galore' --paired --clip_R1 10 --clip_R2 10 --retain_unpaired -r1 35 -r2 35 '/mnt/D/Liaw/WGS-Bi/SRR8100584_1.fastq' '/mnt/D/Liaw/WGS-Bi/SRR8100584_2.fastq'```

Step 8 - Bismark
```
#indexing
/home/user/Downloads/Bismark/bismark_genome_preparation --path_to_bowtie '/usr/local/bin' --verbose '/mnt/D/Liaw/WGS-Bi/arabidopsis'

#alignment
**first, run PE bismark with directional, --unmapped option and default setting
default (paired-end):
'/home/user/Downloads/Bismark/bismark' --path_to_bowtie '/usr/local/bin' --samtools_path /usr/local/bin --bowtie2 --unmapped -o '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' -N 0 -L 20 -R 2 -D 15 --score_min L,0,-0.2 --nucleotide_coverage --genome /mnt/D/Liaw/WGS-Bi/arabidopsis -1 '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_1.fastq' -2 '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_2.fastq' 1>SRR8100584_Bismark_0.2_stringent.log 2>SRR8100584_Bismark_0.2_stringent.err

**second, run SE bismark with directional option and default setting using unmapped read 1
'/home/user/Downloads/Bismark/bismark' --path_to_bowtie '/usr/local/bin' --samtools_path /usr/local/bin --bowtie2 -o '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' -N 0 -L 20 -R 2 -D 15 --score_min L,0,-0.2 --nucleotide_coverage --genome /mnt/D/Liaw/WGS-Bi/arabidopsis --single_end /mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1.fastq_unmapped_reads_1.fq.gz 1>SRR8100584_unmapped_R1_Bismark_0.2_stringent.log 2>SRR8100584_unmapped_R1_Bismark_0.2_stringent.err

**third, run SE bismark with directional, -pbat options and default setting using unmapped read 2
'/home/user/Downloads/Bismark/bismark' --path_to_bowtie '/usr/local/bin' --samtools_path /usr/local/bin --bowtie2 -o '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' -pbat -N 0 -L 20 -R 2 -D 15 --score_min L,0,-0.2 --nucleotide_coverage --genome /mnt/D/Liaw/WGS-Bi/arabidopsis --single_end /mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100591_2.fastq_unmapped_reads_2.fq.gz 1>SRR8100584_unmapped_R2_Bismark_0.2_stringent.log 2>SRR8100584_unmapped_R2_Bismark_0.2_stringent.err

#deduplication
**for PE bismark
'/home/user/Downloads/Bismark/deduplicate_bismark' --paired --output_dir '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' --bam --samtools_path /usr/local/bin '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.bam'

**for SE bismark of unmapped read 1/2
'/home/user/Downloads/Bismark/deduplicate_bismark' --single --output_dir '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' --bam --samtools_path /usr/local/bin '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1.fastq_unmapped_reads_1_bismark_bt2_SE.bam'

'/home/user/Downloads/Bismark/deduplicate_bismark' --single --output_dir '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' --bam --samtools_path /usr/local/bin '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100591_2.fastq_unmapped_reads_2_bismark_bt2_SE.bam'

#bismark_methylation_extractor
'/home/user/Downloads/Bismark/bismark_methylation_extractor' --paired-end --gzip --bedGraph --cytosine_report --CX --multicore 4 --comprehensive --genome_folder /mnt/D/Liaw/WGS-Bi/arabidopsis '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.bam'

'/home/user/Downloads/Bismark-master/bismark_methylation_extractor' --bedGraph '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.bam'

'/home/user/Downloads/Bismark-master/bismark_methylation_extractor' --gzip '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100591_1_bismark_bt2_pe.deduplicated.bam'

**for PE bismark, extract methylation using --no_overlap
'/home/user/Downloads/Bismark/bismark_methylation_extractor' --paired-end --no_overlap --gzip --bedGraph --cytosine_report --CX --multicore 4 --comprehensive --genome_folder /mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10 '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.bam'

**for SE bismark of unmapped read 1/2, extract methylation using default setting
'/home/user/Downloads/Bismark/bismark_methylation_extractor' --single-end --gzip --bedGraph --cytosine_report --CX --multicore 4 --comprehensive --genome_folder /mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10 '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1.fastq_unmapped_reads_1_bismark_bt2_SE.deduplicated.bam'

'/home/user/Downloads/Bismark/bismark_methylation_extractor' --single-end --gzip --bedGraph --cytosine_report --CX --multicore 4 --comprehensive --genome_folder /mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10 '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_2.fastq_unmapped_reads_2_bismark_bt2_SE.deduplicated.bam'

**combined CpG_context output files
cat CpG_context_SRR8100584*.txt > CpG_context_SRR8100584_merged.txt

**combined CHG_context output files
cat CHG_context_SRR8100584*.txt > CHG_context_SRR8100584_merged.txt

**combined CHH_context output files
cat CHH_context_SRR8100584*.txt > CHH_context_SRR8100584_merged.txt

#bismark2bedGraph
/home/user/Downloads/Bismark/bismark2bedGraph -o CpG_context_SRR8100584_merged.bedGraph '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/CpG_context_SRR8100584_merged.txt'

/home/user/Downloads/Bismark/bismark2bedGraph --CX -o CHG_context_SRR8100584_merged.bedGraph '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/CHG_context_SRR8100584_merged.txt'

/home/user/Downloads/Bismark/bismark2bedGraph --CX -o CHH_context_SRR8100591_merged.bedGraph '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/CHHcontext_SRR8100584_merged.txt'

*output is .bedGraph and .cov

#alignment
default (single-end):
'/home/user/Downloads/Bismark/bismark' --path_to_bowtie '/usr/local/bin' --samtools_path /usr/local/bin --bowtie2 --phred33-quals --unmapped -o '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' -N 0 -L 20 -R 2 -D 15 --score_min L,0,-0.2 --non_directional --nucleotide_coverage --genome /mnt/D/Liaw/WGS-Bi/arabidopsis --single_end '/home/user/Downloads/Bismark/test_data.fastq' 1>Bismark_0.2_stringent.log 2>Bismark_0.2_stringent.err

alternative 1 (single-end, 1 mismatch allowed & fast/less sensitive alignment)
'/home/user/Downloads/Bismark/bismark' --path_to_bowtie '/usr/local/bin' --samtools_path '/usr/local/bin' --bowtie2 --phred33-quals --unmapped -o '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' -N 1 -L 32 --score_min L,0,-0.2 --non_directional --nucleotide_coverage --genome /mnt/D/Liaw/WGS-Bi/arabidopsis --single_end '/home/user/Downloads/Bismark/test_data.fastq' 1>Bismark_0.2_stringent.log 2>Bismark_0.2_stringent.err

alternative 1 (paired-end; 1 mismatch allowed & fast/less sensitive alignment)
'/home/user/Downloads/Bismark/bismark' --path_to_bowtie '/usr/local/bin' --samtools_path '/usr/local/bin' --bowtie2 --phred33-quals --unmapped -o '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' -N 1 -L 32 --score_min L,0,-0.2 --non_directional --nucleotide_coverage --genome /mnt/D/Liaw/WGS-Bi/arabidopsis -1 '/home/user/Downloads/bwa-meth-master/example/t_R1.fastq.gz' -2 '/home/user/Downloads/bwa-meth-master/example/t_R2.fastq.gz' 1>Bismark_0.2_stringent.log 2>Bismark_0.2_stringent.err

Other options:
'/home/user/Downloads/Bismark/bismark' --path_to_bowtie '/usr/local/bin' --samtools_path /usr/local/bin --bowtie2 --phred33-quals --unmapped -o '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' -N 1 -L 20 -R 2 -D 30 --score_min L,20,-1 --non_directional --nucleotide_coverage --genome /mnt/D/Liaw/WGS-Bi/arabidopsis -1 '/home/user/Downloads/bwa-meth-master/example/t_R1.fastq.gz' -2 '/home/user/Downloads/bwa-meth-master/example/t_R2.fastq.gz' 1>Bismark_D30_L,20,-1_N1.log 2>Bismark_D30_L,20,-1_N1.err

'/home/user/Downloads/Bismark/bismark' --path_to_bowtie '/usr/local/bin' --samtools_path /usr/local/bin --bowtie2 --phred33-quals --unmapped -o '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' -N 0 -L 20 -R 2 -D 30 --score_min L,20,-1 --non_directional --nucleotide_coverage --genome /mnt/D/Liaw/WGS-Bi/arabidopsis -1 '/home/user/Downloads/bwa-meth-master/example/t_R1.fastq.gz' -2 '/home/user/Downloads/bwa-meth-master/example/t_R2.fastq.gz' 1>Bismark_D30_L,20,-1.log 2>Bismark_D30_L,20,-1.err

'/home/user/Downloads/Bismark/bismark' --path_to_bowtie '/usr/local/bin' --samtools_path /usr/local/bin --bowtie2 --phred33-quals --unmapped -o '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' -N 0 -L 20 -R 2 -D 30 --score_min L,0,-1.5 --non_directional --nucleotide_coverage --genome /mnt/D/Liaw/WGS-Bi/arabidopsis -1 '/home/user/Downloads/bwa-meth-master/example/t_R1.fastq.gz' -2 '/home/user/Downloads/bwa-meth-master/example/t_R2.fastq.gz' 1>Bismark_D30_L,0,-1.5.log 2>Bismark_D30_L,0,-1.5.err

'/home/user/Downloads/Bismark/bismark' --path_to_bowtie '/usr/local/bin' --samtools_path /usr/local/bin --bowtie2 --phred33-quals --unmapped -o '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' -N 0 -L 20 -R 2 -D 30 --score_min L,0,-0.6 --non_directional --nucleotide_coverage --genome /mnt/D/Liaw/WGS-Bi/arabidopsis -1 '/home/user/Downloads/bwa-meth-master/example/t_R1.fastq.gz' -2 '/home/user/Downloads/bwa-meth-master/example/t_R2.fastq.gz' 1>Bismark_D30_L,0,-0.6.log 2>Bismark_D30_L,20,-0.6.err

#bam2nuc
'/home/user/Downloads/Bismark/bam2nuc' --dir '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10' --samtools_path /usr/local/bin/samtools-1.3 --genomic_composition_only --genome_folder '/mnt/D/Liaw/WGS-Bi/arabidopsis' ./SRR8100584_1_bismark_bt2_pe.deduplicated.bam

#bismark2report
'/home/user/Downloads/Bismark/bismark2report' --alignment_report ./SRR8100584_1_bismark_bt2_pe_report.txt --splitting_report ./SRR8100584_1_bismark_bt2_pe.deduplicated_splitting_report.txt --mbias_report ./SRR8100584_1_bismark_bt2_pe.deduplicated.M-bias.txt

'/home/user/Downloads/Bismark/bismark2report' --alignment_report ./SRR8100584_1_bismark_bt2_pe_report.txt --dir .

#coverage2cytosine
'/home/user/Downloads/Bismark/coverage2cytosine' -o ./SRR8100584_1_bismark_bt2_pe.deduplicated --dir . --genome_folder '/mnt/D/Liaw/WGS-Bi/arabidopsis/' --CX --gc --ff ./SRR8100584_1_bismark_bt2_pe.deduplicated.bismark.cov.gz

*output is SRR8100584_1_bismark_bt2_pe.deduplicated.CX_report.txt

#extract CG data from result of coverage2cytosine
awk '{if($6 == "CG") print}' '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.CX_report.txt' > '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.CG_report.txt'

#extract CHG data from result of coverage2cytosine
awk '{if($6 == "CHG") print}' '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.CX_report.txt' > '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.CHG_report.txt'

#extract CHH data from result of coverage2cytosine
awk '{if($6 == "CHH") print}' '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.CX_report.txt' > '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.CHH_report.txt'

#count number of CG in column 6
awk '$6 == "CG" {count++} END {print count}' /mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.CpG_report.txt

#Summary report (for more than 1 bam files)
'/home/user/Downloads/Bismark/bismark2summary' --title test_data_summary_report --basename ./SRR8100584_1_bismark_bt2_summary_report ./SRR8100584_1_bismark_bt2.deduplicated.bam ./SRR8100583_1_bismark_bt2.deduplicated.bam

#analyse differential methylation using TEA, only for TAIR10
java -jar /home/user/Downloads/EpiMolas.jar '/mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.CX_report.txt' '/mnt/D/Liaw/WGS-Bi/arabidopsis/Arabidopsis_thaliana.TAIR10.42.gtf' > /mnt/D/Liaw/WGS-Bi/arabidopsis/SRR8100584_TAIR10/SRR8100584_1_bismark_bt2_pe.deduplicated.CX_report.mtable

upload mtable online to http://tea.iis.sinica.edu.tw/tea/molas.html

#Extract Reads-Pairs Aligned Concordantly Exactly 1 Time
samtools view -hf 0x2 alignments.bam | grep -v "XS:i:" > filtered.alignments.sam
/home/user/Downloads/Bismark/foo.py filtered.alignments.sam > filtered.alignments.bam

samtools view -hf 0x2 alignments.bam | grep -v "XS:i:" | /home/user/Downloads/Bismark/foo.py > filtered.alignments.bam
```

Step 9 - Deeptools
```
## first need to convert bedgraph to bigwig file (bedgraph need to be sorted first)
## sort bedgraph 
# file could not be opened
sort -k1,1 -k2,2n '/mnt/D/Liaw/HN00107334/analysis/bismark/TAIR10/trimmed/CpG_meth.bedgraph.gz' > CpG_meth.bedgraph.sorted.gz

# use bedtools sort
bedtools sort -i '/mnt/D/Liaw/HN00107334/analysis/bismark/TAIR10/trimmed/CpG_meth.bedgraph.gz' > CpG_meth.bedgraph.bedtools.sorted.gz

## add Chr in front
awk '{ print "Chr" $0 }' CpG_meth.bedgraph.bedtools.sorted > CpG_meth_sorted_Chrprefix.bedgraph

# bedgraph to bigwig (from bismark methylation extractor)
/mnt/D/Ubuntu/Downloads/bedGraphToBigWig '/mnt/D/Liaw/HN00107334/analysis/bismark/TAIR10/trimmed/CpG_meth.bedgraph.bedtools.sorted'  '/mnt/D/Liaw/ref/TAIR10/TAIR10_chromsize.bed'  h300d_TAIR10_CpG.chr.bw

## use deeptools to compute matrix (bed file need to be regions with genes or repeats, so the program can binned the meth onto that regions)
# -b -a means upstream downstread (I used 2000, means 2 kb in this case)
# -S is your bigwig file
computeMatrix scale-regions -R '/home/user/Desktop/Liaw/ref/TAIR10/TAIR10.bed' -S h300d_TAIR10_CpG.bw -b 2000 -a 2000 --binSize 500 -o h300d_TAIR10_CpG_500bp_bin.matrix; plotProfile -m h300d_TAIR10_CpG_500bp_bin.matrix -o h300d_TAIR10_CpG_500bp_bin.png

--binSize 100 (default: 10)

# the matrix is used to plot the graph
plotProfile -m h300d_allC_context.matrix --perGroup --samplesLabel CG CHG CHH -o h300d_TAIR_allC.png

```




# 2. Differential Methylation Analysis
Step 1 - Download SRA
Step 2 - FastQC
Step 3 - Bismark and convert to bedgraph using bismark2bedGraph
Step 4 - methylKitR to identify differentially methylated region (DMR) by Fisher exact test
Step 5 - Export DMR as bedgraphs and analyzed for repeats using Bedtools intersect
extract_genes.sh
```
#!/bin/bash
# to intersect bed region with gff to extract genes and check overlap with microarray

# while getopts "sb:" opt; do
#    case $opt in
#        s) use_sra=1 ;;
#        b) branch="$OPTARG" ;;
#        *) echo "Usage: $0 [-s] [-b <branch_name>]" && exit 1
#    esac
# done
# shift $(($OPTIND - 1))

region=hypo.promoter
sample=SRR578938

# compare against a ( -a should be gff file)
bedtools intersect -wo -a '/mnt/D/Liaw/ref/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff' -b  '/mnt/D/Liaw/WGBS/methylKit/SRR578938/hypo.promoter.SRR578938.bed' > ${region}_genes_$sample.txt

#${region}_${sample}.txt

# gene ID at column 9
awk '{print $9}' ${region}_genes_$sample.txt | awk -F ';' '{print $3}' | awk -F ',' '{print $1}' | awk -F '=' '{ print $1 "\n" $2 }' | awk -F: '{print $2}' | sort -u | grep AT > Diff25.${region}_genesID_$sample.txt

rm ${region}_genes_$sample.txt

# -w word-regexp (match whole word) 
# -F fixed-string (PATTERN is a set of newline-separated strings)
# -f file (obtain PATTERN from FILE)
grep -wFf '/mnt/D/Liaw/WGBS/DMR_genes/bedtools_intersect_testing/microarray_genes_sorted.txt' Diff25.${region}_genesID_$sample.txt   > WGBSvsMicroarray_${region}_${sample}_sorted.txt

wc -l WGBSvsMicroarray_${region}_${sample}_sorted.txt

# depends which type of chromosome naming used
# '/mnt/D/Liaw/ref/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff'
# '/mnt/D/Liaw/ref/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic_ChrNotNC.gff'
```

bedtools_intersect.sh
```
#find overlap in reads and reference and extract genelist (through bedtools interesect)

# -wo	Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported. Restricted by -f and -r.
bedtools intersect -wo -f 0.5 -a '/mnt/D/Liaw/HN00102140/bam/Hybridcell-01.bed' -b /mnt/D/Liaw/WGS-Bi/arabidopsis/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff > genesinhybrid-01.txt  

# -wa	Write the original entry in A for each overlap.
bedtools intersect -wb -f 0.5 -abam /mnt/D/Liaw/HN00102140/bam/Hybridcell-01_rmdup.bam -b /mnt/D/Liaw/WGS-Bi/arabidopsis/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff -bed > genesinhybrid-01-frombamgff.txt 
 
samtools view -H '/mnt/D/Liaw/HN00102140/bam/genesinhybrid-01.txt' 

awk '{print $1, $15, $16}' '/mnt/D/Liaw/HN00102140/bam/genesinhybrid-01.txt' > genelist.txt

grep gene genelist.txt | sort -u -t > filteronlygenes.txt 
sort -u -t filteronlygenes.txt > filteronlygenes.txt
sort -u filteronlygenes.txt > filteronlygenes.txt
grep GeneID genelist.txt | sort -u > filteronlygenes.txt 

bedtools intersect -wo -wb -f 0.95 -a '/mnt/D/Liaw/HN00102140/bam/Hybridcell-01.bed' -b /mnt/D/Liaw/WGS-Bi/arabidopsis/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff > genesinhybrid-01.txt  

bedtools intersect -wo -f 0.95 -a '/mnt/D/Liaw/HN00102140/bam/Hybridcell-01.bed' -b /mnt/D/Liaw/WGS-Bi/arabidopsis/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff > genesinhybrid-01.txt  

################################################################################

# -u unique
# -v Reporting the absence of any overlapping features
# -f Requiring a minimal overlap fraction

bedtools intersect -wo -f 0.5 -a '/mnt/D/Liaw/WGBS/methylKit/hypo_region25.bed' -b '/mnt/D/Liaw/ref/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff' > hypo_genes_WGBS.txt

awk '{print $18}' hypo_genes_WGBS.txt | awk -F ';' '{print $3}' | awk -F ',' '{print $1}' > hypo_genes.txt

awk -F: '{ print $1 "\n" $2 }' WGBS_hypo_genes_sorted.txt | sort -u > WGBS_hypo_genes_sorted_removedTAIRAraportwithawk.txt

# -w word-regexp (match whole word) 
# -F fixed-string (PATTERN is a set of newline-separated strings)
# -f file (obtain PATTERN from FILE)
# -o --only-matching show only the part of a line matching PATTERN
grep -wFf microarray_genes_sorted.txt WGBS_hypo_genes_sorted_removedTAIRAraportwithawk.txt > WGBSvsMicroarray_sorted.txt

# -v to reverse -- only output unmatched -- check which genes in microarray that is not present in WGBS hypo region
grep -vf WGBSvsMicroarray_sorted.txt microarray_genes_sorted.txt > microarrayNotInWGBShypo.txt

```



