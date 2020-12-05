# DGE pipeline for coral RNAseq data 2020

# Step 1: Trim Data

#Example
java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 4 [path to file 1] [path to file 2] [sample]_1.trimmed.fq.gz [sample]_1un.trimmed.fq.gz [sample]_2.trimmed.fq.gz [sample]_2un.trimmed.fq.gz ILLUMINACLIP:/opt/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25

#Trim using loop
for infile in *_1.fq.gz
do
   base=$(basename ${infile} _1.fq.gz)
   java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 4 ${infile} ${base}_2.fq.gz ${base}_1.trimmed.fq.gz ${base}_1un.trimmed.fq.gz ${base}_2.trimmed.fq.gz ${base}_2un.trimmed.fq.gz ILLUMINACLIP:/opt/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25

done >> output.txt 2>&1


#to check basenames

for infile in *_1.fq.gz
do
    base=$(basename ${infile} _1.fq.gz)
    echo ${base}
done >> output.txt 2>&1


#move trimmed files to QC folder
#delete copies of raw data

#now trimmomatic with new adapter file 
#since adapter contanination still apparent, seen with fast qc

#example
nohup java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 4 /data/corals/QCdata_full/rawdata/CoPtE2S_CKDL200153112-1B_H7WTWBBXX_L4_1.fq.gz /data/corals/QCdata_full/rawdata/CoPtE2S_CKDL200153112-1B_H7WTWBBXX_L4_2.fq.gz CoPtE2S_1.trimmed.fq.gz CoPtE2S_1un.trimmed.fq.gz CoPtE2S_2.trimmed.fq.gz CoPtE2S_2un.trimmed.fq.gz ILLUMINACLIP:/opt/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 >> output.txt

#new loop trim 
for infile in *_1.fq.gz
do
   base=$(basename ${infile} _1.fq.gz)
   java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 4 ${infile} ${base}_2.fq.gz ${base}_1.trimmed.fq.gz ${base}_1un.trimmed.fq.gz ${base}_2.trimmed.fq.gz ${base}_2un.trimmed.fq.gz ILLUMINACLIP:/opt/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 >> output.txt
done

#third time trimming, to remove polyT tails with new adapter file
#this adapter file has poly T tail sequence and reverse complement of adapter sequences added
#loop
for infile in *_1.fq.gz
do
   base=$(basename ${infile} _1.fq.gz)
   java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 4 ${infile} ${base}_2.fq.gz ${base}_1.trimmed.fq.gz ${base}_1un.trimmed.fq.gz ${base}_2.trimmed.fq.gz ${base}_2un.trimmed.fq.gz ILLUMINACLIP:/data/corals/QCdata_full/adapters_polyt.fasta:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 >> output.txt
done >> output.txt 2>&1

# Mapping

#first attempt at mapping to Barshis 2015 transcriptome
#used round 3 trimmed files 
#very low alignment, deleted output files from this mapping 
nohup perl /opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts 33496_Ahyacinthus_CoralContigs.fasta --seqType fq --samples_file coral_samples.txt --est_method RSEM --aln_method bowtie2 --prep_reference --output_dir /data/corals/QCdata_full/mapping > RSEM_out 2>&1 &

#second attempt at mapping
#used round 3 trimmed files
#used Palumbi transcriptome from csv file reformatted into fasta file 
nohup perl /opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Hyacinthus_genes_palumbi.fasta --seqType fq --samples_file coral_samples.txt --est_method RSEM --aln_method bowtie2 --prep_reference --output_dir RSEM_mapping > RSEM_out 2>&1 &

#third attempt, same as above but reformat reference transcriptome
#deleted output from attempt 2
nohup perl /opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Hyacinthus_genes_palumbi2.fasta --seqType fq --samples_file coral_samples.txt --est_method RSEM --aln_method bowtie2 --prep_reference --output_dir RSEM_mapping > RSEM_out 2>&1 &

#code to make break in header to make cleaner fasta file from csv 
#single digits before .m1 have no break between header and seq.
tr 'm1.[0-9].m1' 'm1.[0-9].m1\n' < Hyacinthus_genes_palumbi2.fasta > Hyacinthus_genes_palumbi3.fasta

#forth attempt at mapping using newly formatted ref transcriptome
nohup perl /opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Hyacinthus_genes_fastafile2.fasta --seqType fq --samples_file coral_samples.txt --est_method RSEM --aln_method bowtie2 --prep_reference --output_dir RSEM_mapping > RSEM_out 2>&1 &

#run again on clean fasta 
#remove carriage returns from fasta 
sed 's/\r//' Hyacinthus_genes_edited_09012020.fasta > Hyacinthus_genes_edited_09012020_nocarr.fasta

#fifth attempt at mapping using newly formatted ref transcriptome 
#removed carriage returns which may have been the problem
nohup perl /opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Hyacinthus_genes_edited_09012020_nocarr.fasta --seqType fq --samples_file coral_samples.txt --est_method RSEM --aln_method bowtie2 --prep_reference --output_dir RSEM_mapping > RSEM_out 2>&1 &


#sixth mapping attempt 
#made sure no duplicates (uniq -d filename.fasta) and removed carriage returns on ref transcriptome
#may have been duplicates in the last ref
nohup perl /opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Hyacinthus_genes_edited_09012020.fasta --seqType fq --samples_file coral_samples.txt --est_method RSEM --aln_method bowtie2 --prep_reference --output_dir RSEM_mapping > RSEM_out 2>&1 &

#try to convert txt file to fasta using some awk thing online
awk '{ if ( $0 ~ /^gi/ ) {gsub(" ","_",$0); print ">"$0 } else { print } }' Hyacinthus_genes_edited_fasta.txt > Hyacinthus_genes_edited_fasta.fasta

#seventh attempt at mapping using newly formatted ref transcriptome 
nohup perl /opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Hyacinthus_genes_edited_fasta.fasta --seqType fq --samples_file coral_samples.txt --est_method RSEM --aln_method bowtie2 --prep_reference --output_dir RSEM_mapping > RSEM_out 2>&1 &

#sed to reformat transcriptome
#bom might have been the problem, remove BOM
sed '1s/^\xEF\xBB\xBF//' < Hyacinthus_genes_edited_09012020_nocarr.fasta > Hyacinthus_genes_edited_09012020_nocarr_noBOM.fasta

#eighth attempt at mapping using newly formatted ref transcriptome
#no BOM
#worked!!!!!!! FINALLY 
#alignment rates much better 
nohup perl /opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Hyacinthus_genes_edited_09012020_nocarr_noBOM.fasta --seqType fq --samples_file coral_samples.txt --est_method RSEM --aln_method bowtie2 --prep_reference --output_dir RSEM_mapping > RSEM_out 2>&1 &


# Counts Matrix

nohup perl /opt/trinityrnaseq/util/abundance_estimates_to_matrix.pl --est_method RSEM  --gene_trans_map none --name_sample_by_basedir CoPt_E_rep1/RSEM.genes.results CoPt_H_rep1/RSEM.genes.results CoPt_E_rep2/RSEM.genes.results CoPt_H_rep2/RSEM.genes.results CoPt_E_rep3/RSEM.genes.results CoPt_H_rep3/RSEM.genes.results CoPt_E_rep4/RSEM.genes.results CoPt_H_rep4/RSEM.genes.results Falu_E_rep1/RSEM.genes.results Falu_H_rep1/RSEM.genes.results Falu_E_rep4/RSEM.genes.results Falu_H_rep4/RSEM.genes.results Falu_E_rep6/RSEM.genes.results Falu_H_rep6/RSEM.genes.results Ftele_E_rep3/RSEM.genes.results Ftele_H_rep3/RSEM.genes.results Ftele_E_rep5/RSEM.genes.results Ftele_H_rep5/RSEM.genes.results Ftele_E_rep6/RSEM.genes.results Ftele_H_rep6/RSEM.genes.results > RSEMae_out 2>&1 &

# Compare reps

#cor matrix
perl /opt/trinityrnaseq/Analysis/DifferentialExpression/PtR --matrix RSEM.isoform.counts.matrix --min_rowSums 10 -s coral_samples.txt --log2 --CPM --sample_cor_matrix

#tried this but no fast cluster package
#copy cor matrix to treebeard and samples txt file
scp naug7321@khaleesi.csumb.edu:/data/corals/QCdata_full/RSEM.isoform.counts.matrix /home/naug7321/RSEM.isoform.counts.matrix
scp naug7321@khaleesi.csumb.edu:/data/corals/QCdata_full/coral_samples.txt home/naug7321/coral_samples.txt

#now try again - on treebeard, works!
perl /opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix RSEM.isoform.counts.matrix --min_rowSums 10 -s coral_samples.txt --log2 --CPM --sample_cor_matrix

#PCA
#also doesn't work on khaleesi, now try treebeard 
perl /opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix RSEM.isoform.counts.matrix -s coral_samples.txt --min_rowSums 10 --log2 --CPM --center_rows --prin_comp 3 

#copy these output files to khaleesi 
scp naug7321@treebeard.csumb.edu:/home/naug7321/RSEM.isoform.counts.matrix.minRow10.CPM.log2* /data/corals/QCdata_full

# Run DGE
perl /opt/trinityrnaseq/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix RSEM.isoform.counts.matrix --method edgeR
#error - lack of replicates

#need to create replicates txt file
#replicates.txt file shows relationship between treatment and replicates
scp /Users/Melissa/Desktop/CSUMB/Thesis/Data\ analysis/Bioinformatics/DGE_RNAseq_Aug2020/replicates.txt naug7321@khaleesi.csumb.edu:/data/corals/QCdata_full

#run DGE again with reps file
perl /opt/trinityrnaseq/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix RSEM.isoform.counts.matrix --method edgeR --samples_file replicates.txt

#DGE not working, need to go back to counts matrix bc two samples missing heheh oops
#deleted all old files

#counts matrix redo
nohup perl /opt/trinityrnaseq/util/abundance_estimates_to_matrix.pl --est_method RSEM  --gene_trans_map none --name_sample_by_basedir CoPt_E_rep1/RSEM.genes.results CoPt_H_rep1/RSEM.genes.results CoPt_E_rep2/RSEM.genes.results CoPt_H_rep2/RSEM.genes.results CoPt_E_rep3/RSEM.genes.results CoPt_H_rep3/RSEM.genes.results CoPt_E_rep4/RSEM.genes.results CoPt_H_rep4/RSEM.genes.results Falu_E_rep1/RSEM.genes.results Falu_H_rep1/RSEM.genes.results Falu_E_rep3/RSEM.genes.results Falu_E_rep4/RSEM.genes.results Falu_H_rep4/RSEM.genes.results Falu_E_rep6/RSEM.genes.results Falu_H_rep6/RSEM.genes.results Ftele_H_rep1/RSEM.genes.results Ftele_E_rep3/RSEM.genes.results Ftele_H_rep3/RSEM.genes.results Ftele_E_rep5/RSEM.genes.results Ftele_H_rep5/RSEM.genes.results Ftele_E_rep6/RSEM.genes.results Ftele_H_rep6/RSEM.genes.results > RSEMae_out 2>&1 &

#copy counts matrix to treebeard
scp naug7321@khaleesi.csumb.edu:/data/corals/QCdata_full/RSEM.isoform.counts.matrix /home/naug7321/RSEM.isoform.counts.matrix

#now run cor matrix and PCR on treebeard
perl /opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix RSEM.isoform.counts.matrix --min_rowSums 10 -s coral_samples.txt --log2 --CPM --sample_cor_matrix
perl /opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix RSEM.isoform.counts.matrix -s coral_samples.txt --min_rowSums 10 --log2 --CPM --center_rows --prin_comp 3 

#copy these output files to khaleesi 
scp naug7321@treebeard.csumb.edu:/home/naug7321/RSEM.isoform.counts.matrix.minRow10.CPM.log2* /data/corals/QCdata_full

#run DGE on khaleesi
perl /opt/trinityrnaseq/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix RSEM.isoform.counts.matrix --method edgeR --samples_file replicates.txt

#cd into edgeR folder, then run this 
#can run with diff p values and diff fold changes (-c 2 is 4 fold change, -c 1 is 2 fold change) 
perl /opt/trinityrnaseq/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix /data/corals/QCdata_full/RSEM.isoform.TMM.EXPR.matrix -P 0.001 -C 2 --samples /data/corals/QCdata_full/replicates.txt

#need fast cluster, move edgeR folder to treebeard and run there 
scp -r naug7321@khaleesi.csumb.edu:/data/corals/QCdata_full/edgeR.22083.dir /home/naug7321/edgeR.22083.dir
scp naug7321@khaleesi.csumb.edu:/data/corals/QCdata_full/replicates.txt /home/naug7321/
scp naug7321@khaleesi.csumb.edu:/data/corals/QCdata_full/RSEM.isoform.TMM.EXPR.matrix /home/naug7321/

#now run in edgeR folder in treebeard
#first with r = 0.001
perl /opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix /home/naug7321/RSEM.isoform.TMM.EXPR.matrix -P 0.001 -C 2 --samples /home/naug7321/replicates.txt

#then with p = 0.05
perl /opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix /home/naug7321/RSEM.isoform.TMM.EXPR.matrix -P 0.05 -C 2 --samples /home/naug7321/replicates.txt

#done!!
#move edgeR folder back to khaleesi 
scp -r naug7321@treebeard.csumb.edu:/home/naug7321/edgeR.22083.dir /data/corals/QCdata_full/edgeR.22083.dir_DGE_ALL



























# RNAseq_DGE
