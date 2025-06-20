# Amplicon 18s Pipeline
# Author: Rebecca Key

# Tutorials: 
## https://john-quensen.com/tutorials/processing-18s-sequences-with-qiime2-and-dada2/
## https://otagoedna.github.io/getting_started_with_qiime2/first_workflow.html
## https://rpubs.com/katerinka/qiime2-tutorial

# This script processes 18s amplicon sequence variants (ASVs) derived from SCOPE ocean samples. Sequences were amplified from the 566F_1200R primer pair (634bp long) and aligned to the most up-to-date pr2-mothur database 

#-----------------------------------------------------------------------------------

# Step 1. Install qiime2-2022.8 in your system's conda environment
# https://forum.qiime2.org/t/qiime-2-2022-8-is-now-available/23994
conda create -n qiime2-community-2022.8-BETA \
  -c conda-forge -c bioconda \
  -c https://packages.qiime2.org/qiime2/2022.8/passed/community/ \
  qiime2-community

# Activate your environment before performing any commands in your terminal
conda activate qiime2-2022.11

#-----------------------------------------------------------------------------------

# Step 2. Create a Qiime-friendly sample file; view basic sample stats. What you need: forward seq, reverse seq, and manifest text file. Your main directory should have the following sub-directories: prok, euk, qza, qzv, and fastqc. Your main directory should 4 files: 2 .txt files and 2 .sh files. 
qiime tools import \
  	--type 'SampleData[SequencesWithQuality]' \
  	--input-path "manifest_g2_18s.txt" \
  	--output-path "./qza/g2_18s_rawSeqs.qza" \
  	--input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
	--i-data "./qza/g2_18s_rawSeqs.qza" \
	--o-visualization "./qzv/g2_18s_rawSeqs.qzv"  
	# NOTE: .qzv are visual files; .qza are artifact files
	# View .qzv files at https://view.qiime2.org/
	# Drag n drop 

# Step 2a. Generate FASTQC files from couple of samples (optional). In addition to qiime2's method for visualizing quality, we also want another application for extra quality assurance. Choose a R1 and a R2 fastq for these steps
chmod 770 ./euk/*gz # optional command, this will reset to full permissions 
cp ./euk/..._R1.fastq.gz ./fastqc
cp ./euk/..._R2.fastq.gz ./fastqc
gzip -d ./fastqc/*gz 
fastqc -t 6 ./fastqc/*fastq

# Step 3. Remove adapters (primer, barcode information); recheck sample stats
qiime cutadapt trim-single \
	--i-demultiplexed-sequences "./qza/g2_18s_rawSeqs.qza"  --p-cores 20 \
	--p-front CAGCAGCCGCGGTAATTCC ACACTGACGACATGGTTCTACA AATGATACGGCGACCACCGAGATCT \
	--p-adapter AGACCAAGTCTCTGC CAAGCAGAAGACGGCATACGAGAT CCCGTGTTGAGTCAAATTAAGC \
	--p-discard-untrimmed  --p-no-indels --p-error-rate 0.2  --o-trimmed-sequences "./qza/g2_18s_trimmed.qza"

qiime demux summarize \
	--i-data "./qza/g2_18s_trimmed.qza" \
	--o-visualization "./qzv/g2_18s_trimmed.qzv"
	
# Step 3a. Generate FASTQC files post-trim (optional). First, break down the trim qza. Enter into the directory, move the files to your fastqc directory in main 
qiime tools export \
  --input-path "./qza/g1_16s_trimmed.qza" \
  --output-path "./qza/trim_export"

cp "./qza/trim_export/<filename>.fastq" ./fastqc
cp "./qza/trim_export/<filename>.fastq" ./fastqc
fastqc -t 6 *fastq

# Step 4. Install Figaro, Acquire truncation lengths. https://github.com/Zymo-Research/figaro
# We use Figaro to remove human subjectiveness (increase reproducibility), remove low quality sequence before merging F and R sequences, and maximize read retention  

# Step 4a. Install Figaro in conda environment
wget http://john-quensen.com/wp-content/uploads/2020/03/figaro.yml
conda env create -n figaro -f figaro.yml
wget https://github.com/Zymo-Research/figaro/archive/master.zip
unzip master.zip
rm master.zip
cd figaro-master/figaro
chmod 755 *.py # Lifts permissions
wget https://github.com/jfq3/Miscellaneous-scripts/raw/master/test_figaro.sh

# Step 4b. Acquire best truncation lengths. Move two fastaq into the Figaro sub-directory. One F seq and one R seq will do. Unzip, rename (optional); run the script; trim from the file what you need; take positions from the highest retention percent score. Note: Watch the max error rates on your highest retention score. Ideally, you want to choose quality over quantity. Don't use truncation lengths from the highest percent score if one of the values is >=3! 
dir="./figaro-master/figaro/test_figaro"
cp ./euk/F566Euk_R1200Euk-68487_GACTCATGCT.R1.fastq.gz ${dir}
cp ./euk/F566Euk_R1200Euk-68487_GACTCATGCT.R2.fastq.gz ${dir}
gzip -d ${dir}/*.fastq.gz

mv ${dir}/F566Euk_R1200Euk-68487_GACTCATGCT.R1.fastq ${dir}/sam1_16s_R1.fastq
mv ${dir}/F566Euk_R1200Euk-68487_GACTCATGCT.R2.fastq ${dir}/sam1_16s_R2.fastq
cd ./figaro-master/figaro
python figaro.py -i ./test_figaro/ -o ./test_figaro/ -f 10 -r 10 -a 253 -F zymo
cd ./test_figaro
less trimParameters.json 
q
cd ../../..

## Step 5. Merge your forward and reverse sequences together using dada. Note: If you chose truncation lengths with >3 max error rate, every read w/ >3 nt errors will be thrown out by the command's default filter (--p-e-max 2). I recommend "quality over quantity". Also, the 18s primer pair spans ~600bp, we can get away with ~200 from the forward read and still make good taxonomy classifications
qiime dada2 denoise-single \
	--i-demultiplexed-seqs "./qza/g2_18s_trimmed.qza" \
	--p-trunc-len 220 --p-max-ee 2 \
	--p-n-threads 20 \
	--o-representative-sequences "./qza/g2_18s_reprSeqs.qza" \
	--o-table "./qza/g2_18s_reprSeqs_countTable.qza" \
	--o-denoising-stats "./qza/g2_18s_denoiseStats.qza"

# Step 5a. Obtain your count matrix and denoising stats for visuals
qiime feature-table tabulate-seqs \
	--i-data "./qza/g2_18s_reprSeqs.qza" \
	--o-visualization "./qzv/g2_18s_reprSeqs.qzv"

qiime metadata tabulate \
	--m-input-file "./qza/g2_18s_denoiseStats.qza" \
	--o-visualization "./qzv/g2_18s_denoiseStats.qzv"
	
# Step 6. Taxonomy Alignments. You are now ready to align your representatives to the taxonomy database. Note: This step uses a lot of computer resources. To be quick, I recommend 8 Milan CPUs, 400GB memory on HiperGator as provided below. OR, you can test 1 CPU on your system, determine how much memory is needed (that won't crash your computer, ie. deplete both MEM and SWP memory allocation), then optimize for your system 

# Make the Hipergator script... nano and paste what is below. You will need two inputs within the working directory: the classifier and the representatives
nano

#!/bin/bash
#SBATCH --job-name=outFiles/taxclass.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --account=b.durham
#SBATCH --qos=b.durham-b
#SBATCH -t 95:00:00
#SBATCH --partition=hpg-milan
#SBATCH --output=outFiles/taxclass.%A.out
#SBATCH --array=1-1

date;hostname;pwd
 
module load qiime2/2022.8

dir="g2/qza/" # Change to fit your pathway structure

qiime feature-classifier classify-sklearn \
  --p-n-jobs 8 \
  --i-classifier "./${dir}/silvaClassifier_nr99_515F.806R.qza" \
  --i-reads "./${dir}/g2_16s_reprSeqs.qza" \
  --o-classification "./${dir}/g2_18s_reprSeqs_taxAligned.qza"
  
# Exit out, sbatch <filename>.sh
# Get out of HiperGator, push files to respected folders in your computer 
#----------------------------------------------------------------------------

# Step 7. Check your classifications 
qiime metadata tabulate \
  --m-input-file "./qza/g2_18s_reprSeqs_taxAligned.qza" \
  --o-visualization "./qzv/g2_18s_reprSeqs_taxAligned.qzv"

# Note: I recommend going through all your qzv and downloading data frames/results. There's important stuff in there!

# Step 8. Construct a phylogenetic tree from the alignments
qiime phylogeny align-to-tree-mafft-raxml \
  --p-n-threads 11 \
  --i-sequences "./qza/g2_18s_reprSeqs.qza" \
  --o-alignment "./qza/g2_18s_reprSeqs_align.qza" \
  --o-masked-alignment "./qza/g2_18s_reprSeqs_alignMask.qza" \
  --o-tree "./trees/g2_18s_treeUnRoot.qza" \
  --o-rooted-tree "./trees/g2_18s_treeRoot.qza"

# Step 9 (optional). Export qza files to files for analysis in the R-environment.

# Step 9a. Make a tree file. Make sure you have a folder named phyloseq in your working directory!
qiime tools export \
  --input-path "./trees/g2_18s_treeUnRoot.qza" \
  --output-path "./trees"
mv ./trees/tree.nwk ./trees/g2_18s_raxmlUnroot.nwk

qiime tools export \
  --input-path "./trees/g2_18s_treeRoot.qza" \
  --output-path "./trees"
mv ./trees/tree.nwk ./trees/g2_18s_raxmlRoot.nwk

## Step 9a. Generate a biom for ASV counts, then convert it to a tsv file
qiime tools export \
  --input-path "./qza/g2_18s_reprSeqs_countTable.qza" \
  --output-path "./phyloseq"

biom convert \
  -i "./phyloseq/feature-table.biom" \
  -o "./phyloseq/g2_18s_countTable.tsv" \
  --to-tsv

sed -i '1d' ./phyloseq/g2_18s_countTable.tsv
sed -i 's/#OTU ID//' ./phyloseq/g2_18s_countTable.tsv

## Step 9b. Export representative sequences. Convert to txt file
qiime tools export \
  --input-path "./qza/g2_18s_reprSeqs.qza"  \
  --output-path "./phyloseq"
awk 'BEGIN{RS=">"}{print "#"$1"\t"$2;}' ./phyloseq/dna-sequences.fasta | tail -n+2 > ./phyloseq/g2_18s_dnaSeqs.txt

## Step 9c. Export Taxonomy
qiime tools export \
  --input-path "./qza/g2_18s_reprSeqs_taxAligned.qza"  \
  --output-path "./phyloseq"
mv ./phyloseq/taxonomy.tsv ./phyloseq/g2_18s_taxonomy.tsv

# END 






