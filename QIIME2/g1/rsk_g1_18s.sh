# SCOPE Amplicon 18s Pipeline
# Author: Rebecca Key
# Last updated: July 2023

# Tutorials that helped me understand: 
## https://john-quensen.com/tutorials/processing-18s-sequences-with-qiime2-and-dada2/
## https://otagoedna.github.io/getting_started_with_qiime2/first_workflow.html
## https://rpubs.com/katerinka/qiime2-tutorial

# Processing script for SCOPE ocean samples' 18s ASVs from 566F_1200R primers (634bp). Aligned to the current pr2 database.

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# STEP 1. CHECKS AND BALANCES
# Set up the most updated qiime2 (I used v2022.11 during my time). Activate in your environment.
conda activate qiime2-2022.11

# Define cruise
cruise="g1"

# Navigate to '${cruise}' folder before running commands.
# Check required content in the '${cruise}' folder. It should have 8 folders: prok, euk, qza, qzv, fastqc, phyloseq, trees, figaro-master. # And 4 files: two manifests, two .sh files 

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# STEP 2. CREATE A QIIME-FRIENDLY SAMPLE FILE
qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "manifest_${cruise}_18s.txt" \
    --output-path "./qza/${cruise}_18s_rawSeqs.qza" \
    --input-format SingleEndFastqManifestPhred33V2

# View basic sample stats.
qiime demux summarize \
    --i-data "./qza/${cruise}_18s_rawSeqs.qza" \
    --o-visualization "./qzv/${cruise}_18s_rawSeqs.qzv"  
    # NOTE: .qzv are visual files and .qza are artifact files
    # View .qzv files at https://view.qiime2.org/ . Drag 'n drop.

# STEP 2A. GENERATE OPTIONAL FASTQC FILES FOR SELECTED SAMPLES. 
# Use an additional quality assurance tool alongside qiime2. select an r1 and an r2 fastq for this process.
chmod 770 ./euk/*gz # reset to full permissions 
cp ./euk/..._R1.fastq.gz ./fastqc
cp ./euk/..._R2.fastq.gz ./fastqc
gzip -d ./fastqc/*gz 
fastqc -t 6 ./fastqc/*fastq

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# STEP 3. TRIM ADAPTERS (PRIMER, BARCODE) AND RECHECK SAMPLE STATS.
# Adjust core usage based on your computer's capacity.
qiime cutadapt trim-single \
    --i-demultiplexed-sequences "./qza/${cruise}_18s_rawSeqs.qza"  --p-cores 20 \
    --p-front CAGCAGCCGCGGTAATTCC ACACTGACGACATGGTTCTACA AATGATACGGCGACCACCGAGATCT \
    --p-adapter AGACCAAGTCTCTGC CAAGCAGAAGACGGCATACGAGAT CCCGTGTTGAGTCAAATTAAGC \
    --p-discard-untrimmed  --p-no-indels --p-error-rate 0.2  --o-trimmed-sequences "./qza/${cruise}_18s_trimmed.qza"

qiime demux summarize \
    --i-data "./qza/${cruise}_18s_trimmed.qza" \
    --o-visualization "./qzv/${cruise}_18s_trimmed.qzv"
    
# STEP 3A. VERIFY TRIMS WITH OPTIONAL POST-TRIM FASTQC FILES.
# Extract trim qza contents, then move files to your main fastqc directory.
qiime tools export \
  --input-path "./qza/${cruise}_16s_trimmed.qza" \
  --output-path "./qza/trim_export"

cp "./qza/trim_export/<filename>.fastq" ./fastqc
cp "./qza/trim_export/<filename>.fastq" ./fastqc
fastqc -t 6 *fastq

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# STEP 4. INSTALL AND USE FIGARO TO DETERMINE TRUNCATION LENGTHS. GitHub: https://github.com/Zymo-Research/figaro
# Figaro ensures objectivity, eliminates low-quality sequences before merging, and maximizes read retention.

# STEP 4A. INSTALL FIGARO
wget http://john-quensen.com/wp-content/uploads/2020/03/figaro.yml
conda env create -n figaro -f figaro.yml
wget https://github.com/Zymo-Research/figaro/archive/master.zip
unzip master.zip
rm master.zip
cd figaro-master/figaro
chmod 755 *.py # Lifts permissions
wget https://github.com/jfq3/Miscellaneous-scripts/raw/master/test_figaro.sh

# STEP 4B. OBTAIN OPTIMAL TRUNCATION LENGTHS USING FIGARO.
# Place an F and an R fastq into Figaro's sub-directory. Extract, possibly rename, and execute the script.
# Extract desired positions from the highest retention percent score, prioritizing quality over quantity.
# Caution: Avoid truncation lengths with values >=3 in max error rates for highest retention score.
dir="./figaro-master/figaro/test_figaro"
cp ./euk/F566Euk_R1200Euk-63683_CACGAAGAGC.R1.fastq.gz ${dir}
cp ./euk/F566Euk_R1200Euk-63683_CACGAAGAGC.R2.fastq.gz ${dir}
gzip -d ${dir}/*.fastq.gz

mv ${dir}/F566Euk_R1200Euk-63683_CACGAAGAGC.R1.fastq ${dir}/sam1_16s_R1.fastq
mv ${dir}/F566Euk_R1200Euk-63683_CACGAAGAGC.R2.fastq ${dir}/sam1_16s_R2.fastq
cd ./figaro-master/figaro
python figaro.py -i ./test_figaro/ -o ./test_figaro/ -f 10 -r 10 -a 253 -F zymo
cd ./test_figaro
less trimParameters.json 
q
cd ../../..

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# STEP 5. MERGE FORWARD AND REVERSE SEQUENCES USING DADA.
# Caution: If truncation lengths have >3 max error rate, default filter (--p-e-max 2) discards reads w/ >3 nt errors.
# Prioritize quality. Primer pair spans ~600bp; ~220 from forward read retains good taxonomy info.
qiime dada2 denoise-single \
    --i-demultiplexed-seqs "./qza/${cruise}_18s_trimmed.qza" \
    --p-trunc-len 220 --p-max-ee 2 \
    --p-n-threads 20 \
    --o-representative-sequences "./qza/${cruise}_18s_reprSeqs.qza" \
    --o-table "./qza/${cruise}_18s_reprSeqs_countTable.qza" \
    --o-denoising-stats "./qza/${cruise}_18s_denoiseStats.qza"

# STEP 5A. OBTAIN YOUR COUNT MATRIX AND DENOISING STATS FOR VISUALS
qiime feature-table tabulate-seqs \
    --i-data "./qza/${cruise}_18s_reprSeqs.qza" \
    --o-visualization "./qzv/${cruise}_18s_reprSeqs.qzv"

qiime metadata tabulate \
    --m-input-file "./qza/${cruise}_18s_denoiseStats.qza" \
    --o-visualization "./qzv/${cruise}_18s_denoiseStats.qzv"

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# STEP 6. PERFORM TAXONOMY ALIGNMENTS.
# Resource-intensive step. Suggested: 8 Milan CPUs, 400GB memory on HiperGator.
# Alternatively, test with 1 CPU and determine memory use. Then optimize based on system capacity.

# I recommend Hipergator (HP) to access hella resources. You will need the classifier.qza and the reprSeqs.qza
# Enter into HP. Prepare the following .sh script...

# Create a new file using the 'nano' text editor and paste the following content...
nano

#!/bin/bash
#SBATCH --job-name=outFiles/taxclass_euk.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=300GB
#SBATCH --account=b.durham
#SBATCH --qos=b.durham-b
#SBATCH -t 72:00:00
#SBATCH --partition=hpg-milan
#SBATCH --output=outFiles/taxclass_euk.%A.out
#SBATCH --array=1-1

date;hostname;pwd
module load conda 
conda activate qiime2-2022.11

db="./db/mothurPR2/pr2Classifier_v5_566F.1200R.qza"

echo "${cruise}-tax align"
qiime feature-classifier classify-sklearn \
  --p-n-jobs 8 \
  --i-classifier ${db} \
  --i-reads "./${cruise}/${cruise}_18s_reprSeqs.qza" \
  --o-classification "./results/${cruise}_18S_reprSeqs_taxAligned_pr2v5.qza"
  
echo "g2-tax align"
qiime feature-classifier classify-sklearn \
  --p-n-jobs 8 \
  --i-classifier ${db} \
  --i-reads "./g2/g2_18s_reprSeqs.qza" \
  --o-classification "./results/g2_18S_reprSeqs_taxAligned_pr2v5.qza"
  
echo "g3-tax align"
qiime feature-classifier classify-sklearn \
  --p-n-jobs 8 \
  --i-classifier ${db} \
  --i-reads "./g3/g3_18s_reprSeqs.qza" \
  --o-classification "./results/g3_18S_reprSeqs_taxAligned_pr2v5.qza"
  
# Exit out. Then sbatch file to run it.
# Transfer files from HiperGator to your local machine via FileZilla or Globus.

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# STEP 7. CHECK YOUR CLASSIFICATIONS   
  qiime metadata tabulate \
  --m-input-file "./qza/${cruise}_18S_reprSeqs_taxAligned_pr2v5_220.qza" \
  --o-visualization "./qzv/${cruise}_18s_reprSeqs_taxAligned_pr2v5_220.qzv"

# Important: Review all qzv files and download data frames/results.
# Valuable information is contained in these files! It's good record keeping and some of them you absolutely need.

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# STEP 8. MAKE FILES FOR PHYLOSEQ ANALYSIS IN THE R-ENVIRONMENT
# Generate ASV count biom, then convert to tsv. This is your ASV count table
qiime tools export \
  --input-path "./qza/${cruise}_18s_reprSeqs_countTable.qza" \
  --output-path "./phyloseq"

biom convert \
  -i "./phyloseq/feature-table.biom" \
  -o "./phyloseq/${cruise}_18s_countTable.tsv" \
  --to-tsv

sed -i '1d' ./phyloseq/${cruise}_18s_countTable.tsv
sed -i 's/#OTU ID//' ./phyloseq/${cruise}_18s_countTable.tsv

# Export representative sequences and convert to txt. These are your dna sequences for alignments
qiime tools export \
  --input-path "./qza/${cruise}_18S_reprSeqs_taxAligned_pr2Mar.qza"  \
  --output-path "./phyloseq"
awk 'BEGIN{RS=">"}{print "#"$1"\t"$2;}' ./phyloseq/dna-sequences.fasta | tail -n+2 > ./phyloseq/${cruise}_18s_dnaSeqs.txt

# Export Taxonomy. This is your taxonomy file
qiime tools export \
  --input-path "./g1-18S-base-artifact-taxonomy.qza"  \
  --output-path "."
mv ./phyloseq/taxonomy.tsv ./phyloseq/${cruise}_18s_taxonomy_pr2mar.tsv

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# STEP 9. BUILD RAXML PHYLOGENETIC TREE FROM ALIGNMENTS.
# Raxml process is time-intensive due to iterations and choosing the best tree
# Suggested: Use HP to get more cores

# Do you need a raxml of everything or just a subset you're interest in? 

# STEP 9 - NO SUBSET.
# Caution: Expect extended processing time! Be patient and think: Do I really need this?
qiime phylogeny align-to-tree-mafft-raxml \
  --p-n-threads 22 \
  --i-sequences "./qza/${cruise}_18s_reprSeqs.qza" \
  --o-alignment "./qza/${cruise}_18s_reprSeqs_align.qza" \
  --o-masked-alignment "./qza/${cruise}_18s_reprSeqs_alignMask.qza" \
  --o-tree "./trees/${cruise}_18s_treeUnRoot.qza" \
  --o-rooted-tree "./trees/${cruise}_18s_treeRoot.qza"

# STEP 9 - SUBSET. 'Focus on phytoplankton populations' example. 
# Convert your tsv taxonomy file to txt
mv "./phyloseq/${cruise}_18s_taxonomy_pr2mar.tsv" "./phyloseq/${cruise}_18s_taxonomy_pr2mar.txt"

# Update the tsv file. Maintain backups for safety. Trust me! 📂🔒 
qiime tools export \
  --input-path "./qza/${cruise}_18S_reprSeqs_taxAligned_pr2Mar.qza"  \
  --output-path "./phyloseq"
mv ./phyloseq/taxonomy.tsv ./phyloseq/${cruise}_18s_taxonomy_pr2mar.tsv

## Convert to qzas for subset
qiime tools import \
  --input-path "./phyloseq/feature-table.biom" \
  --type 'FeatureTable[Frequency]' 
  --output-path "./phyloseq/feature-table.qza"
  
qiime tools import \
  --input-path "./phyloseq/${cruise}_18s_taxonomy_pr2mar.txt" \
  --type 'FeatureData[Taxonomy]' \
  --output-path "./phyloseq/${cruise}_18s_taxonomy_pr2mar.qza"

# Apply filter: Preserve entries with "[phyto selection]" in the string.
qiime taxa filter-table \
  --i-table "./phyloseq/feature-table.qza" \
  --i-taxonomy "./phyloseq/${cruise}_18s_taxonomy_pr2mar.qza" \
  --p-include "Archaeplastida","Stramenopiles","Alveolata","Haptophyta","Cyanobacteria" \
  --o-filtered-table "./phyloseq/${cruise}_asvTable_phytos.qza"

## Extract representative sequences from this subset.
qiime feature-table filter-seqs \
  --i-data "./qza/${cruise}_18s_reprSeqs.qza" \
  --i-table "./phyloseq/${cruise}_asvTable_phytos.qza" \
  --o-filtered-data "./phyloseq/${cruise}_asvTable_phytos_sequences.qza"
  
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# STEP 10. RUN RAXML AND MAKE .NWK FILES FOR GGTREE ANALYSIS IN THE R-ENVIRONMENT
# Again, I recommend HP (so much faster!) 

# Example script to nano/save......
 
#!/bin/bash
#SBATCH --job-name=outFiles/raxml.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --account=b.durham
#SBATCH --qos=b.durham-b
#SBATCH -t 70:00:00
#SBATCH --partition=hpg-milan
#SBATCH --output=outFiles/raxml.%A.out
#SBATCH --array=1-1

date;hostname;pwd

ml conda
conda activate qiime2-2022.11

cruise="${cruise}"
echo "${cruise} Raxml phyto"
date
qiime phylogeny align-to-tree-mafft-raxml \
  --i-sequences "./${cruise}/${cruise}_asvTable_phytos_sequences.qza" \
  --o-alignment "./trees/${cruise}_18s_asvTable_align_phytos.qza" \
  --o-masked-alignment "./trees/${cruise}_18s_asvTable_alignMask_phytos.qza" \
  --o-tree "./trees/${cruise}_18s_treeUnRoot_phytos.qza" \
  --o-rooted-tree "./trees/${cruise}_18s_treeRoot_phytos.qza"
# Transfer results to your local computer folders for the next step.
  
# STEP 10A. EXPORT QZA FILES TO FILES FOR ANALYSIS IN THE R-ENVIRONMENT.
## Make useable tree files as R-studio input. Make sure you have a folder named trees in your working directory!
qiime tools export \
  --input-path "./trees/${cruise}_18s_treeUnRoot_strams.qza" \
  --output-path "./trees"
mv ./trees/tree.nwk ./trees/${cruise}_18s_raxmlUnroot_strams.nwk

qiime tools export \
  --input-path "./trees/${cruise}_18s_treeRoot_strams.qza" \
  --output-path "./trees"
mv ./trees/tree.nwk ./trees/${cruise}_18s_raxmlRoot_strams.nwk

















