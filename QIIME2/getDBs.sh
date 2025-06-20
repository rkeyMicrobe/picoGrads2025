# Database Classifiers in Qiime2
# Author: Rebecca Key (rebeccakey@ufl.edu)

# Tutorials that helped me understand: 
## https://otagoedna.github.io/getting_started_with_qiime2/first_workflow.html
## https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494

# This script creates amplicon-specific taxonomy classifiers. These classifiers can then be used to align 16s, 18s representative sequences in your qiime2 pipeline. 
# Please run these before you get to Step 6 in each 16s, 18s script (sh) file. These files are located in each working cruise directory. 

# First, activate the qiime2 environment before performing any commands in terminal
conda activate qiime2-2022.11

# make sure you have the rescript plug-in. This'll need to be done if you never used it on your computer!
# https://github.com/bokulich-lab/RESCRIPt/


#------------------------------------------------------------ 
#------------------------------------------------------------ 
#------------------------------------------------------------ 

# Silva138.1 --- Prepare a current Silva database (clustering-99%-SSU) specific for the V4.515F_V4.806 primer pair. We will be using the Rescript plug-in
# Step 1. Retrieve your database. Covert rna to dna sequences
qiime rescript get-silva-data \
    --p-version '138.1' \
    --p-target 'SSURef_NR99' \
    --p-include-species-labels \
    --o-silva-sequences "../db/silva138.1/silva138.1_nr99_rna.qza" \
    --o-silva-taxonomy "../db/silva138.1/silva138.1_nr99_tax.qza"  

qiime rescript reverse-transcribe \
	--i-rna-sequences "../db/silva138.1/silva138.1_nr99_rna.qza" \
	--o-dna-sequences "../db/silva138.1/silva138.1_nr99_dna.qza"

# Step 2. Clean db sequences. Remove sequences with >5 ambiguous bases and >8 homopolymer regions 
qiime rescript cull-seqs \
	--i-sequences "../db/silva138.1/silva138.1_nr99_dna.qza" \
	--o-clean-sequences "../db/silva138.1/silva138.1_nr99_dna_clean.qza"

# Step 3. Filter db sequences. Silva contains archae, bacteria, and eukaryotes. We don't want to ignore archae and euk ASVs even though these samples are 16s amplicon sequencing variants
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences "../db/silva138.1/silva138.1_nr99_dna_clean.qza"\
    --i-taxonomy "../db/silva138.1/silva138.1_nr99_tax.qza"  \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs "../db/silva138.1/silva138.1_nr99_dna_filter.qza" \
    --o-discarded-seqs "../db/silva138.1/silva138.1_nr99_dna_discard.qza" 
    # We will push discards into their own file

# Step 4. Dereplicate db sequences
qiime rescript dereplicate \
    --i-sequences "../db/silva138.1/silva138.1_nr99_dna_filter.qza"  \
    --i-taxa "../db/silva138.1/silva138.1_nr99_tax.qza" \
    --p-rank-handles 'silva' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "../db/silva138.1/silva138.1_nr99_dna_derep1.qza" \
    --o-dereplicated-taxa "../db/silva138.1/silva138.1_nr99_tax_derep1.qza"

# Step 5. From your sequences, pull out primer-specific reads. Dereplicate again
qiime feature-classifier extract-reads \
    --i-sequences "../db/silva138.1/silva138.1_nr99_dna_derep1.qza" \
    --p-f-primer GTGYCAGCMGCCGCGGTAA \
    --p-r-primer GGACTACNVGGGTWTCTAAT \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads "../db/silva138.1/silva138.1_nr99_515F.806R_seqs.qza"

qiime rescript dereplicate \
    --i-sequences "../db/silva138.1/silva138.1_nr99_515F.806R_seqs.qza" \
    --i-taxa "../db/silva138.1/silva138.1_nr99_tax_derep1.qza" \
    --p-rank-handles 'silva' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "../db/silva138.1/silva138.1_nr99_515F.806R_seqs_derep.qza" \
    --o-dereplicated-taxa  "../db/silva138.1/silva138.1_nr99_515F.806R_tax_derep.qza"

# Step 6. Create your final product before aligning to your representatives... the classifier file! 
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "../db/silva138.1/silva138.1_nr99_515F.806R_seqs_derep.qza" \
    --i-reference-taxonomy "../db/silva138.1/silva138.1_nr99_515F.806R_tax_derep.qza" \
    --o-classifier "../db/silva138.1/silvaClassifier_nr99_515F.806R.qza"
    
#------------------------------------------------------------ 
#------------------------------------------------------------ 
#------------------------------------------------------------ 

# MothurPR2 v5.0.0 --- Prepare a protist-specific database for the 66F_1200R primer pair (634bp long). We will be again using the Rescript plug-in. Grab the files here: https://github.com/pr2database/pr2database/releases

# Step 1. Create qza files for qiime Rescript
gzip -d ../db/mothurPR2/*.gz
    
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "./db/mothurPR2/pr2_version_5.0.0_SSU_mothur.fasta" \
    --output-path "./db/mothurPR2/pr2_v5.0.0_dna.qza"

cd ./db/mothurPR2
nano pr2_version_4.14.0_SSU_mothur.tax # added two header values: 'Feature ID', 'Taxon'. Saved over
cd ../..
   
qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-path "./db/mothurPR2/pr2_version_5.0.0_SSU_mothur.tax" \
    --output-path "./db/mothurPR2/pr2_v5.0.0_tax.qza"
    #--input-format HeaderlessTSVTaxonomyFormat
    
# Step 2. Clean db sequences
# Remove sequences with >5 ambiguous bases and >8 homopolymer regions 
qiime rescript cull-seqs \
    --i-sequences "./db/mothurPR2/pr2_v5.0.0_dna.qza" \
    --o-clean-sequences "./db/mothurPR2/pr2_v5.0.0_dna_clean.qza"

# Step 3. Filter db sequences
# Silva contains archae, bacteria, and eukaryotes. Naturally, your seawater ASVs will contain all 
# We don't want to ignore archae and special euk ASVs even though we are doing bacteria
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences "./db/mothurPR2/pr2_v5.0.0_dna_clean.qza"\
    --i-taxonomy "./db/mothurPR2/pr2_v5.0.0_tax.qza" \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs "./db/mothurPR2/pr2_v5.0.0_dna_filter.qza" \
    --o-discarded-seqs "./db/mothurPR2/pr2_v5.0.0_dna_discarded.qza"
    # We will push discards into their own file

# Step 4. Dereplicate db sequences
qiime rescript dereplicate \
    --i-sequences "./db/mothurPR2/pr2_v5.0.0_dna_filter.qza"  \
    --i-taxa "./db/mothurPR2/pr2_v5.0.0_tax.qza"  \
    --p-rank-handles 'disable' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "./db/mothurPR2/pr2_v5.0.0_dna_derep1.qza" \
    --o-dereplicated-taxa "./db/mothurPR2/pr2_v5.0.0_tax_derep1.qza"

# Step 5. From your sequences, pull out primer-specific reads. Dereplicate again
qiime feature-classifier extract-reads \
    --i-sequences "./db/mothurPR2/pr2_v5.0.0_dna_derep1.qza" \
    --p-f-primer CAGCAGCCGCGGTAATTCC \
    --p-r-primer CCCGTGTTGAGTCAAATTAAGC \
    --p-n-jobs 16 \
    --p-read-orientation 'forward' \
    --o-reads "./db/mothurPR2/pr2_v5.0.0_566F.1200R_seqs.qza"

qiime rescript dereplicate \
    --i-sequences "./db/mothurPR2/pr2_v5.0.0_566F.1200R_seqs.qza" \
    --i-taxa "./db/mothurPR2/pr2_v5.0.0_tax_derep1.qza" \
    --p-rank-handles 'disable' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "./db/mothurPR2/pr2_v5.0.0_566F.1200R_seqs_derep.qza" \
    --o-dereplicated-taxa  "./db/mothurPR2/pr2_v5.0.0_566F.1200R_tax_derep.qza"

# Step 6. Create your final product before aligning to your representatives... the classifier file! 
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "./db/mothurPR2/pr2_v5.0.0_566F.1200R_seqs_derep.qza" \
    --i-reference-taxonomy "./db/mothurPR2/pr2_v5.0.0_566F.1200R_tax_derep.qza" \
    --o-classifier "./db/mothurPR2/pr2Classifier_566F.1200R.qza"
    
#------------------------------------------------------------ 
#------------------------------------------------------------ 
#------------------------------------------------------------ 

# Marferret Diatom Update, V1.0.0 -- Sacha Coesel 
# This code takes updates taxonomy classifications for Bacilliophyta in the PR2 artifact table above.
# Here, we generate a Marferret db then merge to the PR2 v5.0.0. Updated classifications will overwrite outdated PR2 sequences

# Sacha's original fasta contains lower case DNA seqs. The file was modified... 
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' marferret_v1_diatom_18S.fasta > marferret_v1_diatom_18S_rsk.fasta
# Refer to 'marfarret.taxFile_make.R' for further modifications. Final input before import is '...rsk2.fasta' and '..rsk.tax' files

# Import info
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "./db/marferret/marferret_v1_diatom_18S_rsk2.fasta" \
    --output-path "./db/marferret/marferret_v1_dna.qza"

cd ./db/marferret
nano marferret_v1_diatom_18S_rsk.tax # added two header values: 'Feature ID', 'Taxon'. Saved over
cd ../..

qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-path "./db/marferret/marferret_v1_diatom_18S_rsk.tax" \
    --output-path "./db/marferret/marferret_v1_tax.qza"

# Step 2. Clean db sequences
# Remove sequences with >5 ambiguous bases and >8 homopolymer regions 
qiime rescript cull-seqs \
    --i-sequences "./db/marferret/marferret_v1_dna.qza" \
    --o-clean-sequences "./db/marferret/marferret_v1_dna_clean.qza"

# Step 3. Filter db sequences
# Silva contains archae, bacteria, and eukaryotes. Naturally, your seawater ASVs will contain all 
# We don't want to ignore archae and special euk ASVs even though we are doing bacteria
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences "./db/marferret/marferret_v1_dna_clean.qza"\
    --i-taxonomy "./db/marferret/marferret_v1_tax.qza" \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs "./db/marferret/marferret_v1_dna_filter.qza" \
    --o-discarded-seqs "./db/marferret/marferret_v1_dna_discarded.qza"
    # We will push discards into their own file

# Step 4. Dereplicate db sequences
qiime rescript dereplicate \
    --i-sequences "./db/marferret/marferret_v1_dna_filter.qza"  \
    --i-taxa "./db/marferret/marferret_v1_tax.qza"   \
    --p-rank-handles 'disable' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "./db/marferret/marferret_v1_dna_derep1.qza" \
    --o-dereplicated-taxa "./db/marferret/marferret_v1_tax_derep1.qza"

# Step 5. From your sequences, pull out primer-specific reads. Dereplicate again
qiime feature-classifier extract-reads \
    --i-sequences "./db/marferret/marferret_v1_dna_derep1.qza" \
    --p-f-primer CAGCAGCCGCGGTAATTCC \
    --p-r-primer CCCGTGTTGAGTCAAATTAAGC \
    --p-n-jobs 16 \
    --p-read-orientation 'forward' \
    --o-reads "./db/marferret/marferret_v1_566F.1200R_seqs.qza"

qiime rescript dereplicate \
    --i-sequences "./db/marferret/marferret_v1_566F.1200R_seqs.qza" \
    --i-taxa "./db/marferret/marferret_v1_tax_derep1.qza" \
    --p-rank-handles 'disable' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "./db/marferret/marferret_v1_566F.1200R_seqs_derep.qza" \
    --o-dereplicated-taxa  "./db/marferret/marferret_v1_566F.1200R_tax_derep.qza"

# Step 6. Create your final product before aligning to your representatives... the classifier file! 
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "./db/marferret/marferret_v1_566F.1200R_seqs_derep.qza" \
    --i-reference-taxonomy "./db/marferret/marferret_v1_566F.1200R_tax_derep.qza" \
    --o-classifier "./db/marferret/marferretClassifier_566F.1200R.qza"
    
# Step 7. Merge reads and taxonomies from two databases
qiime feature-table merge-seqs \
  --i-data "./db/mothurPR2/pr2_v5.0.0_566F.1200R_seqs_derep.qza" \
  --i-data "./db/marferret/marferret_v1_566F.1200R_seqs_derep.qza" \
  --o-merged-data "./db/marferret/merged-reads.qza"
  
qiime feature-table merge-taxa \
  --i-data "./db/mothurPR2/pr2_v5.0.0_566F.1200R_tax_derep.qza" \
  --i-data "./db/marferret/marferret_v1_566F.1200R_tax_derep.qza" \
  --o-merged-data "./db/marferret/merged-taxonomy.qza"
  
# Step 8. Train a new classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "./db/marferret/merged-reads.qza" \
  --i-reference-taxonomy "./db/marferret/merged-taxonomy.qza" \
  --o-classifier "./db/marferret/PR2MarClassifier_566F.1200R.qza"

#------------------------------------------------------------ 
#------------------------------------------------------------ 
#------------------------------------------------------------ 
# NCBI --- Prepare a current NCBI database, specific for our primer pair.

# Step 1. Choose and call a current 18S database from NCBI. TaxIDs were derived from the NCBI Taxonomy Page and capture main groups. (I recommend performing this task during working hours, I found I got plugin errors due to bad HTTP requests) After getting each chunk, merge all tax and seqs into master qza files and convert cDNA to DNA

nano 

# copy, paste all below
taxID_list="554915 2686027 554296 1401294 2608240 3027 2611352 38254 2608109 \
            2611341 33154 2763 2698737 33090 2683617 42452 61964"
            

for id in $taxID_list
    do 
	command="qiime rescript get-ncbi-data \
        --p-query 'txid$id[ORGN] AND (18S[TITLE] OR small ribosomal subunit[TITLE] OR \
                                      SSU[TITLE) OR 18S ribosomal rna[TITLE] or \
                                      18S ribosomal RNA[TITLE] OR 18S rRNA biogenesis protein[TITLE])' \
        --p-n-jobs 8 \
        --o-sequences './db/ncbi/ncbi_18S_rna_txid$id.qza' \
        --o-taxonomy './db/ncbi/ncbi_18S_tax_txid$id.qza'"
                echo "Processing : $id"
                echo $command
                eval $command
    done;
# save as getDB_ncbi.sh

bash getDB_ncbi.sh
 	
ncbi_18S_tax_txid1401294.qza 
ncbi_18S_tax_txid1401294.qza


dir="./db/ncbi/ncbi_18S_"    
qiime feature-table merge-taxa \
	--i-data ${dir}tax_txid1401294.qza \
	         ${dir}tax_txid2608109.qza \
	         ${dir}tax_txid2608240.qza \
	         ${dir}tax_txid2611341.qza \
	         ${dir}tax_txid2611352.qza \
	         ${dir}tax_txid2683617.qza \
	         ${dir}tax_txid2686027.qza \
	         ${dir}tax_txid2698737.qza \
	         ${dir}tax_txid2763.qza \
	         ${dir}tax_txid3027.qza \
	         ${dir}tax_txid38254.qza \
	         ${dir}tax_txid42452.qza \
	         ${dir}tax_txid554915.qza \
	         ${dir}tax_txid33154.qza \
	         ${dir}tax_txid33090.qza \
	         ${dir}tax_txid554296.qza \
	         ${dir}tax_txid61964.qza \
	--o-merged-data ${dir}tax.qza 
	
qiime feature-table merge-seqs \
	--i-data ${dir}rna_txid1401294.qza \
	         ${dir}rna_txid2608109.qza \
	         ${dir}rna_txid2608240.qza \
	         ${dir}rna_txid2611341.qza \
	         ${dir}rna_txid2611352.qza \
	         ${dir}rna_txid2683617.qza \
	         ${dir}rna_txid2686027.qza \
	         ${dir}rna_txid2698737.qza \
	--o-merged-data ${dir}rna.qza 

qiime feature-table merge-seqs \
	--i-data ${dir}rna_txid2763.qza \
	         ${dir}rna_txid3027.qza \
	         ${dir}rna_txid38254.qza \
	         ${dir}rna_txid554296.qza \
	         ${dir}rna_txid554915.qza \
	         ${dir}rna_txid33154.qza \
	         ${dir}rna_txid33090.qza \
	         ${dir}rna_txid554296.qza \
	         ${dir}rna_txid61964.qza \
	         ${dir}rna.qza \
	--o-merged-data ${dir}rna.qza         
	

#qiime rescript reverse-transcribe \
#	--i-rna-sequences "${dir}rna.qza" \
#	--o-dna-sequences "${dir}dna.qza" ---------- already DNA seqs!

# Step 2. Clean db sequences. Remove sequences with >5 ambiguous bases and >8 homopolymer regions 
qiime rescript cull-seqs \
	--i-sequences "${dir}rna.qza" \
	--o-clean-sequences "${dir}dna_clean.qza"

# Step 3. Filter db sequences. I've included Archaea just in case here. I don't expect them to be present but wish to include
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences "${dir}dna_clean.qza" \
    --i-taxonomy "${dir}tax.qza"  \
    --p-labels Archaea Eukaryota \
    --p-min-lens 900 1400 \
    --o-filtered-seqs "${dir}dna_filter.qza" \
    --o-discarded-seqs "${dir}dna_discard.qza"

# Step 4. Dereplicate db sequences
qiime rescript dereplicate \
    --i-sequences "${dir}dna_filter.qza"  \
    --i-taxa "${dir}tax.qza" \
    --p-rank-handles 'disable' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "${dir}dna_derep1.qza" \
    --o-dereplicated-taxa "${dir}tax_derep1.qza"

# Step 5. From your sequences, pull out primer-specific reads. Dereplicate again
qiime feature-classifier extract-reads \
    --i-sequences "${dir}dna_derep1.qza" \
    --p-f-primer CAGCAGCCGCGGTAATTCC \
    --p-r-primer CCCGTGTTGAGTCAAATTAAGC \
    --p-n-jobs 16 \
    --p-read-orientation 'forward' \
    --o-reads "${dir}dna_566F.1200R_seqs.qza"

qiime rescript dereplicate \
    --i-sequences "${dir}dna_566F.1200R_seqs.qza" \
    --i-taxa "${dir}tax_derep1.qza" \
    --p-rank-handles 'disable' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "${dir}dna_566F.1200R_seqs_derep.qza" \
    --o-dereplicated-taxa  "${dir}tax_derep.qza"

# Step 6. Create your final product before aligning to your representatives... the classifier file! 
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "${dir}dna_566F.1200R_seqs_derep.qza" \
    --i-reference-taxonomy "${dir}tax_derep.qza" \
    --o-classifier "./db/ncbi/ncbiClassifier_18S_566F.1200R.qza"

#------------------------------------------------------------ 
#------------------------------------------------------------ 
#------------------------------------------------------------ 

# CyanoSeq_silva138.1 v1.1.0 --- Prepare a protist-specific database for the 66F_1200R primer pair (634bp long). We will be again using the Rescript plug-in. Grab the files here: https://github.com/pr2database/pr2database/releases

















