This README.txt file was generated on 2022-02-23 by Julio Cesar Cortes


GENERAL INFORMATION

1. Title of Dataset: Battenberg reference (hg38)

2. Author Information

	A. Principal Investigator Contact Information

		Name: Prof. David C. Wedge

		Institution: University of Manchester

		Address: The University of Manchester - Manchester Cancer Research Centre (MCRC)

		Email: david.wedge@manchester.ac.uk

	B. Alternate Contact Information

		Name: Julio Cesar Cortes Rios

		Institution: University of Manchester

		Address: The University of Manchester - Wolfson Molecular Imaging Centre (WMIC)

		Email: juliocesar.cortesrios@manchester.ac.uk

3. General Description: Reference files for hg38 to be used with the Battenberg package as available in: https://github.com/Wedge-lab/battenberg/


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: None

2. Links to publications that cite or use the data: 
  The Life history of 21 Breast Cancers, Nik-Zainal et al, Cell, 2012; PMID 22608083; https://www.cell.com/fulltext/S0092-8674(12)00527-2 

3. Links to other publicly accessible locations of the data: None

4. Links/relationships to ancillary data sets: 
  Previous incomplete hg38 reference set: https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52
  Previous hg19 reference set: https://ora.ox.ac.uk/objects/uuid:2c1fec09-a504-49ab-9ce9-3f17bac531bc


5. Was data derived from another source? yes

	A. If yes, list source(s): 
           https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download
           http://faculty.washington.edu/browning/beagle/beagle.html#download
           https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download


6. Recommended citation for this dataset: 
  The Life history of 21 Breast Cancers, Nik-Zainal et al, Cell, 2012; PMID 22608083; https://www.cell.com/fulltext/S0092-8674(12)00527-2


DATA & FILE OVERVIEW

1. File List: 

	- 1000G_loci_hg38_chr.zip - ZIP file containing the 24 files with the 1000 Genomes SNP alleles data (with "chr" prefix) and 48 files with the 1000 Genomes SNP loci data (with "chr" prefix) for each chromosome (including soft/hard link for chr23) - 276 Mb
	- 1000G_loci_hg38.zip - ZIP file containing the 24 files with the 1000 Genomes SNP alleles data (without "chr" prefix) and 48 files with the 1000 Genomes SNP loci data (without "chr" prefix) for each chromosome (including soft/hard link for chr23) - 275 Mb
	- beagle_chr.zip - ZIP file containing the JAR file for Beagle 5.3, the 23 files for the phasing information for Beagle and for each chromosome (with "chr" prefix), and the 25 files with the mapping information for Beagle (including PAR/nonPAR regions for chrX) - 7.6 Gb
	- beagle.zip - ZIP file containing the JAR file for Beagle 5.3, the 23 files for the phasing reference information for Beagle and for each chromosome (without "chr" prefix), and the 25 files with the plink mapping information for Beagle (including PAR/nonPAR regions for chrX) - 7.6 Gb
	- GC_correction_hg38_chr.zip  - ZIP file containing the 24 files compressed files, one for each chromosome (with "chr" prefix and including a hard link for chr23), with GC correction content - 1.2 Gb
	- GC_correction_hg38.zip - ZIP file containing the 24 files compressed files, one for each chromosome (without "chr" prefix and including a soft link for chr23), with GC correction content - 1.3 Gb
	- imputation_chr.zip - ZIP file containing the 25 files (including PAR/nonPAR regions for chrX) with the mapping reference information for Impute and for each chromosome (with "chr" prefix) - 54 Mb
	- imputation.zip - ZIP file containing the 25 files (including PAR/nonPAR regions for chrX) with the mapping reference information for Impute and for each chromosome (without "chr" prefix) - 54 Mb
	- probloci_chr.zip - ZIP file containing one file with the problem loci data that contains SNP loci that should be filtered out (with "chr" prefix) - 29 Mb
	- probloci.zip - ZIP file containing one file with the problem loci data that contains SNP loci that should be filtered out (without "chr" prefix) - 29 Mb
	- README.txt - This file, with information about the dataset and instructions on how to use it - 7.5 Kb
	- sha1sums.txt - Text file containing the SHA-1 hash values for checksum verification (integrity) of the 12 ZIP files contained in this dataset - 771 b
	- shapeit2_chr.zip - ZIP file containing the 25 files, one for each chromosome (with "chr" prefix and including PAR/nonPAR regions for chrX), with the haplotype estimation provided by SHAPEIT, and 25 files with the associated legends, required by Impute - 6.4 Gb
    - shapeit2.zip - ZIP file containing the 25 files, one for each chromosome (without "chr" prefix and including PAR/nonPAR regions for chrX), with the haplotype estimation provided by SHAPEIT, and 25 files with the associated legends, required by Impute - 6.4 Gb

2. Relationship between files, if important: The files with suffix "_chr.zip" contain the same information as the other zip files but adding the string "chr" to the label of the chromosomes

3. Additional related data collected that was not included in the current data package: None

4. Are there multiple versions of the dataset? yes

	A. If yes, name of file(s) that was updated: beagle.zip, *_chr.zip, README.txt

		i. Why was the file updated? beagle.zip contains now the required *.bref3 and *.map to invoke Beagle for the phasing stage. *_chr.zip files were included to allow the execution of Battenberg with a chromosome nomenclature that prefixes the chromosome label with the string "chr". The README.txt file was updated to include more information about the dataset and clearer instructions on how to use it.

		ii. When was the file updated? beagle.zip and *_chr.zip files were updated on 23/02/2022. The README.txt file was updated on 23/02/2022


INSTRUCTIONS TO USE THE REFERENCE FILES

1. A dedicated directory shoud be created to download there all the files in this dataset, e.g., reference_hg38

2. Download the *.zip files to the dedicated directory created in step 1.

3. Download the sha1sums.txt file to the same folder and execute a checksum, to verify the integrity of the files downloaded, with the following command (Linux): sha1sum -c sha1sums.txt 

4. Once verified, all the zip files must be uncompressed (using UNZIP or other tool that supports the ZIP file format) into the directory created in step 1.

5. Once step 4. is completed the zip files can be deleted from the dedicated folder

6. The impute_info.txt file must be edited, replacing the string "${REF_PATH}" with the full path to the files contained in the imputation and shapeit2 folders

7. If needed, create the following symbolic links to allow Battenberg to use chr23 instead of chrX (from the main dedicated folder):
        
 - ln -s 1kg.phase3.v5a_GRCh38nounref_allele_index_chrX.txt 1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_allele_index_chr23.txt
 - ln -s 1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chrX.txt 1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr23.txt
 - ln -s 1kg.phase3.v5a_GRCh38nounref_loci_chrX.txt 1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_chr23.txt
 - ln -s 1000G_GC_chrX.txt.gz GC_correction_hg38/1000G_GC_chr23.txt.gz
 - ln -s 1000G_RT_chrX.txt.gz RT_correction_hg38/1000G_RT_chr23.txt.gz
 - ln -s chrX.1kg.phase3.v5a.b37.bref3 beagle/chr23.1kg.phase3.v5a.b37.bref3
 - ln -s chrX.1kg.phase3.v5a_GRCh38nounref.vcf.gz beagle/chr23.1kg.phase3.v5a_GRCh38nounref.vcf.gz
 - ln -s plink.chrX.GRCh38.map beagle/plink.chr23.GRCh38.map
 - ln -s plink.chrX_par1.GRCh38.map beagle/plink.chr23_par1.GRCh38.map
 - ln -s plink.chrX_par2.GRCh38.map beagle/plink.chr23_par2.GRCh38.map
