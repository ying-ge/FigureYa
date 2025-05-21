This folder contains information about the results of the ANOVA-based identification of differentially
regulated processes with age, differentially expressed genes with age and the disease alignment analysis. 

Information on differential process activity with age across species is contained in the "speciesComparison" subdirectory.
The name of each ontology is followed by the time points
that are compared: 1 - young, 2 - mature_1, 3 - mature_2, 4 - old_1 and 5 - old_2. The analysis between young
and the two old age groups is contained in the folders "<ontology>-1-4+5".

In each of these folders you find up to eleven files that contain results for the individual comparisons as follows:
"all.csv" - age-associated transcriptional changes across all species
"blood.csv" - age-associated transcriptional changes across all blood data
"brain.csv" - age-associated transcriptional changes across all brain data
"liver.csv" - age-associated transcriptional changes across all liver data
"skin.csv" - age-associated transcriptional changes across all skin data
"drerio.csv" - age-associated transcriptional changes for Danio rerio
"hsapiens.csv" - age-associated transcriptional changes for human
"mmusculus.csv" - age-associated transcriptional changes for mouse
"nfurzeri.csv" - age-associated transcriptional changes for Nothobranchius furzeri
"fibroblasts.csv" - age-associated transcriptional changes for fibroblast cell cultures
"resample.csv" - lists for each process how often it was found in 1000 repetitions where samples were randomly removed
                 until only a single sample for each individual remained

The folder "SpeciesComparison-All-Ages" contains results for the comparison across all age groups for the individual 
ontologies.
 
Data within the files are delimited by "$". The files are organized as follows:
1st column - process ID (only for gene ontology, otherwise process name)
2nd column - process name
3rd/4th column - direction of the age-associated change
5th column - factors that remained after model-reduction in the ANOVA
6th column - p-value of the ANOVA of the age-associated transcriptional change for this process (not corrected for multiple testing)
8th column - FDR-adjusted p-value
9th column (only in "all.csv") - whether there is a regulation across species following the criteria in the 
                                 main manuscript (significant regulation in at least one fish and one mammal)

The same analysis for a different scheme of deriving process activity whereby
gene expression values are rank normalized for each tissue and species before
process activity is determined is contained in the folders ending with 
"gene-normalization".

The folder "DEGs" contains lists of differentially expressed genes for the ageing analysis.
Differentially expressed genes for each tissue between the different age states following the same naming conventions
as above ('1'=young, ...) are provided in the corresponding sub-folders. Also the differentially expressed genes 
derived from differentially expressed processes are contained in this folder. 

The folder "DiseaseAlignment" contains information about AMDA scores for each disease and 
ageing data set (the score themselves and their p-values) and DAC scorese for each ageing-regulated process and disease. 
Moreover, the folder contains normalized fold-changes between cases and controls on the gene as well as
process level. "DiseaseAlignment-MiddleAge" contains the corresponding information for the respective analysis
during middle age.
