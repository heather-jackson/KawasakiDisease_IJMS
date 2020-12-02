# KawasakiDisease_IJMS
Code to accompany the IJMS Kawasaki Disease patient stratification paper analysis.

Functions are found in <i>functions.R</i>
<i>analysis.R</i> contains code to read in data files and execute functions to produce the main figures. 

The code and the data presented here are for the transcriptomic arm of the study. The code can easily be modified for proteomic data as all steps of the analytical pipeline are identical, except for the inclusion of immune cell proportions. 

<i>discovery_transcriptomics.RData</i> contains the transcriptomic expression data for the discovery dataset and the accompanying metadata. 
<i>validation_transcriptomics.RData</i> contains the transcriptomic expression data for the validation dataset and the accompanying metadata. 

Please note, these are normalised versions of these datasets. The normalisation processes are described in the method section of the paper. In brief, the data was normalised using RSN normalisation from the Lumi R package (Du et al., 2008). 

The raw versions of the datasets can be found at their accompanying GEO entry. Furthermore, GeneSpring normalised versions of the datasets are also found there. 

<b>Discovery transcriptomic dataset:</b> GSE73461

<b>Validation transcriptomic dataset:</b> GSE73462, GSE73463

<i>References:</i> 

Du, P., Kibbe, W.A., Lin, S.M. (2008). “lumi: a pipeline for processing Illumina microarray.” Bioinformatics.
