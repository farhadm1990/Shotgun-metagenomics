# An R workflow for wholegenome sequence (WGS) data analysis
Unlike 16S rRNA amplicons, shotgun metagenomics targets all DNA present in the sample, e.g. colon. This means your samples will contain DNA from bacteria, host, archeae, and DNA-virum. Therefore, in the first step the host DNA must be removed if it is not of your interest. After decontamination, short reads will be assembled to form Metagnomics Assembled Genomes (MAGs) or contigs. For taxonomic annotations MAGs will be binned based on neucleotide identity (NI) threshold and will be blasted against the database. 
All these steps were done using [ATLAS](https://github.com/metagenome-atlas/atlas) Snakmake workflow and the resultant was analysed as demostrated in this R markdown. 
