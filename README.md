# An R workflow for wholegenome sequence (WGS) data analysis
Unlike 16S rRNA amplicons, shotgun metagenomics targets all DNA present in the sample, e.g. colon. This means your samples will contain DNA from bacteria, host, archeae, and DNA-virum. Therefore, in the first step the host DNA must be removed if it is not of your interest. After decontamination, short reads will be assembled to form Metagnomics Assembled Genomes (MAGs) or contigs. For taxonomic annotations MAGs were binned based on neucleotide identity (NI) threshold and will be blasted against the database. 
All these steps were done using [ATLAS](https://github.com/metagenome-atlas/atlas) Snakmake workflow and the resultant was analysed as demostrated in this [R markdown](https://github.com/farhadm1990/Shotgun-metagenomics/blob/main/Megagenomics%20analysis%20of%20DSS%20and%20B.%20pilosicoli.Rmd). 

## 2.2. Removing singletones based on abundance
```{r}
#A function to find singletones. You need to be careful about this step!
out.ASV = function(phyloseq, threshold =1, binwidth = 0.01) {
  
#Loading necessary pkgs      
  pacman::p_load(glue, tidyverse, reshape2, ggrepel, S4Vectors) # nolint
#This function requires phyloseq, tidyverse and glue packages to be loaded. 
    if (sum(colSums(otu_table(phyloseq)))/ncol(otu_table(phyloseq)) == 100 ) {#making the relative abundance table
                    rel_abund = as(t(otu_table(phyloseq)), "matrix")
    } else if (sum(colSums(otu_table(phyloseq)))/ncol(otu_table(phyloseq)) == 1) {
                    rel_abund = as(t(otu_table(phyloseq)), "matrix")
                    } else {
                    rel_abund = as(t(apply(otu_table(phyloseq), 
                    ifelse(taxa_are_rows(phyloseq), 1,2), 
                    function(x) x/sum(x))), "matrix")  
                    } 
                      
                      
                      names.single = apply(rel_abund, 1, function(x){ifelse(x == threshold, TRUE, ifelse(x == sum(x),
                      TRUE, FALSE))}) %>% reshape2::melt() %>% filter(value == TRUE) %>% dplyr::select(2) %>%
                      pull   %>% as.vector()
                      
                        
                        if (length(names.single) == 0 ) {
                        print(glue("WOW! {length(names.single)} singletones detected in this dataset"))
                        qplot.noSing = qplot(rel_abund, geom = "histogram", binwidth = binwidth, 
                        show.legend = F, main = "Frequency count of relative abundance, no singletones detected") +
                        xlab ("Relative abundance in samples") + ylab("Frequency") + theme_bw()
                            
                        
                        return(structure(list(qplot.noSing)))
                            
                        } else { 
                             
                       single.ASV = rel_abund[rownames(rel_abund) %in% names.single,]
                       single.ASV[single.ASV == 0] <- NA # A separate dataset for annotation of singletones on the barplot
                            
                       qplot.withSing = qplot(rel_abund, geom = "histogram", binwidth = binwidth, 
                       main = "Frequency count of relative abundance with singletones") +
                       geom_bar(aes(single.ASV), fill = "red",  color = NA, width = binwidth)+
                       xlab ("Relative abundance in samples") + ylab("Frequency") + 
                       geom_label_repel(aes(x = 1, y =length(rel_abund)/5), 
                       label.padding =  unit(0.55, "lines"), 
                       label = glue("{length(names.single)}\n Singletones"), color = "black") + theme_bw()
                            
                       qplot.rmSing = qplot(rel_abund[!rownames(rel_abund) %in% names.single, ], geom = "histogram",
                       binwidth = binwidth, main = "Frequency count of relative abundance without singletones") +
                       xlab ("Relative abundance in samples") + ylab("Frequency")+ theme_bw()
                            
                       print(glue('Oh no..! {length(names.single)} singletones detected in the dataset'))
                       return(structure(list(qplot.withSing, qplot.rmSing, unlist(names.single))) )
                    
                        }                        
    
                             
        }
                        
single.test = out.ASV(phyloseq = pst.count, threshold = 1, binwidth = 0.1)
#singletones = single.test[[3]] #here you can extract the names of the singletones
single.test[[1]]#to show the plot with singletones
#single.test[[2]]#to show the plot without singletones
#Now you can remove the singletones from your pst file as follows:
#pst.no.single = subset_taxa(ps, !taxa_names(ps)%in% singletones)
#ps = pst.no.single
rm(single.test)
```

![phylogenetic tree of species](https://github.com/farhadm1990/Shotgun-metagenomics/blob/main/Tree.species.circular.jpeg)
