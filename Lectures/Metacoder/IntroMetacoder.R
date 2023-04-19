###########METACODER SCRIPT FOR BIOMERG9905

# Install metacoder and other packages

#install.packages("metacoder")
#install.packages("taxa")
#install.packages("ggplot2")
#install.packages("readr")
#install.packages("dplyr")
#install.packages("tidyverse")

#you can also install the development version for the newest functions
#install.packages("devtools")
#devtools::install_github("grunwaldlab/metacoder")

#Read libraries
library(metacoder)
library(ggplot2)
library(taxa)
library(readr)
library(dplyr)
library(tidyverse)



# Read your files, containing OTU table and Metadata
#make sure that the taxonomy is attached as a column in your OTU table

otu_data<-read_tsv("rarotutable_small.txt")
print(otu_data)
View(otu_data)

sample_data<- read_tsv("sampledata.txt")
                    
print(sample_data)

#check that the taxonomy is attached at the end
head(otu_data$taxonomy, 10)

                 

#Parsing taxa (this sorts your taxa using the package taxa into a 
#taxmap object)
obj <- parse_tax_data(otu_data,
                      class_cols = "taxonomy",
                      class_sep = ";", # What each taxon is separated by
                      class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                      class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))

#View your taxmap object:

print(obj)

#Check if it returns ok taxon names (e.g. "Fungi", "Basidiomycota")
head(taxon_names(obj)) 

#rename the tax_data
names(obj$data) <- "read_abundance"
#to view the entire object
print(obj)


#to make heat trees we need to calculate taxon abundance
obj$data$tax_abund <- calc_taxon_abund(obj, "read_abundance",
                                       cols = sample_data$SampleID)

#For calculating total read abundance per taxa - this gives you a list of the
#sum of all taxa across all samples
obj$data$tax_abund$total <- rowSums(obj$data$tax_abund[, -1]) 
# -1 = taxon_id column

#The dataset may contain uwanted names/taxa, these can be filtered out
#For removing names=unidentified

taxa_to_remove <- c("unidentified")
obj <-filter_taxa(obj, taxon_names %in% taxa_to_remove, invert = TRUE, subtaxa = FALSE, supertaxa = FALSE)

#Intitial tree, no layout specifications (this may take a while to plot):
heat_tree(obj, node_label = taxon_names, node_size = total, 
          node_color = total)


#Heat tree with more appropriate settings for getting an overview:
#Here we have set taxon_ranks="o", which is order level, 
#we clean up the names, and we use 
#the relative read abundance (total) to make the edges and nodes
obj %>%
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree(node_label = taxon_names,
            node_size = total,
            node_color = total,
            node_color_range = c("gray",  "#85CEBF","#018571", "#A78959", "#A6611A"),
            #node_color_trans= "ln area",
            layout = "da", initial_layout = "re")

#if you want to change layout, check 
#?heat_tree()


######## Plotting heat trees for different sample types:

#In this dataset we have samples from roots and soil.
#Let's plot these individually:

view(sample_data) 
sample_data$Sample
#Sample types are Bistorta, Potentilla and Soil

#Calculate read abundances for the different sample types 
#Here you have to specify groups= 

obj$data$type_abund <- calc_taxon_abund(obj, "read_abundance",
                                        cols = sample_data$SampleID,
                                        groups = sample_data$Sample)

#Make a heat tree with only taxa present in Bistorta vivipara roots
# filter_taxa is used to sort out only taxa present in Bistorta with at least 1 read:

obj %>%
  filter_taxa(Bistorta>0) %>% #if you wish to sort out very rare OTUs (i.e. less than 50 reads), you can set (Bistorta>50)
  # taxa:: needed because of phyloseq::filter_taxa
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% 
  # remove "odd" taxa
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE) %>% 
  # to change taxonomic resolutions, 
  #set taxon_ranks to "k","p","c","o","g","s"
  heat_tree(node_label = taxon_names,
            node_size = Bistorta, 
            node_color = Bistorta, 
            node_color_range= quantative_palette(),
            #node_color_range=diverging_palette(),
            #node_color_range = c("gray",  "pink",  "blue"),
            #node_color_trans= "area",
            layout = "da", initial_layout = "re", 
            title = "Taxa in Bistorta vivipara")


#For different colour scemes, it might be a good idea to have a look a RcolorBrewer
#install.packages("RColorBrewer")
library(RColorBrewer)
display.brewer.all()
Greens<-brewer.pal(4,"Greens")
Greens
Greens<-as.list(brewer.pal(4,"Greens"))
Greens
#Plotting the taxa present with more than 50 reads in Potentilla erecta

obj %>%
  filter_taxa(Potentilla>50) %>% # if you wish to sort out very rare OTUs (i.e. less than 50 reads)
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "g", supertaxa = TRUE) %>% # to change taxonomic resolutions, 
  #set taxon_ranks to "k","p","c","o","g","s"
  heat_tree(node_label = taxon_names,
            node_size = Potentilla, 
            node_color = Potentilla, 
            node_color_range=Greens,
            layout = "da", initial_layout = "re", 
            title = "Taxa in Potentilla erecta")


obj %>%
  filter_taxa(Soil>50) %>% #Filtering out everything below 50 reads
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE) %>% # to change taxonomic resolutions, 
  #set taxon_ranks to "k","p","c","o","g","s"
  heat_tree(node_label = taxon_names,
            node_size = Soil, 
            node_color = Soil, 
            node_color_range = c("gray",  "#85CEBF","#018571", "#A78959", "#A6611A"),
            #node_color_trans= "area",
            layout = "da", initial_layout = "re", 
            title = "Taxa in soil")


##############Compare two groups######################
#Here we compare taxa present in Plant roots and soil using Wilcoxon Rank Sum 
#test and the function compare_groups(). 
#Columns should be specified to SampleID and groups must be specified.
#(Ideally the number of samples that are compared should be equal, but there are 
#twice as many plantroots as soil samples in this reduced dataset)

obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols = sample_data$SampleID,
                                      groups = sample_data$Type)

View(obj$data$diff_table)
range(obj$data$diff_table$log2_median_ratio)

#We then need to correct for multiple comparison, in this case we use "false discovery rate"
#but others can be specified, see:
#?p.adjust()

obj <- mutate_obs(obj, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

#then set the non-significant p-values to zero to aid visualization

obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0

#Plotting the tree. Change taxon_ranks=="", for "c" (class),"f" (family), "o" (order), "g"(genus)  "s" (species)
set.seed(2)
obj %>%
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$")) %>%
  heat_tree(node_label = cleaned_names,
            node_size = total, # number of OTUs
            node_color = log2_median_ratio, # difference between groups
            node_color_interval = c(-3, 3), # symmetric interval to make zero (midpoint) grey (next step)
            node_color_range = c("#D55E00", "gray", "#009E73"), # should use diverging colors
            node_size_axis_label = "read abundance",
            node_color_axis_label = "Log 2 ratio of median counts",
            layout = "da", initial_layout = "re", # good layout for large trees
            title = "Plant roots vs soil samples")




########Comparing more than two treatments###################

#Here we compare three different sample types: soil, Potentilla ercta roots and Bistorta vivipara
#roots
#(This is the same tax_abund calculation as earlier in the script, 
#somehow it gets removed and has to be redone after the previous 
#comparison.)
obj$data$tax_abund <- calc_taxon_abund(obj, "read_abundance",
                                       cols = sample_data$SampleID)
obj$data$tax_abund$total <- rowSums(obj$data$tax_abund[, -1])


#Make a diff_table using Wilcoxon Rank Sum test, specifiying groups = 
#This may give you warnings, if there are some taxa that don't overlap between groups. You can ignore this


obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols = sample_data$SampleID,
                                      groups = sample_data$Sample)



#Look at the diff_table
#View(obj$data$diff_table)


#We then need to correct for multiple comparisons and 
#set non-significant differences to zero as above

obj <- mutate_obs(obj, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0

view(obj$data$diff_table)

#Plotting the tree. 
#Change taxon_ranks=="", for "c" (class),"f" (famili), "o" (order), "g"(genus)  "s" (species)
obj %>%
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,
                   node_size = total, # read abundance
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "read abundance",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)


#Let's try on genus level:
 
obj %>%
  filter_taxa(taxon_ranks == "g", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,
                   node_size = total, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "read abundance",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)



################ Downstream analyses #######################
###Only to show you some examples on how the data you now have can easily be used 
###in downstream analyses#####



###rarefaction curves and species accumulation curves####

library(vegan)
library(agricolae)

obj$data$taxon_abund <- calc_taxon_abund(obj, data = "tax_abund", cols= sample_data$SampleID)
obj$data$taxon_abund$total <- rowSums(obj$data$taxon_abund[, -1]) # -1 = taxon_id column

#Rarefaction curves - you have to transpose the data using t()
rarecurve(t(obj$data$taxon_abund[, sample_data$SampleID]),
          col = "blue", cex = 1,label=T )
#Species accumulation curve
spec<-specaccum(t(obj$data$taxon_abund[, sample_data$SampleID]),permutations = 999)
plot(spec, ci.col="lightgrey",xlab="# of Samples", ylab="# of OTUs")


###NMDS Ordination###
#Again the obj data has to be transposed, since vegan expects OTUs as coloumns
y<-(t(obj$data$taxon_abund[, sample_data$SampleID]))

#metaMDS with standard settings, but three dimensions (k=3)
ordination<-metaMDS(y, k=3)

plot(ordination) #initial plot

mds_data<-as.data.frame(ordination$points) #making a dataframe for ggplot2
mds_data$SampleID<-rownames(mds_data) #making a column for SampleID
mds_data<-left_join(mds_data,sample_data, by="SampleID") #Joining mds scores and sampledata by SampleID
#look at the data structure
str(mds_data) ##The numeric values are characters!


#library(dplyr) #if you haven't loaded this already
mds_data[,c(10:23)] <- mds_data[,c(10:23)] %>% mutate_if(is.character, as.numeric)
str(mds_data) #now the numeric values are numeric


###Fitting numerical environmental variables to the ordination
env<-envfit(mds_data[,c(1:2)],mds_data[,c(10:23)], 999, na.rm=T)
env.scores <- as.data.frame(scores(env, display = "vectors"))
env.scores <- cbind(env.scores, names = rownames(env.scores))
env.scores<- cbind(env.scores, pval = env$vectors$pvals)
sig.env.scrs <- subset(env.scores, pval<=0.05)

#####plotting the ordination diagram, coloured by sample type
library(ggplot2)
ordination_plot<-ggplot(mds_data, aes(x = MDS1, y = MDS2)) +
  geom_point(mapping=aes(shape=Sample, color=Sample), size=4,stroke=1.5) +
  scale_color_manual(values=c( "lightgreen","darkgreen","#A6611A"))+
  scale_shape_manual(values=c(17,8,19 ))+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(y= "NMDS2", x = "NMDS1")

ordination_plot

#Adding the significant environmental variables fitted using the env.fit() function
ordination_plot+
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=MDS1, y=0, yend=MDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3)+
  geom_text(data = sig.env.scrs, aes(x = MDS1*1.1, y = MDS2*1.1, label = names), size = 5,
            nudge_x = 0.001, nudge_y = 0.001, colour = "grey5")

