#The objective of my project is to build two classifiers that can diferentiate between two mitochondial protein coding genes within the Order Urodela, or Salamanders. The genes I have chosen, COX1 and CytB, are common targets for enviroemntal DNA studies (Deiner et al. 2017). If I wanted to create a duplexed assay specific to the COX1 and CytB genes for salamanders, perhaps to do general amphibian monitoring, I would want to develop primers that were specific to these genes. Since these genes are both mitochondrial protein coding genes they are likely quite similar, and therfore I need to be careful to chose a primer sequence which will be specific to only one gene and not overlap between the two. Using the classifers I built I can test any primer sequence candiadate to make sure it has good specificity to the gene I'm trying to target. 

#I chose to do my project on Salamanders because they are considered an indicator species (Davic and Welsh 2004) which means we can monitor the health of an ecosystem by monitoring Salamander populations, and a generalized assay could be of use as a rough edged way to do enviromental health monitoring. I also chose them because they are awesome and I did my undergrad thesis on Salamanders, so they have a special place in my heart :)
#First step is loading the required packages
library(rentrez)
#install.packages("seqinr")
library(seqinr)
#BiocManager::install("Biostrings")
library(Biostrings)
#install.packages("tidyverse")
library(tidyverse)
#For this project I want to compare two mitochondrial genes (cytB and COX1) for the Order Urudela (Salamanders).
#First I'm going to search the NCBI Database nuccore for cytb sequences for the Oder Urodela that are between 600 and 1000 bp in length. This lets me make sure I'm not getting whole geneome or partial genes in my search. 
cytb_Urodela_Search <- entrez_search(db = "nuccore", term = "(Urodela[ORGN] AND CytB[Gene] AND 600:1000[SLEN])")
#Now I want to know how many "hits" my search got. This number will be important when I pull my data into a FASTA file. 
cytb_Urodela_Search$count
#Now I'm searching the same NCBI database, nuccore, for COX1 gene sequences for the Order Urodela. However I do not have a sequence length range in my search, like I did in the cytb sequence call. This is because I could not find a sequence range for the COX1 gene so none was used in the search call. Obviously this creates a imbalance in the "refinemnt" of the two sets of sequences I'm searching for. Further down in the code I will use other methods refine my COX1 data without knowing the length of the gene.
COX1_Urodela_Search <- entrez_search(db = "nuccore", term = "(Urodela[ORGN] AND COX1[Gene])")
#Now I want to know how many "hits" my search got. This number will be important when I fetch my data into a fasta file. 
COX1_Urodela_Search$count
#I want to get all the hits from my searches but since there are thousands of sequences I can't fetch the sequences directly to a FASTA file. I need to create a web history object that holds my call on the NCBI servers. 
cytb_Urodela_Search <- entrez_search(db = "nuccore", term = "Urodela[ORGN] AND CytB[Gene] AND 600:1000[SLEN]", retmax = 4410, use_history = T)

cytb_Urodela_Search$web_history #Storing call on server

#Now I'm going to use the function Sally wrote to fetch my web history object in FASTA format and write it to a file. 
source("Entrez_Functions.R")

#Now I'm going to write a series of FASTA files that each contain 100 sequences from my NCBI search call. 

FetchFastaFiles(searchTerm = "Urodela[ORGN] AND CytB[Gene] AND 600:1000[SLEN]", seqsPerFile = 100, fastaFileName = "Urodela_cytb")
# Now I'm taking all the FASTA files in my working directory that have Urodela_cytb in the name and combine them into a data frame. 
dfUrodela_Cytb <- MergeFastaFiles(filePattern = "Urodela_cytb*")
#Quick check to make sure it worked
view(dfUrodela_Cytb)
#Now I'm repeating the process with my COX1 data call
COX1_Urodela_Search <- entrez_search(db = "nuccore", term =  "Urodela[ORGN] AND COX1[Gene]", retmax = 2302, use_history = T)

COX1_Urodela_Search
COX1_Urodela_Search$web_history 
source("Entrez_Functions.R")

FetchFastaFiles(searchTerm = "Urodela[ORGN] AND COX1[Gene]", seqsPerFile = 100, fastaFileName = "Urodela_COX1")

dfUrodela_COX1 <- MergeFastaFiles(filePattern = "Urodela_COX1*")

view(dfUrodela_COX1)
#Now we need to load the packages needed for Random Forest 

library(randomForest)

#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("Biostrings")
library(Biostrings)
#I want to take a look at the distribution of sequence lengths for my two data frames, so I'm making a histogram of the sequence lengths using the fuction nchar which counts the number of letters (nucleotides) in each sequence. 
hist(nchar(dfUrodela_COX1$Sequence), xlab = "Sequence Length", ylab = "
     Frequency", main = "Frequency Histogram of COX1 Sequence Lengths for Order Urodela")
#Adjusting axis do daa doesnt go beyond the axis
axis(1,at=seq(0,25000,by=5000))
#I can see there are some sequences that are much bigger then the rest (>5000bp) these are likely whole mitochondrial genomes.

# Now I will make the same chart for the CytB data.
hist(nchar(dfUrodela_Cytb$Sequence) ,xlab = "Sequence Length", ylab = "
     Frequency", main = "Frequency Histogram of cytB Sequence Lengths for Order Urodela")
#As expected the cytB has a much more narrow distribution of sequence lengths. 
#I want to improve the quality of the sequences in my data frames and remove poor quality sequences that I don't want in my analysis.
#First I'm creating a new column called sequences2 that I will be editing the sequences in. This is so I have a redundant copy of the original unedited sequences.
#Then I will be removing any gaps or N's from the beggining of the sequences in the Sequences2 Column
#^[-N]+ can be transletd to "multiple instances(+) of gaps and unknown nucleotids(-N) at the beggining of the sequence(^).
#In line three Im again removing instances of gaps and N's, but this time from the end($) of the sequence. 
#In line four I'm continuing the edit the Sequence2 column, this time I'm removing all gaps from the sequences. 
#In line 5 I'm taking the sequences in the Sequences2 column and removing any that still #have greater or equal to 5% of their sequence made up of N's. 
dffiltered_COX1 <- dfUrodela_COX1 %>%
  mutate(Sequences2 = str_remove(Sequence, "^[-N]+")) %>%
  mutate(Sequences2 = str_remove(Sequences2, "[-N]+$")) %>%
  mutate(Sequences2 = str_remove_all(Sequences2, "-+")) %>%
  filter(str_count(Sequences2, "N") <= (0.05 * str_count(Sequence)))
#Check to see that some rows were removed 
count(dffiltered_COX1)
# Now I'm repeating the same preocess for the cytb data
dffiltered_cytb <- dfUrodela_Cytb %>%
  mutate(Sequences2 = str_remove(Sequence, "^[-N]+")) %>%
  mutate(Sequences2 = str_remove(Sequences2, "[-N]+$")) %>%
  mutate(Sequences2 = str_remove_all(Sequences2, "-+")) %>%
  filter(str_count(Sequences2, "N") <= (0.05 * str_count(Sequence)))
count(dffiltered_cytb)
#Now I need a way to refine my COX1 data to remove whole genomes and partial gene fragments.
#To do this I'm only going to keep the 50% of sequences at the centre of the distribution. #First I'll need to identify the first and third quartile to set the cut off points.
q1 <- quantile(nchar(dffiltered_COX1$Sequences2), probs = 0.25, na.rm = TRUE)
q1

q3 <- quantile(nchar(dffiltered_COX1$Sequences2), probs = 0.75, na.rm = TRUE)
q3
#Setting q1 and q3 as variables makes the following lines of code much easier.
#What I am doing now is creating a new dataframe for COX1 where all Sequences (from the Sequences2 column) whose length falls outside the first and third quartile threshold are excluded. 
dfCOX1_Urodela <- dffiltered_COX1 %>%
  filter(((str_count(Sequences2) >= q1 & str_count(Sequences2) <= q3))) 
#Now I'm going to change the name of the cytb data to be more congruent with the COX1 data frame 
dfcytb_Urodela <- dffiltered_cytb
#For my supervised classification problem I want to use k-mer frequencies, and to calculate k-mer frequencies I need to use the Biostrings package. Therefore I need to conver the sequences in Sequences2 to DNAStringset format so that I can use Biostrings functions on them. 
dfCOX1_Urodela <- as.data.frame(dfCOX1_Urodela)
dfCOX1_Urodela$Sequences2 <- DNAStringSet(dfCOX1_Urodela$Sequences2)
#Checking to see that it worked. 
class(dfCOX1_Urodela$Sequences2)
#Repeating the same process for the cytb data.
dfcytb_Urodela <- as.data.frame(dfcytb_Urodela)
dfcytb_Urodela$Sequences2 <- DNAStringSet(dfcytb_Urodela$Sequences2)
class(dfcytb_Urodela$Sequences2)
#I'm now going to use the letterFrequency function to calculate the frequency of individual nucleotides in each of the sequences in Sequences2, I'm then using Cbind to merge this new information with my original data frame. So now each row contains the sequence and nucleotide prequencies. 
dfcytb_Urodela <- cbind(dfcytb_Urodela, as.data.frame(letterFrequency(dfcytb_Urodela$Sequences2, letters = c("A", "C","G", "T"))))

#Now I'm going to calculate nucleotide proportions. This will be useful for when I calculate dinucleotide frequency as by using nucleotide proportions I can account for the variability in sequence length when calculating dinucleotide frequencies. 
dfcytb_Urodela$Aprop <- (dfcytb_Urodela$A) / (dfcytb_Urodela$A + dfcytb_Urodela$T + dfcytb_Urodela$C + dfcytb_Urodela$G)

dfcytb_Urodela$Tprop <- (dfcytb_Urodela$T) / (dfcytb_Urodela$A + dfcytb_Urodela$T + dfcytb_Urodela$C + dfcytb_Urodela$G)

dfcytb_Urodela$Gprop <- (dfcytb_Urodela$G) / (dfcytb_Urodela$A + dfcytb_Urodela$T + dfcytb_Urodela$C + dfcytb_Urodela$G)

# Once ATG proportions there is no need to calculate C because it is now implicit since the other nucelotide propotions are known. PropC = 1-(Aprop + Gprop + Tprop). 

# Now I'm using the function dinucleotideFrequency to caluclate the dinucleotide frequencies and adding these values to my data frame using cbind.
dfcytb_Urodela <- cbind(dfcytb_Urodela, as.data.frame(dinucleotideFrequency(dfcytb_Urodela$Sequences2, as.prob = TRUE)))

#Check data frame to see if it worked 
view(dfcytb_Urodela)


#Same for COX1
dfCOX1_Urodela <- cbind(dfCOX1_Urodela, as.data.frame(letterFrequency(dfCOX1_Urodela$Sequences2, letters = c("A", "C","G", "T"))))


dfCOX1_Urodela$Aprop <- (dfCOX1_Urodela$A) / (dfCOX1_Urodela$A + dfCOX1_Urodela$T + dfCOX1_Urodela$C + dfCOX1_Urodela$G)

dfCOX1_Urodela$Tprop <- (dfCOX1_Urodela$T) / (dfCOX1_Urodela$A + dfCOX1_Urodela$T + dfCOX1_Urodela$C + dfCOX1_Urodela$G)

dfCOX1_Urodela$Gprop <- (dfCOX1_Urodela$G) / (dfCOX1_Urodela$A + dfCOX1_Urodela$T + dfCOX1_Urodela$C + dfCOX1_Urodela$G)

dfCOX1_Urodela <- cbind(dfCOX1_Urodela, as.data.frame(dinucleotideFrequency(dfCOX1_Urodela$Sequences2, as.prob = TRUE)))

view(dfCOX1_Urodela)

#Now that I've calculated the k-mer frequencies I can convert the Sequences2 column from DNAStringset back to character data. 
dfcytb_Urodela$Sequences2 <- as.character(dfcytb_Urodela$Sequences2)

dfCOX1_Urodela$Sequences2 <- as.character(dfCOX1_Urodela$Sequences2)
#Lets look at the number of sequences in each data frame
count(dfcytb_Urodela)
count(dfCOX1_Urodela)
#There is a severe class imbalance, to even this out I'm going to randomly select 1300 rows from the cytb sequences so that the data frames are equal size. 
set.seed(62)
dfcytb <- dfcytb_Urodela %>%
  sample_n(1300)
#Check that it worked
count(dfcytb)
#Now I'm going to make a new column in each database that holds the classification value. 
Gene <- c("cytB")

dfcytb$Gene_Name <- Gene
#Check that it worked
view(dfcytb)
#Same thing for COX1.
Gene2 <- c("COX1")

dfCOX1_Urodela$Gene_Name <- Gene2
#Check that it worked
view(dfCOX1_Urodela)
#Now I'm going to split my data into training and validation data sets. First I'm going to set the seed so that my analysis is reproducible. 
set.seed(315)
#Im going to do a 80/20 training-validation split. So first I will randomly sample 80% of the rows and create a training data frame. 
dfcytb_Training <- dfcytb %>%
  sample_frac(0.8)
view(dfcytb_Training)
#Now I'm repeating this process with the COX1 data frame. 
set.seed(81)

dfCOX1_Training <- dfCOX1_Urodela %>%
  sample_frac(0.8)
view(dfCOX1_Training)
#Now I am using rbind to combine my training data frames together 
dfTraining <- rbind(dfcytb_Training,dfCOX1_Training)
view(dfTraining)
#Now I'm combining the original gene data frames.
dfUrodela_Genes <-rbind(dfCOX1_Urodela, dfcytb)
view(dfUrodela_Genes)
#Now I'm going to create my validation data set. To do this I'm going to take the data frame of the of the combined COX1 and CytB data (pre training data selection) and filter out all the rows that are present in my training data frame. This will produce a data frame that is 20% of the original Urodela Genes data set, and is distinct from the training data set (!Title %in% dfTraining$Title).


set.seed(13)
dfValidation <- dfUrodela_Genes %>%
  filter(!Title %in% dfTraining$Title)
#No need for sample frac because it will automatically be 20%.
view(dfValidation)
#Now I'm going to build a Random forest classifier using the data sets. I'm going to start with dimer frequencies. I'm setting the tree number to 840 because this is 40% of the rows in my training data frame. 

gene_classifier_dimer <- randomForest::randomForest(x = dfTraining[,11:26], y = as.factor(dfTraining$Gene_Name), ntree = 840, importance = TRUE)

#Let's look at the results.
gene_classifier_dimer
#Now I'm going to take a look at some key features of my Random forest classifier. 
#Now I'm looking at the relative importance of each of the dimers. 
gene_classifier_dimer$importance

#Now I'm looking how many times each row was left out of bag and predicted, this will give me and indication as to whether I used the appropriate number of trees. 
gene_classifier_dimer$oob.times

#Checking the error rate 
gene_classifier_dimer$err.rate
# We can see that each row was estimate many times and the error rate is very low so we don't need to adjust the number of trees. 
#confusion matrix. Checking the false positives and negatives for each of the genes.  
gene_classifier_dimer$confusion

#Seeing how many votes each gene glass got for every row. 
gene_classifier_dimer$votes
#Now I'm using the validation data frame to test my Random Forest Classifier on unseen data. I'm the prediction form the classifier will be entered into the final row of the data frame. 
predictValidation <- predict(gene_classifier_dimer, dfValidation[, c(27,11:26)])

#Let's have a look at the data format here. There are the predictions.
predictValidation
class(predictValidation)
length(predictValidation)

#Creating a confusion matrix to look at the false negatives and positives from my classifier. 
table(observed = dfValidation$Gene_Name, predicted = predictValidation)
# Now I'm going to analyze my data using a different classification Algortithm, Classification and Regression Trees or CART. However I think CART is trademarked so the R package is called Recursive Partitioning and Regression Trees or rpart (https://www.rdocumentation.org/packages/rpart/versions/4.1-15/topics/rpart)
#install.packages("rpart")

library(rpart) 
#To build my CART classifier I'm going to use the function mytree on my training data frame I made previously. I'm setting the predicted variable as the Gene Name (CytB or COX1) and I'm using the Dimer frequencies as the independent variable. I'm using method "class" because y (Gene Name) is a factor. I'm using the default for everything else. 
mytree <- rpart(
  Gene_Name ~ AA + AC + AG + AT + CA + CC + CG + CT + GA + GC + GG + GT + TA + TC + TG + TT, 
  data = dfTraining, 
  method = "class")

#Try # plot mytree
#install.packages("rattle")
#install.packages("rpart.plot")
#install.packages("RColorBrewer")
library(rattle)
library(rpart.plot)
library(RColorBrewer)

# Now I'm going to plot my tree, the descision tree shows me what variable i.e dimers were important in determining between the two genes, as well as what frequency value for each dimer is significant in differentiating the two genes. 
fancyRpartPlot(mytree, caption = NULL)


#Now I'm going to go through the process of pruning my tree to prevent over fitting and reduce the complexity of the classifier. 
mytree$variable.importance
printcp(mytree)
#To prune my tree I'm going to set the cp value to just above the cp value of the tree with the smallest non zero xerror rate. In the case of my tree the cp was 0.01 so I will prune my tree with a cp value of 0.011. https://rdrr.io/cran/rpart/f/inst/doc/longintro.pdf
mytree <- prune(mytree, cp = 0.011)
fancyRpartPlot(mytree, caption = "Descision Tree for determining COX1 or Cytb based on dinucleotide frequency in the Order Urodela")
#Now I'm going to add two columns to my data frame, predicted gene name and predicted gene probability. Then I'm going to populate these columns using mytree to predict gene type, and then the probability based on the descision tree that this prediction is true. 
dfValidation$Predicted_gene_name <- predict(mytree, newdata = dfValidation, type = "class")
dfValidation$Predicted_gene_prob<- predict(mytree, newdata = dfValidation, type = "prob")
#Now I'm going to look at the rate of false negatives and positives for my classifier. 
table(observed = dfValidation$Gene_Name, predicted = dfValidation$Predicted_gene_name)
#Conclusion 
#Both my RandomForest and CART classifier worked very well using only dinucleotide frequencies, however my Random Forest classifier was the more accurate of the two. There are a few key caveats to my study. The first is the refinemnt in of the COX1 Gene sequences, because I didn't know what nucleotide length this gene was I was hoping that partial sequences and whole genomes would be excluded by removing the first and fourth quartile of data, and this assumption may have impacted the quality of my classifier to accuraltey identify the COX1 gene by including non COX1 sequences or partial sequences. For future research on this topic I would need to do more research into the lenght of the COX1 gene. Another caveat to my project is that I don't provide justification for using an 80/20 Train-Test split other then it is commonly used for ML problems. If I were to expand this project I would like to explore more rigorous methods of determining the training test split such as the method describe in Guyon (1997).  
#Academic References 
#Deiner, K, Bik, HM, Mächler, E, et al. Environmental DNA metabarcoding: Transforming how we survey animal and plant communities. Mol Ecol. 2017; 26: 5872- 5895. https://doi.org/10.1111/mec.14350
#Davic RD, Welsh HH Jr. 2004. On the ecological roles of salamanders. Annu Rev Ecol Evol Syst. 35(1):405-434.
#Guyon, Isabelle. "A Scaling Law for the Validation-Set Training-Set Size Ratio." (1997).

#I also found these blog posts to be helpful 
##https://towardsdatascience.com/decision-trees-d07e0f420175
##https://www.gormanalysis.com/blog/decision-trees-in-r-using-rpart/