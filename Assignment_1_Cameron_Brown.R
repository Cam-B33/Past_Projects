#In recent decades there has been a push in the conservation community to increase the amount of Marine area protected to preserve marine biodiversity (CBD, 2010; IUCN 2014). Countries have in turn significantly increased the amount of Marine area they are protecting. In 2013-2014 the USA, New Caledonia, Australia and South Africa increased their Marine Protected areas by a combined 3.5 million square kilometres (Boonzaier and Pauly, 2016). For this assignment I wanted to explore the impact of Marine protected areas on the diversity of a specific taxonomic group, the Phylum Porifera (Sponges).  This taxonomic group was suitable for my research because of the Porifera's global distribution (Renard et al. 2013), which would hopefully allow me to compare as many data points as possible. My hypothesis is that the proportion of a countries territorial waters that are protected will have a significant positive correlation with the diversity (Shannon Index) of Porifera specimens collected from that country. 
#For my data on percent marine protected area per country I used a dataset from the world bank, available at https://data.worldbank.org/indicator/ER.MRN.PTMR.ZS 
#First I'll load all the packages I will need for my analysis
library(readr)
library(tidyverse)
library(vegan)
library(dpylr)
library(ggplot)
# Now I'll load my data directly from BOLD into our R workspace as a dataframe
dfBOLD_Assignment_1 <- read_tsv(file = "http://www.boldsystems.org/index.php/API_Public/combined?taxon=Porifera&format=tsv")
#For consistency in future analyses we'll save the BOLD data as a TSV file for use in subsequent sessions 
write_tsv(dfBOLD_Assignment_1, "Porifera_BOLD_data.tsv")
# If I had already done the last step in a previous I would use this line to load the BOLD data into R as a data frame. 
dfBOLD_Assignment_1 <- read_tsv(file = "Porifera_BOLD_data.tsv")
# Some basic housekeeping to explore my data set 
class(dfBOLD_Assignment_1)
names(dfBOLD_Assignment_1)
dim(dfBOLD_Assignment_1)
summary(dfBOLD_Assignment_1)
# For my analysis I'm using country and BIN data, so I want to organize (group) my data set by country and BIN, and I want to get a speciment count per country and per BIN
dfWorld.Count.by.BIN <- dfBOLD_Assignment_1 %>%
  group_by(bin_uri,country) %>%
  count(bin_uri,country)
# Since I'm using the package Vegan I will need to get my data into the community data object format, using the function Pivot Wider. This will organize my data so that my rows are countries and columns are BIN's. Rows will now contain all the specimens found in a certain country, organized by BIN. This format will allow me to calculate diversity by country. 
data <-pivot_wider(data = dfWorld.Count.by.BIN, names_from  = bin_uri, values_from = n)
view(data)
# I'm now going to get rid of the Country = NA row since I can't use it for my analysis.
data.no.na.country <-data[-7,]
# Quick check to see that it worked
view(data.no.na.country)
# Now I want to get my data ready so that I can calculate Porifera diversity by country (row). This was the trickiest part of this assignment for me, see acknowledgments for more details.
# I need to calculate the diversity of each row in my data frame, but to do so I need to remove the first column, which has the country names, because the diversity function can't calculate data of charachter class type. However, just removing the first column would result in a NULL row names error (R Documentation). Therefore I'm assigning each rows row name based upon the data found in the first column of each row i.e country
rownames(data.no.na.country)<- data.no.na.country$country
# Now I can remove the first column.
data.no.na.country <- data.no.na.country[,-1]
# Quick check to see that it worked
view(data.no.na.country)
# Currently all cells of the data frame that don't have a specimen count have value NA.This is the same problem as before. Since the diversity function can't calculate the charachter class type I need to change those NA's to 0's. 
data.no.na.country[is.na(data.no.na.country)] <- 0
# Quick check to see that it worked
view(data.no.na.country)
# Now I want to calculate the Porifera diversity, using the Shannon Index, per country. I then want to create a data frame that only includes the diversity value of each row and that preserves the order of the rows, so that the country name can be easily added back in. 
data.diversity <- as.data.frame(diversity(data.no.na.country, index = "shannon")) 

#I'm now creating a new column in data.diversity called country, and I'm populating it with country names form my 'data' data frame. Minus the 7th row because that was the Country = NA row. 
data.diversity$country <- data$country[-7]
# Quick check to see that it worked
view(data.diversity)
# There are a couple countires that have a Shannon diversity value of zero so I'm going to remove them from my data frame. 
df.diversity.no.zero <- data.diversity %>% 
  filter(data.diversity$`diversity(data.no.na.country, index = "shannon")` >0)
view(df.diversity.no.zero)
#I'm now going to pull my protected marine area by country data from the World Bank into my workspace as a data frame
# It should be noted that the I was having alot of issues with subsetting the original data in R, so for the sake of time I made some format changes in excel before pulling the data into R. 
Protected_Marine_Area_by_Country_Formatted <- read_csv("Protected Marine Area by Country Formatted.csv")
View(Protected_Marine_Area_by_Country_Formatted)
#I'm subsetting the data frame so that it only includes rows with a country name, and only includes the column with the most recent data (2018)
dfPMA_Subset <- Protected_Marine_Area_by_Country_Formatted[5:269, c(1,4)]
#Quick check to see that it worked
view(dfPMA_Subset)
# Now I want to remove rows that don't have a value for 'Percent territorial waters protected'
#Currently that column name is ...4, terrible name, but I will change that soon.
dfPMA_Subset_rm.na <- dfPMA_Subset %>% 
  filter(!is.na(...4))
#I'm now going to find what countries in the World Bank Dataset are only present in df.diversity.no.zero data frame.
#First I will make a list of the country names from my df.diversity.no.zero data frame
Country_list <-unique(df.diversity.no.zero$country)
# I know want to see whcih rows from my protected marine area data frame match the country names in Country_list
match(dfPMA_Subset_rm.na$`Data Source`, Country_list)
# I'm now going to subset my protected marine area data frame so that it only includes rows that match the country names in Country_List
dfPMA.in.BOLD <-dfPMA_Subset_rm.na[dfPMA_Subset_rm.na$`Data Source`%in%Country_list,]
#Now because my dfPMA.in.BOLD data frame has less rows then the df.divesity.no.zero data frame, I now need to do the same subsetting process but in reverse. 
Country_list_PMA <-unique(dfPMA.in.BOLD$'Data Source')
match(df.diversity.no.zero$country, Country_list_PMA)
df.diversity.no.zero <-df.diversity.no.zero[df.diversity.no.zero$country%in%Country_list_PMA,]
# Since I want to combine the two data frames I wan the country data to be in the same order for both. Thefore I'm going to order each data frame alphabetically by country
df.diversity.no.zero<- df.diversity.no.zero[order(df.diversity.no.zero$country),]
dfPMA.in.BOLD<- dfPMA.in.BOLD[order(dfPMA.in.BOLD$'Data Source'),]
#Fixing the terrible column names so that the data frames can be joined by the country name 
colnames(dfPMA.in.BOLD) <- c('country','pma')
#Now i want to combine the data frames using the country column, which is identical for both. 

df.diversity.to.pma <- inner_join(df.diversity.no.zero, dfPMA.in.BOLD, by = "country")
# The column name for diversity index is problematic since R thinks Im trying to use the diversity function when Im subsetting that column. 
colnames(df.diversity.to.pma) <- c('Shannon_Index', 'country','pma')
# Now Im going to take a look at the data frame I'll use for my analysis
view(df.diversity.to.pma)
# Now Im going to make a histogram showing the distribution of values for % protected marine area by number of countires 
ggplot(data = df.diversity.to.pma) +
  geom_histogram(mapping = aes(x = pma)) +
  labs(title = "Histogram of Protected Marine Area", x = " Protected Marine Area by Percent of Territorial Waters", y = "Number of Countries")
# All the data looks reasonable, except for Slovenia. The data is saying that Slovenia protects 213% of it's territorial waters. This seems highly improbable and potentially an error, so I'm going to exclude it from my dataset. 
df.diversity.to.pma<- df.diversity.to.pma[-38,]
# Now Im going to make a histogram showing the distribution of values for the Shannon diversity index by number of countires. 
ggplot(data = df.diversity.to.pma) +
  geom_histogram(mapping = aes(x = Shannon_Index)) +
  labs(title = "Histogram of Shannon Diversity Index", x = " Shannon Diversity Index", y = "Number of Countries")
# Now I'm going to create a scatter plot so I can visualize the relationship between the two variables. Based on my hypothesis I will be using % protected marine are as the independent variable and Shannon index as the dependent variable. 
ggplot(data = df.diversity.to.pma) +
  geom_point(mapping = aes(x = pma, y = Shannon_Index)) +
  labs(title = "Protected Marine Area vs Diversity of Order Poriferia", x = "Percent Protected Marine Area", y = "Poriferia Diversity")
# Now I want to run a linear regression to test my hypothesis. As stated above % protected marine are is the independent variable and the  Shannon index score is the dependent variable. 
pma.spong_diversity.lm <- lm(Shannon_Index ~ pma, data = df.diversity.to.pma)
summary(pma.spong_diversity.lm)
# We can see from the R squared value that the protected marine area has a negligble influence on sponge diversity by country, so I can reject my hypothesis. 
#Results and Discussion
#Based on the results of my analysis I can reject my hypothesis, as the linear regression shows that percent protected marine area does not influence Porifera diversity (very small R squared value). There are many potential causes for this. For example sponge diversity is known to differentiate with latitude (Ruzicka and Gleason, 2008) which is a factor I did not take into account with my analysis. While I did not use any species accumulation curves in my analysis, it would be reasonable to expect that not all countries are equally well sampled, which could also impact the results of my analysis. Another potential factor is that protected marine areas do not all have the same restrictions in place (De Santo, 2013), so treating all areas as equal is a major assumption.
#If I were to further develop this research question into a larger project I would need to treat these potential co factors with more serious consideration. This would mean developing a list of potential hypothesis (forest hypothesis methodology) and use an appropriate statistical test (ANOVA, multiple regression, etc) to compare the hypothesis to identify which cofactors, if any, significantly impact diversity by country. 
#Acknowledgments 
#Thankyou to Abi Yogasekaram for helping me trouble shoot the error messages I was getting when I was attempting to apply the Shannon diversity index function to each country. 
#Below are some of the error messages I was getting
#Error in UseMethod("rowwise") :
#no applicable method for 'rowwise' applied to an object of class "NULL"  
#Error in diversity(x, index = shannon) : input data must be numeric
#I thought the errors I was getting were because the diversity function wasn't being applied by row, or that it wasn't working because I needed to create a new column to populate with the diversity values. Abi pointed out to me (which in retrospect is quite obvious from the error messages) that the errors were mainly occurring because I was attempting to apply diversity to Null and Character data types when the diversity function only works on numeric data. Once she pointed that out the solution of changing NA values to 0 and removing the first column (country) became self evident. The main lesson I learned from this is to not assume why the error messages are occurring, but to read them and directly address what the error message is saying. My assumptions were only overcomplicating the problem because I started trying to use functions such as rowwise() and mutate() which were unnecessary.  She also suggested I check the documentation for the match function in r when I was trying to match up my two data frames. Previously I had been trying to filter countries by Country_List (filter(country == Country_List)), but this won't work because R thinks I'm trying to find cells that contain the entire list as is value. 
#References 