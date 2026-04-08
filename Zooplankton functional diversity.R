

# General script to obtain the Zooplankton Funciontal Diversity (FD) from the Ebro´s reservoirs

#The mainly goal is the data analysis for the paper and for the presentation of the AIL congress.  
#Anyway, this script is also useful to obtain the FD of any other group of organisms, 
#as long as you have the trait data and the abundance data.

#Set the working directory

#check in which directory you are now
getwd()

#Set the working directory to the folder where you have the data
# removed by privacy
#setwd("your_path_here"


#Load the libraries ####

library(dplyr)
library(FD)
library(missForest)
library(party)
library(caret)
library(permimp)
library(usdm)


# Load the data ####

#To perform the analysis, you need two data frames: 
#one with the abundance data and another with the trait data.

#Load the abundance data 
#this df must have the species as columns and sites as rows 

abundance_raw <- read.csv("matrix_abundancias_t.csv", row.names = 1)

View(abundance_raw)

#Load the trait data

fd_traits <- read.csv("Rasgosfuncionales_allspecies.csv", header = T, row.names = 1)

str(fd_traits)

#Data manipulation ####

#To calculate the functional diversity, we need to use the dbFD function from the FD package.
#This function requires the abundance data and the trait data as input.

#Calculate the functional diversity
fd_results <- dbFD(fd_traits, abundance_raw)

# We receive an error text, it turns out that in the abundance dataset, this is beacuse
# he names of the species in the columns are separated with a dot(.), for example,
#Filinia.longiseta while in traits df they are fine, without the point.

head(names(abundance_raw))
head(row.names(fd_traits))


#to solve this, we can use the gsub function to replace the dots with spaces.
names(abundance_raw) <- gsub("\\.", " ", names(abundance_raw))

#check
head(names(abundance_raw))

#Now when we try to do again the fd calculation, but it has another problem with the names

fd_results <- dbFD(fd_traits, abundance_raw)

#The problems are that in the traits df we have some points or dots (e.g, cyclop sp.)
# then we remove the dots and spaces in both datasets to make sure are equal 

names(abundance_raw) <- gsub("\\.", "", names(abundance_raw))
names(abundance_raw) <- gsub(" ", "", names(abundance_raw))

row.names(fd_traits) <- gsub("\\.", "", row.names(fd_traits))
row.names(fd_traits) <- gsub(" ", "", row.names(fd_traits))

#now all names are like Acanthocyclopsamericanus without spaces or dots. 

#check both datasets to see if the names are now identical

identical(row.names(fd_traits), colnames(abundance_raw))
#true

#Almost forget to transform the NAs into zeros 

abundance_raw[is.na(abundance_raw)] <- 0

#The columns in the trait df must be numeric or as factor
#We have some columns as character, we need to transform them into factors

fd_traits <- fd_traits %>%
  mutate_if(is.character, as.factor)

str(fd_traits)

#new try
fd_results <- dbFD(fd_traits, abundance_raw)
#success! now we have our results 


#There is a problem with fric. It cannont be calculated if in one reservoir are less 
# than three different functional species. So first we need to check if in all 
#reservoirs we have at least three species.

#obtaining the reservoirs with less than five species just to check 
abundance_raw %>%
  mutate(Embalse = rownames(.),
         Richness = rowSums(. > 0)) %>%
  filter(Richness < 5) %>%
  rowwise() %>%
  mutate(Especies = paste(colnames(abundance_raw)[c_across(-c(Embalse, Richness)) > 0],
                          collapse = ", ")) %>%
  select(Embalse, Richness, Especies) -> few_species_reservoirs

#Anyway the results of FD in those reservoirs are NA, so when do it any other analyses 
#they will be removed.

#Saving the results

#Una ves que obtuvimos datos correctos, nos da una gran cantidad de información, asi que, que hacemos con esto?. Pues guardamos los indices más importantes como un nuevo dataframes para poder verlos y compararlos màs agusto.

fd_indices <- cbind(fd_results$nbsp,
                    fd_results$FRic,
                    fd_results$FEve,
                    fd_results$FDiv,
                    fd_results$FDis, 
                    fd_results$RaoQ)

colnames(fd_indices) <- c("NumbSpecies", "FRic", "FEve", "FDiv", "FDis", "Rao")

#Now we have a new df with the functional diversity indices for each reservoir. 
#We can use this df to perform further analyses, such as comparing the fd between reservoirs, 
#or correlating it with other environmental variables.

#Saving as df
fd_indices <- as.data.frame(fd_indices)

head(fd_indices)

# TIP: if running the dbFD() function appears a message like this:
# "QH6044 qhull option error: cannot open file "vert.txt" something is wrong with the
# R session. So to fix it do the magic informatic solution. Close and open again the programm
#likely it will solve it. At least it solves to me now.


# environmental variables ####

#load the environmental data
envi_raw <- read.csv("environmental_data.csv", header = T, row.names = 1)

str(envi_raw)

#Data manipulation

#transforming characters into factors 

envi_raw <- envi_raw %>%
  mutate_if(is.character, as.factor)

str(envi_raw)


# there are some NAs in the df. So to have a complete df we will try to obtain the 
# data using a random forest technique to complete the missing data.

library(missForest)

envi_missforest <- missForest(envi_raw)

#now we have an object with some lists. 
#extracting the df from the list from the previous output 

envi_complete <- envi_missforest$ximp 

View(envi_complete)

#great now we have the environmental data. 

#Matching fd indices and environmental data ##### 

#Checking the envi data there are more reservoirs in the FD data than in the environmental 
#so we need to check which are the reservoirs with no environmental data and indices data

#First we will check the names of the reservoirs in both datasets
rownames(fd_indices)
rownames(envi_complete)

#Now we will check which reservoirs are in the fd_indices but not in the envi_complete
setdiff(rownames(fd_indices), rownames(envi_complete))

#Now we will check which reservoirs are in the envi_complete but not in the fd_indices
setdiff(rownames(envi_complete), rownames(fd_indices))

#I found that in fd indices the reservoir code of gallipuen is GALL2019 with double LL,
#it has to be only with one L, so we will change it in the fd_indices df

rownames(fd_indices)[rownames(fd_indices) == "GALL2019"] <- "GAL2019"

#Also i noted that GUI2010, ULL2023 were not in the envi, so I checked the original files and 
#probably was my mistake not included. Now it is included in the raw file. 
#Were like 10-12 reservoirs in total, probably because i remove them previously 
#for a normality test and the n I forget to include them again.
#now all the reservoirs are in both datasets, so should be 366 samples each. 

#The names of MAR2013-2016 were changed in excel due to the program used to change it 
# to march-2013. To solve it they were changed to MAG2012-2016. So we need to restore the original names

rownames(envi_complete)[rownames(envi_complete) == "MAG2012"] <- "MAR2012"
rownames(envi_complete)[rownames(envi_complete) == "MAG2013"] <- "MAR2013"
rownames(envi_complete)[rownames(envi_complete) == "MAG2014"] <- "MAR2014"
rownames(envi_complete)[rownames(envi_complete) == "MAG2015"] <- "MAR2015"
rownames(envi_complete)[rownames(envi_complete) == "MAG2016"] <- "MAR2016"

setdiff(rownames(fd_indices), rownames(envi_complete))
#no results, so it means now it correct and we can continue with the analysis.


#Alpha diversity index #####

#Now we will obtain the alpha diversity index 
#this will be done through the vegan package, so we need to load it first

library(vegan)

#now we obtain the different indices using the diversity function from the vegan package.

 data.frame(
  Shannon = diversity(abundance_raw, index = "shannon"),
  Simpson = diversity(abundance_raw, index = "simpson"),
  Richness = rowSums(abundance_raw > 0)) %>% 
  mutate(Pielous = Shannon/log(Richness)) %>% 
  #now we bind the indices with the fd_indices to have all in one df by reservoir
  cbind(fd_indices) -> fd_indices

 
 setdiff(rownames(envi_complete), rownames(fd_indices))
 #we do this again since we modified one name
 rownames(fd_indices)[rownames(fd_indices) == "GALL2019"] <- "GAL2019"
 
 setdiff(rownames(envi_complete), rownames(fd_indices))
 
 #now finally we can combine the fd_indices with the envi_complete to have all the data in one df.
 # and start perfoming some data analysis 
 
#clean and complete dataframe #### 
 
#combining the fd_indices with the envi_complete by rownames 
 
 #Only run it once since if you run it again, it will create a new column 
 #with duplicate colnames and then you will have problems to merge it again.
merge(fd_indices, envi_complete, by = "row.names", all.x = TRUE) %>% 
  #now we will remove the row.names column and set it as rownames again
  tibble::column_to_rownames(var = "Row.names") -> fd_envi_df

#we save this new df as a csv file to have it for future analyses and to not 
#to re-run all the script again in a future.

write.csv(fd_envi_df, "complete_dataset.csv")

#Removing those reservoirs with two or less species which fd indices are NA

fd_envi_df %>% 
  filter(NumbSpecies > 2) -> clean_fd_envi



#Normalized environmental df for PCA and RDA ####

clean_fd_envi %>% 
  #removing all alpha and fd indices and categorical variables
  select(-c(Shannon, Simpson, Richness, Pielous, NumbSpecies,
            FRic, FEve, FDiv, FDis, Rao,
            Trophic.state, TSI.total, Ecological.Potential,
            TS, TSI.PT, TSI.SD, TSI.Chla, Location, Geology, Type)) %>% 
  #removing pH since it is already on scale
  select(-pH) %>% 
  log1p() %>%
  #now pasting all previous removed columns
  cbind(clean_fd_envi %>% 
          select(c(pH, Shannon, Simpson, Richness, Pielous, NumbSpecies,
                   FRic, FEve, FDiv, FDis, Rao,
                   Trophic.state, TSI.total, Ecological.Potential,
                   TS, TSI.PT, TSI.SD, TSI.Chla, Location, Geology,
                   Type)))  -> normalized_fd_envi



#VIF analysis for collinearity #####


library(usdm)

#First we remove all the categorical and TSI indices 

vifstep(normalized_fd_envi[,-c(34:43)], th = 10)

#result: "7 variables from the 37 input variables have collinearity problem: 
# Richness Shannon FDis Simpson Photic.Zone Volume TN 

# NumbSpecies and Richness is the same variable, duplicate, removed.
#Shannon, simpson and Pielous were colinear, removed the last two. 
#Rao with FDis, removed Rao.
#Nitrite with TN and Nitrate, removed nitrite and TN
#Photic zone with Secchi. Removed photic zone.
# Volumen with volume max, removed volume.max

normalized_fd_envi %>% 
  select(-c(NumbSpecies, Simpson, Pielous, Rao, Nitrite, TN,
            Photic.Zone, volume.max)) -> normalized_fd_envi

#checking again
vifstep(normalized_fd_envi[,-c(26:35)], th = 10)

#No variable from the 25 input variables has collinearity problem. 


normalized_fd_envi %>% 
  mutate(Type = as.factor(Type )) -> normalized_fd_envi


normalized_fd_envi %>% 
  mutate(Trophic.state = factor(Trophic.state, 
                                   levels = c("Oligotrophic",
                                              "Mesotrophic",
                                              "Eutrophic",
                                              "Hypereutrophic") )) -> normalized_fd_envi


# PCA #####

#Now that we have the normalized data we can perfomed a PCA


library(factoextra)
library(FactoMineR)



#PCA function 

PCA_fd_envi <- PCA(normalized_fd_envi[-c(26:35)],
                   graph = T,
                   scale.unit = T)

#checking resuls
summary(PCA_fd_envi)

# 21.12% and 11.67 of variance explained by the two first axes

fviz_eig(PCA_fd_envi, addlabels = TRUE, ylim = c(0, 25))



#grafica biplot

fviz_pca_biplot(PCA_fd_envi,
                repel = T,
                geom = "point",
                pointsize = 6,
                col.var = "black",
                labelsize = 7,
                arrowsize = 1.2,
                habillage = normalized_fd_envi$Trophic.state,
                addEllipses = F,
                mean.point = FALSE,
                title = "")+
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 15, face = "bold", color = "black")) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-6, 6)) +
  xlab("PC1 (21.12%)") +
  ylab("PC2 (11.67%)") 
  




#ggplot option

library(ggplot2)
library(dplyr)
library(ggrepel)

# individuals coordinates
sites <- as.data.frame(PCA_fd_envi$ind$coord)
sites$Trophic.state <- normalized_fd_envi$Trophic.state
sites$Type <- normalized_fd_envi$Type

# Variable coordinates
vars <- as.data.frame(PCA_fd_envi$var$coord)
vars$variable <- rownames(vars)

# explained variance
eig <- PCA_fd_envi$eig
pc1 <- round(eig[1, 2], 2)
pc2 <- round(eig[2, 2], 2)

# Arrow escalated to fit the plot 
mult <- min(
  (max(sites$Dim.1) - min(sites$Dim.1)) / (max(vars$Dim.1) - min(vars$Dim.1)),
  (max(sites$Dim.2) - min(sites$Dim.2)) / (max(vars$Dim.2) - min(vars$Dim.2))
) * 0.75

vars$Dim.1 <- vars$Dim.1 * mult
vars$Dim.2 <- vars$Dim.2 * mult

# PCA graphic
 
#by Trophic state

ggplot() +
  geom_point(
    data = sites,
    aes(x = Dim.1, y = Dim.2,
        color = Trophic.state,
        shape = Type),
    size = 5,
    alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  geom_segment(
    data = vars,
    aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 1, color = "black"
  ) +
  geom_text_repel(
    data = vars,
    aes(x = Dim.1, y = Dim.2, label = variable),
    size = 7,
    color = "black",
    segment.color = "grey60",
    max.overlaps = Inf
  ) +
  labs(
    x = paste0("PC1 (", pc1, "%)"),
    y = paste0("PC2 (", pc2, "%)"),
    color = "Trophic state",
    shape = "Type"
  ) +
  coord_cartesian(xlim = c(-10, 10),
                  ylim = c(-6, 4.5)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    text = element_text(size = 17.5, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)) +
  guides(color = guide_legend(nrow = 1, byrow = T), # only in one row
         shape = guide_legend(nrow = 1, byrow = T)) +
   scale_color_manual(values = c(
     "Oligotrophic" = "dodgerblue3",
     "Mesotrophic" = "gold3",
     "Eutrophic" = "darkolivegreen",
     "Hypereutrophic" = "firebrick"
   )) +
  scale_shape_manual(values = c(
    "1"  = 15,
    "7"  = 16,
    "9"  = 17,
    "10" = 18,
    "11" = 19,
    "12" = 20,
    "13" = 8
  ))

#PCA fig saved in 1300*790 size


#By only type

ggplot() +
  geom_point(
    data = sites,
    aes(x = Dim.1, y = Dim.2,
        color = Type),
    size = 5,
    alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  geom_segment(
    data = vars,
    aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 1, color = "black"
  ) +
  geom_text_repel(
    data = vars,
    aes(x = Dim.1, y = Dim.2, label = variable),
    size = 7,
    color = "black",
    segment.color = "grey60",
    max.overlaps = Inf
  ) +
  labs(
    x = paste0("PC1 (", pc1, "%)"),
    y = paste0("PC2 (", pc2, "%)"),
    color = "Type",
  ) +
  coord_cartesian(xlim = c(-10, 10),
                  ylim = c(-6, 4.5)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    text = element_text(size = 17.5, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)) +
  guides(color = guide_legend(nrow = 1, byrow = T))


#PCA fig saved in 1300*790 size


#random forest approach  ####


#Random forest to determine the most important variables that explain the functional diversity indices.
#and acting as a screening tool to select the important variables for a GAM or GAMMs models
 
library(party)
library(caret)
library(permimp)

#the function cforest from the party package to perform the random forest analysis.
#the function permimp from the permimp package to obtain the important variables 
#we will use the normalized df from the PCA


set.seed(1951) #for reproducibility

#FRic  ####

#First we need to train the model to check if we need to change the mtry value,
#which is the number of variables randomly sampled as candidates at each split.

train(FRic ~ ., data = normalized_fd_envi[-c(20,21,23:35)], method = "cforest",
      trControl = trainControl(method = "cv", number = 5),
      tuneGrid = expand.grid(mtry = 1:10))
#the best mtry value is 6, so we will update the model.

#creating the best model 

cf.fric <- cforest(FRic ~ .,
                   data = normalized_fd_envi[-c(20,21,23:35)], #without the another indices and unnecesary data
                   controls = cforest_unbiased(mtry = 6))

#checking statistics of the model
cforestStats(cf.fric)


#obtaining the important variables using the permimp function from the permimp package.
FRic.perimp <- permimp(cf.fric,
                    conditional = TRUE,
                    progressBar = FALSE)

#plotting the important variables from all the df

plot(FRic.perimp,
     type = "bar",
     horizontal = TRUE,
     main = "FRic \n Conditional Permutation Importance (threshold = 0.95)")


#obtaining the 10 most important variables 

FRic.perimp$values %>% 
  sort(decreasing = TRUE) %>%
  head(10) -> FRic_top10

#Results as tibble
FRic_top10 %>% 
  tibble::enframe(name = "Variable", value = "Importance")


#plotting

#marginals beacase the names of the variables are long and we want to see them all,
#we will increase the margins of the plot to have more space for the names.

#modify as needed

par(mar = c(5, #down
            12, #left
            5, #upper
            5 #right
            ))

barplot(
  rev(FRic_top10), #reversing the order to have the most important at the begging
     horiz = TRUE,
     las = 1,
  cex.names = 1.2,
     main = "Top 10 variables explaining FRic",
     sub = "Conditional Permutation Importance (threshold = 0.95) \n",
     col = "steelblue")


#metafiles saved as 1200*787 size



#FEve  ####

#training the model to check if we need to change the mtry value,
#which is the number of variables randomly sampled as candidates at each split.

train(FEve ~ ., data = normalized_fd_envi[-c(20:22,24:35)], method = "cforest",
      trControl = trainControl(method = "cv", number = 5),
      tuneGrid = expand.grid(mtry = 1:10))
#the best mtry value is 5, so we will keep it in the model.


cf.FEve <- cforest(FEve ~ .,
                   data = normalized_fd_envi[-c(20:22,24:35)], #without the another indices
                   controls = cforest_unbiased(mtry = 5))

#checking statistics of the model
cforestStats(cf.FEve) 

#obtaining the important variables using the permimp function from the permimp package.
FEve.perimp <- permimp(cf.FEve,
                       conditional = TRUE,
                       progressBar = FALSE)

#plotting the important variables

plot(FEve.perimp,
     type = "bar",
     horizontal = TRUE,
     main = "FEve \n Conditional Permutation Importance (threshold = 0.95) \n")


#obtaining the 10 most important variables 

FEve.perimp$values %>% 
  sort(decreasing = TRUE) %>%
  head(10) -> FEve_top10


FEve_top10 %>% 
  tibble::enframe(name = "Variable", value = "Importance")

#plotting

#marginals because the names of the variables are long and we want to see them all,
#we will increase the margins of the plot to have more space for the names.

#modify as needed

par(mar = c(5, #down
            13, #left
            5, #upper
            5 #right
))

barplot(
  rev(FEve_top10), #reversing the order to have the most important at the begging
  horiz = TRUE,
  las = 1,
  cex.names = 1.2,
  main = "Top 10 variables explaining FEve",
  sub = "Conditional Permutation Importance (threshold = 0.95)\n",
  col = "steelblue")

#Fig saved as 1200*787 size


#FDis  ####

#training the model to check if we need to change the mtry value,
#which is the number of variables randomly sampled as candidates at each split.

train(FDis ~ ., data = normalized_fd_envi[-c(20:24,26:35)], method = "cforest",
      trControl = trainControl(method = "cv", number = 5),
      tuneGrid = expand.grid(mtry = 1:10))

#the best mtry value is 7, so we will add it into the model.

#creating the best model 

cf.FDis <- cforest(FDis ~ .,
                   data = normalized_fd_envi[-c(20:24,26:35)], #without the another indices
                   controls = cforest_unbiased(mtry = 7))

#checking statistics of the model
cforestStats(cf.FDis)


#obtaining the important variables using the permimp function from the permimp package.
FDis.perimp <- permimp(cf.FDis,
                       conditional = TRUE,
                       progressBar = FALSE)

#plotting the important variables

plot(FDis.perimp,
     type = "bar",
     horizontal = TRUE,
     main = "FDis \n Conditional Permutation Importance (threshold = 0.95) \n")


#obtaining the 10 most important variables 

FDis.perimp$values %>% 
  sort(decreasing = TRUE) %>%
  head(10) -> FDis_top10


#Results as tibble 

FDis_top10 %>% 
  tibble::enframe(name = "Variable", value = "Importance")



#plotting

#marginals because the names of the variables are long and we want to see them all,
#we will increase the margins of the plot to have more space for the names.

#modify as needed

par(mar = c(5, #down
            13, #left
            5, #upper
            5 #right
))

barplot(
  rev(FDis_top10), #reversing the order to have the most important at the begging
  horiz = TRUE,
  las = 1,
  cex.names = 1.2,
  main = "Top 10 variables explaining FDis",
  sub = "Conditional Permutation Importance (threshold = 0.95) \n",
  col = "steelblue")





#FDiv  ####

#training the model to check if we need to change the mtry value,
#which is the number of variables randomly sampled as candidates at each split.

train(FDiv ~ ., data = normalized_fd_envi[-c(20:23,25:35)], method = "cforest",
      trControl = trainControl(method = "cv", number = 5),
      tuneGrid = expand.grid(mtry = 1:10))
#the best mtry value is 4, so we will keep it in the model.


cf.FDiv <- cforest(FDiv ~ .,
                   data = normalized_fd_envi[-c(20:23,25:35)], #without the another indices
                   controls = cforest_unbiased(mtry = 4))

#checking statistics of the model
cforestStats(cf.FDiv)


#obtaining the important variables using the permimp function from the permimp package.
FDiv.perimp <- permimp(cf.FDiv,
                       conditional = TRUE,
                       progressBar = FALSE)

#plotting the important variables

plot(FDiv.perimp,
     type = "bar",
     horizontal = TRUE,
     main = "FDiv \n Conditional Permutation Importance (threshold = 0.95)\n")


#obtaining the 10 most important variables 

FDiv.perimp$values %>% 
  sort(decreasing = TRUE) %>%
  head(10) -> FDiv_top10

#Results as tibble
FDiv_top10 %>% 
  tibble::enframe(name = "Variable", value = "Importance")

#plotting

#marginals because the names of the variables are long and we want to see them all,
#we will increase the margins of the plot to have more space for the names.

#modify as needed

par(mar = c(5, #down
            13, #left
            5, #upper
            5 #right
))

barplot(
  rev(FDiv_top10), #reversing the order to have the most important at the begging
  horiz = TRUE,
  las = 1,
  cex.names = 1.2,
  main = "Top 10 variables explaining FDiv",
  sub = "Conditional Permutation Importance (threshold = 0.95)\n",
  col = "steelblue")


# fig saved as 1200*787 size 



#Richness  ####

#training the model to check if we need to change the mtry value,
#which is the number of variables randomly sampled as candidates at each split.

train(Richness ~ ., data = normalized_fd_envi[-c(20,22:35)], method = "cforest",
      trControl = trainControl(method = "cv", number = 5),
      tuneGrid = expand.grid(mtry = 1:10))
#the best mtry value is 10, so we will keep it in the model.


cf.richness <- cforest(Richness ~ .,
                   data = normalized_fd_envi[-c(20,22:35)], #without the another indices
                   controls = cforest_unbiased(mtry = 10))

#checking statistics of the model
cforestStats(cf.richness)


#obtaining the important variables using the permimp function from the permimp package.
richness.perimp <- permimp(cf.richness,
                       conditional = TRUE,
                       progressBar = FALSE)

#plotting the important variables

plot(richness.perimp,
     type = "bar",
     horizontal = TRUE,
     main = "Richness \n Conditional Permutation Importance (threshold = 0.95)\n")


#obtaining the 10 most important variables 

richness.perimp$values %>% 
  sort(decreasing = TRUE) %>%
  head(10) -> richness_top10

#Results as tibble
richness_top10 %>% 
  tibble::enframe(name = "Variable", value = "Importance")

#plotting

#marginals because the names of the variables are long and we want to see them all,
#we will increase the margins of the plot to have more space for the names.

#modify as needed

par(mar = c(5, #down
            13, #left
            5, #upper
            5 #right
))

barplot(
  rev(richness_top10), #reversing the order to have the most important at the begging
  horiz = TRUE,
  las = 1,
  cex.names = 1.2,
  main = "Top 10 variables explaining richness",
  sub = "Conditional Permutation Importance (threshold = 0.95)\n",
  col = "steelblue")


# fig saved as 1200*787 size 


#Shannon  ####

#training the model to check if we need to change the mtry value,
#which is the number of variables randomly sampled as candidates at each split.

train(Shannon ~ ., data = normalized_fd_envi[-c(21:35)], method = "cforest",
      trControl = trainControl(method = "cv", number = 5),
      tuneGrid = expand.grid(mtry = 1:10))
#the best mtry value is 3, so we will keep it in the model.


cf.shannon <- cforest(Shannon ~ .,
                       data = normalized_fd_envi[-c(21:35)], #without the another indices
                       controls = cforest_unbiased(mtry = 3))

#checking statistics of the model
cforestStats(cf.shannon)


#obtaining the important variables using the permimp function from the permimp package.
shannon.perimp <- permimp(cf.shannon,
                           conditional = TRUE,
                           progressBar = FALSE)

#plotting the important variables

plot(shannon.perimp,
     type = "bar",
     horizontal = TRUE,
     main = "Shannon \n Conditional Permutation Importance (threshold = 0.95)\n")


#obtaining the 10 most important variables 

shannon.perimp$values %>% 
  sort(decreasing = TRUE) %>%
  head(10) -> shannon_top10

#Results as tibble
shannon_top10 %>% 
  tibble::enframe(name = "Variable", value = "Importance")

#plotting

#marginals because the names of the variables are long and we want to see them all,
#we will increase the margins of the plot to have more space for the names.

#modify as needed

par(mar = c(5, #down
            13, #left
            5, #upper
            5 #right
))

barplot(
  rev(shannon_top10), #reversing the order to have the most important at the begging
  horiz = TRUE,
  las = 1,
  cex.names = 1.2,
  main = "Top 10 variables explaining shannon index",
  sub = "Conditional Permutation Importance (threshold = 0.95)\n",
  col = "steelblue")


# fig saved as 1200*787 size 



#results as tibble from  FRic_top10, FEve_top10, FDis_top10 and FDiv_top10 together
FRic_top10 %>% 
  tibble::enframe(name = "FRic", value = "Importance") %>%  
  bind_cols(FEve_top10 %>% 
              tibble::enframe(name = "FEve", value = "Importance")) %>% 
  bind_cols(FDis_top10 %>% 
              tibble::enframe(name = "FDis", value = "Importance")) %>% 
  bind_cols(FDiv_top10 %>% 
              tibble::enframe(name = "FDiv", value = "Importance")) %>% 
  bind_cols(richness_top10 %>% 
              tibble::enframe(name = "Richness", value = "Importance")) %>%
  bind_cols(shannon_top10 %>% 
              tibble::enframe(name = "Shannon", value = "Importance")
            ) -> results_important_variables_indices

View(results_important_variables_indices)



# results from cforeststats for all the models together in a tibble

cforestStats(cf.fric) %>%
  tibble::enframe(name = "FRic", value = "Value") %>% 
  bind_cols(cforestStats(cf.FEve) %>% 
              tibble::enframe(name = "FEve", value = "Value")) %>% 
  bind_cols(cforestStats(cf.FDis) %>% 
              tibble::enframe(name = "FDis", value = "Value")) %>% 
  bind_cols(cforestStats(cf.FDiv) %>% 
              tibble::enframe(name = "FDiv", value = "Value")) %>% 
  bind_cols(cforestStats(cf.richness) %>% 
              tibble::enframe(name = "Richness", value = "Value")) %>%
  bind_cols(cforestStats(cf.shannon) %>% 
              tibble::enframe(name = "Shannon", value = "Value")
            ) -> results_cforest_stats_indices

View(results_cforest_stats_indices)


#Here I stop for today, continue with GAM models


