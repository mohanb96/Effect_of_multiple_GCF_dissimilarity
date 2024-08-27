################################################
###. Part 1 Calculating Dissimilarity Index. ###
################################################


# 1. Loading R packages

library(tidyverse)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(splitstackshape)
library(dplyr)

# 2. Import data

df<-read.csv("/Data/singal_factor.csv")
treatments = c("AntiB","Cu","Drought","Fung","Herb","Insect","Li","MP","Nit","PFAS","Salin","Surf")
responses = c("actv_ace","actv_cello","actv_gluco","actv_phos","Decom","PH","WSA")
n_iter =10000
# 3. Clustering of single factors

# 3.1 Calculate mean effect sizes of single factors

# The clustering of 12 global change factors is based on the experimental measurement of the effects of single factors on soil, each factor has 8 replicates. 

response_ES_bs =list()
for (i_response in responses) {
  response_ES_bs[[i_response]] = Bootstrap_ES_rep(i_response)
}
mean_sig_ES=data.frame(matrix(data = NA, nrow = 12, ncol = 7))
colnames(mean_sig_ES)=responses
rownames(mean_sig_ES)=treatments

for (i_response in responses) {
  for (treatment in treatments) {
    mean_sig_ES[treatment,i_response]=response_ES_bs[[i_response]][[treatment]]
  }
}

write.csv(mean_sig_ES, "/mean_sig_ES.csv")

mean_sig_ES=read.csv("/Data/mean_sig_ES.csv", header = TRUE)
rownames(mean_sig_ES)<-mean_sig_ES$X
mean_sig_ES=select(mean_sig_ES, 2:8)

# 3.2 Normalization by response standard deviation
standardized_sig_facts=mean_sig_ES%>%
  mutate(PH = PH/sd(mean_sig_ES[,"PH"]),Decom = Decom/sd(mean_sig_ES[,"Decom"]),WSA = WSA/sd(mean_sig_ES[,"WSA"]), actv_ace= actv_ace/sd(mean_sig_ES[,"actv_ace"]), actv_cello= actv_cello/sd(mean_sig_ES[,"actv_cello"]), actv_gluco= actv_gluco/sd(mean_sig_ES[,"actv_gluco"]), actv_phos= actv_phos/sd(mean_sig_ES[,"actv_phos"]),)

standardized_sig_facts=standardized_sig_facts%>%
  mutate(id = treatments)

write.csv(standardized_sig_facts, "/Data/standardized_sig_facts.csv")

#3.3 Building factors traits based on expert opinion (qualitative)
#The factor traits are used from this paper: https://onlinelibrary.wiley.com/doi/10.1111/gcb.15577

df01 <- readxl::read_excel("/Data/Rillig_etal_2021_Classifying human influences.xlsx", sheet = "Traits")

Lithium <-data.frame(matrix(data = NA, ncol = 30, nrow = 1))
Lithium[1,] <-c("Lithium","0","1","0","0","0","0","0","1","0","0","0","1","0","0","1","0","0","0","0","1","1","0","0","0","0","0","0","0","0")
colnames(Lithium) <- colnames(df01)
df01 <- rbind(df01,Lithium)

gP <- c(5,7,8,9,12,13) # plant-related and animal-related variables to be eliminated from the tanglegram analysis 
# because in the soil experiment there was no plant and animal included
Effect_direct<- c(10,11)
Effect_mode <- c(6,16,17,18,19,20)

df_traits_12factors <- df01[c(2,3,6,11,13,14,17,18,19,20,21,31),]

write.csv(df_traits_12factors, "/Data/opinion_based_12factor_traits.csv")

#3.4 Hierarchical clustering the 12 factors

#The original code of hierarchical clustering analyse is written by *Masahiro Ryo*. More information can be found from the appendix of [*Rillig, Ryo, & Lehmann (2021) Classifying human influences on terrestrial ecosystems in Global Change Biology*](https://masahiroryo.github.io/Classifying-human-influences/Rillig_etal_2020_Appendix_Final.html).

#The hierarchical clustering of single factors is based on the euclidean distances. 
df_traits_12factors$Factor <- c("Nit","Drought","MP","Surf","Cu","Salin","PFAS","Fung","AntiB","Herb","Insect","Li")
rownames(df_traits_12factors) <-df_traits_12factors$Factor

hc_traits <- df_traits_12factors %>%
  dist(method = "euclidean") %>%
  hclust(method = 'average')%>%
  as.dendrogram()

hc_experiment<-standardized_sig_facts %>%
  dist(method = 'euclidean')%>%
  hclust(method = 'average')%>%
  as.dendrogram()


library(reshape2)

hmap<-melt(standardized_sig_facts, id = 'id')
P=ggplot(hmap, aes(id, variable, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#649599",
                       high = "#D7AD46",
                       mid = "white",
                       guide = "colorbar")
P


#The two dendrograms are now compared as a tanglegram.

dendlist(hc_experiment, hc_traits)%>%
  untangle(method = "step1side")%>%
  tanglegram(lwd = 0.5,
             columns_width = c(1,0.3,1),
             highlight_distinct_edges = F,
             highlight_branches_lwd = F,
             margin_inner = 15)

# Pre-define of three functions:

#- 1)dendro_data_k
#- 2)set_labels_params
#- 3)plot_ggdendro

# 4. Building distance matrix 

#Our distance matrix of 12 global change factors is built by the experimental data of single factors. The distance of each two factors is calculated by bray-curtis distance.

# 4.1 loading packages

library(vegan)
library(ape)
library(hagis)
library(tidyr)
library(stringr)
library(ecodist)

# 4.2 Reordering factors

#Reorder factors by the label numbers from [master list](https://docs.google.com/spreadsheets/d/1hQJQGpTI8dOpMIeIG_e77gVqqPud2nY0/edit#gid=1499916030).

sig_fac<-read.csv("/Data/standardized_sig_facts.csv")
fac_order<-c(5,2,12,8,9,6,3,10,4,1,11,7) #replace factor names with the label numbers based on the sheet in master list.
sig_fac<-sig_fac%>%
  mutate(order = fac_order)%>%
  arrange(order)%>%
  select(2:8)

# 4.3 Calculating Euclidean distance.

distances<-vegdist(sig_fac,method = "euclidean")

# 4.4 Visualization of Euclidean distance by PCoA.

Eucli_dis_pcoa <- pcoa(distances)
#calculate the percentage of variation that each principal coordinate accounts for.
Axis1_percent<-Eucli_dis_pcoa$values$Relative_eig[1]*100
Axis2_percent<-Eucli_dis_pcoa$values$Relative_eig[2]*100
#Selecting the first two dimension that account for the most variation in the data.
Eucli_dis_pcoa_df<-data.frame(pcoa1 = Eucli_dis_pcoa$vectors[,1],
                             pcoa2 = Eucli_dis_pcoa$vectors[,2])
#attaching metadata
Fac_nr<-c(1,2,3,4,5,6,7,8,9,10,11,12)
Fac_Name<-c("PFAS","Copper","Lithium","Nitrogen","Antibiotic","Insecticide","Surfactants","Fungicide","Herbicide","Microplastic","Salinity","Drought")

Eucli_dis_pcoa_df<-Eucli_dis_pcoa_df%>%
  mutate(Nr = Fac_nr)%>%
  mutate(GCFs = Fac_Name)

#Plotting.
Eucli_dis_plot<- ggplot(data = Eucli_dis_pcoa_df, aes(x=pcoa1, y=pcoa2, color = GCFs))+
  geom_point(aes(size = 5))+
  xlab(paste("PCOA1 - ", round(Axis1_percent, 2), "%", sep = " "))+
  ylab(paste("PCOA2 - ", round(Axis2_percent, 2), "%", sep = " "))+
  theme_bw()
Eucli_dis_plot

# 5. Calculating dissimilarities of multi-factor combinations

# Now we are calculating the dissimilarities of each multi-factor combinations based on the single factor distance matrix.

# Extract data of 2, 5 and 8 factors combinations

alldata<-read.csv("/Data/alldata.csv")
f2<-alldata%>%
  filter(Factor.Num==2)
f5<-alldata%>%
  filter(Factor.Num==5)
f8<-alldata%>%
  filter(Factor.Num==8)

# 5.1 Two-factors combination

# For the two-factor group, the dissimilaritiy of each combination equals to the eucliean distance. 
# Transfer Tube.label into two columns of factor numbers

f2_mut<-f2%>%
  select(5,7:13)%>%
  separate(Tube.label,c("Nr","2f", "lable1", "lable2"))
f2_mut<-f2_mut%>%
  mutate(Label = paste(f2_mut$lable1, f2_mut$lable2))%>%
  select(5:12)%>%
  mutate(f1 = str_extract(Label, "(1)"),
         f2 = str_extract(Label, "(2)"),
         f3 = str_extract(Label, "(3)"),
         f4 = str_extract(Label, "(4)"),
         f5 = str_extract(Label, "(5)"),
         f6 = str_extract(Label, "(6)"),
         f7 = str_extract(Label, "(7)"),
         f8 = str_extract(Label, "(8)"),
         f9 = str_extract(Label, "(9)"),
         f10 = str_extract(Label, "[M]"),
         f11 = str_extract(Label, "[S]"),
         f12 = str_extract(Label, "[D]"),
  )
f2_mut$f10<-str_replace(f2_mut$f10,"[M]","10")
f2_mut$f11<-str_replace(f2_mut$f11,"[S]","11")
f2_mut$f12<-str_replace(f2_mut$f12,"[D]","12")

f2_mut<-f2_mut%>%
  unite("label",9:20,na.rm = TRUE, remove = TRUE)%>%
  select(1:7,9)%>%
  separate(label,c("label_1","label_2"))

# Attaching dissimilarity values

df<-as.matrix(distances)

label_1<-f2_mut$label_1
label_2<-f2_mut$label_2
dissim_2f = vector()

for (i in 1:nrow(f2_mut)) {
  dissim_2f[i]<-df[label_1[i],label_2[i]]
}
f2_mut<-f2_mut%>%
  mutate(dissim = dissim_2f)

write.csv(f2_mut,"/Data/f2_dissimilarity.csv")

# 5.2 Five-factors combination

# For five factor group, the dissimilarity of each five-factor combination is calculated based on the arithmetic mean value of the distances of all factor permutations. 

f5_mut<-f5%>%
  select(5,7:13)%>%
  separate(Tube.label,c("Nr","5f", "lable1", "lable2"))
f5_mut<-f5_mut%>%
  mutate(Label = paste(f5_mut$lable1, f5_mut$lable2))%>%
  select(5:12)%>%
  mutate(f1 = str_extract(Label, "(1)"),
         f2 = str_extract(Label, "(2)"),
         f3 = str_extract(Label, "(3)"),
         f4 = str_extract(Label, "(4)"),
         f5 = str_extract(Label, "(5)"),
         f6 = str_extract(Label, "(6)"),
         f7 = str_extract(Label, "(7)"),
         f8 = str_extract(Label, "(8)"),
         f9 = str_extract(Label, "(9)"),
         f10 = str_extract(Label, "[M]"),
         f11 = str_extract(Label, "[S]"),
         f12 = str_extract(Label, "[D]"),
  )
f5_mut$f10<-str_replace(f5_mut$f10,"[M]","10")
f5_mut$f11<-str_replace(f5_mut$f11,"[S]","11")
f5_mut$f12<-str_replace(f5_mut$f12,"[D]","12")
f5_mut<-f5_mut%>%
  unite("label",9:20,na.rm = TRUE, remove = TRUE)%>%
  select(1:7,9)%>%
  separate(label,c("label_1","label_2","label_3","label_4","label_5"))

# Attaching dissimilarity values

df<-as.matrix(distances)

a5<-f5_mut$label_1
b5<-f5_mut$label_2
c5<-f5_mut$label_3
d5<-f5_mut$label_4
e5<-f5_mut$label_5

X_1 = vector()
X_2 = vector()
X_3 = vector()
X_4 = vector()
X_5 = vector()
X_6 = vector()
X_7 = vector()
X_8 = vector()
X_9 = vector()
X_10 = vector()
Y = vector()

for (i in 1:nrow(f5_mut)) {
  X_1[i]<-df[a5[i],b5[i]]
  X_2[i]<-df[a5[i],c5[i]]
  X_3[i]<-df[a5[i],d5[i]]
  X_4[i]<-df[a5[i],e5[i]]
  X_5[i]<-df[b5[i],c5[i]]
  X_6[i]<-df[b5[i],d5[i]]
  X_7[i]<-df[b5[i],e5[i]]
  X_8[i]<-df[c5[i],e5[i]]
  X_9[i]<-df[c5[i],d5[i]]
  X_10[i]<-df[d5[i],e5[i]]
  Y[i]<-(X_1[i]+X_2[i]+X_3[i]+X_4[i]+X_5[i]+X_6[i]+X_7[i]+X_8[i]+X_9[i]+X_10[i]) #/10
}

f5_mut<-f5_mut%>%
  mutate(dissim = Y)

write.csv(f5_mut,"/Effect_of_multiple_GCF_dissimilarity/Data/f5_dissimilarity.csv")

# 5.3 Eight-factors combination

# For eight factor group, the dissimilarity of each eight-factor combination is calculated based on the arithmetic mean value of the distances of all factor permutations. 

f8_mut<-f8%>%
  select(5,7:13)%>%
  separate(Tube.label,c("Nr","5f", "lable1", "lable2"))
f8_mut<-f8_mut%>%
  mutate(Label = paste(f8_mut$lable1, f8_mut$lable2))%>%
  select(5:12)%>%
  mutate(f1 = str_extract(Label, "(1)"),
         f2 = str_extract(Label, "(2)"),
         f3 = str_extract(Label, "(3)"),
         f4 = str_extract(Label, "(4)"),
         f5 = str_extract(Label, "(5)"),
         f6 = str_extract(Label, "(6)"),
         f7 = str_extract(Label, "(7)"),
         f8 = str_extract(Label, "(8)"),
         f9 = str_extract(Label, "(9)"),
         f10 = str_extract(Label, "[M]"),
         f11 = str_extract(Label, "[S]"),
         f12 = str_extract(Label, "[D]"),
  )
f8_mut$f10<-str_replace(f8_mut$f10,"[M]","10")
f8_mut$f11<-str_replace(f8_mut$f11,"[S]","11")
f8_mut$f12<-str_replace(f8_mut$f12,"[D]","12")
f8_mut<-f8_mut%>%
  unite("label",9:20,na.rm = TRUE, remove = TRUE)%>%
  select(1:7,9)%>%
  separate(label,c("label_1","label_2","label_3","label_4","label_5","label_6","label_7","label_8"))

# Attaching dissimilarity values

df<-as.matrix(distances)

a8<-f8_mut$label_1
b8<-f8_mut$label_2
c8<-f8_mut$label_3
d8<-f8_mut$label_4
e8<-f8_mut$label_5
f8<-f8_mut$label_6
g8<-f8_mut$label_7
h8<-f8_mut$label_8

N<-1:50
X_1 = vector()
X_2 = vector()
X_3 = vector()
X_4 = vector()
X_5 = vector()
X_6 = vector()
X_7 = vector()
X_8 = vector()
X_9 = vector()
X_10 = vector()
X_11 = vector()
X_12 = vector()
X_13 = vector()
X_14 = vector()
X_15 = vector()
X_16 = vector()
X_17 = vector()
X_18 = vector()
X_19 = vector()
X_20 = vector()
X_21 = vector()
X_22 = vector()
X_23 = vector()
X_24 = vector()
X_25 = vector()
X_26 = vector()
X_27 = vector()
X_28 = vector()
Y = vector()
for (i in N) {
  X_1[i]<-df[a8[i],b8[i]]
  X_2[i]<-df[a8[i],c8[i]]
  X_3[i]<-df[a8[i],d8[i]]
  X_4[i]<-df[a8[i],e8[i]]
  X_5[i]<-df[a8[i],f8[i]]
  X_6[i]<-df[a8[i],g8[i]]
  X_7[i]<-df[a8[i],h8[i]]
  X_8[i]<-df[b8[i],c8[i]]
  X_9[i]<-df[b8[i],d8[i]]
  X_10[i]<-df[b8[i],e8[i]]
  X_11[i]<-df[b8[i],f8[i]]
  X_12[i]<-df[b8[i],g8[i]]
  X_13[i]<-df[b8[i],h8[i]]
  X_14[i]<-df[c8[i],d8[i]]
  X_15[i]<-df[c8[i],e8[i]]
  X_16[i]<-df[c8[i],f8[i]]
  X_17[i]<-df[c8[i],g8[i]]
  X_18[i]<-df[c8[i],h8[i]]
  X_19[i]<-df[d8[i],e8[i]]
  X_20[i]<-df[d8[i],f8[i]]
  X_21[i]<-df[d8[i],g8[i]]
  X_22[i]<-df[d8[i],h8[i]]
  X_23[i]<-df[e8[i],f8[i]]
  X_24[i]<-df[e8[i],g8[i]]
  X_25[i]<-df[e8[i],h8[i]]
  X_26[i]<-df[f8[i],g8[i]]
  X_27[i]<-df[f8[i],h8[i]]
  X_28[i]<-df[g8[i],h8[i]]
  Y[i]<-(X_1[i]+X_2[i]+X_3[i]+X_4[i]+X_5[i]+X_6[i]+X_7[i]+X_8[i]+X_9[i]+X_10[i]+X_11[i]+X_12[i]+X_13[i]+X_14[i]+X_15[i]+X_16[i]+X_17[i]+X_18[i]+X_19[i]+X_20[i]+X_21[i]+X_22[i]+X_23[i]+X_24[i]+X_25[i]+X_26[i]+X_27[i]+X_28[i]) #/28
}

f8_mut<-f8_mut%>%
  mutate(dissim = Y)
write.csv(f8_mut,"/Effect_of_multiple_GCF_dissimilarity/Data/f8_dissimilarity.csv")
