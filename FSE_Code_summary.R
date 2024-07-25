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

df<-read.csv("/Effect_of_multiple_GCF_dissimilarity/Data/singal_factor.csv")
treatments = c("AntiB","Cu","Drought","Fung","Herb","Insect","Li","MP","Nit","PFAS","Salin","Surf")
responses = c("actv_ace","actv_cello","actv_gluco","actv_phos","Decom","PH","WSA")

# 3. Clustering of single factors

# 3.1 Calculate mean effect sizes of single factors

# The clustering of 12 global change factors is based on the experimental measurement of the effects of single factors on soil, each factor has 8 replicates. 

Bootstrap_ES_rep = function(response, data=df, target=treatments, n_perm = n_iter){
  resampled = list()
  for (treatment in target) {
    bs = numeric(0)
    population_TR =data[data["Factor.Spec"]==treatment,response]
    population_CK =data[data["Factor.Spec"]=="Control",response]
    size_CK = length(population_CK)
    size_TR = length(population_TR)
    
    for (id in 1:n_perm) {
      k_CK = mean(sample(population_CK, size_CK, replace = T), na.rm = TRUE)
      k_TR = mean(sample(population_TR, size_TR, replace = T), na.rm = TRUE)
      bs = append(bs, k_TR - k_CK)
    }
    resampled[[treatment]] = mean(bs, na.rm = TRUE)
  }
  return(resampled)
}


n_iter =10000
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


mean_sig_ES=read.csv("/mean_sig_ES.csv", header = TRUE)
rownames(mean_sig_ES)<-mean_sig_ES$X
mean_sig_ES=select(mean_sig_ES, 2:8)

# 3.2 Normalization by response standard deviation
standardized_sig_facts=mean_sig_ES%>%
  mutate(PH = PH/sd(mean_sig_ES[,"PH"]),Decom = Decom/sd(mean_sig_ES[,"Decom"]),WSA = WSA/sd(mean_sig_ES[,"WSA"]), actv_ace= actv_ace/sd(mean_sig_ES[,"actv_ace"]), actv_cello= actv_cello/sd(mean_sig_ES[,"actv_cello"]), actv_gluco= actv_gluco/sd(mean_sig_ES[,"actv_gluco"]), actv_phos= actv_phos/sd(mean_sig_ES[,"actv_phos"]),)

standardized_sig_facts=standardized_sig_facts%>%
  mutate(id = treatments)

write.csv(standardized_sig_facts, "/Effect_of_multiple_GCF_dissimilarity/Data/standardized_sig_facts.csv")

#3.3 Building factors traits based on expert opinion (qualitative)
#The factor traits are used from this paper: https://onlinelibrary.wiley.com/doi/10.1111/gcb.15577

df01 <- readxl::read_excel("/Users/mohanbi/Documents/GitHub/FSE-factor-similarity-experiment-/Rillig_etal_2021_Classifying human influences.xlsx", sheet = "Traits")

Lithium <-data.frame(matrix(data = NA, ncol = 30, nrow = 1))
Lithium[1,] <-c("Lithium","0","1","0","0","0","0","0","1","0","0","0","1","0","0","1","0","0","0","0","1","1","0","0","0","0","0","0","0","0")
colnames(Lithium) <- colnames(df01)
df01 <- rbind(df01,Lithium)

gP <- c(5,7,8,9,12,13) # plant-related and animal-related variables to be eliminated from the tanglegram analysis 
# because in the soil experiment there was no plant and animal included
Effect_direct<- c(10,11)
Effect_mode <- c(6,16,17,18,19,20)

df_traits_12factors <- df01[c(2,3,6,11,13,14,17,18,19,20,21,31),]


write.csv(df_traits_12factors, "/Users/mohanbi/Documents/GitHub/FSE-factor-similarity-experiment-/opinion_based_12factor_traits.csv")

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


dendro_data_k <- function(hc, k) {
  
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  
  hcdata
}

set_labels_params <- function(nbLabels,
                              direction = c("tb", "bt", "lr", "rl"),
                              fan       = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}


plot_ggdendro <- function(hcdata,
                          direction   = c("lr", "rl", "tb", "bt"),
                          fan         = FALSE,
                          scale.color = NULL,
                          branch.size = 1,
                          label.size  = 3,
                          nudge.label = 0.01,
                          expand.y    = 0.1) {
  
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)
  
  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  branch.size)
  
  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }
  
# labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle
  
  p <- p +
    geom_text(data        =  label(hcdata),
              aes(x       =  x,
                  y       =  y,
                  label   =  label,
                  colour  =  factor(clust),
                  angle   =  angle),
              vjust       =  labelParams$vjust,
              hjust       =  labelParams$hjust,
              nudge_y     =  ymax * nudge.label,
              size        =  label.size,
              show.legend =  FALSE)
  
  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }
  
  ylim <- -round(ymax * expand.y, 1)
  p    <- p + expand_limits(y = ylim)
  
  p
}

# visualization of the dendrogram


p <- plot_ggdendro(hc_experiment %>%
                     dendro_data_k(3),
                   direction ="lr",
                   label.size = 3.5,
                   branch.size = 0.1,
                   expand.y =3)
p<- p + theme_void()
p

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

sig_fac<-read.csv("/Effect_of_multiple_GCF_dissimilarity/Data/standardized_sig_facts.csv")
fac_order<-c(5,2,12,8,9,6,3,10,4,1,11,7) #replace factor names with the label numbers based on the sheet in master list.
sig_fac<-sig_fac%>%
  mutate(order = fac_order)%>%
  arrange(order)%>%
  select(2:8)

# 4.3 Calculating bray-curtis distance.

distances<-vegdist(sig_fac_1,method = "euclidean")

# 4.4 Visualization of Euclidean distance by PCoA.

bray_dis_pcoa <- pcoa(bray_dis)
#calculate the percentage of variation that each principal coordinate accounts for.
Axis1_percent<-bray_dis_pcoa$values$Relative_eig[1]*100
Axis2_percent<-bray_dis_pcoa$values$Relative_eig[2]*100
#Selecting the first two dimension that account for the most variation in the data.
bray_dis_pcoa_df<-data.frame(pcoa1 = bray_dis_pcoa$vectors[,1],
                             pcoa2 = bray_dis_pcoa$vectors[,2])
#attaching metadata
Fac_nr<-c(1,2,3,4,5,6,7,8,9,10,11,12)
Fac_cluster<-c("Cluster_1","Cluster_2","Cluster_2","Cluster_2","Cluster_2","Cluster_1","Cluster_1","Cluster_2","Cluster_1","Cluster_1","Cluster_3","Cluster_4") #the clustering is based on euclidean distance
Fac_Name<-c("PFAS","Copper","Lithium","Nitrogen","Antibiotic","Insecticide","Surfactants","Fungicide","Herbicide","Microplastic","Salinity","Drought")
#Adding metadata to the factors
#The metadata classification is based on the GCB paper "Classifying human influences on terrestrial ecosystems"
Nature_of_factor<-c("Chemical","Chemical","Chemical","Chemical","Chemical","Chemical","Chemical","Chemical","Chemical","Chemical&physical","Chemical","Chemical")
Effect_mechanism<-c("toxicant","toxicant","toxicant","Resource","toxicant","toxicant","physicalmodifier","toxicant","toxicant","physicalmodifier","osmotic","osmotic")
bray_dis_pcoa_df<-bray_dis_pcoa_df%>%
  mutate(Nr = Fac_nr)%>%
  mutate(GCFs = Fac_Name)%>%
  mutate(cluster = Fac_cluster)%>%
  mutate(Nature_of_factor = Nature_of_factor)%>%
  mutate(Effect_mechanism = Effect_mechanism)
#Plotting.
bray_dis_plot<- ggplot(data = bray_dis_pcoa_df, aes(x=pcoa1, y=pcoa2,color = GCFs))+
  geom_point(aes(size = 5))+
  xlab(paste("PCOA1 - ", round(Axis1_percent, 2), "%", sep = " "))+
  ylab(paste("PCOA2 - ", round(Axis2_percent, 2), "%", sep = " "))+
  theme_bw()

# 5. Calculating dissimilarities of multi-factor combinations

# Now we are calculating the dissimilarities of each multi-factor combinations based on the single factor distance matrix.

# Extract data of 2, 5 and 8 factors combinations

alldata<-read.csv("/Effect_of_multiple_GCF_dissimilarity/Data/alldata.csv")
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

df<-as.matrix(bray_dis)

a<-f2_mut$label_1
b<-f2_mut$label_2

N<-1:47
c = vector()
for (i in N) {
  c[i]<-df[a[i],b[i]]
}
f2_mut<-f2_mut%>%
  mutate(dissim = c)
write.csv(f2_mut,"/Effect_of_multiple_GCF_dissimilarity/Data/f2_dissimilarity.csv")

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

df<-as.matrix(bray_dis)

a5<-f5_mut$label_1
b5<-f5_mut$label_2
c5<-f5_mut$label_3
d5<-f5_mut$label_4
e5<-f5_mut$label_5

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
Y = vector()
for (i in N) {
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

df<-as.matrix(bray_dis)

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


################################################
###.         Part 2 Data preparetion         ###
################################################

# 1. Import data from 2, 5, 8 factor combination groups.

library(tidyverse)
library(ggplot2)

f2=read.csv("/Effect_of_multiple_GCF_dissimilarity/Data/f2_dissimilarity.csv")
f5=read.csv("/Effect_of_multiple_GCF_dissimilarity/Data/f5_dissimilarity.csv")
f8=read.csv("/Effect_of_multiple_GCF_dissimilarity/Data/f8_dissimilarity.csv")

# 2. Transfer factor label number into 0/1 matrix

# 2.1 Transfer factor label number into 0/1 matrix of 2 factor group.

c=c(rep(0,47))
f2=mutate(f2, P=c,C=c,L=c,N=c,A=c,I=c,SU=c,FU=c,H=c,M=c,S=c,D=c)
FA=1:12
N=1:47
Lab=1:2
L=vector()
for (j in Lab) {
  L[j]=f2[j+8]
  for (i in N) {
    for (k in FA) {
      f2[[k+11]][i]=ifelse(L[[j]][i]==k,1,f2[[k+11]][i])
    }
  }
}
r=c(rep(2,47))
f2=f2%>%
  select(2:8,11:23)%>%
  mutate(remark=r,Lv=r)

# 2.2 Transfer factor label number into 0/1 matrix of 5 factor group.

c=c(rep(0,50))
f5=mutate(f5, P=c,C=c,L=c,N=c,A=c,I=c,SU=c,FU=c,H=c,M=c,S=c,D=c)
Lab=1:5
N=1:50
FA=1:12
L=vector()
for (j in Lab) {
  L[j]=f5[j+8]
  for (i in N) {
    for (k in FA) {
      f5[[k+14]][i]=ifelse(L[[j]][i]==k,1,f5[[k+14]][i])
    }
  }
}
r=c(rep(5,50))
f5=f5%>%
  select(2:8,14:26)%>%
  mutate(remark=r,Lv=r)


# 2.3 Transfer factor label number into 0/1 matrix of 8 factor group.

c=c(rep(0,50))
f8=mutate(f8, P=c,C=c,L=c,N=c,A=c,I=c,SU=c,FU=c,H=c,M=c,S=c,D=c)
Lab=1:8
N=1:50
L=vector()
for (j in Lab) {
  L[j]=f8[j+8]
  for (i in N) {
    for (k in FA) {
      f8[[k+17]][i]=ifelse(L[[j]][i]==k,1,f8[[k+17]][i])
    }
  }
}
r=c(rep(8,50))
f8=f8%>%
  select(2:8,17:29)%>%
  mutate(remark=r,Lv=r)

# 3. Normalize dissimilarity data in each group

library(caret)

f2_dis<-as.data.frame(f2$dissim)
colnames(f2_dis)<-c("dissim")
pro_f2<-preProcess(f2_dis, method = c("range"))
f2_dis<-predict(pro_f2, f2_dis)
c<-vector()
c<-f2_dis$dissim
f2<-f2%>%
  mutate(dissim = c)

f5_dis<-as.data.frame(f5$dissim)
colnames(f5_dis)<-c("dissim")
pro_f5<-preProcess(f5_dis, method = c("range"))
f5_dis<-predict(pro_f5, f5_dis)
c<-vector()
c<-f5_dis$dissim
f5<-f5%>%
  mutate(dissim = c)

f8_dis<-as.data.frame(f8$dissim)
colnames(f8_dis)<-c("dissim")
pro_f8<-preProcess(f8_dis, method = c("range"))
f8_dis<-predict(pro_f8, f8_dis)
c<-vector()
c<-f8_dis$dissim
f8<-f8%>%
  mutate(dissim = c)

# 4. Import single factor data

sf=read.csv("/Effect_of_multiple_GCF_dissimilarity/Data/singal_factor.csv")
c=c(rep(0,123))
r=c(rep(1,123))
sf=sf%>%
  mutate(dissim=c,P=c,C=c,L=c,N=c,A=c,I=c,SU=c,FU=c,H=c,M=c,S=c,D=c,remark=r,Lv=r)

# 5. Single factor data, as the group of factor level "1".

N=1:123
for (i in N) {
  sf$P[i]=ifelse(sf$Factor.Spec[i]=="PFAS",1,0)
  sf$C[i]=ifelse(sf$Factor.Spec[i]=="Cu",1,0)
  sf$L[i]=ifelse(sf$Factor.Spec[i]=="Li",1,0)
  sf$N[i]=ifelse(sf$Factor.Spec[i]=="Nit",1,0)
  sf$A[i]=ifelse(sf$Factor.Spec[i]=="AntiB",1,0)
  sf$I[i]=ifelse(sf$Factor.Spec[i]=="Insect",1,0)
  sf$SU[i]=ifelse(sf$Factor.Spec[i]=="Surf",1,0)
  sf$FU[i]=ifelse(sf$Factor.Spec[i]=="Fung",1,0)
  sf$H[i]=ifelse(sf$Factor.Spec[i]=="Herb",1,0)
  sf$M[i]=ifelse(sf$Factor.Spec[i]=="MP",1,0)
  sf$S[i]=ifelse(sf$Factor.Spec[i]=="Salin",1,0)
  sf$D[i]=ifelse(sf$Factor.Spec[i]=="Drought",1,0)
}

sf1=sf
sf1=sf1%>%
  filter(sf1$Factor.Spec=="Control"|sf1$Factor.Spec=="Water Control", .preserve = TRUE)
sf2=setdiff(sf,sf1)
sf2=sf2%>%
  select(6:27)

# Repet single factor data, with specific factor labels.

for (i in N) {
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="Control","CT",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="Water Control","WC",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="PFAS","P",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="Cu","C",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="Li","L",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="Nit","N",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="AntiB","A",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="Insect","I",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="Surf","SU",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="Fung","FU",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="Herb","H",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="MP","M",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="Salin","S",sf$remark[i])
  sf$remark[i]=ifelse(sf$Factor.Spec[i]=="Drought","D",sf$remark[i])
  sf$Lv[i]=ifelse(sf$Factor.Spec[i]=="Water Control"|sf$Factor.Spec[i]=="Control",0,1)
  
}
sf=sf%>%
  select(6:27)

# 6. combine data

trans_data=rbind(sf,sf2,f2,f5,f8)
write.csv(trans_data,"/Effect_of_multiple_GCF_dissimilarity/Data/CombinedData.csv")

#####################################################
###. Part 3 Random forest & Correlation analysis  ###
#####################################################

#############################
###     Packages          ###
#############################

library(ggplot2)
library(ggepi)
library(ggridges)
library(ggpubr)
library(patchwork)
library(party)
library(caret)
library(dplyr)
library(randomForest)

#############################
### Import data           ###
#############################

df=read.csv("/Effect_of_multiple_GCF_dissimilarity/Data/CombinedData.csv")
df=df[df["remark"]!="WC",]
df$remark=factor(df$remark, levels = unique(df$remark))
levels = c("1", "2", "5", "8")
stressors = c("P","C","L","N","A","I","SU","FU","H","M","S","D")
responses = c("actv_ace","actv_cello","actv_gluco","actv_phos","Decom","PH","WSA")
treatment = as.vector(unique(df$remark))
n_iter = 1000 # final calculation 100,000

#############################
### Loading functions     ###
#############################
#Estimating mean and its 95%confidence interval
BootStrap_mean = function(response, data=df, target = treatment, n_perm = n_iter){
  summary = list()
  
  for(treatment in target){
    bs = numeric(0)
    if(treatment=="1") population = data[data$remark%in%stressors, response]
    if(treatment!="1") population = data[data$remark==treatment, response]
    size = length(population)-sum(is.na(population))
    
    for(id in c(1:n_perm)){
      k = mean(sample(population, size, replace = T), na.rm = TRUE)
      bs = append(bs, k)
    }
    summary[[treatment]] = c(quantile(bs, .025,na.rm = TRUE), mean(bs,na.rm = TRUE), quantile(bs, .975,na.rm = TRUE))
    names(summary[[treatment]]) = c("2.5%", "mean", "97.5%")
  }
  summary = t(data.frame(summary))
  summary = data.frame("target" = target, summary); row.names(summary) = c()
  return(summary)
}

BootStrap_ES_rep = function(response, data=df, target = treatment, n_perm = n_iter){
  resampled = list()
  
  population_CT = data[data$remark=="CT", response]
  
  for(treatment in target){
    bs = numeric(0)
    if(treatment=="1") population_TR = data[data$remark%in%stressors, response]
    if(treatment!="1") population_TR = data[data$remark==treatment, response]
    size_CT = length(population_CT)-sum(is.na(population_CT))
    size_TR = length(population_TR)-sum(is.na(population_TR))
    
    for(id in c(1:n_perm)){
      k_CT = mean(sample(population_CT, size_CT, replace = T), na.rm = TRUE)
      k_TR = mean(sample(population_TR, size_TR, replace = T), na.rm = TRUE)
      bs = append(bs, k_TR - k_CT)
    }
    resampled[[treatment]] = bs
  }
  resampled[["CT"]] = rep(0, n_perm)
  return(resampled)
}

BootStrap_ES_summary = function(data){
  summary = list()
  p = 0
  summary[["CT"]] = c(0,0,0,1)
  target = names(data)
  
  for(treatment in target[-1]){
    bs = data[[treatment]]
    p = length(which(bs>0))/length(bs)
    p = min(p, 1-p)
    summary[[treatment]] = c(quantile(bs, .025,na.rm = TRUE), mean(bs,na.rm = TRUE), quantile(bs, .975,na.rm = TRUE), p)
  }
  summary = t(data.frame(summary))
  colnames(summary) = c("2.5%", "mean", "97.5%", "p_value")
  summary = data.frame(target, summary); row.names(summary) = c()
  
  return(summary)
}

Null_distribution_rep = function(response, data=df, n_perm=n_iter){
  
  output = list()
  for(Lv in levels){
    
    resampled = list()
    
    # Checking which stressor combinations were jointly tested
    if(Lv=="1") combination = data[data$remark%in%stressors,c(10:21)]
    if(Lv!="1") combination = data[data["remark"]==Lv,c(10:21)]
    Level = sum(combination[1,]) 
    
    
    # Null distributions can be taken based on three different assumptions
    for(type in c("Additive", "Multiplicative", "Dominative")){
      
      population_CT = df[df$remark=="CT", response]
      size_CT = length(population_CT)-sum(is.na(population_CT)) ##subtract NA value
      
      # For each combination, bootstrap resampling is conducted
      for(j in c(1:nrow(combination))){
        bs = numeric(0)
        selected_stressors = stressors[which(combination[j,]==1)]
        sub_n_perm = ceiling(n_perm/nrow(combination))*5 
        
        # bootstrap resampling
        for(id in c(1:sub_n_perm)){
          each_effect = numeric(0)
          k_CT = mean(sample(population_CT, size_CT, replace = T),na.rm = TRUE) 
          
          for(treatment in selected_stressors){
            population_TR = df[df$remark==treatment, response]
            size_TR = length(population_TR)-sum(is.na(population_TR))
            k_TR = mean(sample(population_TR, size_TR, replace = T),na.rm = TRUE)
            
            # ES estimate depending on the type of null hypotheses
            if(type=="Additive")       each_effect = append(each_effect, (k_TR - k_CT))
            if(type=="Multiplicative") each_effect = append(each_effect, (k_TR - k_CT)/k_CT)
            if(type=="Dominative")      each_effect = append(each_effect, (k_TR - k_CT))
          }
          
          # Calculating an expected ES after collecting the ESs of all relevant single stressors
          if(type=="Additive")       joint_effect = sum(each_effect)
          if(type=="Multiplicative"){
            z = 1
            for(m in c(1:Level)) z = z * (1 + each_effect[m])
            joint_effect = (z - 1)*k_CT
          }
          if(type=="Dominative")      joint_effect = each_effect[which(max(abs(each_effect))==abs(each_effect))]
          
          bs = append(bs, joint_effect)
        }
        resampled[[type]][[j]] = bs
      }
      
    }
    output[[Lv]] = resampled
  }  
  return(output)
}



Null_distribution_rep_transform = function(data){
  output = list()
  for(Lv in levels){
    for(type in c("Additive", "Multiplicative", "Dominative")){
      output[[Lv]][[type]] = sample(unlist(data[[Lv]][[type]]), n_iter, replace=F)
    }
  }
  return(output)
}

NHST_summary = function(null_data, Actual_data){
  output = list()
  for(Lv in levels){
    summary = list()
    summary[["Actual"]] = c(quantile(Actual_data[[Lv]], .025,na.rm = TRUE), mean(Actual_data[[Lv]],na.rm = TRUE), quantile(Actual_data[[Lv]], .975,na.rm = TRUE), 1)
    p = 0
    assumptions = c("Additive", "Multiplicative", "Dominative")
    
    for(i_assumption in assumptions){
      bs   = (Actual_data[[Lv]] - null_data[[Lv]][[i_assumption]])
      p = length(which(bs>0))/length(bs)
      p = min(p, 1-p)
      summary[[i_assumption]] = c(quantile(null_data[[Lv]][[i_assumption]], .025,na.rm = TRUE), mean(null_data[[Lv]][[i_assumption]],na.rm = TRUE), quantile(null_data[[Lv]][[i_assumption]], .975,na.rm = TRUE), p)
    }
    summary = t(data.frame(summary))
    colnames(summary) = c("2.5%", "mean", "97.5%", "p_value")
    summary = data.frame(ES = c("Actual", "Additive","Multiplicative","Dominative"), summary); row.names(summary) = c()
    
    output[[Lv]] = summary
  }
  
  return(output)
}

NHST_summary_transform = function(data){
  output = list()
  for(i in 1:4){
    summary = rbind(data[["1"]][i, 2:4], data[["2"]][i, 2:4], data[["5"]][i, 2:4],
                    data[["8"]][i, 2:4])
    summary = cbind(levels, summary)
    colnames(summary) = c("Lv", "Low", "Mean", "High")
    output[[c("Actual", "Additive", "Multiplicative", "Dominative")[i]]] = summary
  }
  return(output)
}

Expected_response_for_each = function(data){
  output = numeric(0)
  for(type in c("Additive", "Multiplicative", "Dominative")){
    tmp = numeric(0)
    for(Lv in levels){
      n_len = length(data[[Lv]][[type]])
      for(i in 1:n_len){
        tmp = append(tmp, mean(data[[Lv]][[type]][[i]])+response_mean[1,"mean"])
      }
    }
    output = cbind(output,tmp)
  }
  colnames(output)= c("R1", "R2", "R3")
  return(output)
}


postResample <- function(pred, obs)
{
  isNA <- is.na(pred)
  pred <- pred[!isNA]
  obs <- obs[!isNA]
  if (!is.factor(obs) && is.numeric(obs))
  {
    if(length(obs) + length(pred) == 0)
    {
      out <- rep(NA, 3)
    } else {
      if(length(unique(pred)) < 2 || length(unique(obs)) < 2)
      {
        resamplCor <- NA
      } else {
        resamplCor <- try(cor(pred, obs, use = "pairwise.complete.obs"), silent = TRUE)
        if (inherits(resamplCor, "try-error")) resamplCor <- NA
      }
      mse <- mean((pred - obs)^2)
      mae <- mean(abs(pred - obs))
      out <- c(sqrt(mse), resamplCor^2, mae)
    }
    names(out) <- c("RMSE", "Rsquared", "MAE")
  } else {
    if(length(obs) + length(pred) == 0)
    {
      out <- rep(NA, 2)
    } else {
      pred <- factor(pred, levels = levels(obs))
      requireNamespaceQuietStop("e1071")
      out <- unlist(e1071::classAgreement(table(obs, pred)))[c("diag", "kappa")]
    }
    names(out) <- c("Accuracy", "Kappa")
  }
  if(any(is.nan(out))) out[is.nan(out)] <- NA
  out
}

###################################
###     Null model analysis     ###
###################################

response_mean_all = list()
response_ES_all   = list()
joint_ES_null_all = list()
ES_for_each_all   = list()
g_rawdata_all     = list()
g_ms_all          = list()
Pred_respon_for_each_all = list()
###Correlation figures

g_2fac_cor = list()
g_5fac_cor = list()
g_8fac_cor = list()


for (i_response in responses) {
  
  #bootstrap mean
  response_mean = BootStrap_mean(i_response)
  response_mean_all[[i_response]] = response_mean
  #effect size
  response_ES_bs  = BootStrap_ES_rep(i_response)
  response_ES     = BootStrap_ES_summary(response_ES_bs)
  response_ES_all[[i_response]]   = response_ES
  #Null model
  Null_ES_bs0     = Null_distribution_rep(i_response)
  Null_ES_bs      = Null_distribution_rep_transform(Null_ES_bs0)
  
  joint_ES_null   = NHST_summary(Null_ES_bs, response_ES_bs)
  joint_ES_null_all[[i_response]] = bind_rows(joint_ES_null, .id = "column_label")
  
  ES_plot         = NHST_summary_transform(joint_ES_null)
  ES_for_each     = Expected_response_for_each(Null_ES_bs0)
  Pred_respon_for_each   = Expected_response_for_each(Null_ES_bs0)
  ES_for_each_all[[i_response]]   = ES_for_each
  Pred_respon_for_each_all[[i_response]] = Pred_respon_for_each
  
  g_rawdata_all[[i_response]] =  local({
    i_response = i_response
    response_mean = response_mean
    Mycolor_TR=c("#616161A0","#616161A0","#616161A0","#616161A0","#616161A0","#616161A0","#616161A0","#616161A0","#616161A0","#616161A0","#616161A0","#616161A0","#616161A0","#616161A0","#495B46A0","#A9796AA0","#BAA43FA0")
    Mycolor=c("#616161","#616161","#616161","#616161","#616161","#616161","#616161","#616161","#616161","#616161","#616161","#616161","#616161","#616161","#495B46","#A9796A","#BAA43F")
    ggplot() +
      theme_bw()+
      theme(legend.position = 'none', axis.title.x=element_blank())+
      xlab(i_response) + 
      coord_flip() +
      stat_density_ridges(data=df, aes_string(x = i_response, y = "remark",color= "remark", fill= "remark"),
                          geom = "density_ridges_gradient",
                          rel_min_height = 0.01, 
                          jittered_points = TRUE,
                          position = position_points_jitter(height = .2, yoffset = .15),
                          point_size = 1, point_alpha = .3, 
                          scale = 0.5) +
      scale_fill_manual(values  = Mycolor_TR)+
      geom_estci(data=response_mean, aes(x = mean, y = target, xmin=X2.5., xmax=X97.5., 
                                         xintercept=response_mean[1,"mean"], color =target), center.linecolour = "black",
                 size=0.6, ci.linesize = 0.5, position=position_nudge(y = -0.15))+
      scale_color_manual(values  = Mycolor)+
      theme(panel.background = element_rect(fill = "#4D728530", color = "white"),
            panel.grid.major = element_line(color = "white", size = 0.1),
            panel.grid.minor = element_line(color = "white", size = 0.5))
    
  })
  
  g_2fac_cor[[i_response]] =  local({
    
    ggplot(df[df[,"Lv"]=="2",], aes_string(x = "dissim", y = i_response))+
      geom_point(color = "#495B46", alpha = .7, size = 2)+
      geom_smooth(formula = 'y ~ x', method = lm, fill = "#495B46", color = "#495B46", alpha = .2)+
      theme_bw()+
      stat_cor(method = "spearman",label.x = 0)+
      theme(legend.position = 'none', axis.title.x=element_blank(), axis.title.y =element_blank())+
      theme(panel.background = element_rect(fill = "#4D728530", color = "white"),
            panel.grid.major = element_line(color = "white", size = 0.1),
            panel.grid.minor = element_line(color = "white", size = 0.5))
    
  })
  
  g_5fac_cor[[i_response]] =  local({
    
    ggplot(df[df[,"Lv"]=="5",], aes_string(x = "dissim", y = i_response))+
      geom_point(color = "#A9796A", alpha = .7, size = 2)+
      geom_smooth(formula = 'y ~ x', method = lm, fill = "#A9796A", color = "#A9796A", alpha = .2)+
      theme_bw()+
      stat_cor(method = "spearman",label.x = 0)+
      theme(legend.position = 'none', axis.title.x=element_blank(), axis.title.y =element_blank())+
      theme(panel.background = element_rect(fill = "#4D728530", color = "white"),
            panel.grid.major = element_line(color = "white", size = 0.1),
            panel.grid.minor = element_line(color = "white", size = 0.5))
    
  })
  
  g_8fac_cor[[i_response]] =  local({
    
    ggplot(df[df[,"Lv"]=="8",], aes_string(x = "dissim", y = i_response))+
      geom_point(color = "#BAA43F", alpha = .7, size = 2)+
      geom_smooth(formula = 'y ~ x', method = lm, fill = "#BAA43F", color = "#BAA43F", alpha = .2)+
      theme_bw()+
      stat_cor(method = "spearman",label.x = 0)+
      theme(legend.position = 'none', axis.title.x=element_blank(), axis.title.y =element_blank())+
      theme(panel.background = element_rect(fill = "#4D728530", color = "white"),
            panel.grid.major = element_line(color = "white", size = 0.1),
            panel.grid.minor = element_line(color = "white", size = 0.5))
    
  })
  
}


####################################
###     Random forest analysis   ###
####################################
df<- na.roughfix(df) 
df.rf = df[df[, "remark"] %in% levels,]  

lv_list  = unique(df.rf[,"Lv"])
id_lv1   = which(df.rf[,"Lv"]==1)
id_lvh   = which(df.rf[,"Lv"]>1)

n_data   = nrow(df.rf)
n_lv1    = sum(df.rf[,"Lv"]==1)
n_lvh    = sum(df.rf[,"Lv"]!=1)
n_eachlv = 50 # 50 for balanced resampling. Our experiment has 50 replicates for high levels.
n_tree   = 100
n_iter2  = 100

rf.r2 = data.frame(matrix(NA, ncol=3, nrow=7*n_iter2*length(responses)))
rf.r2[,1] = rep(responses,each=7*n_iter2)
rf.r2[,2] = rep(c("Null","Lv","Dis","Null+Lv","Null+Dis","Null+Lv+Dis","All"), n_iter2*length(responses))
rf.r2[,2] = factor(rf.r2[,2], levels=c("Null","Lv","Dis","Null+Lv","Null+Dis","Null+Lv+Dis","All"))
colnames(rf.r2) = c("Response", "Model", "R2")

rf.pred = data.frame(matrix(NA, ncol=5, nrow=n_iter2*length(responses)))
rf.pred[,1] = rep(responses,each=n_iter2)
colnames(rf.pred) = c("Response", lv_list)

rf.vimp = data.frame(matrix(NA, ncol=18, nrow=n_iter2*length(responses)))
rf.vimp[,1] = rep(responses,each=n_iter2) 
colnames(rf.vimp) = c("Response", "Lv", "R1","R2","R3","dissim", stressors)

rf.prediction.all = list()

j = 0
for(i_response in responses){
  
  df.rf.tmp = cbind(df.rf, Pred_respon_for_each_all[[i_response]])
  #define the formulas for three random forest models
  eval(parse(text=(paste("fml      = formula(",i_response,"~", paste(c("Lv","R1","R2","R3","dissim", stressors), collapse=" + "),")", sep=""))))
  eval(parse(text=(paste("fml.Null = formula(",i_response,"~", paste(c("R1","R2","R3"), collapse=" + "),")", sep=""))))
  eval(parse(text=(paste("fml.Lv   = formula(",i_response,"~Lv)", sep=""))))
  eval(parse(text=(paste("fml.Dis   = formula(",i_response,"~dissim)", sep=""))))
  eval(parse(text=(paste("fml.NullLv      = formula(",i_response,"~", paste(c("R1","R2","R3","Lv"), collapse=" + "),")", sep=""))))
  eval(parse(text=(paste("fml.NullDis      = formula(",i_response,"~", paste(c("R1","R2","R3","dissim"), collapse=" + "),")", sep=""))))
  eval(parse(text=(paste("fml.NullLvDis      = formula(",i_response,"~", paste(c("R1","R2","R3","Lv","dissim"), collapse=" + "),")", sep=""))))
  
  for(i in 1:n_iter2){
    j = j + 1
    # bootstrap resampling
    set.seed(j)
    # take 50 sample from Lv1, and take 150 sample from the other levels for balanced resampling
    rid = c(sample(id_lv1, n_eachlv, replace=T), sample(id_lvh, n_lvh, replace=T))
    rdf = df.rf.tmp[rid,]
    rf_model      = tryCatch({
      cforest(fml,      data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))},
      error = function(e) {cforest(fml.Lv,   data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))})
    rf_model.Null   = cforest(fml.Null,   data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
    rf_model.Lv   = cforest(fml.Lv,   data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
    rf_model.Dis   = cforest(fml.Dis,   data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
    rf_model.NullLv   = cforest(fml.NullLv,   data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
    rf_model.NullDis   = cforest(fml.NullDis,   data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
    rf_model.NullLvDis   = cforest(fml.NullLvDis,   data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
    
    # shuffling corresponding variables for evaluating R2 reduction
    # rdf.test.LvID = rdf
    # rdf.test.LvID[,colnames(ES_for_each_all[[i_response]])] = apply(rdf.test.LvID[,colnames(ES_for_each_all[[i_response]])],2,sample)  
    
    # evaluating fitting performance
    
    rf_prdct.Null   = tryCatch({predict(rf_model.Null, OOB=T)}, error=function(e){rep(0,length(rid))})
    rf_prdct.Lv   = tryCatch({predict(rf_model.Lv, OOB=T)}, error=function(e){rep(0,length(rid))})
    rf_prdct.Dis   = tryCatch({predict(rf_model.Dis, OOB=T)}, error=function(e){rep(0,length(rid))})
    rf_prdct.NullLv = tryCatch({predict(rf_model.NullLv, OOB=T)}, error=function(e){rep(0,length(rid))})
    rf_prdct.NullDis = tryCatch({predict(rf_model.NullDis, OOB=T)}, error=function(e){rep(0,length(rid))})
    rf_prdct.NullLvDis = tryCatch({predict(rf_model.NullLvDis, OOB=T)}, error=function(e){rep(0,length(rid))})
    rf_prdct.All  = tryCatch({predict(rf_model, OOB=T)}, error=function(e){rep(0,length(rid))})
    
    rf.r2[7*j-6,3]    = postResample(rf_prdct.Null,rdf[,i_response])[2]
    rf.r2[7*j-5,3]    = postResample(rf_prdct.Lv,rdf[,i_response])[2]
    rf.r2[7*j-4,3]    = postResample(rf_prdct.Dis,rdf[,i_response])[2]
    rf.r2[7*j-3,3]    = postResample(rf_prdct.NullLv,rdf[,i_response])[2]
    rf.r2[7*j-2,3]    = postResample(rf_prdct.NullDis,rdf[,i_response])[2]
    rf.r2[7*j-1,3]    = postResample(rf_prdct.NullLvDis,rdf[,i_response])[2]
    rf.r2[7*j-0,3]    = postResample(rf_prdct.All,rdf[,i_response])[2]
    
    # variable importance
    set.seed(j)
    tmp.vimp = numeric(0)
    for(itmp in 1:5) tmp.vimp = rbind(tmp.vimp,varimp(rf_model))
    tmp.vimp = apply(tmp.vimp,2,mean)
    tmp.vimp = tmp.vimp/sum(tmp.vimp)
    rf.vimp[j, 2:(length(tmp.vimp)+1)] = tmp.vimp*rf.r2[2*j-0,3]*100
    
    # fitting curve
    tmp.curve = c()
    for(i_lv in lv_list){
      tmp.curve = append(tmp.curve, mean(rf_prdct.All[rdf[,"Lv"]==i_lv]))
    }
    rf.pred[j,2:5]  =  tmp.curve
  }
  
  rf_model      = cforest(fml,      data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  rf_model.Null   = cforest(fml.Null,   data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  rf_model.Lv   = cforest(fml.Lv,   data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  rf_model.Dis   = cforest(fml.Dis,   data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  rf_model.NullLv   = cforest(fml.NullLv,   data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  rf_model.NullDis   = cforest(fml.NullDis,   data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  rf_model.NullLvDis = cforest(fml.NullLvDis, data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  
  rf.prediction.all[[i_response]] = data.frame(
    Model     = rep(c("Null","Lv","Dis","Null+Lv","Null+Dis","Null+Lv+Dis","All"), each = nrow(df.rf.tmp)),
    Predicted = c(predict(rf_model.Null, OOB=F),predict(rf_model.Lv, OOB=F),predict(rf_model.Dis, OOB=F),predict(rf_model.NullLv, OOB=F),predict(rf_model.NullDis, OOB=F),predict(rf_model.NullLvDis, OOB=F),predict(rf_model, OOB=F)),
    Observed  = rep(df.rf.tmp[,i_response], 7))
  
}
rf.r2.summary = data.frame(matrix(NA,ncol=5,nrow=7*length(responses)))
colnames(rf.r2.summary) = c("Response","Model", "CI.low", "Mean", "CI.high")
rf.r2.summary[,1] = rep(responses,each=7)
rf.r2.summary[,2] = rep(c("Null","Lv","Dis","Null+Lv","Null+Dis","Null+Lv+Dis","All"),length(responses))
rf.r2.summary[,2] = factor(rf.r2.summary[,2],levels=c("Null","Lv","Dis","Null+Lv","Null+Dis","Null+Lv+Dis","All"))


rf.pred.summary = data.frame(matrix(NA,ncol=5,nrow=4*length(responses)))
colnames(rf.pred.summary) = c("Response","Lv","CI.low", "Mean", "CI.high")
rf.pred.summary[,1] = rep(responses,each=length(lv_list))
rf.pred.summary[,2] = rep(lv_list, length(responses))

rf.vimp.summary = data.frame(matrix(NA,ncol=5,nrow=17*length(responses)))
colnames(rf.vimp.summary) = c("Response","Variable","CI.low", "Mean", "CI.high")
rf.vimp.summary[,1] = rep(responses,each=17)
rf.vimp.summary[,2] = rep(colnames(rf.vimp)[2:18], length(responses))
rf.vimp.summary[,2] = factor(rf.vimp.summary[,2], levels = unique(rf.vimp.summary[,2]))

j  = 0
for(i_response in responses){
  jj = 0
  jjj = 0
  rf.r2.summary[7*j+1,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="Null"),   3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[7*j+2,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="Lv"),   3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[7*j+3,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="Dis"),3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[7*j+4,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="Null+Lv"),  3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[7*j+5,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="Null+Dis"),   3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[7*j+6,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="Null+Lv+Dis"),3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[7*j+7,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="All"),  3],c(.025,.50,.975), na.rm=T)
  
  for(i_lv in lv_list){
    jj = jj + 1
    rf.pred.summary[4*j+jj,3:5]     = quantile(rf.pred[which(rf.pred[,1]==i_response),as.character(i_lv)],c(.025,.50,.975), na.rm=T)
  }#change 5-4*j+jj,3:5]
  
  for(i_var in colnames(rf.vimp)[2:18]){
    jjj = jjj + 1
    rf.vimp.summary[17*j+jjj,3:5]     = quantile(rf.vimp[which(rf.vimp[,1]==i_response),as.character(i_var)],c(.025,.50,.975), na.rm=T)
  }
  
  j = j + 1
  
}

#write.csv(rf.r2.summary,"/rf.r2.summary.csv")
######################################
### Plotting random forest figures ###
######################################

g_rf_all     = list()
g_vimp_all   = list()
g_r2_all     = list()
g_r2_four     = list()
g_correl_all = list()
level_order_four=c("Null","Null+Lv","Null+Lv+Dis","All")
rf.r2_four = rf.r2[rf.r2[,2]%in%level_order_four,]
rf.r2.summary_four = rf.r2.summary[rf.r2.summary[,2]%in%level_order_four,]

for(i in 1:(length(responses))){
  
  g_vimp_all[[i]] = local({
    i = i
    ggplot(data=rf.vimp.summary[rf.vimp.summary[,1]==responses[i],],aes(x=factor(Variable, level = c("Lv","R1","R2","R3","dissim","P","C","L","N","A","I","SU","FU","H","M","S","D")), y=Mean))+
      geom_bar(stat="identity", fill="#999999", alpha=0.5) +
      xlab("Variability explained [%]")+
      geom_errorbar(aes(ymax=CI.high, ymin=CI.low), width=.2, position=position_dodge(width=0.0)) +
      theme_bw()+
      theme(legend.position = 'none', axis.title.x=element_blank(),axis.title.y=element_blank())
  })
  
  g_r2_all[[i]] = local({
    i = i
    rf.r2 = rf.r2
    rf.r2.summary = rf.r2.summary 
    level_order=c("Null","Lv","Dis","Null+Lv","Null+Dis","Null+Lv+Dis","All")
    ggplot(data=rf.r2[rf.r2[,1]==responses[i],],aes(x=factor(Model,level = level_order),y=R2, fill=Model))+
      geom_violin( color="#00000000",alpha=.7,position=position_dodge(width=0.3),trim=F)+
      geom_pointrange(data=rf.r2.summary[rf.r2.summary[,1]==responses[i],], aes(y=Mean, ymax=CI.high, ymin=CI.low,color=Model),fatten = 4, position=position_dodge(width=0.3)) +
      scale_fill_manual(values  = c("#95A5A6","#95A5A6","#95A5A6","#7F8D9D","#7F8D9D","#34495E","#BDC3C7"))+ 
      scale_color_manual(values = c("#95A5A6","#95A5A6","#95A5A6","#7F8D9D","#7F8D9D","#34495E","#BDC3C7"))+ 
      theme_bw() + ylim(c(0,1.0))+
      theme(axis.title.y=element_text(size=8))+
      ylab('Variability explained (R2%)')+
      theme(legend.position = 'none', axis.title.x=element_blank())+
      theme(panel.background = element_rect(fill = "#367DAD30", color = "white"),
            panel.grid.major = element_line(color = "white", size = 1),
            panel.grid.minor = element_line(color = "white", size = 1))
  })
  
  g_r2_four[[i]] = local({
    i = i
    rf.r2_four = rf.r2_four
    rf.r2.summary = rf.r2.summary 
    level_order_four=c("Null","Null+Lv","Null+Lv+Dis","All")
    ggplot(data=rf.r2_four[rf.r2_four[,1]==responses[i],],aes(x=factor(Model,level = level_order_four),y=R2, fill=Model))+
      geom_violin( color="#00000000",alpha=.7,position=position_dodge(width=0.3),trim=F)+
      geom_pointrange(data=rf.r2.summary_four[rf.r2.summary_four[,1]==responses[i],], aes(y=Mean, ymax=CI.high, ymin=CI.low,color=Model),fatten = 4, position=position_dodge(width=0.3)) +
      scale_fill_manual(values  = c("#95A5A6","#7F8D9D","#34495E","#BDC3C7"))+ 
      scale_color_manual(values = c("#95A5A6","#7F8D9D","#34495E","#BDC3C7"))+ 
      theme_bw() + ylim(c(0,1.0))+
      theme(axis.title.y=element_text(size=8))+
      ylab('Variability explained (R2%)')+
      theme(legend.position = 'none', axis.title.x=element_blank())+
      theme(panel.background = element_rect(fill = "#367DAD30", color = "white"),
            panel.grid.major = element_line(color = "white", size = 1),
            panel.grid.minor = element_line(color = "white", size = 1))
  })
  
  g_correl_all[[i]] = local({
    i = i
    rf.prediction.all = rf.prediction.all
    level_order=c("Null","Lv","Dis","Null+Lv","Null+Dis","Null+Lv+Dis","All")
    ggplot(data = rf.prediction.all[[i]], aes(y=Predicted, x=Observed, group=factor(Model,level = level_order), color=factor(Model,level = level_order)))+
      theme_bw()+
      ylab(responses[i]) + 
      geom_abline(slope=1, intercept=0)+
      geom_point()+
      geom_smooth(formula = 'y ~ x',method="lm", fullrange=F,size= 1)+
      xlim(range(rf.prediction.all[[i]][,2:3]))+
      ylim(range(rf.prediction.all[[i]][,2:3]))+
      theme(legend.position = 'none', axis.title.x=element_blank(),axis.title.y=element_blank()) +
      scale_color_manual(values = c("#34988870","#34988870","#34988870","#7F800070","#7F800070","#95000070","#34495E70"))
  }) 
  
}

#############################################
###     Permutation based Random Forest   ###
#############################################

###########################################################################
#
# To define the two functions for permutation-based random forest
#   1. myvarimp
#   2. RF_Hapfelmeier
#
#   INPUT: none
#   OUTPUT: none (This script is called by the main script: "SML_application.R")
#
###########################################################################

package.list = c("party", "caret","foreach","doSNOW")
tmp.install = which(lapply(package.list, require, character.only = TRUE)==FALSE)
if(length(tmp.install)>0) install.packages(package.list[tmp.install])
lapply(package.list, require, character.only = TRUE)

###-------------------------------------------------------------------------
### 1. myvarimp
###    function for calculating variable importance measure (modified by Masahiro Ryo)
###-------------------------------------------------------------------------
myvarimp = function(object, mincriterion = 0, conditional = FALSE, 
                    pre1.0_0 = conditional, varID) {
  response = object@responses
  input = object@data@get("input")
  xnames = colnames(input)
  inp = initVariableFrame(input, trafo = NULL)
  y = object@responses@variables[[1]]
  if (length(response@variables) != 1) stop("cannot compute variable importance measure")
  CLASS = all(response@is_nominal)
  ORDERED = all(response@is_ordinal)
  if (CLASS) {
    error = function(x, oob) mean((levels(y)[sapply(x, which.max)] != y)[oob])
  }else {
    if (ORDERED) {
      error = function(x, oob) mean((sapply(x, which.max) != y)[oob])
    } else {
      error = function(x, oob) mean((unlist(x) - y)[oob]^2)}
  }
  perror = rep(0, length(object@ensemble))
  for (b in 1:length(object@ensemble)) {
    tree = object@ensemble[[b]]
    w = object@weights[[b]]
    w[w == 0] = 1
    oob = object@weights[[b]] == 0
    p = .Call("R_predict", tree, inp, mincriterion, -1L, PACKAGE = "party")
    eoob = error(p, oob)
    j = varID
    p = .Call("R_predict", tree, inp, mincriterion, as.integer(j), PACKAGE = "party")
    perror[(b - 1)] = (error(p,oob) - eoob)
  }
  return(MeanDecreaseAccuracy = mean(perror))
}
environment(myvarimp) = environment(varimp)

###-------------------------------------------------------------------------
### 2. Permutation-based random forest function (modified by Masahiro Ryo)
###   
###-------------------------------------------------------------------------
RF_permutation = function(formula, data, nperm = 500, ntree = 100, ncore=4, alpha=0.05)
  # formula: object of class "formula".
  # data: data frame containing the variables.
  # nperm: Number of permutation steps used for the permutation test.
  # ntree: Number of trees in the Random Forest.
  # ncore: Number of cores used for parallel computing
{
  x.names = all.vars(formula)[-1]
  y.names = all.vars(formula)[1]
  terms. = terms.formula(formula)
  x.formula = attr(terms., "term.labels")
  y.formula = as.character(attr(terms., "variables"))[2]
  mtry = ceiling(sqrt(length(x.formula)))
  dat = subset(data, select = c(y.names, x.names))
  forest = party::cforest(formula, data = dat, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
  obs.varimp = varimp(forest)
  perm.mat = matrix(NA, ncol = length(x.names), nrow = nperm, dimnames = list(1:nperm, x.names))
  
  cl=makeCluster(ncore) #change to your number of CPU cores
  registerDoSNOW(cl)
  
  for (j in x.names) {
    cat("\r", "Processing variable ", which(j == x.names), " of ", length(x.names)); flush.console()
    perm.dat = dat
    perm.mat[, j] = unlist(foreach(i = 1:nperm, .packages='party',.export="myvarimp") %dopar% {
      perm.dat[, j] = sample(perm.dat[, j]);#variable j is permuted to destroy its relation to response
      myvarimp(party::cforest(formula, data = perm.dat, controls = cforest_unbiased(mtry = mtry, ntree = ntree)), varID = which(x.names == j))
    })
  }
  stopCluster(cl)
  p.vals = sapply(x.names, function(x) sum(perm.mat[, x] >= obs.varimp[which(x == x.names)]) / nperm) #Assess the p-value for each variable by means of the empirical distributions and the original importance measures.
  p.vals.bonf = p.vals * length(p.vals)
  
  if (any(p.vals < alpha)) {
    selection = names(p.vals)[which(p.vals < alpha)]
    mtry = ceiling(sqrt(length(selection)))
    forest = cforest(as.formula(paste(y.formula, ".", sep = " ~ ")), data = subset(dat, select = c(y.names, selection)),
                     controls = cforest_unbiased(mtry = mtry, ntree = ntree))
    p = p.vals[which(p.vals < alpha)]
  }
  if (any(p.vals.bonf < alpha)) {             
    selection.bonf = names(p.vals.bonf)[which(p.vals.bonf < alpha)]                         
    mtry = ceiling(sqrt(length(selection.bonf)))
    forest.bonf = party::cforest(as.formula(paste(y.formula, ".", sep = " ~ ")), data = subset(dat, select = c(y.names, selection.bonf)),
                                 controls = cforest_unbiased(mtry = mtry, ntree = ntree))
    p.bonf = p.vals.bonf[which(p.vals.bonf < alpha)]
    varimp.bonf = varimp(forest.bonf)
    accuracy.fitting = postResample(predict(forest.bonf), subset(dat, select = y.names))
    accuracy.validation = caret:::cforestStats(forest.bonf)
    varimp.R2 = accuracy.fitting[2]*varimp.bonf/sum(varimp.bonf)
    residual = subset(dat, select = y.names) - predict(forest.bonf)
  }
  if (!any(p.vals < alpha)) {
    selection = c(); forest = c(); p = c()
  }
  if (!any(p.vals.bonf < alpha)) {
    selection.bonf = c(); forest.bonf = c(); p.bonf = c(); varimp.bonf = c();
    accuracy.fitting = c(); accuracy.validation = c(); varimp.R2 = c(); residual = c()
  }
  Y = as.numeric(as.character(dat[,y.names]))
  oob.error = ifelse(length(selection) != 0, mean((Y - as.numeric(as.character(predict(forest, OOB = T))))^2), mean((Y - ifelse(all(Y %in% 0:1), round(mean(Y)), mean(Y)))^2))
  oob.error.bonf = ifelse(length(selection.bonf) != 0, mean((Y - as.numeric(as.character(predict(forest.bonf, OOB = T))))^2), mean((Y - ifelse(all(Y %in% 0:1), round(mean(Y)), mean(Y)))^2))
  cat("\n", "\r"); flush.console()
  return(list("selection" = selection, "forest" = forest, "oob.error" = oob.error, "p.values" = p.vals,
              "selection.bonf" = selection.bonf, "forest.bonf" = forest.bonf, "oob.error.bonf" = oob.error.bonf, "p.values.bonf" = p.vals.bonf,
              "varimp" =obs.varimp, "varimp.selection"=varimp.bonf, "fitting" = accuracy.fitting, "validation" = accuracy.validation, "varimp.R2"=varimp.R2, "residual" =residual))
}

library(mlr)
library(party)
####
variable.importance_all = list()
model.Hapfelmeier.RF = list()
pd = list()
j = 0

for (i_response in responses) {
  df.rf.tmp = cbind(df.rf, Pred_respon_for_each_all[[i_response]])
  seed = j
  #formula = as.formula(paste(i_response, paste(c("R1","R2","R3","Lv","dissim"), collapse =" + "), sep=" ~ "))
  formula = as.formula(paste(i_response, paste(c("R1","R2","R3","Lv","dissim"), collapse =" + "), sep=" ~ "))
  
  set.seed(seed)
  model.Hapfelmeier.RF[[i_response]] = RF_permutation(formula, data=df.rf.tmp, nperm=1000, ntree=100, ncore=4)
  variable.importance = cbind(model.Hapfelmeier.RF[[i_response]]$varimp, model.Hapfelmeier.RF[[i_response]]$p.values.bonf)
  variable.importance_all[[i_response]] = variable.importance
  
  j=j+1
  
}

#add adjusted p value based on "BH" method.
p_importance = list()

for (i_response in responses) {
  variable.importance = cbind(c("R1","R2","R3","Lv","dissim"),model.Hapfelmeier.RF[[i_response]]$varimp, model.Hapfelmeier.RF[[i_response]]$p.values, p.adjust(model.Hapfelmeier.RF[[i_response]]$p.values,method = "BH"))
  variable.importance =data.frame(variable.importance)
  variable.importance$X2 = as.numeric(variable.importance$X2)
  variable.importance$X3 = as.numeric(variable.importance$X3)
  variable.importance$X4 = as.numeric(variable.importance$X4)
  variable.importance_all[[i_response]] = variable.importance
  
  p_importance[[i_response]] = ggplot(data = variable.importance, aes(x=X1, y=X2))+
    geom_segment( aes(xend=X1, yend=0)) +
    geom_point( size=2, color="orange") +
    theme_bw()
}



g_r2_all[[1]]+p_importance[["actv_ace"]]+g_r2_all[[5]]+p_importance[["Decom"]]+
  g_r2_all[[2]]+p_importance[["actv_cello"]]+g_r2_all[[6]]+p_importance[["PH"]]+
  g_r2_all[[3]]+p_importance[["actv_gluco"]]+g_r2_all[[7]]+p_importance[["WSA"]]+
  g_r2_all[[4]]+p_importance[["actv_phos"]]+
  plot_layout(ncol=4, widths=c(1,1,1,1))
  



###########################
###    Combine plots    ###
###########################

g_rawdata_all[[1]]+g_2fac_cor[[1]]+g_5fac_cor[[1]]+g_8fac_cor[[1]]+g_r2_four[[1]]+
  g_rawdata_all[[2]]+g_2fac_cor[[2]]+g_5fac_cor[[2]]+g_8fac_cor[[2]]+g_r2_four[[2]]+
  g_rawdata_all[[3]]+g_2fac_cor[[3]]+g_5fac_cor[[3]]+g_8fac_cor[[3]]+g_r2_four[[3]]+
  g_rawdata_all[[4]]+g_2fac_cor[[4]]+g_5fac_cor[[4]]+g_8fac_cor[[4]]+g_r2_four[[4]]+
  plot_layout(ncol=5, widths=c(2,1,1,1,1))

g_rawdata_all[[5]]+g_2fac_cor[[5]]+g_5fac_cor[[5]]+g_8fac_cor[[5]]+g_r2_four[[5]]+
  g_rawdata_all[[6]]+g_2fac_cor[[6]]+g_5fac_cor[[6]]+g_8fac_cor[[6]]+g_r2_four[[6]]+
  g_rawdata_all[[7]]+g_2fac_cor[[7]]+g_5fac_cor[[7]]+g_8fac_cor[[7]]+g_r2_four[[7]]+
  plot_spacer()+
  plot_layout(ncol=5, widths=c(2,1,1,1,1))

g_r2_all[[1]]+g_correl_all[[1]]+g_r2_all[[5]]+g_correl_all[[5]]+
  g_r2_all[[2]]+g_correl_all[[2]]+g_r2_all[[6]]+g_correl_all[[6]]+
  g_r2_all[[3]]+g_correl_all[[3]]+g_r2_all[[7]]+g_correl_all[[7]]+
  g_r2_all[[4]]+g_correl_all[[4]]+
  plot_layout(ncol=4, widths=c(1,1,1,1))

g_r2_all[[1]]+g_r2_all[[5]]+
  g_r2_all[[2]]+g_r2_all[[6]]+
  g_r2_all[[3]]+g_r2_all[[7]]+
  g_r2_all[[4]]+
  plot_layout(ncol=2, widths=c(1,1))

############################################################################################
### Hierarchical linear model framework (alternative for random forest model analysis)   ###
############################################################################################

### 1.Decomposition

h_levels = levels[-1]
df.multiple_decom = df.rf.decom[df.rf.decom[,"Lv"]%in%h_levels,]


m_R_decom <- lm(Decom ~ R1+R2+R3,
                 data = df.multiple_decom)
m_lv_decom <- lm(Decom ~ R1+R2+R3+Lv,
                  data = df.multiple_decom)
m_L_decom <- lm(Decom ~ Lv,
                 data = df.multiple_decom)
m_D_decom <- lm(Decom ~ dissim,
                 data = df.multiple_decom)
m_dissim_decom <- lm(Decom ~ R1+R2+R3+dissim,
                  data = df.multiple_decom)

m_all_decom <- lm(Decom ~ R1+R2+R3+Lv+dissim,
                 data = df.multiple_decom)
m_all_sig_decom<- lm(Decom ~ R1+R2+R3+Lv+dissim+P+C+L+N+A+I+SU+FU+H+M+S+D,
                     data = df.multiple_decom)

AIC(m_R_decom,m_L_decom,m_D_decom, m_lv_decom, m_dissim_decom, m_all_decom, m_all_sig_decom)
logLik(m_R_decom)
logLik(m_L_decom)
logLik(m_D_decom)
logLik(m_lv_decom)
logLik(m_dissim_decom)
logLik(m_all_decom)
logLik(m_all_sig_decom)

anova(m_R_decom,m_lv_decom)
anova(m_R_decom,m_dissim_decom)
anova(m_lv_decom,m_all_decom)
anova(m_dissim_decom,m_all_decom)
anova(m_all_decom,m_all_sig_decom)

### 2.Water stable aggregate

df.rf.wsa = cbind(df.rf, Pred_respon_for_each_all[["WSA"]])
df.multiple_WSA = df.rf.wsa[df.rf.wsa[,"Lv"]%in%h_levels,]

m_R_WSA <- lm(WSA ~ R1+R2+R3,
             data = df.multiple_WSA) #M1
m_L_WSA <- lm(WSA ~ Lv,
             data = df.multiple_WSA) #M2
m_D_WSA <- lm(WSA ~ dissim,
             data = df.multiple_WSA) #M3
m_lv_WSA <- lm(WSA ~ R1+R2+R3+Lv,
              data = df.multiple_WSA) #M4
m_dissim_WSA <- lm(WSA ~ R1+R2+R3+dissim,
                  data = df.multiple_WSA) #M5
m_all_WSA <- lm(WSA ~ R1+R2+R3+Lv+dissim,
               data = df.multiple_WSA) #M6
m_all_sig_WSA<- lm(WSA ~ R1+R2+R3+Lv+dissim+P+C+L+N+A+I+SU+FU+H+M+S+D,
                  data = df.multiple_WSA) #M7

AIC(m_R_WSA,m_L_WSA,m_D_WSA, m_lv_WSA, m_dissim_WSA, m_all_WSA, m_all_sig_WSA)
logLik(m_R_WSA)
logLik(m_L_WSA)
logLik(m_D_WSA)
logLik(m_lv_WSA)
logLik(m_dissim_WSA)
logLik(m_all_WSA)
logLik(m_all_sig_WSA)

anova(m_R_WSA,m_lv_WSA)
anova(m_R_WSA,m_dissim_WSA)
anova(m_lv_WSA,m_all_WSA)
anova(m_dissim_WSA,m_all_WSA)
anova(m_all_WSA,m_all_sig_WSA)


### 3.Soil pH
df.rf.ph = cbind(df.rf, Pred_respon_for_each_all[["PH"]])
df.multiple_PH = df.rf.ph[df.rf.ph[,"Lv"]%in%h_levels,]

m_R_PH <- lm(PH ~ R1+R2+R3,
                data = df.multiple_PH) #M1
m_L_PH <- lm(PH ~ Lv,
                data = df.multiple_PH) #M2
m_D_PH <- lm(PH ~ dissim,
                data = df.multiple_PH) #M3
m_lv_PH <- lm(PH ~ R1+R2+R3+Lv,
              data = df.multiple_PH) #M4
m_dissim_PH <- lm(PH ~ R1+R2+R3+dissim,
                     data = df.multiple_PH) #M5
m_all_PH <- lm(PH ~ R1+R2+R3+Lv+dissim,
                  data = df.multiple_PH) #M6
m_all_sig_PH<- lm(PH ~ R1+R2+R3+Lv+dissim+P+C+L+N+A+I+SU+FU+H+M+S+D,
                     data = df.multiple_PH) #M7

AIC(m_R_PH,m_L_PH,m_D_PH, m_lv_PH, m_dissim_PH, m_all_PH, m_all_sig_PH)
logLik(m_R_PH)
logLik(m_L_PH)
logLik(m_D_PH)
logLik(m_lv_PH)
logLik(m_dissim_PH)
logLik(m_all_PH)
logLik(m_all_sig_PH)

anova(m_R_PH,m_lv_PH)
anova(m_R_PH,m_dissim_PH)
anova(m_lv_PH,m_all_PH)
anova(m_dissim_PH,m_all_PH)
anova(m_all_PH,m_all_sig_PH)

### 3.Ace activity
df.rf.ace = cbind(df.rf, Pred_respon_for_each_all[["actv_ace"]])
df.multiple_actv_ace = df.rf.ace[df.rf.ace[,"Lv"]%in%h_levels,]

m_R_actv_ace <- lm(actv_ace ~ R1+R2+R3,
             data = df.multiple_actv_ace) #M1
m_L_actv_ace <- lm(actv_ace ~ Lv,
             data = df.multiple_actv_ace) #M2
m_D_actv_ace <- lm(actv_ace ~ dissim,
             data = df.multiple_actv_ace) #M3
m_lv_actv_ace <- lm(actv_ace ~ R1+R2+R3+Lv,
              data = df.multiple_actv_ace) #M4
m_dissim_actv_ace <- lm(actv_ace ~ R1+R2+R3+dissim,
                  data = df.multiple_actv_ace) #M5
m_all_actv_ace <- lm(actv_ace ~ R1+R2+R3+Lv+dissim,
               data = df.multiple_actv_ace) #M6
m_all_sig_actv_ace<- lm(actv_ace ~ R1+R2+R3+Lv+dissim+P+C+L+N+A+I+SU+FU+H+M+S+D,
                  data = df.multiple_actv_ace) #M7

AIC(m_R_actv_ace,m_L_actv_ace,m_D_actv_ace, m_lv_actv_ace, m_dissim_actv_ace, m_all_actv_ace, m_all_sig_actv_ace)
logLik(m_R_actv_ace)
logLik(m_L_actv_ace)
logLik(m_D_actv_ace)
logLik(m_lv_actv_ace)
logLik(m_dissim_actv_ace)
logLik(m_all_actv_ace)
logLik(m_all_sig_actv_ace)

anova(m_R_actv_ace,m_lv_actv_ace)
anova(m_R_actv_ace,m_dissim_actv_ace)
anova(m_lv_actv_ace,m_all_actv_ace)
anova(m_dissim_actv_ace,m_all_actv_ace)
anova(m_all_actv_ace,m_all_sig_actv_ace)

### 3.Cellu activity
df.rf.cello = cbind(df.rf, Pred_respon_for_each_all[["actv_cello"]])
df.multiple_actv_cello = df.rf.cello[df.rf.cello[,"Lv"]%in%h_levels,]

m_R_actv_cello <- lm(actv_cello ~ R1+R2+R3,
                   data = df.multiple_actv_cello) #M1
m_L_actv_cello <- lm(actv_cello ~ Lv,
                   data = df.multiple_actv_cello) #M2
m_D_actv_cello <- lm(actv_cello ~ dissim,
                   data = df.multiple_actv_cello) #M3
m_lv_actv_cello <- lm(actv_cello ~ R1+R2+R3+Lv,
                    data = df.multiple_actv_cello) #M4
m_dissim_actv_cello <- lm(actv_cello ~ R1+R2+R3+dissim,
                        data = df.multiple_actv_cello) #M5
m_all_actv_cello <- lm(actv_cello ~ R1+R2+R3+Lv+dissim,
                     data = df.multiple_actv_cello) #M6
m_all_sig_actv_cello<- lm(actv_cello ~ R1+R2+R3+Lv+dissim+P+C+L+N+A+I+SU+FU+H+M+S+D,
                        data = df.multiple_actv_cello) #M7

AIC(m_R_actv_cello,m_L_actv_cello,m_D_actv_cello, m_lv_actv_cello, m_dissim_actv_cello, m_all_actv_cello, m_all_sig_actv_cello)
logLik(m_R_actv_cello)
logLik(m_L_actv_cello)
logLik(m_D_actv_cello)
logLik(m_lv_actv_cello)
logLik(m_dissim_actv_cello)
logLik(m_all_actv_cello)
logLik(m_all_sig_actv_cello)

anova(m_R_actv_cello,m_lv_actv_cello)
anova(m_R_actv_cello,m_dissim_actv_cello)
anova(m_lv_actv_cello,m_all_actv_cello)
anova(m_dissim_actv_cello,m_all_actv_cello)
anova(m_all_actv_cello,m_all_sig_actv_cello)

### 3.Gluco activity
df.rf.gluco = cbind(df.rf, Pred_respon_for_each_all[["actv_gluco"]])
df.multiple_actv_gluco = df.rf.gluco[df.rf.gluco[,"Lv"]%in%h_levels,]

m_R_actv_gluco <- lm(actv_gluco ~ R1+R2+R3,
                     data = df.multiple_actv_gluco) #M1
m_L_actv_gluco <- lm(actv_gluco ~ Lv,
                     data = df.multiple_actv_gluco) #M2
m_D_actv_gluco <- lm(actv_gluco ~ dissim,
                     data = df.multiple_actv_gluco) #M3
m_lv_actv_gluco <- lm(actv_gluco ~ R1+R2+R3+Lv,
                      data = df.multiple_actv_gluco) #M4
m_dissim_actv_gluco <- lm(actv_gluco ~ R1+R2+R3+dissim,
                          data = df.multiple_actv_gluco) #M5
m_all_actv_gluco <- lm(actv_gluco ~ R1+R2+R3+Lv+dissim,
                       data = df.multiple_actv_gluco) #M6
m_all_sig_actv_gluco<- lm(actv_gluco ~ R1+R2+R3+Lv+dissim+P+C+L+N+A+I+SU+FU+H+M+S+D,
                          data = df.multiple_actv_gluco) #M7

AIC(m_R_actv_gluco,m_L_actv_gluco,m_D_actv_gluco, m_lv_actv_gluco, m_dissim_actv_gluco, m_all_actv_gluco, m_all_sig_actv_gluco)
logLik(m_R_actv_gluco)
logLik(m_L_actv_gluco)
logLik(m_D_actv_gluco)
logLik(m_lv_actv_gluco)
logLik(m_dissim_actv_gluco)
logLik(m_all_actv_gluco)
logLik(m_all_sig_actv_gluco)

anova(m_R_actv_gluco,m_lv_actv_gluco)
anova(m_R_actv_gluco,m_dissim_actv_gluco)
anova(m_lv_actv_gluco,m_all_actv_gluco)
anova(m_dissim_actv_gluco,m_all_actv_gluco)
anova(m_all_actv_gluco,m_all_sig_actv_gluco)

### 3.Phos activity
df.rf.phos = cbind(df.rf, Pred_respon_for_each_all[["actv_phos"]])
df.multiple_actv_phos = df.rf.phos[df.rf.phos[,"Lv"]%in%h_levels,]

m_R_actv_phos <- lm(actv_phos ~ R1+R2+R3,
                     data = df.multiple_actv_phos) #M1
m_L_actv_phos <- lm(actv_phos ~ Lv,
                     data = df.multiple_actv_phos) #M2
m_D_actv_phos <- lm(actv_phos ~ dissim,
                     data = df.multiple_actv_phos) #M3
m_lv_actv_phos <- lm(actv_phos ~ R1+R2+R3+Lv,
                      data = df.multiple_actv_phos) #M4
m_dissim_actv_phos <- lm(actv_phos ~ R1+R2+R3+dissim,
                          data = df.multiple_actv_phos) #M5
m_all_actv_phos <- lm(actv_phos ~ R1+R2+R3+Lv+dissim,
                       data = df.multiple_actv_phos) #M6
m_all_sig_actv_phos<- lm(actv_phos ~ R1+R2+R3+Lv+dissim+P+C+L+N+A+I+SU+FU+H+M+S+D,
                          data = df.multiple_actv_phos) #M7

AIC(m_R_actv_phos,m_L_actv_phos,m_D_actv_phos, m_lv_actv_phos, m_dissim_actv_phos, m_all_actv_phos, m_all_sig_actv_phos)
logLik(m_R_actv_phos)
logLik(m_L_actv_phos)
logLik(m_D_actv_phos)
logLik(m_lv_actv_phos)
logLik(m_dissim_actv_phos)
logLik(m_all_actv_phos)
logLik(m_all_sig_actv_phos)

anova(m_R_actv_phos,m_lv_actv_phos)
anova(m_R_actv_phos,m_dissim_actv_phos)
anova(m_lv_actv_phos,m_all_actv_phos)
anova(m_dissim_actv_phos,m_all_actv_phos)
anova(m_all_actv_phos,m_all_sig_actv_phos)

#####################################################
###.    Part 4 Factor net interaction analysis    ###
#####################################################

# 1. Import data

library(ggplot2)
library(ggridges)
library(patchwork)
library(party)
library(caret)
library(dplyr)
library(randomForest)
library(dabestr)

df=read.csv("/Effect_of_multiple_GCF_dissimilarity/Data/CombinedData.csv")
df=df[df["remark"]!="WC",]
df$remark=factor(df$remark, levels = unique(df$remark))
levels = c( "2", "5", "8")
stressors = c("P","C","L","N","A","I","SU","FU","H","M","S","D")
responses = c("actv_ace","actv_cello","actv_gluco","actv_phos","Decom","PH","WSA")
df<- na.roughfix(df) 

# 2. Loading functions

# Functions for calculating null model predictions for specific factor combinations
# Tutorial: https://mohanb96.github.io/Nullmodels.html

NullModel = function(object, selected_factors=vector(), n_perm = 500){
  
  output = list()
  CT = object[["summary"]][["mean"]][[1]]
  input=object[["data"]]
  population_CT= input[input[,as.character(object[["x"]][[2]])] == object[["result"]][["control_group"]][1], object[["result"]][["variable"]][1]]
  size_CT=length(population_CT)
  
  for (type in c("additive","multiplicative", "dominative")) {
    bs = numeric(0)
    for (id in c(1:n_perm)) {
      each_effect = numeric(0)
      k_CT = mean(sample(population_CT, size_CT, replace = T))
      for (treatment in selected_factors) {
        population_TR = input[input[,as.character(object[["x"]][[2]])] == object[["result"]][["test_group"]][which(object[["result"]][["test_group"]]==treatment)], object[["result"]][["variable"]][1]]
        size_TR = length(population_TR)
        k_TR = mean(sample(population_TR,size_TR,replace = T))
        
        # ES estimate depending on the type of null hypothesis
        if(type == "additive")        each_effect = append(each_effect,(k_TR - k_CT))
        if(type == "multiplicative")  each_effect = append(each_effect, (k_TR-k_CT)/k_CT)
        if(type == "dominative")      each_effect = append(each_effect,(k_TR - k_CT))
      }
      
      if(type == "additive") {
        joint_effect = sum(each_effect)
        pre_response = joint_effect+CT
      }
      
      if(type=="multiplicative"){
        z = 1
        for(m in c(1:length(selected_factors))) {
          z = z * (1 + each_effect[m])
          joint_effect = (z - 1)*k_CT
        }
        pre_response =joint_effect+CT
      }
      
      if(type=="dominative")  {
        joint_effect = each_effect[which(max(abs(each_effect))==abs(each_effect))]
        pre_response =joint_effect+CT
      }
      
      bs = append(bs, pre_response)
    }
    output[[type]] = bs
  }
  return(output)
}

NullModel_summary = function(null_data, actual_data = vector()){
  output = list()
  
  if(length(actual_data)<=2){
    output[["actual"]] = c(mean(actual_data),mean(actual_data),mean(actual_data),1)
    for (type in c("additive","multiplicative", "dominative")) {  
      bs = rep(mean(actual_data),length(null_data[[type]])) - null_data[[type]]
      p = length(which(bs>0))/length(bs)
      p = min(p, 1-p)
      output[[type]] = c(quantile(null_data[[type]], .025), mean(null_data[[type]]), quantile(null_data[[type]], .975), p)
    }
  }
  
  if(length(actual_data)>=3){
    #re-sampling of actual data
    size_actul=length(actual_data)
    bs_actual = numeric(0)
    for (i in c(1:length(null_data[[1]]))) {
      k_actual = mean(sample(actual_data, size_actul, replace = T))
      bs_actual = append(bs_actual, k_actual)
    }
    
    output[["actual"]] = c(quantile(bs_actual, .025), mean(bs_actual), quantile(bs_actual, .975), 1)
    
    for (type in c("additive","multiplicative", "dominative")) { 
      bs = bs_actual - null_data[[type]]
      p = length(which(bs>0))/length(bs)
      p = min(p, 1-p)
      output[[type]] = c(quantile(null_data[[type]], .025), mean(null_data[[type]]), quantile(null_data[[type]], .975), p)
    }
  }
  
  output = t(data.frame(output))
  colnames(output) = c("X2.5%", "mean", "X97.5%", "p_value")
  output = data.frame(ES = c("actual", "additive","multiplicative","dominative"), output)
  row.names(output) = c()
  
  return(output)
}

# 3. Effect size calculation

ES.plot_ace<-df%>%
  dabest(remark, actv_ace,
         idx = c("CT","I","A","P","SU","H","S","M","D", "C", "FU", "L", "N","1","2","5","8"),
         paired = FALSE)
ES.plot.meandiff_ace<-mean_diff(ES.plot_ace,reps = 500)

ES.plot_cel<-df%>%
  dabest(remark, actv_cello,
         idx = c("CT","I","A","P","SU","H","S","M","D", "C", "FU", "L", "N","1","2","5","8"),
         paired = FALSE)
ES.plot.meandiff_cel<-mean_diff(ES.plot_cel,reps = 500)

ES.plot_glu<-df%>%
  dabest(remark, actv_gluco,
         idx = c("CT","I","A","P","SU","H","S","M","D", "C", "FU", "L", "N","1","2","5","8"),
         paired = FALSE)
ES.plot.meandiff_glu<-mean_diff(ES.plot_glu,reps = 500)

ES.plot_pho<-df%>%
  dabest(remark, actv_phos,
         idx = c("CT","I","A","P","SU","H","S","M","D", "C", "FU", "L", "N","1","2","5","8"),
         paired = FALSE)
ES.plot.meandiff_pho<-mean_diff(ES.plot_pho,reps = 500)

ES.plot_dec<-df%>%
  dabest(remark, Decom,
         idx = c("CT","I","A","P","SU","H","S","M","D", "C", "FU", "L", "N","1","2","5","8"),
         paired = FALSE)
ES.plot.meandiff_dec<-mean_diff(ES.plot_dec,reps = 500)

ES.plot_ph<-df%>%
  dabest(remark, PH,
         idx = c("CT","I","A","P","SU","H","S","M","D", "C", "FU", "L", "N","1","2","5","8"),
         paired = FALSE)
ES.plot.meandiff_ph<-mean_diff(ES.plot_ph,reps = 500)

ES.plot_wsa<-df%>%
  dabest(remark, WSA,
         idx = c("CT","I","A","P","SU","H","S","M","D", "C", "FU", "L", "N","1","2","5","8"),
         paired = FALSE)
ES.plot.meandiff_wsa<-mean_diff(ES.plot_wsa,reps = 500)

plot(ES.plot.meandiff_wsa)


# 4. Main analysis (Null model)

ES.plot.meandiff = list()

ES.plot.meandiff[["actv_ace"]]=ES.plot.meandiff_ace
ES.plot.meandiff[["actv_cello"]]=ES.plot.meandiff_cel
ES.plot.meandiff[["actv_gluco"]]=ES.plot.meandiff_glu
ES.plot.meandiff[["actv_phos"]]=ES.plot.meandiff_pho
ES.plot.meandiff[["Decom"]]=ES.plot.meandiff_dec
ES.plot.meandiff[["PH"]]=ES.plot.meandiff_ph
ES.plot.meandiff[["WSA"]]=ES.plot.meandiff_wsa

df_H=df[df[,"remark"]%in%levels,]
df_CK=df[df[,"remark"]=="CT",]
n_perm = 1000

Diviation_additive = data.frame(matrix(data = NA, nrow = 147, ncol = 6))
colnames(Diviation_additive)=c("Lv","diviation","P","dissim","actual_ES","null_ES")

Diviation_multiplicative = data.frame(matrix(data = NA, nrow = 147, ncol = 6))
colnames(Diviation_multiplicative)=c("Lv","diviation","P","dissim","actual_ES","null_ES")

Diviation_dominative = data.frame(matrix(data = NA, nrow = 147, ncol = 6))
colnames(Diviation_dominative)=c("Lv","diviation","P","dissim","actual_ES","null_ES")

Diviation_3_models = list()


for (i_response in responses) {
  
  for (treatment_i in 1:nrow(df_H)) {
    combination_i=c()
    combination_i=stressors[which(df_H[treatment_i,10:21]==1)]
    Null_modle_i_Treatment = NullModel(ES.plot.meandiff[[i_response]], selected_factors = combination_i,n_perm = n_perm)
    Null_modle_summary_i_Treatment = NullModel_summary(Null_modle_i_Treatment, actual_data = df_H[treatment_i,i_response]) #df_H[treatment_i,6] refere to decomposition rate
    
    Diviation_additive[treatment_i,2]=(Null_modle_summary_i_Treatment[1,3]-Null_modle_summary_i_Treatment[2,3])/Null_modle_summary_i_Treatment[2,3]  #changed
    Diviation_additive[treatment_i,3]=Null_modle_summary_i_Treatment[2,5]
    Diviation_additive[treatment_i,4]=df_H[treatment_i,9]
    Diviation_additive[treatment_i,1]=df_H[treatment_i,23]
    Diviation_additive[treatment_i,5]=Null_modle_summary_i_Treatment[1,3]-mean(df_CK[,i_response])
    Diviation_additive[treatment_i,6]=Null_modle_summary_i_Treatment[2,3]-mean(df_CK[,i_response])
    
    Diviation_multiplicative[treatment_i,2]=(Null_modle_summary_i_Treatment[1,3]-Null_modle_summary_i_Treatment[3,3])/Null_modle_summary_i_Treatment[3,3]
    Diviation_multiplicative[treatment_i,3]=Null_modle_summary_i_Treatment[3,5]
    Diviation_multiplicative[treatment_i,4]=df_H[treatment_i,9]
    Diviation_multiplicative[treatment_i,1]=df_H[treatment_i,23]
    Diviation_multiplicative[treatment_i,5]=Null_modle_summary_i_Treatment[1,3]-mean(df_CK[,i_response])
    Diviation_multiplicative[treatment_i,6]=Null_modle_summary_i_Treatment[3,3]-mean(df_CK[,i_response])
    
    Diviation_dominative[treatment_i,2]=(Null_modle_summary_i_Treatment[1,3]-Null_modle_summary_i_Treatment[4,3])/Null_modle_summary_i_Treatment[4,3]
    Diviation_dominative[treatment_i,3]=Null_modle_summary_i_Treatment[4,5]
    Diviation_dominative[treatment_i,4]=df_H[treatment_i,9]
    Diviation_dominative[treatment_i,1]=df_H[treatment_i,23]
    Diviation_dominative[treatment_i,5]=Null_modle_summary_i_Treatment[1,3]-mean(df_CK[,i_response])
    Diviation_dominative[treatment_i,6]=Null_modle_summary_i_Treatment[4,3]-mean(df_CK[,i_response])
  }
  
  Diviation_3_models[[i_response]][["additive"]]=Diviation_additive
  Diviation_3_models[[i_response]][["multiplicative"]]=Diviation_multiplicative
  Diviation_3_models[[i_response]][["dominative"]]=Diviation_dominative
  
}

# 5.1 Identifying interaction type


positive_responses = c("actv_ace","actv_cello","actv_gluco","actv_phos","PH")
negativ_responses = c("Decom", "WSA")
null_models = c("additive", "multiplicative", "dominative")

for (i_response in negativ_responses) {
  for (n_model in null_models) {
    for (i in 1:147) {
      if(Diviation_3_models[[i_response]][[n_model]][["diviation"]][i]>0 ) Diviation_3_models[[i_response]][[n_model]][["interaction_type"]][i] = "Antagonistic"
      if(Diviation_3_models[[i_response]][[n_model]][["diviation"]][i]<0 ) Diviation_3_models[[i_response]][[n_model]][["interaction_type"]][i] = "Synergistic"
      if(abs(Diviation_3_models[[i_response]][[n_model]][["diviation"]][i])<0.0 ) Diviation_3_models[[i_response]][[n_model]][i,"interaction_type"] = n_model ### threshold adjustment
      if(Diviation_3_models[[i_response]][[n_model]][i,"P"]>=0.05) Diviation_3_models[[i_response]][[n_model]][i,"interaction_type"] = n_model
    }
  }
}

for (i_response in positive_responses) {
  for (n_model in null_models) {
    for (i in 1:147) {
      if(Diviation_3_models[[i_response]][[n_model]][["diviation"]][i]>0 ) Diviation_3_models[[i_response]][[n_model]][["interaction_type"]][i] = "Synergistic"
      if(Diviation_3_models[[i_response]][[n_model]][["diviation"]][i]<0 ) Diviation_3_models[[i_response]][[n_model]][["interaction_type"]][i] = "Antagonistic"
      if(abs(Diviation_3_models[[i_response]][[n_model]][["diviation"]][i])<0.0 ) Diviation_3_models[[i_response]][[n_model]][i,"interaction_type"] = n_model ### threshold adjustment
      if(Diviation_3_models[[i_response]][[n_model]][i,"P"]>=0.05) Diviation_3_models[[i_response]][[n_model]][i,"interaction_type"] = n_model
    }
  }
}

# Ploting (point)

p_diviation = list()
p_diviation_lv = list()
p_diviation_all = list()

for (i_response in responses) {
  
  p_diviation[[i_response]][["additive"]] = ggplot(Diviation_3_models[[i_response]][["additive"]], aes(x=dissim, y=diviation))+
    geom_point(size =2, aes(color = interaction_type))+
    scale_color_manual(values = c("Antagonistic"="#5499C780","additive"="#90949780","Synergistic"="#F1666780"))+
    geom_smooth(method = lm, aes(color = "#999999"))+
    stat_cor(method = "spearman",label.x = 0)+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
    ylim(-1,20)
  
  p_diviation[[i_response]][["multiplicative"]] = ggplot(Diviation_3_models[[i_response]][["multiplicative"]], aes(x=dissim, y=diviation))+
    geom_point(size =2, aes(color = interaction_type))+
    scale_color_manual(values = c("Antagonistic"="#5499C780","multiplicative"="#90949780","Synergistic"="#F1666780"))+
    geom_smooth(method = lm, aes(color = "#999999"))+
    stat_cor(method = "spearman",label.x = 0)+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  ylim(-5,5)
  
  p_diviation[[i_response]][["dominative"]] = ggplot(Diviation_3_models[[i_response]][["dominative"]], aes(x=dissim, y=diviation))+
    geom_point(size =2, aes(color = interaction_type))+
    scale_color_manual(values = c("Antagonistic"="#5499C780","dominative"="#90949780","Synergistic"="#F1666780"))+
    geom_smooth(method = lm, aes(color = "#999999"))+
    stat_cor(method = "spearman",label.x = 0)+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  ylim(-10,10)
}

for (i_response in responses) {
  p_diviation_lv[[i_response]][["additive"]] = 
    ggplot()+
    geom_jitter(data = Diviation_3_models[[i_response]][["additive"]], aes(x=Lv, y=diviation, color = interaction_type),width = 0.5, size = 1)+
    geom_boxplot(data = Diviation_3_models[[i_response]][["additive"]], aes( x=Lv, y=diviation, group = Lv), width = 1, fill = "#00000000", color = "#72665880")+
    scale_x_continuous(breaks = c(2,5,8))+
    scale_color_manual(values = c("Antagonistic"="#5499C780","additive"="#90949780","Synergistic"="#F1666780"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    stat_compare_means(method = "t.test")+
    xlab(i_response)
  
  p_diviation_lv[[i_response]][["multiplicative"]] = 
    ggplot()+
    geom_jitter(data = Diviation_3_models[[i_response]][["multiplicative"]], aes(x=Lv, y=diviation, color = interaction_type),width = 0.5, size = 1)+
    geom_boxplot(data = Diviation_3_models[[i_response]][["multiplicative"]], aes( x=Lv, y=diviation, group = Lv), width = 1, fill = "#00000000", color = "#72665880")+
    scale_x_continuous(breaks = c(2,5,8))+
    scale_color_manual(values = c("Antagonistic"="#5499C780","multiplicative"="#90949780","Synergistic"="#F1666780"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  
  p_diviation_lv[[i_response]][["dominative"]] = 
    ggplot()+
    geom_jitter(data = Diviation_3_models[[i_response]][["dominative"]], aes(x=Lv, y=diviation, color = interaction_type),width = 0.5, size = 1)+
    geom_boxplot(data = Diviation_3_models[[i_response]][["dominative"]], aes( x=Lv, y=diviation, group = Lv), width = 1, fill = "#00000000", color = "#72665880")+
    scale_x_continuous(breaks = c(2,5,8))+
    scale_color_manual(values = c("Antagonistic"="#5499C780","dominative"="#90949780","Synergistic"="#F1666780"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
}

Diviation_models_all = list()

for (i_response in responses) {
  for (n_model in null_models) {
    Diviation_3_models[[i_response]][[n_model]]["model_type"] = n_model
  }
  Diviation_models_all[[i_response]] = rbind(Diviation_3_models[[i_response]][["additive"]],Diviation_3_models[[i_response]][["multiplicative"]],Diviation_3_models[[i_response]][["dominative"]])
}

for (i_response in responses) {
  p_diviation_all[[i_response]] = local({
    i_response = i_response
    Diviation_models_all = Diviation_models_all
    Diviation_models_all[[i_response]] = data.frame(Diviation_models_all[[i_response]])%>%
      arrange(model_type)%>%
      mutate(model_type = factor(model_type, levels = c("additive","multiplicative","dominative")))
    
    ggplot(data = Diviation_models_all[[i_response]], aes(x = model_type, y=diviation))+
    geom_jitter(aes(color = interaction_type), width = 0.2, size = 1)+
    scale_color_manual(values = c("Antagonistic"="#5499C780","additive"="#90949780","multiplicative"="#90949780","dominative"="#90949780","Synergistic"="#F1666780"))+
    geom_boxplot(aes(x = model_type, y=diviation, group = model_type), width = 0.5, fill = "#00000000", color = "#72665880")+
    scale_x_discrete("model_type", labels = c("additive","multiplicative","dominative"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  })
}

p_diviation[["Decom"]][["additive"]]+p_diviation[["Decom"]][["multiplicative"]]+p_diviation[["Decom"]][["dominative"]]+
  p_diviation[["PH"]][["additive"]]+p_diviation[["PH"]][["multiplicative"]]+p_diviation[["PH"]][["dominative"]]+
  p_diviation[["WSA"]][["additive"]]+p_diviation[["WSA"]][["multiplicative"]]+p_diviation[["WSA"]][["dominative"]]+
  plot_spacer()+
  plot_layout(ncol=3, widths=c(2,2,2))

p_diviation[["actv_ace"]][["additive"]]+p_diviation[["actv_ace"]][["multiplicative"]]+p_diviation[["actv_ace"]][["dominative"]]+
  p_diviation[["actv_cello"]][["additive"]]+p_diviation[["actv_cello"]][["multiplicative"]]+p_diviation[["actv_cello"]][["dominative"]]+
  p_diviation[["actv_gluco"]][["additive"]]+p_diviation[["actv_gluco"]][["multiplicative"]]+p_diviation[["actv_gluco"]][["dominative"]]+
  p_diviation[["actv_phos"]][["additive"]]+p_diviation[["actv_phos"]][["multiplicative"]]+p_diviation[["actv_phos"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))

p_diviation_lv[["Decom"]][["additive"]]+p_diviation_lv[["Decom"]][["multiplicative"]]+p_diviation_lv[["Decom"]][["dominative"]]+
  p_diviation_lv[["PH"]][["additive"]]+p_diviation_lv[["PH"]][["multiplicative"]]+p_diviation_lv[["PH"]][["dominative"]]+
  p_diviation_lv[["WSA"]][["additive"]]+p_diviation_lv[["WSA"]][["multiplicative"]]+p_diviation_lv[["WSA"]][["dominative"]]+
  plot_spacer()+
  plot_layout(ncol=3, widths=c(2,2,2))

p_diviation_lv[["actv_ace"]][["additive"]]+p_diviation_lv[["actv_ace"]][["multiplicative"]]+p_diviation_lv[["actv_ace"]][["dominative"]]+
  p_diviation_lv[["actv_cello"]][["additive"]]+p_diviation_lv[["actv_cello"]][["multiplicative"]]+p_diviation_lv[["actv_cello"]][["dominative"]]+
  p_diviation_lv[["actv_gluco"]][["additive"]]+p_diviation_lv[["actv_gluco"]][["multiplicative"]]+p_diviation_lv[["actv_gluco"]][["dominative"]]+
  p_diviation_lv[["actv_phos"]][["additive"]]+p_diviation_lv[["actv_phos"]][["multiplicative"]]+p_diviation_lv[["actv_phos"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))

p_lv_lm[["Decom"]][["additive"]]+p_lv_lm[["Decom"]][["multiplicative"]]+p_lv_lm[["Decom"]][["dominative"]]+
  p_lv_lm[["PH"]][["additive"]]+p_lv_lm[["PH"]][["multiplicative"]]+p_lv_lm[["PH"]][["dominative"]]+
  p_lv_lm[["WSA"]][["additive"]]+p_lv_lm[["WSA"]][["multiplicative"]]+p_lv_lm[["WSA"]][["dominative"]]+
  plot_spacer()+
  plot_layout(ncol=3, widths=c(2,2,2))

p_lv_lm[["actv_ace"]][["additive"]]+p_lv_lm[["actv_ace"]][["multiplicative"]]+p_lv_lm[["actv_ace"]][["dominative"]]+
  p_lv_lm[["actv_cello"]][["additive"]]+p_lv_lm[["actv_cello"]][["multiplicative"]]+p_lv_lm[["actv_cello"]][["dominative"]]+
  p_lv_lm[["actv_gluco"]][["additive"]]+p_lv_lm[["actv_gluco"]][["multiplicative"]]+p_lv_lm[["actv_gluco"]][["dominative"]]+
  p_lv_lm[["actv_phos"]][["additive"]]+p_lv_lm[["actv_phos"]][["multiplicative"]]+p_lv_lm[["actv_phos"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))

# 5.2 Summarize proportions of occurence of factor net interaction types by three null model assumptions

# 5.2.1 Summarize proportion of occurence of factor net interaction types by factor level

InterType_percentage_lv = list()

InterType_lv_i = data.frame(matrix(data = NA, nrow = 9, ncol = 5))
colnames(InterType_lv_i)=c("Lv","InteractionType","N","n","P")
InterType_lv_i$Lv=c(rep(2,3),rep(5,3),rep(8,3))
InterType_lv_i$InteractionType=c(rep(c("Synergistic","Null","Antagonistic"),3))

for (i_response in responses) {
  for (n_model in null_models) {
    j=1
    for (Lv in levels) {
      InterType_lv_i[j,"n"] = length(which(Diviation_3_models[[i_response]][[n_model]][Diviation_3_models[[i_response]][[n_model]][["Lv"]]==Lv,"interaction_type"]=="Synergistic"))
      InterType_lv_i[j+2,"n"] = length(which(Diviation_3_models[[i_response]][[n_model]][Diviation_3_models[[i_response]][[n_model]][["Lv"]]==Lv,"interaction_type"]=="Antagonistic"))
      
      InterType_lv_i[j,"N"] =nrow(Diviation_3_models[[i_response]][[n_model]][Diviation_3_models[[i_response]][[n_model]][["Lv"]]==Lv,])
      InterType_lv_i[j+1,"N"]=nrow(Diviation_3_models[[i_response]][[n_model]][Diviation_3_models[[i_response]][[n_model]][["Lv"]]==Lv,])
      InterType_lv_i[j+2,"N"]=nrow(Diviation_3_models[[i_response]][[n_model]][Diviation_3_models[[i_response]][[n_model]][["Lv"]]==Lv,])
      
      InterType_lv_i[j+1,"n"] = InterType_lv_i[j,"N"]-InterType_lv_i[j,"n"]-InterType_lv_i[j+2,"n"]
      
      InterType_lv_i[j,"P"]=InterType_lv_i[j,"n"]/InterType_lv_i[j,"N"]
      InterType_lv_i[j+1,"P"]=InterType_lv_i[j+1,"n"]/InterType_lv_i[j+1,"N"]
      InterType_lv_i[j+2,"P"]=InterType_lv_i[j+2,"n"]/InterType_lv_i[j+2,"N"]
      
      j=j+3
    }
    InterType_percentage_lv[[i_response]][[n_model]]=InterType_lv_i 
  }
}

# Ploting percentage barplot

p_percent_level = list()

for (i_response in responses) {
  p_percent_level[[i_response]][["additive"]]=ggplot(InterType_percentage_lv[[i_response]][["additive"]], aes(fill = InteractionType, y=P, x=Lv))+
    geom_bar(position = "fill", stat = "identity")+
    scale_fill_manual(values = c("Antagonistic"="#5499C7","Null"="#909497","Synergistic"="#F16667"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  
  p_percent_level[[i_response]][["multiplicative"]]=ggplot(InterType_percentage_lv[[i_response]][["multiplicative"]], aes(fill = InteractionType, y=P, x=Lv))+
    geom_bar(position = "fill", stat = "identity")+
    scale_fill_manual(values = c("Antagonistic"="#5499C7","Null"="#909497","Synergistic"="#F16667"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  
  p_percent_level[[i_response]][["dominative"]]=ggplot(InterType_percentage_lv[[i_response]][["dominative"]], aes(fill = InteractionType, y=P, x=Lv))+
    geom_bar(position = "fill", stat = "identity")+
    scale_fill_manual(values = c("Antagonistic"="#5499C7","Null"="#909497","Synergistic"="#F16667"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
}

p_percent_level[["Decom"]][["additive"]]+p_percent_level[["Decom"]][["multiplicative"]]+p_percent_level[["Decom"]][["dominative"]]+
  p_percent_level[["PH"]][["additive"]]+p_percent_level[["PH"]][["multiplicative"]]+p_percent_level[["PH"]][["dominative"]]+
  p_percent_level[["WSA"]][["additive"]]+p_percent_level[["WSA"]][["multiplicative"]]+p_percent_level[["WSA"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))

p_percent_level[["actv_ace"]][["additive"]]+p_percent_level[["actv_ace"]][["multiplicative"]]+p_percent_level[["actv_ace"]][["dominative"]]+
  p_percent_level[["actv_cello"]][["additive"]]+p_percent_level[["actv_cello"]][["multiplicative"]]+p_percent_level[["actv_cello"]][["dominative"]]+
  p_percent_level[["actv_gluco"]][["additive"]]+p_percent_level[["actv_gluco"]][["multiplicative"]]+p_percent_level[["actv_gluco"]][["dominative"]]+
  p_percent_level[["actv_phos"]][["additive"]]+p_percent_level[["actv_phos"]][["multiplicative"]]+p_percent_level[["actv_phos"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))

# 5.2.2 Summarize proportion of occurence of factor net interaction types by factor dissimilarity

InterType_percentage_dis= list()

InterType_dis_i = data.frame(matrix(data = NA, nrow = 12, ncol = 6))
colnames(InterType_dis_i)=c("dissim_range","Lv","InteractionType","N","n","P")
InterType_dis_i$dissim_range=c(rep("[0,0.25)",3),rep("[0.25,0.5)",3),rep("[0.5,0.75)",3),rep("[0.75,1]",3))
InterType_dis_i$Lv=c(rep(1,3),rep(2,3),rep(3,3),rep(4,3))
InterType_dis_i$InteractionType=c(rep(c("Synergistic","Null","Antagonistic"),4))

id=list()
id[[1]]=which(Diviation_3_models[["Decom"]][["additive"]][["dissim"]]<0.25)
id[[2]]=which(Diviation_3_models[["Decom"]][["additive"]][["dissim"]]>=0.25 & Diviation_3_models[["Decom"]][["additive"]][["dissim"]]<0.5)
id[[3]]=which(Diviation_3_models[["Decom"]][["additive"]][["dissim"]]>=0.5 & Diviation_3_models[["Decom"]][["additive"]][["dissim"]]<0.75)
id[[4]]=which(Diviation_3_models[["Decom"]][["additive"]][["dissim"]]>=0.75)

InterType_dis_i$N=c(rep(length(id[[1]]),3), rep(length(id[[2]]),3), rep(length(id[[3]]),3), rep(length(id[[4]]),3))

for (i_response in responses) {
  for (n_model in null_models) {
    j=1
    for (Lv in 1:4) {
      InterType_dis_i[j,"n"]=length(which(Diviation_3_models[[i_response]][[n_model]][id[[Lv]],"interaction_type"]=="Synergistic"))
      InterType_dis_i[j+2,"n"]=length(which(Diviation_3_models[[i_response]][[n_model]][id[[Lv]],"interaction_type"]=="Antagonistic"))
      InterType_dis_i[j+1,"n"]=InterType_dis_i[j+1,"N"]-InterType_dis_i[j,"n"]-InterType_dis_i[j+2,"n"]
      
      InterType_dis_i[j,"P"]=InterType_dis_i[j,"n"]/InterType_dis_i[j,"N"]
      InterType_dis_i[j+1,"P"]=InterType_dis_i[j+1,"n"]/InterType_dis_i[j+1,"N"]
      InterType_dis_i[j+2,"P"]=InterType_dis_i[j+2,"n"]/InterType_dis_i[j+2,"N"]
      
      j=j+3  
    }
    
    InterType_percentage_dis[[i_response]][[n_model]]=InterType_dis_i 
  }
}

# Ploting

p_percent_dis = list()

for (i_response in responses) {
  p_percent_dis[[i_response]][["additive"]]=ggplot(InterType_percentage_dis[[i_response]][["additive"]], aes(fill = InteractionType, y=P, x=Lv))+
    geom_bar(position = "fill", stat = "identity")+
    scale_fill_manual(values = c("Antagonistic"="#5499C7","Null"="#909497","Synergistic"="#F16667"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  
  p_percent_dis[[i_response]][["multiplicative"]]=ggplot(InterType_percentage_dis[[i_response]][["multiplicative"]], aes(fill = InteractionType, y=P, x=Lv))+
    geom_bar(position = "fill", stat = "identity")+
    scale_fill_manual(values = c("Antagonistic"="#5499C7","Null"="#909497","Synergistic"="#F16667"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  
  p_percent_dis[[i_response]][["dominative"]]=ggplot(InterType_percentage_dis[[i_response]][["dominative"]], aes(fill = InteractionType, y=P, x=Lv))+
    geom_bar(position = "fill", stat = "identity")+
    scale_fill_manual(values = c("Antagonistic"="#5499C7","Null"="#909497","Synergistic"="#F16667"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
}

p_percent_dis[["Decom"]][["additive"]]+p_percent_dis[["Decom"]][["multiplicative"]]+p_percent_dis[["Decom"]][["dominative"]]+
  p_percent_dis[["PH"]][["additive"]]+p_percent_dis[["PH"]][["multiplicative"]]+p_percent_dis[["PH"]][["dominative"]]+
  p_percent_dis[["WSA"]][["additive"]]+p_percent_dis[["WSA"]][["multiplicative"]]+p_percent_dis[["WSA"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))

p_percent_dis[["actv_ace"]][["additive"]]+p_percent_dis[["actv_ace"]][["multiplicative"]]+p_percent_dis[["actv_ace"]][["dominative"]]+
  p_percent_dis[["actv_cello"]][["additive"]]+p_percent_dis[["actv_cello"]][["multiplicative"]]+p_percent_dis[["actv_cello"]][["dominative"]]+
  p_percent_dis[["actv_gluco"]][["additive"]]+p_percent_dis[["actv_gluco"]][["multiplicative"]]+p_percent_dis[["actv_gluco"]][["dominative"]]+
  p_percent_dis[["actv_phos"]][["additive"]]+p_percent_dis[["actv_phos"]][["multiplicative"]]+p_percent_dis[["actv_phos"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))


############# Deviation t_test for each factor level ##############

t.test(Diviation_3_models[["Decom"]][["additive"]][Diviation_3_models[["Decom"]][["additive"]]["Lv"]=="8","diviation"], y = NULL)
t.test(Diviation_3_models[["Decom"]][["multiplicative"]][Diviation_3_models[["Decom"]][["multiplicative"]]["Lv"]=="8","diviation"], y = NULL)
t.test(Diviation_3_models[["Decom"]][["dominative"]][Diviation_3_models[["Decom"]][["dominative"]]["Lv"]=="8","diviation"], y = NULL)



t.test(Diviation_3_models[["PH"]][["additive"]][Diviation_3_models[["PH"]][["additive"]]["Lv"]=="2","diviation"], y = NULL)
t.test(Diviation_3_models[["PH"]][["multiplicative"]][Diviation_3_models[["PH"]][["multiplicative"]]["Lv"]=="2","diviation"], y = NULL)
t.test(Diviation_3_models[["PH"]][["dominative"]][Diviation_3_models[["PH"]][["dominative"]]["Lv"]=="2","diviation"], y = NULL)



t.test(Diviation_3_models[["WSA"]][["additive"]][Diviation_3_models[["WSA"]][["additive"]]["Lv"]=="5","diviation"], y = NULL)
t.test(Diviation_3_models[["WSA"]][["multiplicative"]][Diviation_3_models[["WSA"]][["multiplicative"]]["Lv"]=="5","diviation"], y = NULL)
t.test(Diviation_3_models[["WSA"]][["dominative"]][Diviation_3_models[["WSA"]][["dominative"]]["Lv"]=="8","diviation"], y = NULL)



t.test(Diviation_3_models[["actv_ace"]][["additive"]][Diviation_3_models[["actv_ace"]][["additive"]]["Lv"]=="8","diviation"], y = NULL)
t.test(Diviation_3_models[["actv_ace"]][["multiplicative"]][Diviation_3_models[["actv_ace"]][["multiplicative"]]["Lv"]=="8","diviation"], y = NULL)
t.test(Diviation_3_models[["actv_ace"]][["dominative"]][Diviation_3_models[["actv_ace"]][["dominative"]]["Lv"]=="8","diviation"], y = NULL)



t.test(Diviation_3_models[["actv_cello"]][["additive"]][Diviation_3_models[["actv_cello"]][["additive"]]["Lv"]=="8","diviation"], y = NULL)
t.test(Diviation_3_models[["actv_cello"]][["multiplicative"]][Diviation_3_models[["actv_cello"]][["multiplicative"]]["Lv"]=="8","diviation"], y = NULL)
t.test(Diviation_3_models[["actv_cello"]][["dominative"]][Diviation_3_models[["actv_cello"]][["dominative"]]["Lv"]=="2","diviation"], y = NULL)


t.test(Diviation_3_models[["actv_gluco"]][["additive"]][Diviation_3_models[["actv_gluco"]][["additive"]]["Lv"]=="5","diviation"], y = NULL)
t.test(Diviation_3_models[["actv_gluco"]][["multiplicative"]][Diviation_3_models[["actv_gluco"]][["multiplicative"]]["Lv"]=="8","diviation"], y = NULL)
t.test(Diviation_3_models[["actv_gluco"]][["dominative"]][Diviation_3_models[["actv_gluco"]][["dominative"]]["Lv"]=="8","diviation"], y = NULL)


t.test(Diviation_3_models[["actv_phos"]][["additive"]][Diviation_3_models[["actv_phos"]][["additive"]]["Lv"]=="8","diviation"], y = NULL)
t.test(Diviation_3_models[["actv_phos"]][["multiplicative"]][Diviation_3_models[["actv_phos"]][["multiplicative"]]["Lv"]=="8","diviation"], y = NULL)
t.test(Diviation_3_models[["actv_phos"]][["dominative"]][Diviation_3_models[["actv_phos"]][["dominative"]]["Lv"]=="2","diviation"], y = NULL)

############# Deviation t_test for each null model ##############

t.test(Diviation_models_all[["Decom"]][Diviation_models_all[["Decom"]]["model_type"] =="additive","diviation"], y = NULL)
t.test(Diviation_models_all[["Decom"]][Diviation_models_all[["Decom"]]["model_type"] =="multiplicative","diviation"], y = NULL)
t.test(Diviation_models_all[["Decom"]][Diviation_models_all[["Decom"]]["model_type"] =="dominative","diviation"], y = NULL)
sd(Diviation_models_all[["Decom"]][Diviation_models_all[["Decom"]]["model_type"] =="additive","diviation"])
sd(Diviation_models_all[["Decom"]][Diviation_models_all[["Decom"]]["model_type"] =="multiplicative","diviation"])
sd(Diviation_models_all[["Decom"]][Diviation_models_all[["Decom"]]["model_type"] =="dominative","diviation"])
sum((Diviation_models_all[["Decom"]][Diviation_models_all[["Decom"]]["model_type"] =="additive","diviation"])^2)
sum((Diviation_models_all[["Decom"]][Diviation_models_all[["Decom"]]["model_type"] =="multiplicative","diviation"])^2)
sum((Diviation_models_all[["Decom"]][Diviation_models_all[["Decom"]]["model_type"] =="dominative","diviation"])^2)


t.test(Diviation_models_all[["PH"]][Diviation_models_all[["PH"]]["model_type"] =="additive","diviation"], y = NULL)
t.test(Diviation_models_all[["PH"]][Diviation_models_all[["PH"]]["model_type"] =="multiplicative","diviation"], y = NULL)
t.test(Diviation_models_all[["PH"]][Diviation_models_all[["PH"]]["model_type"] =="dominative","diviation"], y = NULL)
sd(Diviation_models_all[["PH"]][Diviation_models_all[["PH"]]["model_type"] =="additive","diviation"])
sd(Diviation_models_all[["PH"]][Diviation_models_all[["PH"]]["model_type"] =="multiplicative","diviation"])
sd(Diviation_models_all[["PH"]][Diviation_models_all[["PH"]]["model_type"] =="dominative","diviation"])
sum((Diviation_models_all[["PH"]][Diviation_models_all[["PH"]]["model_type"] =="additive","diviation"])^2)
sum((Diviation_models_all[["PH"]][Diviation_models_all[["PH"]]["model_type"] =="multiplicative","diviation"])^2)
sum((Diviation_models_all[["PH"]][Diviation_models_all[["PH"]]["model_type"] =="dominative","diviation"])^2)

t.test(Diviation_models_all[["WSA"]][Diviation_models_all[["WSA"]]["model_type"] =="additive","diviation"], y = NULL)
t.test(Diviation_models_all[["WSA"]][Diviation_models_all[["WSA"]]["model_type"] =="multiplicative","diviation"], y = NULL)
t.test(Diviation_models_all[["WSA"]][Diviation_models_all[["WSA"]]["model_type"] =="dominative","diviation"], y = NULL)
sum((Diviation_models_all[["WSA"]][Diviation_models_all[["WSA"]]["model_type"] =="additive","diviation"])^2)
sum((Diviation_models_all[["WSA"]][Diviation_models_all[["WSA"]]["model_type"] =="multiplicative","diviation"])^2)
sum((Diviation_models_all[["WSA"]][Diviation_models_all[["WSA"]]["model_type"] =="dominative","diviation"])^2)


t.test(Diviation_models_all[["actv_ace"]][Diviation_models_all[["actv_ace"]]["model_type"] =="additive","diviation"], y = NULL)
t.test(Diviation_models_all[["actv_ace"]][Diviation_models_all[["actv_ace"]]["model_type"] =="multiplicative","diviation"], y = NULL)
t.test(Diviation_models_all[["actv_ace"]][Diviation_models_all[["actv_ace"]]["model_type"] =="dominative","diviation"], y = NULL)
sum((Diviation_models_all[["actv_ace"]][Diviation_models_all[["actv_ace"]]["model_type"] =="additive","diviation"])^2,na.rm = TRUE)
sum((Diviation_models_all[["actv_ace"]][Diviation_models_all[["actv_ace"]]["model_type"] =="multiplicative","diviation"])^2)
sum((Diviation_models_all[["actv_ace"]][Diviation_models_all[["actv_ace"]]["model_type"] =="dominative","diviation"])^2)


t.test(Diviation_models_all[["actv_cello"]][Diviation_models_all[["actv_cello"]]["model_type"] =="additive","diviation"], y = NULL)
t.test(Diviation_models_all[["actv_cello"]][Diviation_models_all[["actv_cello"]]["model_type"] =="multiplicative","diviation"], y = NULL)
t.test(Diviation_models_all[["actv_cello"]][Diviation_models_all[["actv_cello"]]["model_type"] =="dominative","diviation"], y = NULL)
sum((Diviation_models_all[["actv_cello"]][Diviation_models_all[["actv_cello"]]["model_type"] =="additive","diviation"])^2,na.rm = TRUE)
sum((Diviation_models_all[["actv_cello"]][Diviation_models_all[["actv_cello"]]["model_type"] =="multiplicative","diviation"])^2)
sum((Diviation_models_all[["actv_cello"]][Diviation_models_all[["actv_cello"]]["model_type"] =="dominative","diviation"])^2)


t.test(Diviation_models_all[["actv_gluco"]][Diviation_models_all[["actv_gluco"]]["model_type"] =="additive","diviation"], y = NULL)
t.test(Diviation_models_all[["actv_gluco"]][Diviation_models_all[["actv_gluco"]]["model_type"] =="multiplicative","diviation"], y = NULL)
t.test(Diviation_models_all[["actv_gluco"]][Diviation_models_all[["actv_gluco"]]["model_type"] =="dominative","diviation"], y = NULL)
sum((Diviation_models_all[["actv_gluco"]][Diviation_models_all[["actv_gluco"]]["model_type"] =="additive","diviation"])^2,na.rm = TRUE)
sum((Diviation_models_all[["actv_gluco"]][Diviation_models_all[["actv_gluco"]]["model_type"] =="multiplicative","diviation"])^2)
sum((Diviation_models_all[["actv_gluco"]][Diviation_models_all[["actv_gluco"]]["model_type"] =="dominative","diviation"])^2)



t.test(Diviation_models_all[["actv_phos"]][Diviation_models_all[["actv_phos"]]["model_type"] =="additive","diviation"], y = NULL)
t.test(Diviation_models_all[["actv_phos"]][Diviation_models_all[["actv_phos"]]["model_type"] =="multiplicative","diviation"], y = NULL)
t.test(Diviation_models_all[["actv_phos"]][Diviation_models_all[["actv_phos"]]["model_type"] =="dominative","diviation"], y = NULL)
sum((Diviation_models_all[["actv_phos"]][Diviation_models_all[["actv_phos"]]["model_type"] =="additive","diviation"])^2,na.rm = TRUE)
sum((Diviation_models_all[["actv_phos"]][Diviation_models_all[["actv_phos"]]["model_type"] =="multiplicative","diviation"])^2)
sum((Diviation_models_all[["actv_phos"]][Diviation_models_all[["actv_phos"]]["model_type"] =="dominative","diviation"])^2)
