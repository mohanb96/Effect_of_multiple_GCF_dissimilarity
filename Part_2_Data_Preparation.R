################################################
###.         Part 2 Data preparetion         ###
################################################

# 1. Import data from 2, 5, 8 factor combination groups.

library(tidyverse)
library(ggplot2)

f2=read.csv("/Data/f2_dissimilarity.csv")
f5=read.csv("/Data/f5_dissimilarity.csv")
f8=read.csv("/Data/f8_dissimilarity.csv")

# 2. Transfer factor label number into 0/1 matrix

# 2.1 Transfer factor label number into 0/1 matrix of 2 factor group.

columns <- c(rep(0,47))
f2=mutate(f2, P=columns,C=columns,L=columns,N=columns,A=columns,I=columns,SU=columns,FU=columns,H=columns,M=columns,S=columns,D=columns)
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
r <- c(rep(2,47))
f2 <- f2%>%
  select(2:8,11:23)%>%
  mutate(remark=r,Lv=r)

# 2.2 Transfer factor label number into 0/1 matrix of 5 factor group.

columns <- c(rep(0,50))
f5=mutate(f5, P=columns,C=columns,L=columns,N=columns,A=columns,I=columns,SU=columns,FU=columns,H=columns,M=columns,S=columns,D=columns)
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
r <- c(rep(5,50))
f5 <- f5%>%
  select(2:8,14:26)%>%
  mutate(remark=r,Lv=r)


# 2.3 Transfer factor label number into 0/1 matrix of 8 factor group.

columns <- c(rep(0,50))
f8 <- mutate(f8, P=columns,C=columns,L=columns,N=columns,A=columns,I=columns,SU=columns,FU=columns,H=columns,M=columns,S=columns,D=columns)
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
r <- c(rep(8,50))
f8 <- f8%>%
  select(2:8,17:29)%>%
  mutate(remark=r,Lv=r)

# 3. Normalize dissimilarity data in each group

library(caret)

f2_dis <-as.data.frame(f2$dissim)
colnames(f2_dis) <-c("dissim")
pro_f2 <-preProcess(f2_dis, method = c("range"))
f2_dis <-predict(pro_f2, f2_dis)
column <-vector()
column <-f2_dis$dissim
f2<-f2%>%
  mutate(dissim = column)

f5_dis<-as.data.frame(f5$dissim)
colnames(f5_dis)<-c("dissim")
pro_f5<-preProcess(f5_dis, method = c("range"))
f5_dis<-predict(pro_f5, f5_dis)
column <-vector()
column <-f5_dis$dissim
f5<-f5%>%
  mutate(dissim = column)

f8_dis<-as.data.frame(f8$dissim)
colnames(f8_dis)<-c("dissim")
pro_f8<-preProcess(f8_dis, method = c("range"))
f8_dis<-predict(pro_f8, f8_dis)
column <-vector()
column <-f8_dis$dissim
f8<-f8%>%
  mutate(dissim = column)

# 4. Import single factor data

sf=read.csv("/Users/mohanbi/Documents/GitHub/Effect_of_multiple_GCF_dissimilarity/Data/singal_factor.csv")
column <- c(rep(0,123))
r=c(rep(1,123))
sf=sf%>%
  mutate(dissim=column,P=column,C=column,L=column,N=column,A=column,I=column,SU=column,FU=column,H=column,M=column,S=column,D=column,remark=r,Lv=r)

# 5. Single factor data, as the group of factor level "1".

for (i in 1:nrow(sf)) {
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

for (i in 1:nrow(sf)) {
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
