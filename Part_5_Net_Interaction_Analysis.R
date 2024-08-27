#####################################################
###.    Part 5 Factor net interaction analysis    ###
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
# Tutorial: https://mohanb96.github.io/Nullmodels.html

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

Deviation_additive = data.frame(matrix(data = NA, nrow = 147, ncol = 6))
colnames(Deviation_additive)=c("Lv","Deviation","P","dissim","actual_ES","null_ES")

Deviation_multiplicative = data.frame(matrix(data = NA, nrow = 147, ncol = 6))
colnames(Deviation_multiplicative)=c("Lv","Deviation","P","dissim","actual_ES","null_ES")

Deviation_dominative = data.frame(matrix(data = NA, nrow = 147, ncol = 6))
colnames(Deviation_dominative)=c("Lv","Deviation","P","dissim","actual_ES","null_ES")

Deviation_3_models = list()

for (i_response in responses) {
  
  for (treatment_i in 1:nrow(df_H)) {
    combination_i=c()
    combination_i=stressors[which(df_H[treatment_i,10:21]==1)]
    Null_modle_i_Treatment = NullModel(ES.plot.meandiff[[i_response]], selected_factors = combination_i,n_perm = n_perm)
    Null_modle_summary_i_Treatment = NullModel_summary(Null_modle_i_Treatment, actual_data = df_H[treatment_i,i_response]) #df_H[treatment_i,6] refere to decomposition rate
    
    Deviation_additive[treatment_i,2]=(Null_modle_summary_i_Treatment[1,3]-Null_modle_summary_i_Treatment[2,3])/Null_modle_summary_i_Treatment[2,3]  #changed
    Deviation_additive[treatment_i,3]=Null_modle_summary_i_Treatment[2,5]
    Deviation_additive[treatment_i,4]=df_H[treatment_i,9]
    Deviation_additive[treatment_i,1]=df_H[treatment_i,23]
    Deviation_additive[treatment_i,5]=Null_modle_summary_i_Treatment[1,3]-mean(df_CK[,i_response])
    Deviation_additive[treatment_i,6]=Null_modle_summary_i_Treatment[2,3]-mean(df_CK[,i_response])
    
    Deviation_multiplicative[treatment_i,2]=(Null_modle_summary_i_Treatment[1,3]-Null_modle_summary_i_Treatment[3,3])/Null_modle_summary_i_Treatment[3,3]
    Deviation_multiplicative[treatment_i,3]=Null_modle_summary_i_Treatment[3,5]
    Deviation_multiplicative[treatment_i,4]=df_H[treatment_i,9]
    Deviation_multiplicative[treatment_i,1]=df_H[treatment_i,23]
    Deviation_multiplicative[treatment_i,5]=Null_modle_summary_i_Treatment[1,3]-mean(df_CK[,i_response])
    Deviation_multiplicative[treatment_i,6]=Null_modle_summary_i_Treatment[3,3]-mean(df_CK[,i_response])
    
    Deviation_dominative[treatment_i,2]=(Null_modle_summary_i_Treatment[1,3]-Null_modle_summary_i_Treatment[4,3])/Null_modle_summary_i_Treatment[4,3]
    Deviation_dominative[treatment_i,3]=Null_modle_summary_i_Treatment[4,5]
    Deviation_dominative[treatment_i,4]=df_H[treatment_i,9]
    Deviation_dominative[treatment_i,1]=df_H[treatment_i,23]
    Deviation_dominative[treatment_i,5]=Null_modle_summary_i_Treatment[1,3]-mean(df_CK[,i_response])
    Deviation_dominative[treatment_i,6]=Null_modle_summary_i_Treatment[4,3]-mean(df_CK[,i_response])
  }
  
  Deviation_3_models[[i_response]][["additive"]]=Deviation_additive
  Deviation_3_models[[i_response]][["multiplicative"]]=Deviation_multiplicative
  Deviation_3_models[[i_response]][["dominative"]]=Deviation_dominative
  
}

# 5.1 Identifying interaction type

positive_responses = c("actv_ace","actv_cello","actv_gluco","actv_phos","PH")
negativ_responses = c("Decom", "WSA")
null_models = c("additive", "multiplicative", "dominative")

for (i_response in negativ_responses) {
  for (n_model in null_models) {
    for (i in 1:147) {
      if(Deviation_3_models[[i_response]][[n_model]][["Deviation"]][i]>0 ) 
        Deviation_3_models[[i_response]][[n_model]][["interaction_type"]][i] = "Antagonistic"
      if(Deviation_3_models[[i_response]][[n_model]][["Deviation"]][i]<0 ) 
        Deviation_3_models[[i_response]][[n_model]][["interaction_type"]][i] = "Synergistic"
      if(abs(Deviation_3_models[[i_response]][[n_model]][["Deviation"]][i])<0.0 ) 
        Deviation_3_models[[i_response]][[n_model]][i,"interaction_type"] = n_model ### threshold adjustment
      if(Deviation_3_models[[i_response]][[n_model]][i,"P"]>=0.05) 
        Deviation_3_models[[i_response]][[n_model]][i,"interaction_type"] = n_model
    }
  }
}

for (i_response in positive_responses) {
  for (n_model in null_models) {
    for (i in 1:147) {
      if(Deviation_3_models[[i_response]][[n_model]][["Deviation"]][i]>0 ) 
        Deviation_3_models[[i_response]][[n_model]][["interaction_type"]][i] = "Synergistic"
      if(Deviation_3_models[[i_response]][[n_model]][["Deviation"]][i]<0 ) 
        Deviation_3_models[[i_response]][[n_model]][["interaction_type"]][i] = "Antagonistic"
      if(abs(Deviation_3_models[[i_response]][[n_model]][["Deviation"]][i])<0.0 ) 
        Deviation_3_models[[i_response]][[n_model]][i,"interaction_type"] = n_model ### threshold adjustment
      if(Deviation_3_models[[i_response]][[n_model]][i,"P"]>=0.05) 
        Deviation_3_models[[i_response]][[n_model]][i,"interaction_type"] = n_model
    }
  }
}

# Plotting (dot plots)

p_Deviation = list()
p_Deviation_lv = list()
p_Deviation_all = list()

for (i_response in responses) {
  
  p_Deviation[[i_response]][["additive"]] = ggplot(Deviation_3_models[[i_response]][["additive"]], aes(x=dissim, y=Deviation))+
    geom_point(size =2, aes(color = interaction_type))+
    scale_color_manual(values = c("Antagonistic"="#5499C780","additive"="#90949780","Synergistic"="#F1666780"))+
    geom_smooth(method = lm, aes(color = "#999999"))+
    stat_cor(method = "spearman",label.x = 0)+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  ylim(-1,20)
  
  p_Deviation[[i_response]][["multiplicative"]] = ggplot(Deviation_3_models[[i_response]][["multiplicative"]], aes(x=dissim, y=Deviation))+
    geom_point(size =2, aes(color = interaction_type))+
    scale_color_manual(values = c("Antagonistic"="#5499C780","multiplicative"="#90949780","Synergistic"="#F1666780"))+
    geom_smooth(method = lm, aes(color = "#999999"))+
    stat_cor(method = "spearman",label.x = 0)+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  ylim(-5,5)
  
  p_Deviation[[i_response]][["dominative"]] = ggplot(Deviation_3_models[[i_response]][["dominative"]], aes(x=dissim, y=Deviation))+
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
  p_Deviation_lv[[i_response]][["additive"]] = 
    ggplot()+
    geom_jitter(data = Deviation_3_models[[i_response]][["additive"]], aes(x=Lv, y=Deviation, color = interaction_type),width = 0.5, size = 1)+
    geom_boxplot(data = Deviation_3_models[[i_response]][["additive"]], aes( x=Lv, y=Deviation, group = Lv), width = 1, fill = "#00000000", color = "#72665880")+
    scale_x_continuous(breaks = c(2,5,8))+
    scale_color_manual(values = c("Antagonistic"="#5499C780","additive"="#90949780","Synergistic"="#F1666780"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    stat_compare_means(method = "t.test")+
    xlab(i_response)
  
  p_Deviation_lv[[i_response]][["multiplicative"]] = 
    ggplot()+
    geom_jitter(data = Deviation_3_models[[i_response]][["multiplicative"]], aes(x=Lv, y=Deviation, color = interaction_type),width = 0.5, size = 1)+
    geom_boxplot(data = Deviation_3_models[[i_response]][["multiplicative"]], aes( x=Lv, y=Deviation, group = Lv), width = 1, fill = "#00000000", color = "#72665880")+
    scale_x_continuous(breaks = c(2,5,8))+
    scale_color_manual(values = c("Antagonistic"="#5499C780","multiplicative"="#90949780","Synergistic"="#F1666780"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  
  p_Deviation_lv[[i_response]][["dominative"]] = 
    ggplot()+
    geom_jitter(data = Deviation_3_models[[i_response]][["dominative"]], aes(x=Lv, y=Deviation, color = interaction_type),width = 0.5, size = 1)+
    geom_boxplot(data = Deviation_3_models[[i_response]][["dominative"]], aes( x=Lv, y=Deviation, group = Lv), width = 1, fill = "#00000000", color = "#72665880")+
    scale_x_continuous(breaks = c(2,5,8))+
    scale_color_manual(values = c("Antagonistic"="#5499C780","dominative"="#90949780","Synergistic"="#F1666780"))+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
}

Deviation_models_all = list()

for (i_response in responses) {
  for (n_model in null_models) {
    Deviation_3_models[[i_response]][[n_model]]["model_type"] = n_model
  }
  Deviation_models_all[[i_response]] = rbind(Deviation_3_models[[i_response]][["additive"]],Deviation_3_models[[i_response]][["multiplicative"]],Deviation_3_models[[i_response]][["dominative"]])
}

for (i_response in responses) {
  p_Deviation_all[[i_response]] = local({
    i_response = i_response
    Deviation_models_all = Deviation_models_all
    Deviation_models_all[[i_response]] = data.frame(Deviation_models_all[[i_response]])%>%
      arrange(model_type)%>%
      mutate(model_type = factor(model_type, levels = c("additive","multiplicative","dominative")))
    
    ggplot(data = Deviation_models_all[[i_response]], aes(x = model_type, y=Deviation))+
      geom_jitter(aes(color = interaction_type), width = 0.2, size = 1)+
      scale_color_manual(values = c("Antagonistic"="#5499C780","additive"="#90949780","multiplicative"="#90949780","dominative"="#90949780","Synergistic"="#F1666780"))+
      geom_boxplot(aes(x = model_type, y=Deviation, group = model_type), width = 0.5, fill = "#00000000", color = "#72665880")+
      scale_x_discrete("model_type", labels = c("additive","multiplicative","dominative"))+
      theme_bw()+
      theme(legend.position = 'none',axis.title.x = element_blank())+
      xlab(i_response)
  })
}

p_Deviation[["Decom"]][["additive"]]+p_Deviation[["Decom"]][["multiplicative"]]+p_Deviation[["Decom"]][["dominative"]]+
  p_Deviation[["PH"]][["additive"]]+p_Deviation[["PH"]][["multiplicative"]]+p_Deviation[["PH"]][["dominative"]]+
  p_Deviation[["WSA"]][["additive"]]+p_Deviation[["WSA"]][["multiplicative"]]+p_Deviation[["WSA"]][["dominative"]]+
  plot_spacer()+
  plot_layout(ncol=3, widths=c(2,2,2))

p_Deviation[["actv_ace"]][["additive"]]+p_Deviation[["actv_ace"]][["multiplicative"]]+p_Deviation[["actv_ace"]][["dominative"]]+
  p_Deviation[["actv_cello"]][["additive"]]+p_Deviation[["actv_cello"]][["multiplicative"]]+p_Deviation[["actv_cello"]][["dominative"]]+
  p_Deviation[["actv_gluco"]][["additive"]]+p_Deviation[["actv_gluco"]][["multiplicative"]]+p_Deviation[["actv_gluco"]][["dominative"]]+
  p_Deviation[["actv_phos"]][["additive"]]+p_Deviation[["actv_phos"]][["multiplicative"]]+p_Deviation[["actv_phos"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))

p_Deviation_lv[["Decom"]][["additive"]]+p_Deviation_lv[["Decom"]][["multiplicative"]]+p_Deviation_lv[["Decom"]][["dominative"]]+
  p_Deviation_lv[["PH"]][["additive"]]+p_Deviation_lv[["PH"]][["multiplicative"]]+p_Deviation_lv[["PH"]][["dominative"]]+
  p_Deviation_lv[["WSA"]][["additive"]]+p_Deviation_lv[["WSA"]][["multiplicative"]]+p_Deviation_lv[["WSA"]][["dominative"]]+
  plot_spacer()+
  plot_layout(ncol=3, widths=c(2,2,2))

p_Deviation_lv[["actv_ace"]][["additive"]]+p_Deviation_lv[["actv_ace"]][["multiplicative"]]+p_Deviation_lv[["actv_ace"]][["dominative"]]+
  p_Deviation_lv[["actv_cello"]][["additive"]]+p_Deviation_lv[["actv_cello"]][["multiplicative"]]+p_Deviation_lv[["actv_cello"]][["dominative"]]+
  p_Deviation_lv[["actv_gluco"]][["additive"]]+p_Deviation_lv[["actv_gluco"]][["multiplicative"]]+p_Deviation_lv[["actv_gluco"]][["dominative"]]+
  p_Deviation_lv[["actv_phos"]][["additive"]]+p_Deviation_lv[["actv_phos"]][["multiplicative"]]+p_Deviation_lv[["actv_phos"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))
