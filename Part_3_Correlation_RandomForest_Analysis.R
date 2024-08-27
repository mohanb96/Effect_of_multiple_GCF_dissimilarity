#####################################################
###. Part 3 Effect sizes & Correlation analysis  ###
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

##################################################
###    Effect size and Null model analysis     ###
##################################################

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

#########################################################################################
###     Permutation based Random Forest for calculating predictor important values    ###
#########################################################################################

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
