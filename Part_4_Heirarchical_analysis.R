############################################################################################
### Hierarchical linear model framework (alternative for random forest model analysis)   ###
############################################################################################

### 1.Decomposition

df.rf.decom = cbind(df.rf, Pred_respon_for_each_all[["Decom"]])
df.multiple = df.rf.decom[df.rf.decom[,"Lv"]%in%h_levels,]

m_R_decom <- lm(Decom ~ R1+R2+R3,
                data = df.multiple)
m_lv_decom <- lm(Decom ~ R1+R2+R3+Lv,
                 data = df.multiple)
m_L_decom <- lm(Decom ~ Lv,
                data = df.multiple)
m_D_decom <- lm(Decom ~ dissim,
                data = df.multiple)
m_dissim_decom <- lm(Decom ~ R1+R2+R3+dissim,
                     data = df.multiple)

m_all_decom <- lm(Decom ~ R1+R2+R3+Lv+dissim,
                  data = df.multiple)
m_all_sig_decom<- lm(Decom ~ R1+R2+R3+Lv+dissim+P+C+L+N+A+I+SU+FU+H+M+S+D,
                     data = df.multiple)

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
