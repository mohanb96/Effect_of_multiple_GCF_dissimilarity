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
