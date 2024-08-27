# Functions for calculating null model predictions for specific factor combinations
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
