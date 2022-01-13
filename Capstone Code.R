#Loading Required Libraries
library(matrixcalc)
library(matlib)
library(sampling)
library(coxed)

#Primary Loop for Data Creation
for(t in 1:tot){
  
  for(i in 1:ns){
    frame[i,1]=30+7*rnorm(1)
    frame[i,2]=2+.1*frame[i,1]+.5*rnorm(1)
    frame[i,3]=5+15*frame[i,2]+2*rnorm(1)
    frame[i,4]=2*frame[i,3]+9*rnorm(1)
    frame[i,5]=.4*frame[i,4]+10*rnorm(1)
  }
  
  cov_storage[[t]]=cov(frame)
  mean_storage[[t]]=colMeans(frame)
  total_data[[t]]=frame
  
}

#Computation for full population
for(a in 1:tot){
  inv_cov_storage[[a]]=inv(cov_storage[[a]])
}

running_sum_inv=0

for(b in 1:tot){
  running_sum_inv=running_sum_inv+inv_cov_storage[[b]]
}

final_inv=inv(running_sum_inv)

running_sum_2=0

for(c in 1:tot){
  running_sum_2=running_sum_2+(inv_cov_storage[[c]]%*%mean_storage[[c]])
}

mu_hat_est_pop=final_inv%*%running_sum_2

#Function Start
standard_error2=function(plan_type,n){
  set.seed(15)
  #Start Time Measurement
  start_time=Sys.time()
  
  #Variable Declaration
  B=1000
  
  #Storage Creation
  test_mean=list()
  test_cov=list()
  mean_store=list()
  cov_store=list()
  mu_hat_est_store=list()
  inclu_pass=unlist(lapply(mean_storage,'[[',5))
  
  if(plan_type==1){
    #Sample to be used for Bootstrap
    boot_samp=sample(total_data,n,replace=FALSE)
  }
  
  if(plan_type==2){
    if(n==100){
      #Sample to be used for Bootstrap
      total_data_2=sample(total_data,tot,replace=FALSE)
      boot_samp=total_data_2[seq(1,length(total_data_2),length.out=tot/10000)]
    }
    
    if(n==500){
      #Sample to be used for Bootstrap
      total_data_2=sample(total_data,tot,replace=FALSE)
      boot_samp=total_data_2[seq(1,length(total_data_2),length.out=tot/2000)]
    }
    
    if(n==1000){
      #Sample to be used for Bootstrap
      total_data_2=sample(total_data,tot,replace=FALSE)
      boot_samp=total_data_2[seq(1,length(total_data_2),length.out=tot/1000)]
    }
    
    if(n==2500){
      #Sample to be used for Bootstrap
      total_data_2=sample(total_data,tot,replace=FALSE)
      boot_samp=total_data_2[seq(1,length(total_data_2),length.out=tot/400)]
    }
    
    if(n==5000){
      #Sample to be used for Bootstrap
      total_data_2=sample(total_data,tot,replace=FALSE)
      boot_samp=total_data_2[seq(1,length(total_data_2),length.out=tot/200)]
    }
  }
  
  if(plan_type==3){
    #Sample to be used for Bootstrap
    gc()
    memory.limit(999999999)
    inclu_pass=as.numeric(unlist(lapply(mean_storage,'[[',5)))
    probs=inclusionprobabilities(inclu_pass,n)
    pois_index=UPsampford(probs,max_iter = 1000000)
    pois_index=which(pois_index==1)
    pois_index=pois_index[1:n]
    boot_samp=total_data[c(pois_index)]
    gc()
  }
  
  if(plan_type==4){
    #Sample to be used for Bootstrap
    gc()
    memory.limit(999999999)
    inclu_pass=unlist(lapply(mean_storage,'[[',5))
    pivot_index=UPrandompivotal(inclusionprobabilities(inclu_pass,n+1))
    pivot_index=which(pivot_index==1)
    pivot_index=pivot_index[1:n]
    boot_samp=total_data[c(pivot_index)]
    gc()
  }
  
  if(plan_type==5){
    #Sample to be used for Bootstrap
    gc()
    memory.limit(999999999)
    inclu_pass=unlist(lapply(mean_storage,'[[',5))
    brewer_index=UPbrewer(inclusionprobabilities(inclu_pass,n))
    brewer_index=which(brewer_index==1)
    brewer_index=brewer_index[1:n]
    boot_samp=total_data[c(brewer_index)]
    gc()
  }  
  
  if(plan_type==6){
    index=rep(0,tot)
    inclu_pass=unlist(lapply(mean_storage,'[[',5))
    probs=inclusionprobabilities(inclu_pass,n)
    draw2=runif(tot)
    for(i in 1:tot){
      if(draw2[[i]]>probs[[i]]){
        index[[i]]=0
      }
      else{
        index[[i]]=1
      }
    }
    if(sum(index)!=n){
      lambda=nth((draw2/probs),n)
      for(j in 1:tot){
        if(draw2[[j]]>(lambda*probs[[j]])){
          index[[j]]=0
        }
        else{
          index[[j]]=1
        }
      }
    }
    pois_index=which(index==1)
    boot_samp=total_data[c(pois_index)]
  }
  
  #Loop to store mean and co-variances
  for(j in 1:n){
    test_mean[[j]]=colMeans(boot_samp[[j]])
    test_cov[[j]]=cov(boot_samp[[j]])
  }
  
  #Clearing inverse co-variance storage
  inv_cov_storage=list()
  
  #Loop to invert co-variances
  for(a in 1:n){
    inv_cov_storage[[a]]=inv(test_cov[[a]])
  }
  
  #Clearing auxiliary co-variance running sum  
  running_sum_inv=0
  
  #Loop for summation of inverse co-variances  
  for(b in 1:n){
    running_sum_inv=running_sum_inv+inv_cov_storage[[b]]
  }
  
  #Final inverse of the inverse of sums (for least squares procedure)  
  final_inv=inv(running_sum_inv)
  
  #Clearing auxiliary co-variance*mean running sum
  running_sum_2=0
  
  #Loop for sum of the multiplication of inverse co-variances and mean  
  for(c in 1:n){
    running_sum_2=running_sum_2+(inv_cov_storage[[c]]%*%test_mean[[c]])
  }
  
  #Final estimate multiplication  
  test_est=final_inv%*%running_sum_2
  
  #Primary Bootstrap Loop
  for(i in 1:B){
    boot=sample(boot_samp,n,replace=TRUE)
    
    #Loop to store mean and co-variances
    for(j in 1:n){
      mean_store[[j]]=colMeans(boot[[j]])
      cov_store[[j]]=cov(boot[[j]])
    }
    
    #Clearing inverse co-variance storage
    inv_cov_storage=list()
    
    #Loop to invert co-variances
    for(a in 1:n){
      inv_cov_storage[[a]]=inv(cov_store[[a]])
    }
    
    #Clearing auxiliary co-variance running sum
    running_sum_inv=0
    
    #Loop for summation of inverse co-variances  
    for(b in 1:n){
      running_sum_inv=running_sum_inv+inv_cov_storage[[b]]
    }
    
    #Final inverse of the inverse of sums (for least squares procedure)  
    final_inv=inv(running_sum_inv)
    
    #Clearing auxiliary co-variance*mean running sum
    running_sum_2=0
    
    #Loop for sum of the multiplication of inverse co-variances and mean  
    for(c in 1:n){
      running_sum_2=running_sum_2+(inv_cov_storage[[c]]%*%mean_store[[c]])
    }
    
    #Final estimate multiplication and storage  
    mu_hat_est_store[[i]]=final_inv%*%running_sum_2
  }
  
  #Clearing SE Aux Pass
  se_1=rep(0,B)
  se_2=rep(0,B)
  se_3=rep(0,B)
  se_4=rep(0,B)
  se_5=rep(0,B)
  
  #Loop for Sum of Squares
  for(r in 1:B){
    se_1[r]=mu_hat_est_store[[r]][[1]]
    se_2[r]=mu_hat_est_store[[r]][[2]]
    se_3[r]=mu_hat_est_store[[r]][[3]]
    se_4[r]=mu_hat_est_store[[r]][[4]]
    se_5[r]=mu_hat_est_store[[r]][[5]]
  }
  
  #BCA Interval
  interval_1=bca(se_1)
  interval_2=bca(se_2)
  interval_3=bca(se_3)
  interval_4=bca(se_4)
  interval_5=bca(se_5)
  
  #Finding sample standard deviation
  se_1=sd(se_1)/sqrt(n)
  se_2=sd(se_2)/sqrt(n)
  se_3=sd(se_3)/sqrt(n)
  se_4=sd(se_4)/sqrt(n)
  se_5=sd(se_5)/sqrt(n)
  
  #Finish Time Measurement
  finish_time=Sys.time()
  execution_time=finish_time-start_time
  execution_time
  
  #Function return
  my_list=list(se_1,se_2,se_3,se_4,se_5,execution_time,test_est,interval_1,
               interval_2,interval_3,interval_4,interval_5)
  return(my_list)