library(AER)
Data<-Affairs
index<-which(Data$affairs==0) 

# Count left-censored observations
left_censored_count<-length(index)
censoring<-(left_censored_count/601)*100

## first stage regression 
Data$res<-residuals(lm(rating~age, data=Data))

### To Check Endogenity
# Create binary children variable
Affairs$children_b <- ifelse(Affairs$children == "yes", 1, 0)
Affairs$gender_b<- ifelse(Affairs$gender == "male", 1, 0)


### optimize the function using least absolute loss function
ABL<-function(beta,x){
  g_f<-exp(x$yearsmarried*beta[1]+x$gender_b*beta[2]+x$children_b*beta[3]+x$occupation*beta[4]+x$age*beta[5]+x$rating*beta[6])+x$res*beta[7]
  pred<-pmax(0,g_f)
  L<-sum(abs(x$affairs-pred))
  return(L)
}
est_abl<-optim(c(0,0,0,0,0,0,0.1), function(u) ABL(u, Data))$par


##optimize Huber loss
###define absolute loss function 
HBL <- function(beta, x, w = 1.35) {
  # Calculate the residuals
  g_f<-exp(x$yearsmarried*beta[1]+x$gender_b*beta[2]+x$children_b*beta[3]+x$occupation*beta[4]+x$age*beta[5]+x$rating*beta[6])+x$res*beta[7]
  pred<-pmax(0,g_f)
  abs_diff <- abs(x$affairs -pred)
  
  # Apply the Huber loss formula
  loss <- ifelse(
    abs(abs_diff) <= w,
    0.5 * abs_diff^2,
    w * (abs(abs_diff) - 0.5 * w)
  )
  
  # Return the mean loss
  return(sum(loss))
}

est_hbl<-optim(c(0,0,0,0,0,0,0.1), function(u) HBL(u, Data))$par


###optimzie CLCE loss
CLH <- function(beta, x) {
  g_f<-exp(x$yearsmarried*beta[1]+x$gender_b*beta[2]+x$children_b*beta[3]+x$occupation*beta[4]+x$age*beta[5]+x$rating*beta[6])+x$res*beta[7]
  pred<-pmax(0,g_f)
  L<-sum(log(cosh(x$affairs-pred)))
  return(L)
}
est_clh<-optim(c(0,0,0,0,0,0,0.1), function(u) CLH(u, Data))$par

## Bootstrap samples

B<-seq(100,500,100)

mse_1<-c()
mse_2<-c()
mse_3<-c()
mse_4<-c()
mse_5<-c()
mse_6<-c()
mse_7<-c()

for (b in B) 
{
  bootstrap_sample <- function(data) {
    n <- nrow(data)
    indices <- sample(1:n, replace = TRUE)
    return(data[indices, ])
  }
  
  ## Generate B = 1000 bootstrap samples and store in a list.
  
  Boot_samples<- lapply(1:b, function(ii)bootstrap_sample(Data))
  
  
  es_t<-c()
  ## Bootstrap MSE of tobit parameters 
  for (i in 1:b) 
  {
    bs_i<-Boot_samples[[i]]
    
    ###first stage regression
    bs_i$res<-residuals(lm(rating ~ age, data = bs_i))
    
    ###define absolute loss function 
    ABL<-function(beta,x=bs_i){
      g_f<-exp(x$yearsmarried*beta[1]+x$gender_b*beta[2]+x$children_b*beta[3]+x$occupation*beta[4]+x$age*beta[5]+x$rating*beta[6])+x$res*beta[7]
      pred<-pmax(0,g_f)
      L<-sum(abs(bs_i$affairs-pred))
      return(L)
    }
    es_t[[i]]<-optim(c(0,0,0,0,0,0,0.1),ABL)$par
    
  }
  
  ## Bootstrap MSE
  v1<-sapply(1:b, function(ii)es_t[[ii]][1])
  v2<-sapply(1:b, function(ii)es_t[[ii]][2])
  v3<-sapply(1:b, function(ii)es_t[[ii]][3])
  v4<-sapply(1:b, function(ii)es_t[[ii]][4])
  v5<-sapply(1:b, function(ii)es_t[[ii]][5])
  v6<-sapply(1:b, function(ii)es_t[[ii]][6])
  v7<-sapply(1:b, function(ii)es_t[[ii]][7])
  
  
  
  ## MSE 
  mse_1<-append(mse_1,mean((v1 - est_abl[1])^2),after = length(mse_1))
  mse_2<-append(mse_2,mean((v2 - est_abl[2])^2),after = length(mse_2))
  mse_3<-append(mse_3,mean((v2 - est_abl[3])^2),after = length(mse_3))
  mse_4<-append(mse_4,mean((v2 - est_abl[4])^2),after = length(mse_4))
  mse_5<-append(mse_5,mean((v2 - est_abl[5])^2),after = length(mse_5))
  mse_6<-append(mse_6,mean((v2 - est_abl[6])^2),after = length(mse_6))
  mse_7<-append(mse_7,mean((v2 - est_abl[7])^2),after = length(mse_7))
}
mse_1
mse_2
mse_3
mse_4
mse_5
mse_6
mse_7


