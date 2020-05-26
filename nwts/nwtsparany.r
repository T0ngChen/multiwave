library(rjags)
library(doParallel)
library(foreach)
library(survey)

########################
##                    ##
## influence function ##
##                    ## 
########################

inf.fun <- function(fit) {
  dm <- model.matrix(fit)
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values)))
  ## influence function
  infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
  infl
}

inf.prior = function(dm, para.val){
  fit.val = exp(rowSums(dm %*% diag(para.val))) / (1+exp(rowSums(dm %*% diag(para.val))))
  Ihat <- (t(dm) %*% (dm * fit.val * (1 - fit.val)))
  infl = as.matrix(dm * (nwts$relaps - fit.val)) %*% solve(Ihat)
  infl
}


########################
##                    ##
## Neyman allocation  ##
##                    ## 
########################

integer.neyman.w1 = function(n.strata, NS, sample.size, upper){
  nc = max(upper+1)
  s = 1:nc
  arr = matrix(rep(NS, each = nc)/sqrt(s*(s+1)), nrow = n.strata, byrow = T)
  arr.list = as.list(as.data.frame(t(arr)))
  for(i in 1:length(arr.list)){
    arr.list[[i]][upper[i]:nc] = 0
  }
  arr = do.call(cbind, arr.list)[-1,]
  rk <- order(arr, na.last=TRUE, decreasing=TRUE)[1:(sample.size - 2 * n.strata)]
  re.om.zero = table(rk%/%(nc-1) + 1)
  re.zero = rep(2, n.strata)
  re.zero[1:n.strata %in% names(re.om.zero)] = 2 + re.om.zero 
  re.zero
}


integer.neyman.w2 = function(n.strata, NS, sample.size, upper){
  nc = max(upper + 1)
  s = 1:nc
  arr = matrix(rep(NS, each = nc)/sqrt(s*(s+1)), nrow = n.strata, byrow = T)
  arr.list = as.list(as.data.frame(t(arr)))
  for(i in 1:length(arr.list)){
    arr.list[[i]][upper[i]:nc] = 0
  }
  arr = do.call(cbind, arr.list)
  rk <- order(arr, na.last=TRUE, decreasing=TRUE)[1:(sample.size - n.strata)]
  re.om.zero = table(rk%/%nc + 1) + 1
  re.zero = rep(0, n.strata)
  re.zero[1:n.strata %in% names(re.om.zero)] = re.om.zero 
  re.zero[which.max(upper - re.zero)] =  re.zero[which.max(upper - re.zero)] + length(which(re.zero == 0))
  re.zero
}

#########################
##                     ##
## two-wave with prior ##
##                     ##
#########################

bayesian.fit = function(parameters, data, BUGSMODEL){
  # How many burn-in steps?
  burn_in = 1000
  
  # How many proper steps?
  steps = 1500
  
  # Thinning?
  thin = 1
  
  #compilation of the BUGS model, no Gibbs sampling yet
  foo <- jags.model(textConnection(BUGSMODEL),data=data, n.chains=4) 
  
  #burnin samples
  update(foo, burn_in)
  
  #draws n.iter MCMC samples, monitors the parameters specified in variable.names, thins the
  #output by a factor of thin and stores everything in a MCMC.list object
  out <- coda.samples(model=foo, variable.names=parameters, n.iter=steps, thin=thin)
  
  out
}

wave.sample.size = function(dm, para.val){
  infl.prior = inf.prior(dm, para.val)
  infl.wave = infl.prior[,2]
  sd.wave = numeric()
  for(i in 1:length(strata)){
    sd.wave = c(sd.wave, sd(infl.wave[strata[[i]]]))
  }
  sd.wave
}

two.wave.bayes = function(piror_info, prior_pcs, n, sl, parameters, para.outcome, data, n1, strata){
  imphisw1 = rbinom(nrow(data), 1, piror_info[,1][!rownames(piror_info) %in% parameters])
  w1.sd = wave.sample.size(dm = cbind(rep(1,nrow(data)), imphisw1, data$age1, data$age2, data$stage1,
                                      data$tumdiam, data$stage1 * data$tumdiam), 
                           para.val = piror_info[,1][para.outcome])
  
  w1.size = integer.neyman.w1(n.strata = length(strata), NS = w1.sd * n1, sample.size = n, upper = n1)
  
  s.w1 = list()
  s.w2 = list()
  for (i in 1:length(strata)) {
    s.w1[[i]] <- sample(strata[[i]], w1.size[i])
  }
  
  fit2 = bayesian.fit(parameters = c(paste0("alpha", 0:5), paste0("beta", 0:6), "histol"),
                      data=list(histol = data$histol[ifelse(1:nrow(data) %in% unlist(s.w1), T, NA)], instit = nwts$instit, relaps = nwts$relaps, stage2 = nwts$stage2, age1 = nwts$age1, 
                                age2 = nwts$age2, age3 = nwts$age3, stage1 = nwts$stage1, tumdiam = nwts$tumdiam, study = nwts$study, N = nrow(nwts), 
                                malpha0 = piror_info[,1]["alpha0"], palpha0 = prior_pcs["alpha0"], malpha1 = piror_info[,1]["alpha1"], palpha1 = prior_pcs["alpha1"], 
                                malpha2 = piror_info[,1]["alpha2"], palpha2 = prior_pcs["alpha2"], malpha3 = piror_info[,1]["alpha3"], palpha3 = prior_pcs["alpha3"],
                                malpha4 = piror_info[,1]["alpha4"], palpha4 = prior_pcs["alpha4"], malpha5 = piror_info[,1]["alpha5"], palpha5 = prior_pcs["alpha5"], 
                                mbeta0 = piror_info[,1]["beta0"], pbeta0 = prior_pcs["beta0"], mbeta1 = piror_info[,1]["beta1"], pbeta1 = prior_pcs["beta1"], 
                                mbeta2 = piror_info[,1]["beta2"], pbeta2 = prior_pcs["beta2"], mbeta3 = piror_info[,1]["beta3"], pbeta3 = prior_pcs["beta3"],
                                mbeta4 = piror_info[,1]["beta4"], pbeta4 = prior_pcs["beta4"], mbeta5 = piror_info[,1]["beta5"], pbeta5 = prior_pcs["beta5"], 
                                mbeta6 = piror_info[,1]["beta6"], pbeta6 = prior_pcs["beta6"]),
                      BUGSMODEL = BUGSmodel)
  
  w1.post.info = summary(fit2)[[1]][,1:2]
  
  ix = rbinom(nrow(data), 1, w1.post.info[,1][!rownames(w1.post.info) %in% parameters])
  
  w2.sd = wave.sample.size(dm = cbind(rep(1,nrow(data)), ix, data$age1, data$age2, data$stage1,
                                      data$tumdiam, data$stage1 * data$tumdiam), 
                           para.val = w1.post.info[,1][para.outcome])
  
  upper = n1 - sapply(s.w1,length)
  
  w2.size = integer.neyman.w2(n.strata = length(strata), NS = w2.sd  * n1, sample.size = sl - n, upper = upper)
  
  s.w2 = list()
  for (i in 1:length(strata)) {
    if(length(strata[[i]][-match(s.w1[[i]], strata[[i]])]) != 0){
      s.w2[[i]] <- sample(strata[[i]][-match(s.w1[[i]], strata[[i]])], w2.size[i])
    } else if(length(strata[[i]][-match(s.w1[[i]], strata[[i]])]) == 0){
      s.w2[[i]] <- NULL
    }
  }
  p2.ob = mapply(c, s.w1, s.w2, SIMPLIFY=FALSE)
  
  data$in.subsample <- (1:nrow(data)) %in% unlist(p2.ob)
  twophase.w2 <- twophase(id = list(~1, ~1), strata = list(NULL, ~stra), 
                          subset = ~in.subsample, data = data)
  wgt = svyglm(relaps~histol+age1+age2+stage1*tumdiam, family = quasibinomial, design = twophase.w2)
  ## imputation calibration
  impmodel <- svyglm(histol~instit+age3+study*stage2,family=quasibinomial,design=twophase.w2)
  data$impx <- as.vector(predict(impmodel,newdata=data,type="response",se.fit=FALSE))
  phase1model_imp <- glm(relaps~impx+age1+age2+stage1*tumdiam, family=binomial, data=data)
  inffun_imp <- inf.fun(phase1model_imp)
  colnames(inffun_imp)<-paste0("if",1:ncol(inffun_imp))
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, ~stra), 
                              subset = ~in.subsample, data = cbind(data, inffun_imp), method="simple")
  calformula <- make.formula(colnames(inffun_imp))
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking", force = T)
  svyest_imp<-svyglm(relaps~histol+age1+age2+stage1*tumdiam, family=quasibinomial, design=cal_twophase_imp)
  list(wgt = wgt, cal = svyest_imp, p2.ob = p2.ob)
}






##################
##              ## 
##  Mean score  ##         
##              ##
##################


two.wave.ms <- function(strata, prob, w1.n, p2.n, data) {
  s.w1 <- list()
  s.w2 <- list()
  ## sample wave 1
  if(any(round((w1.n * prob)) > n1)){
    w1.sam = round(w1.n * prob)
    w1.sam[which(round((w1.n * prob)) > n1)] = n1[round((w1.n * prob)) > n1]
    w1.sam[which.max(n1)] = w1.n - sum(w1.sam) + w1.sam[which.max(n1)]
  } else {
    w1.sam = round((w1.n * prob))
    w1.sam[which.max(n1)] = w1.n - sum(w1.sam) + w1.sam[which.max(n1)]
  }
  for (i in 1:length(strata)) {
    s.w1[[i]] <- sample(strata[[i]], w1.sam[i])
  }
  data$in.subsample <- (1:nrow(data)) %in% unlist(s.w1)
  
  twophase.w1 <- twophase(id = list(~1, ~1), strata = list(NULL, ~stra), 
                          subset = ~in.subsample, data = data)
  fit.w1 <- glm(relaps~histol+age1+age2+stage1*tumdiam, family = binomial, data = data[data$in.subsample,], weights = 1/twophase.w1$prob)
  
  n2 = xtabs(~stra, data[data$in.subsample,])
  
  
  
  X <- model.matrix(fit.w1)
  pi.hat <- fit.w1$fitted.values
  Ihat <- (t(as.matrix(X)) %*% (as.matrix(X) * 1/twophase.w1$prob * pi.hat * (1 - pi.hat)))/sum(n1)
  wgt<-n1/n2*(n1-n2)
  S <- X * (data[data$in.subsample,]$relaps - pi.hat)
  invI <- solve(Ihat)
  result <- 0
  varsi<-ar<-array(0,dim=c(nrow(Ihat),nrow(Ihat),length(wgt)))
  for(j in 1:length(wgt)) {
    si <- S[data[data$in.subsample,]$stra == j,  ]
    varsi[,,j]<-var(si)
    ar[,,j]<-invI%*%varsi[,,j]%*%invI
  }
  Wzykk=1
  k = 2
  for (i in 1:length(n1)) Wzykk[i]=ar[k,k,i]
  prev<-n1/sqrt(sum(n1))
  numer<-prev*sqrt(Wzykk)
  denom<-sum(prev*sqrt(Wzykk))
  prop<-numer/denom
  
  ## Neyman allocation
  sl <- p2.n - w1.n
  
  ## correcting if the sample size is greater than population
  n3 = numeric()
  for (k in 1:length(strata)){
    n3[k] = length(strata[[k]][-match(s.w1[[k]], strata[[k]])])
  }
  
  id = 1:length(prop)
  index1 = NULL
  ms.ney = prop*sl
  while (any(ms.ney>n3)){
    index1 = c(index1,id[ms.ney > n3])
    numer1 = prev*sqrt(Wzykk)
    denom1<-sum(sqrt(Wzykk[-index1])*prev[-index1])
    prop<-numer1/denom1
    ms.ney[-index1] = prop[-index1]*(sl-sum(n3[index1]))
    ms.ney[index1] = n3[index1]
    #print(prop)
  }
  
  
  for (i in 1:length(strata)) {
    if(length(strata[[i]][-match(s.w1[[i]], strata[[i]])]) != 0){
      s.w2[[i]] <- sample(strata[[i]][-match(s.w1[[i]], strata[[i]])], ms.ney[i])
    } else if(length(strata[[i]][-match(s.w1[[i]], strata[[i]])]) == 0){
      s.w2[[i]] <- NULL
    }
  }
  p2.ob = mapply(c, s.w1, s.w2, SIMPLIFY=FALSE)
  data$insample <- (1:nrow(data)) %in% unlist(p2.ob)
  twophase.w2 <- twophase(id = list(~1, ~1), strata = list(NULL, ~stra), 
                          subset = ~insample, data = data)
  wgt = svyglm(relaps~histol+age1+age2+stage1*tumdiam, family = quasibinomial, design = twophase.w2)
  ## imputation calibration
  impmodel <- svyglm(histol~instit+age3+study*stage2,family=quasibinomial,design=twophase.w2)
  data$impx <- as.vector(predict(impmodel,newdata=data,type="response",se.fit=FALSE))
  phase1model_imp <- glm(relaps~impx+age1+age2+stage1*tumdiam, family=binomial, data=data)
  inffun_imp <- inf.fun(phase1model_imp)
  colnames(inffun_imp)<-paste0("if",1:ncol(inffun_imp))
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, ~stra), 
                              subset = ~insample, data = cbind(data, inffun_imp), method="simple")
  calformula <- make.formula(colnames(inffun_imp))
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking", force = T)
  svyest_imp<-svyglm(relaps~histol+age1+age2+stage1*tumdiam, family=quasibinomial, design=cal_twophase_imp)
  list(wgt = wgt, cal = svyest_imp, p2.ob = p2.ob)
}


############################
##                        ##
## one-wave without prior ##
##                        ##  
############################

one.wave <- function(strata, prob, p2.n, data) {
  s.w1 <- list()
  ## sample wave 1
  if(any(round((p2.n * prob)) > n1)){
    w1.sam = round(p2.n * prob)
    w1.sam[which(round((p2.n * prob)) > n1)] = n1[round((p2.n * prob)) > n1]
    w1.sam[which.max(n1)] = p2.n - sum(w1.sam) + w1.sam[which.max(n1)]
  } else {
    w1.sam = round((p2.n * prob))
    w1.sam[which.max(n1)] = p2.n - sum(w1.sam) + w1.sam[which.max(n1)]
  }
  for (i in 1:length(strata)) {
    s.w1[[i]] <- sample(strata[[i]], w1.sam[i])
  }
  data$insample <- (1:nrow(data)) %in% unlist(s.w1)
  twophase.w2 <- twophase(id = list(~1, ~1), strata = list(NULL, ~stra), 
                          subset = ~insample, data = data)
  wgt = svyglm(relaps~histol+age1+age2+stage1*tumdiam, family = quasibinomial, design = twophase.w2)
  ## imputation calibration
  impmodel <- svyglm(histol~instit+age3+study*stage2,family=quasibinomial,design=twophase.w2)
  data$impx <- as.vector(predict(impmodel,newdata=data,type="response",se.fit=FALSE))
  phase1model_imp <- glm(relaps~impx+age1+age2+stage1*tumdiam, family=binomial, data=data)
  inffun_imp <- inf.fun(phase1model_imp)
  colnames(inffun_imp)<-paste0("if",1:ncol(inffun_imp))
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, ~stra), 
                              subset = ~insample, data = cbind(data, inffun_imp), method="simple")
  calformula <- make.formula(colnames(inffun_imp))
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking", force = T)
  svyest_imp<-svyglm(relaps~histol+age1+age2+stage1*tumdiam, family=quasibinomial, design=cal_twophase_imp)
  list(wgt = wgt, cal = svyest_imp, p2.ob = s.w1)
}




## data
nwts <- read.table("nwts-share.txt", header=TRUE)
nwts$age1 <- with(nwts, pmin(age, 1))
nwts$age2 <- with(nwts, pmax(age, 1))
nwts$age3 <- ifelse(nwts$age > 10, 1, 0)
nwts$stage1 <- ifelse(nwts$stage > 2, 1, 0)
nwts$stage2 <- ifelse(nwts$stage > 3, 1, 0)
nwts$study <- nwts$study - 3

full <- glm(relaps~histol+age1+age2+stage1*tumdiam, family=binomial, data=nwts)
impute <- glm(histol~instit+age3+study*stage2,family=binomial, data=nwts)

## simulation
## define strata
## define strata 
nwts$stra = 1 + 4 * nwts$relaps + 2 * nwts$instit + nwts$study
## population strata size
n1 = xtabs(~stra, nwts)

strata = list()
for (index in 1:length(n1)){
  strata[[index]] = which(nwts$stra == index)  
}


infl <- inf.fun(full)[, 2]

sd.stra1 = numeric()

for(i in 1:length(strata)){
  sd.stra1 = c(sd.stra1, sd(infl[strata[[i]]]))
}

optimal.size = integer.neyman.w2(n.strata = length(strata), NS = n1 * sd.stra1, sample.size = 720, upper = n1)

## sample wave 1 using prior
# Put the BUGS code in a character string
BUGSmodel = "model
{for(i in 1:N) {
  logit(phistol[i]) <- alpha0 + alpha1 * instit[i] +  alpha2 * age3[i] + alpha3 * study[i] + alpha4 * stage2[i] + 
                       alpha5 * study[i] * stage2[i]
  histol[i] ~ dbern(phistol[i])
  logit(prelaps[i]) <- beta0 + beta1 * histol[i] + beta2 * age1[i] + beta3 * age2[i] + beta4 * stage1[i] + beta5 * tumdiam[i] + 
                       beta6 * stage1[i] * tumdiam[i]
  relaps[i] ~ dbern(prelaps[i])     
  }
  
  alpha0 ~ dnorm(malpha0, palpha0)
  alpha1 ~ dnorm(malpha1, palpha1)
  alpha2 ~ dnorm(malpha2, palpha2)
  alpha3 ~ dnorm(malpha3, palpha3)
  alpha4 ~ dnorm(malpha4, palpha4)
  alpha5 ~ dnorm(malpha5, palpha5)
  beta0 ~ dnorm(mbeta0, pbeta0)
  beta1 ~ dnorm(mbeta1, pbeta1)
  beta2 ~ dnorm(mbeta2, pbeta2)
  beta3 ~ dnorm(mbeta3, pbeta3)
  beta4 ~ dnorm(mbeta4, pbeta4)
  beta5 ~ dnorm(mbeta5, pbeta5)
  beta6 ~ dnorm(mbeta6, pbeta6)
  }
"



m1 = bayesian.fit(parameters = c(paste0("alpha", 0:5), paste0("beta", 0:6), "histol"),
                  data=list(histol = rep(NA, nrow(nwts)), instit = nwts$instit, relaps = nwts$relaps, stage2 = nwts$stage2, age1 = nwts$age1, 
                            age2 = nwts$age2, age3 = nwts$age3, stage1 = nwts$stage1, tumdiam = nwts$tumdiam, study = nwts$study, N = nrow(nwts), 
                            malpha0 = coef(impute)[1]-0.158, palpha0 = 10, malpha1 = coef(impute)[2]-0.158, palpha1 = 10, malpha2 = coef(impute)[3]-0.158, palpha2 = 10,
                            malpha3 = coef(impute)[4]-0.158, palpha3 = 10, malpha4 = coef(impute)[5]-0.158, palpha4 = 10, malpha5 = coef(impute)[6]-0.158, palpha5 = 10, 
                            mbeta0 = coef(full)[1]-0.158, pbeta0 = 10, mbeta1 = coef(full)[2]-0.158, pbeta1 = 10, mbeta2 = coef(full)[3]-0.158, pbeta2 = 10,
                            mbeta3 = coef(full)[4]-0.158, pbeta3 = 10, mbeta4 = coef(full)[5]-0.158, pbeta4 = 10, mbeta5 = coef(full)[6]-0.158, pbeta5 = 10,
                            mbeta6 = coef(full)[7]-0.158, pbeta6 = 10),
                  BUGSMODEL = BUGSmodel)



m2 = bayesian.fit(parameters = c(paste0("alpha", 0:5), paste0("beta", 0:6), "histol"),
                  data=list(histol = rep(NA, nrow(nwts)), instit = nwts$instit, relaps = nwts$relaps, stage2 = nwts$stage2, age1 = nwts$age1, 
                            age2 = nwts$age2, age3 = nwts$age3, stage1 = nwts$stage1, tumdiam = nwts$tumdiam, study = nwts$study, N = nrow(nwts), 
                            malpha0 = coef(impute)[1]-0.158, palpha0 = 1, malpha1 = coef(impute)[2]-0.158, palpha1 = 1, malpha2 = coef(impute)[3]-0.158, palpha2 = 1,
                            malpha3 = coef(impute)[4]-0.158, palpha3 = 1, malpha4 = coef(impute)[5]-0.158, palpha4 = 1, malpha5 = coef(impute)[6]-0.158, palpha5 = 1, 
                            mbeta0 = coef(full)[1]-0.158, pbeta0 = 1, mbeta1 = coef(full)[2]-0.158, pbeta1 = 1, mbeta2 = coef(full)[3]-0.158, pbeta2 = 1,
                            mbeta3 = coef(full)[4]-0.158, pbeta3 = 1, mbeta4 = coef(full)[5]-0.158, pbeta4 = 1, mbeta5 = coef(full)[6]-0.158, pbeta5 = 1,
                            mbeta6 = coef(full)[7]-0.158, pbeta6 = 1),
                  BUGSMODEL = BUGSmodel)


m3 = bayesian.fit(parameters = c(paste0("alpha", 0:5), paste0("beta", 0:6), "histol"),
                  data=list(histol = rep(NA, nrow(nwts)), instit = nwts$instit, relaps = nwts$relaps, stage2 = nwts$stage2, age1 = nwts$age1, 
                            age2 = nwts$age2, age3 = nwts$age3, stage1 = nwts$stage1, tumdiam = nwts$tumdiam, study = nwts$study, N = nrow(nwts), 
                            malpha0 = coef(impute)[1]-0.5, palpha0 = 10, malpha1 = coef(impute)[2]-0.5, palpha1 = 10, malpha2 = coef(impute)[3]-0.5, palpha2 = 10,
                            malpha3 = coef(impute)[4]-0.5, palpha3 = 10, malpha4 = coef(impute)[5]-0.5, palpha4 = 10, malpha5 = coef(impute)[6]-0.5, palpha5 = 10, 
                            mbeta0 = coef(full)[1]-0.5, pbeta0 = 10, mbeta1 = coef(full)[2]-0.5, pbeta1 = 10, mbeta2 = coef(full)[3]-0.5, pbeta2 = 10,
                            mbeta3 = coef(full)[4]-0.5, pbeta3 = 10, mbeta4 = coef(full)[5]-0.5, pbeta4 = 10, mbeta5 = coef(full)[6]-0.5, pbeta5 = 10,
                            mbeta6 = coef(full)[7]-0.5, pbeta6 = 10),
                  BUGSMODEL = BUGSmodel)



m4 = bayesian.fit(parameters = c(paste0("alpha", 0:5), paste0("beta", 0:6), "histol"),
                  data=list(histol = rep(NA, nrow(nwts)), instit = nwts$instit, relaps = nwts$relaps, stage2 = nwts$stage2, age1 = nwts$age1, 
                            age2 = nwts$age2, age3 = nwts$age3, stage1 = nwts$stage1, tumdiam = nwts$tumdiam, study = nwts$study, N = nrow(nwts), 
                            malpha0 = coef(impute)[1]-0.5, palpha0 = 1, malpha1 = coef(impute)[2]-0.5, palpha1 = 1, malpha2 = coef(impute)[3]-0.5, palpha2 = 1,
                            malpha3 = coef(impute)[4]-0.5, palpha3 = 1, malpha4 = coef(impute)[5]-0.5, palpha4 = 1, malpha5 = coef(impute)[6]-0.5, palpha5 = 1, 
                            mbeta0 = coef(full)[1]-0.5, pbeta0 = 1, mbeta1 = coef(full)[2]-0.5, pbeta1 = 1, mbeta2 = coef(full)[3]-0.5, pbeta2 = 1,
                            mbeta3 = coef(full)[4]-0.5, pbeta3 = 1, mbeta4 = coef(full)[5]-0.5, pbeta4 = 1, mbeta5 = coef(full)[6]-0.5, pbeta5 = 1,
                            mbeta6 = coef(full)[7]-0.5, pbeta6 = 1),
                  BUGSMODEL = BUGSmodel)

parameters = c(paste0("alpha", 0:5), paste0("beta", 0:6))
para.outcome = paste0("beta", 0:6)
#Obtain summary statistics for each of the monitored parameters
piror_info_s = summary(m1)[[1]][,1:2]
prior_pcs_s = 1 / piror_info_s[,2]^2

piror_info_w = summary(m2)[[1]][,1:2]
prior_pcs_w = 1 / piror_info_w[,2]^2

bad_info_s = summary(m3)[[1]][,1:2]
bad_pcs_s = 1 / bad_info_s[,2]^2

bad_info_w = summary(m4)[[1]][,1:2]
bad_pcs_w = 1 / bad_info_w[,2]^2


cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)

r <- foreach(j=1:5, .packages = c("rjags", "foreach", "survey")) %dopar% {
  tryCatch({
    n = 120 * j
    st.prior = two.wave.bayes(piror_info = piror_info_s, prior_pcs = prior_pcs_s, n = n, sl = 720, parameters = parameters,
                              para.outcome = para.outcome, data = nwts, n1 = n1, strata = strata)
    wg.prior.st <- st.prior$wgt$coefficients - coef(full)
    wg.prior.st.se <- sqrt(diag(st.prior$wgt$cov.unscaled))
    cal.prior.st <- st.prior$cal$coefficients - coef(full)
    cal.prior.st.se <- sqrt(diag(st.prior$cal$cov.unscaled))
    prior.st.prob <- (unlist(lapply(st.prior$p2.ob, length)) - optimal.size)
    
    ##
    flat.prior = two.wave.bayes(piror_info = piror_info_w, prior_pcs = prior_pcs_w, n = n, sl = 720, parameters = parameters,
                                para.outcome = para.outcome, data = nwts, n1 = n1, strata = strata)
    wg.prior.flat <- flat.prior$wgt$coefficients - coef(full)
    wg.prior.flat.se <- sqrt(diag(flat.prior$wgt$cov.unscaled))
    cal.prior.flat <- flat.prior$cal$coefficients - coef(full)
    cal.prior.flat.se <- sqrt(diag(flat.prior$cal$cov.unscaled))
    prior.flat.prob <- (unlist(lapply(flat.prior$p2.ob, length)) - optimal.size)
    
    ##
    bad.st.prior = two.wave.bayes(piror_info = bad_info_s, prior_pcs = bad_pcs_s, n = n, sl = 720, parameters = parameters,
                                  para.outcome = para.outcome, data = nwts, n1 = n1, strata = strata)
    wg.bad.prior.st <- bad.st.prior$wgt$coefficients - coef(full)
    wg.bad.prior.st.se <- sqrt(diag(bad.st.prior$wgt$cov.unscaled))
    cal.bad.prior.st <- bad.st.prior$cal$coefficients - coef(full)
    cal.bad.prior.st.se <- sqrt(diag(bad.st.prior$cal$cov.unscaled))
    prior.bad.st.prob <- (unlist(lapply(bad.st.prior$p2.ob, length)) - optimal.size)
    
    
    ##
    bad.flat.prior = two.wave.bayes(piror_info = bad_info_w, prior_pcs = bad_pcs_w, n = n, sl = 720, parameters = parameters,
                                    para.outcome = para.outcome, data = nwts, n1 = n1, strata = strata)
    wg.bad.prior.flat <- bad.flat.prior$wgt$coefficients - coef(full)
    wg.bad.prior.flat.se <- sqrt(diag(bad.flat.prior$wgt$cov.unscaled))
    cal.bad.prior.flat <- bad.flat.prior$cal$coefficients - coef(full)
    cal.bad.prior.flat.se <- sqrt(diag(bad.flat.prior$cal$cov.unscaled))
    prior.bad.flat.prob <- (unlist(lapply(bad.flat.prior$p2.ob, length)) - optimal.size)
    
    
    ##
#    pps_prob = pps.prob(n1 = n1, w1.n = n)
#    
#    pps <- two.wave.ms(strata = strata, prob = pps_prob, 
#                    w1.n = n, p2.n = 720, data = nwts)
#    wg.pps <- pps$wgt$coefficients - coef(full)
#    wg.pps.se <- sqrt(diag(pps$wgt$cov.unscaled))
#    cal.pps <- pps$cal$coefficients - coef(full)
#    cal.pps.se <- sqrt(diag(pps$cal$cov.unscaled))
#    pps.w2.prob <- (unlist(lapply(pps$p2.ob, length)) - optimal.size)
    
    ##
    bal <- two.wave.ms(strata = strata,  prob = rep(1/8, 8), 
                    w1.n = n, p2.n = 720, data = nwts)
    wg.bal <- bal$wgt$coefficients - coef(full)
    wg.bal.se <- sqrt(diag(bal$wgt$cov.unscaled))
    cal.bal <- bal$cal$coefficients - coef(full)
    cal.bal.se <- sqrt(diag(bal$cal$cov.unscaled))
    bal.w2.prob <- (unlist(lapply(bal$p2.ob, length)) - optimal.size)
    
    ##
    ## one wave 
    ## pps 
#    pps.one = one.wave(strata = strata, prob = pps_prob, p2.n = 720,  data = nwts)
#    wg.pps.one <- pps.one$wgt$coefficients - coef(full)
#    wg.pps.one.se <- sqrt(diag(pps.one$wgt$cov.unscaled))
#    cal.pps.one <- pps.one$cal$coefficients - coef(full)
#    cal.pps.one.se <- sqrt(diag(pps.one$cal$cov.unscaled))
#    pps.w2.prob.one <- (unlist(lapply(pps.one$p2.ob, length)) - optimal.size)
    
    #balanced
    bal.one <- one.wave(strata = strata,  prob = rep(1/8, 8), p2.n = 720, data = nwts)
    wg.bal.one <- bal.one$wgt$coefficients - coef(full)
    wg.bal.one.se <- sqrt(diag(bal.one$wgt$cov.unscaled))
    cal.bal.one <- bal.one$cal$coefficients - coef(full)
    cal.bal.one.se <- sqrt(diag(bal.one$cal$cov.unscaled))
    bal.w2.prob.one <- (unlist(lapply(bal.one$p2.ob, length)) - optimal.size)
    
    ## optimal
    opt.one <- one.wave(strata = strata,  prob = optimal.size/720, p2.n = 720, data = nwts)
    wg.opt.one <- opt.one$wgt$coefficients - coef(full)
    wg.opt.one.se <- sqrt(diag(opt.one$wgt$cov.unscaled))
    cal.opt.one <- opt.one$cal$coefficients - coef(full)
    cal.opt.one.se <- sqrt(diag(opt.one$cal$cov.unscaled))
    
    
    ## result
    list(coef.wgt = cbind(wg.prior.st, wg.prior.flat, wg.bad.prior.st, wg.bad.prior.flat, wg.bal, wg.bal.one, wg.opt.one),
         se.wgt = cbind(wg.prior.st.se, wg.prior.flat.se, wg.bad.prior.st.se, wg.bad.prior.flat.se, wg.bal.se, wg.bal.one.se, wg.opt.one.se),
         coef.cal = cbind(cal.prior.st, cal.prior.flat, cal.bad.prior.st, cal.bad.prior.flat, cal.bal, cal.bal.one, cal.opt.one),
         se.cal = cbind(cal.prior.st.se, cal.prior.flat.se, cal.bad.prior.st.se, cal.bad.prior.flat.se, cal.bal.se, cal.bal.one.se, cal.opt.one.se),
         prob = cbind(prior.st.prob, prior.flat.prob, prior.bad.st.prob, prior.bad.flat.prob, bal.w2.prob, bal.w2.prob.one))
  }, error = function(e) {
    print(e)
  })
} 
parallel::stopCluster(cl)
saveRDS(r, file = paste0("result", sample(1:10e5,1), ".rds"))
