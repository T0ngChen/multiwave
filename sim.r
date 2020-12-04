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
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) / nrow(dm)
  ## influence function
  infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
  infl
}

inf.prior = function(dm, para.val){
  fit.val = exp(rowSums(dm %*% diag(para.val))) / (1+exp(rowSums(dm %*% diag(para.val))))
  Ihat <- (t(dm) %*% (dm * fit.val * (1 - fit.val))) / nrow(dm) 
  infl = as.matrix(dm * (data1$y - fit.val)) %*% solve(Ihat)
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
  wgt = svyglm(y ~ x + z1 + z2, family = quasibinomial, design = twophase.w2)
  impmodel <- svyglm(x ~ a + z1 + z2, family = quasibinomial, design = twophase.w2)
  data$impx <- as.vector(predict(impmodel,newdata=data,type="response",se.fit=FALSE))
  phase1model_imp <- glm(y ~ impx + z1 + z2, family=binomial, data=data)
  inffun_imp <- inf.fun(phase1model_imp)
  colnames(inffun_imp)<-paste0("if",1:ncol(inffun_imp))
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, ~stra), 
                              subset = ~insample, data = cbind(data, inffun_imp), method="simple")
  calformula <- make.formula(colnames(inffun_imp))
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking")
  svyest_imp<-svyglm(y ~ x + z1 + z2, family=quasibinomial, design=cal_twophase_imp)
  list(wgt = wgt, cal = svyest_imp, p2.ob = s.w1)
}


pps.prob = function(n1 = n1, w1.n = n){
  if(any(n1/sum(n1) * w1.n < 3)){
    prob = n1/sum(n1)
    prob[prob<0.05] = 0.05
    prob[which.max(prob)] = prob[which.max(prob)] - (sum(prob) - 1)
  } else {
    prob = n1/sum(n1)
  }
  prob
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
  fit.w1 <- glm(y ~ x + z1 + z2, family = binomial, data = data[data$in.subsample,], weights = 1/twophase.w1$prob)
  
  n2 = xtabs(~stra, data[data$in.subsample,])
  
  X <- model.matrix(fit.w1)
  pi.hat <- fit.w1$fitted.values
  Ihat <- (t(as.matrix(X)) %*% (as.matrix(X) * 1/twophase.w1$prob * pi.hat * (1 - pi.hat)))/sum(n1)
  wgt<-n1/n2*(n1-n2)
  S <- X * (data[data$in.subsample,]$y - pi.hat)
  invI <- solve(Ihat)
  result <- 0
  varsi<-ar<-array(0,dim=c(nrow(Ihat),nrow(Ihat),length(wgt)))
  for(jj in 1:length(wgt)) {
    si <- S[data[data$in.subsample,]$stra == jj,  ]
    varsi[,,jj]<-var(si)
    ar[,,jj]<-invI%*%varsi[,,jj]%*%invI
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
  for (index in 1:length(strata)){
    n3[index] = length(strata[[index]][-match(s.w1[[index]], strata[[index]])])
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
  wgt = svyglm(y ~ x + z1 + z2, family = quasibinomial, design = twophase.w2)
  impmodel <- svyglm(x ~ a + z1 + z2, family = quasibinomial, design = twophase.w2)
  data$impx <- as.vector(predict(impmodel,newdata=data,type="response",se.fit=FALSE))
  phase1model_imp <- glm(y ~ impx + z1 + z2, family=binomial, data=data)
  inffun_imp <- inf.fun(phase1model_imp)
  colnames(inffun_imp)<-paste0("if",1:ncol(inffun_imp))
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, ~stra), 
                              subset = ~insample, data = cbind(data, inffun_imp), method="simple")
  calformula <- make.formula(colnames(inffun_imp))
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking")
  svyest_imp<-svyglm(y ~ x + z1 + z2, family=quasibinomial, design=cal_twophase_imp)
  list(wgt = wgt, cal = svyest_imp, p2.ob = p2.ob)
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
  w1.sd = wave.sample.size(dm = cbind(rep(1,1000), rbinom(1000, 1, piror_info[,1][!rownames(piror_info) %in% parameters]), data$z1, data$z2), 
                           para.val = piror_info[,1][para.outcome])
  
  w1.size = integer.neyman.w1(n.strata = length(strata), NS = w1.sd * n1, sample.size = n, upper = n1)
  
  s.w1 = list()
  s.w2 = list()
  for (i in 1:length(strata)) {
    s.w1[[i]] <- sample(strata[[i]], w1.size[i])
  }
  
  fit2 = bayesian.fit(parameters = c("alpha0", "alpha1", "alpha2", "alpha3", "beta0", "beta1", "beta2", "beta3", "x"),
                      data=list(x = data$x[ifelse(1:nrow(data) %in% unlist(s.w1), T, NA)], a = data$a, z1 = data$z1, z2 = data$z2,
                                y = data$y, N = nrow(data), 
                                malpha0 = piror_info[,1]["alpha0"], palpha0 = prior_pcs["alpha0"], malpha1 = piror_info[,1]["alpha1"], palpha1 = prior_pcs["alpha1"], 
                                malpha2 = piror_info[,1]["alpha2"], palpha2 = prior_pcs["alpha2"], malpha3 = piror_info[,1]["alpha3"], palpha3 = prior_pcs["alpha3"],
                                mbeta0 = piror_info[,1]["beta0"], pbeta0 = prior_pcs["beta0"], mbeta1 = piror_info[,1]["beta1"], pbeta1 = prior_pcs["beta1"], 
                                mbeta2 = piror_info[,1]["beta2"], pbeta2 = prior_pcs["beta2"], mbeta3 = piror_info[,1]["beta3"], pbeta3 = prior_pcs["beta3"]),
                      BUGSMODEL = BUGSmodel)
  
  w1.post.info = summary(fit2)[[1]][,1:2]
  
  ix = rbinom(1000, 1, w1.post.info[,1][!rownames(w1.post.info) %in% parameters])
  
  w2.sd = wave.sample.size(dm = cbind(rep(1,1000), ix, data$z1, data$z2), para.val = w1.post.info[,1][para.outcome])
  
  upper = n1 - sapply(s.w1,length)
  
  w2.size = integer.neyman.w2(n.strata = length(strata), NS = w2.sd * n1, sample.size = sl - n, upper = upper)
  
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
  wgt = svyglm(y ~ x + z1 + z2, family = quasibinomial, design = twophase.w2)
  impmodel <- svyglm(x ~ a + z1 + z2, family = quasibinomial, design = twophase.w2)
  data$impx <- as.vector(predict(impmodel,newdata=data,type="response",se.fit=FALSE))
  phase1model_imp <- glm(y ~ impx + z1 + z2, family=binomial, data=data)
  inffun_imp <- inf.fun(phase1model_imp)
  colnames(inffun_imp)<-paste0("if",1:ncol(inffun_imp))
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, ~stra), 
                              subset = ~in.subsample, data = cbind(data, inffun_imp), method="simple")
  calformula <- make.formula(colnames(inffun_imp))
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking")
  svyest_imp<-svyglm(y ~ x + z1 + z2, family=quasibinomial, design=cal_twophase_imp)
  list(wgt = wgt, cal = svyest_imp, p2.ob = p2.ob)
}

################
##            ##
## simulation ##
##            ##
################

## generate data
data1 <- data.frame(x = rbinom(1000, size = 1, 0.15))
# sensitivity 0.8 specificity 0.8
data1$a = ifelse(data1$x == 1, rbinom(1000,1,0.8), rbinom(1000,1,0.2))
# sensitivity 0.9 specificity 0.9
#data1$a = ifelse(data1$x == 1, rbinom(1000,1,0.9), rbinom(1000,1,0.1))
# sensitivity 0.9 specificity 0.8
#data1$a = ifelse(data1$x == 1, rbinom(1000,1,0.9), rbinom(1000,1,0.2))
# sensitivity 0.8 specificity 0.9
#data1$a = ifelse(data1$x == 1, rbinom(1000,1,0.8), rbinom(1000,1,0.1))

data1$z1 = runif(1000)
data1$z2 = rbinom(1000, size = 1, 0.6)

p_y <- exp(-2 + 1.5 * data1$x + data1$z1 + data1$z2)/(1 + exp(-2 + 1.5 * data1$x + data1$z1 + data1$z2))
data1$y <- rbinom(1000, size = 1, p_y)

tru.val <- c(-2, 1.5, 1, 1)

## simulation
## define strata
## define strata 
data1$stra = 1 + 4 * data1$y + 2 * data1$a + data1$z2
## population strata size
n1 =  xtabs(~stra, data1)

strata = list()
for (index in 1:length(n1)){
  strata[[index]] = which(data1$stra == index)  
}

full = glm(y ~ x + z1 + z2, data = data1, family = binomial)
impute = glm(x ~ a + z1 + z2, data = data1, family = binomial)

infl <- inf.fun(full)[, 2]

sd.stra1 = numeric()

for(ii in 1:length(strata)){
  sd.stra1 = c(sd.stra1, sd(infl[strata[[ii]]]))
}


optimal.size = integer.neyman.w2(n.strata = length(strata), NS = n1 * sd.stra1, sample.size = 300, upper = n1)


## sample wave 1 using prior
# Put the BUGS code in a character string
BUGSmodel = "model
{for(i in 1:N) {
  logit(px[i]) <- alpha0 + alpha1 * a[i] +  alpha2 * z1[i] +  alpha3 * z2[i]
  x[i] ~ dbern(px[i])
  logit(py[i]) <- beta0 + beta1 * x[i] + beta2 * z1[i] + beta3 * z2[i]
  y[i] ~ dbern(py[i])     
  }
  
  alpha0 ~ dnorm(malpha0, palpha0)
  alpha1 ~ dnorm(malpha1, palpha1)
  alpha2 ~ dnorm(malpha2, palpha2)
  alpha3 ~ dnorm(malpha3, palpha3)
  beta0 ~ dnorm(mbeta0, pbeta0)
  beta1 ~ dnorm(mbeta1, pbeta1)
  beta2 ~ dnorm(mbeta2, pbeta2)
  beta3 ~ dnorm(mbeta3, pbeta3)
  }
"
m1 = bayesian.fit(parameters = c("alpha0", "alpha1", "alpha2", "alpha3", "beta0", "beta1", "beta2", "beta3", "x"),
                  data=list(x = rep(NA, nrow(data1)), a = data1$a, z1 = data1$z1, z2 = data1$z2,
                            y = data1$y, N = nrow(data1), malpha0 = coef(impute)[1]-0.1581139, palpha0 = 10,
                            malpha1 = coef(impute)[2]-0.1581139, palpha1 = 10, malpha2 = coef(impute)[3]-0.1581139, palpha2 = 10,
                            malpha3 = coef(impute)[4]-0.1581139, palpha3 = 10,
                            mbeta0 = -2-0.1581139, pbeta0 = 10, mbeta1 = 1.5-0.1581139, pbeta1 = 10, mbeta2 = 1-0.1581139, pbeta2 = 10,
                            mbeta3 = 1-0.1581139, pbeta3 = 10),
                  BUGSMODEL = BUGSmodel)



m2 = bayesian.fit(parameters = c("alpha0", "alpha1", "alpha2", "alpha3", "beta0", "beta1", "beta2", "beta3", "x"),
                  data=list(x = rep(NA, nrow(data1)), a = data1$a, z1 = data1$z1, z2 = data1$z2,
                            y = data1$y, N = nrow(data1), malpha0 = coef(impute)[1]-0.1581139, palpha0 = 1,
                            malpha1 = coef(impute)[2]-0.1581139, palpha1 = 1, malpha2 = coef(impute)[3]-0.1581139, palpha2 = 1,
                            malpha3 = coef(impute)[4]-0.1581139, palpha3 = 1,
                            mbeta0 = -2-0.1581139, pbeta0 = 1, mbeta1 = 1.5-0.1581139, pbeta1 = 1, mbeta2 = 1-0.1581139, pbeta2 = 1,
                            mbeta3 = 1-0.1581139, pbeta3 = 1),
                  BUGSMODEL = BUGSmodel)

m3 = bayesian.fit(parameters = c("alpha0", "alpha1", "alpha2", "alpha3", "beta0", "beta1", "beta2", "beta3", "x"),
                  data=list(x = rep(NA, nrow(data1)), a = data1$a, z1 = data1$z1, z2 = data1$z2,
                            y = data1$y, N = nrow(data1), malpha0 = coef(impute)[1]-0.5, palpha0 = 10,
                            malpha1 = coef(impute)[2]-0.5, palpha1 = 10, malpha2 = coef(impute)[3]-0.5, palpha2 = 10,
                            malpha3 = coef(impute)[4]-0.5, palpha3 = 10,
                            mbeta0 = -2-0.5, pbeta0 = 10, mbeta1 = 1.5-0.5, pbeta1 = 10, mbeta2 = 1-0.5, pbeta2 = 10,
                            mbeta3 = 1-0.5, pbeta3 = 10),
                  BUGSMODEL = BUGSmodel)



m4 = bayesian.fit(parameters = c("alpha0", "alpha1", "alpha2", "alpha3", "beta0", "beta1", "beta2", "beta3", "x"),
                  data=list(x = rep(NA, nrow(data1)), a = data1$a, z1 = data1$z1, z2 = data1$z2,
                            y = data1$y, N = nrow(data1), malpha0 = coef(impute)[1]-0.5, palpha0 = 1,
                            malpha1 = coef(impute)[2]-0.5, palpha1 = 1, malpha2 = coef(impute)[3]-0.5, palpha2 = 1,
                            malpha3 = coef(impute)[4]-0.5, palpha3 = 1,
                            mbeta0 = -2-0.5, pbeta0 = 1, mbeta1 = 1.5-0.5, pbeta1 = 1, mbeta2 = 1-0.5, pbeta2 = 1,
                            mbeta3 = 1-0.5, pbeta3 = 1),
                  BUGSMODEL = BUGSmodel)

parameters = c("alpha0", "alpha1", "alpha2", "alpha3", "beta0", "beta1", "beta2", "beta3")
para.outcome = c("beta0", "beta1", "beta2", "beta3")
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
    n = 50 * j
    st.prior = two.wave.bayes(piror_info = piror_info_s, prior_pcs = prior_pcs_s, n = n, sl = 300, parameters = parameters,
                              para.outcome = para.outcome, data = data1, n1 = n1, strata = strata)
    wg.prior.st <- st.prior$wgt$coefficients - tru.val
    wg.prior.st.se <- sqrt(diag(st.prior$wgt$cov.unscaled))
    cal.prior.st <- st.prior$cal$coefficients - tru.val
    cal.prior.st.se <- sqrt(diag(st.prior$cal$cov.unscaled))
    prior.st.prob <- (unlist(lapply(st.prior$p2.ob, length)) - optimal.size)
    
    ##
    flat.prior = two.wave.bayes(piror_info = piror_info_w, prior_pcs = prior_pcs_w, n = n, sl = 300, parameters = parameters,
                                para.outcome = para.outcome, data = data1, n1 = n1, strata = strata)
    wg.prior.flat <- flat.prior$wgt$coefficients - tru.val
    wg.prior.flat.se <- sqrt(diag(flat.prior$wgt$cov.unscaled))
    cal.prior.flat <- flat.prior$cal$coefficients - tru.val
    cal.prior.flat.se <- sqrt(diag(flat.prior$cal$cov.unscaled))
    prior.flat.prob <- (unlist(lapply(flat.prior$p2.ob, length)) - optimal.size)
    
    ##
    bad.st.prior = two.wave.bayes(piror_info = bad_info_s, prior_pcs = bad_pcs_s, n = n, sl = 300, parameters = parameters,
                                  para.outcome = para.outcome, data = data1, n1 = n1, strata = strata)
    wg.bad.prior.st <- bad.st.prior$wgt$coefficients - tru.val
    wg.bad.prior.st.se <- sqrt(diag(bad.st.prior$wgt$cov.unscaled))
    cal.bad.prior.st <- bad.st.prior$cal$coefficients - tru.val
    cal.bad.prior.st.se <- sqrt(diag(bad.st.prior$cal$cov.unscaled))
    prior.bad.st.prob <- (unlist(lapply(bad.st.prior$p2.ob, length)) - optimal.size)
    
    
    ##
    bad.flat.prior = two.wave.bayes(piror_info = bad_info_w, prior_pcs = bad_pcs_w, n = n, sl = 300, parameters = parameters,
                                    para.outcome = para.outcome, data = data1, n1 = n1, strata = strata)
    wg.bad.prior.flat <- bad.flat.prior$wgt$coefficients - tru.val
    wg.bad.prior.flat.se <- sqrt(diag(bad.flat.prior$wgt$cov.unscaled))
    cal.bad.prior.flat <- bad.flat.prior$cal$coefficients - tru.val
    cal.bad.prior.flat.se <- sqrt(diag(bad.flat.prior$cal$cov.unscaled))
    prior.bad.flat.prob <- (unlist(lapply(bad.flat.prior$p2.ob, length)) - optimal.size)
    
    
    ##
    pps_prob = pps.prob(n1 = n1, w1.n = n)

    pps <- two.wave.ms(strata = strata, prob = pps_prob, 
                    w1.n = n, p2.n = 300, data = data1)
    wg.pps <- pps$wgt$coefficients - tru.val
    wg.pps.se <- sqrt(diag(pps$wgt$cov.unscaled))
    cal.pps <- pps$cal$coefficients - tru.val
    cal.pps.se <- sqrt(diag(pps$cal$cov.unscaled))
    pps.w2.prob <- (unlist(lapply(pps$p2.ob, length)) - optimal.size)
    
    ##
    bal <- two.wave.ms(strata = strata,  prob = rep(1/8, 8), 
                    w1.n = n, p2.n = 300, data = data1)
    wg.bal <- bal$wgt$coefficients - tru.val
    wg.bal.se <- sqrt(diag(bal$wgt$cov.unscaled))
    cal.bal <- bal$cal$coefficients - tru.val
    cal.bal.se <- sqrt(diag(bal$cal$cov.unscaled))
    bal.w2.prob <- (unlist(lapply(bal$p2.ob, length)) - optimal.size)
    
    ## one wave 
    ## pps 
    pps.one = one.wave(strata = strata, prob = pps_prob, p2.n = 300,  data = data1)
    wg.pps.one <- pps.one$wgt$coefficients - tru.val
    wg.pps.one.se <- sqrt(diag(pps.one$wgt$cov.unscaled))
    cal.pps.one <- pps.one$cal$coefficients - tru.val
    cal.pps.one.se <- sqrt(diag(pps.one$cal$cov.unscaled))
    pps.w2.prob.one <- (unlist(lapply(pps.one$p2.ob, length)) - optimal.size)
    
    #balanced
    bal.one <- one.wave(strata = strata,  prob = rep(1/8, 8), p2.n = 300, data = data1)
    wg.bal.one <- bal.one$wgt$coefficients - tru.val
    wg.bal.one.se <- sqrt(diag(bal.one$wgt$cov.unscaled))
    cal.bal.one <- bal.one$cal$coefficients - tru.val
    cal.bal.one.se <- sqrt(diag(bal.one$cal$cov.unscaled))
    bal.w2.prob.one <- (unlist(lapply(bal.one$p2.ob, length)) - optimal.size)
    
    
    
    #opt
    opt.one <- one.wave(strata = strata,  prob = optimal.size/300, p2.n = 300, data = data1)
    wg.opt.one <- opt.one$wgt$coefficients - tru.val
    wg.opt.one.se <- sqrt(diag(opt.one$wgt$cov.unscaled))
    cal.opt.one <- opt.one$cal$coefficients - tru.val
    cal.opt.one.se <- sqrt(diag(opt.one$cal$cov.unscaled))
    opt.w2.prob.one <-  optimal.size
    
    ## result
    list(coef.wgt = cbind(wg.prior.st, wg.prior.flat, wg.bad.prior.st, wg.bad.prior.flat, wg.pps, wg.bal, wg.pps.one, wg.bal.one, wg.opt.one),
         se.wgt = cbind(wg.prior.st.se, wg.prior.flat.se, wg.bad.prior.st.se, wg.bad.prior.flat.se, wg.pps.se, wg.bal.se, wg.pps.one.se, wg.bal.one.se, wg.opt.one.se),
         coef.cal = cbind(cal.prior.st, cal.prior.flat, cal.bad.prior.st, cal.bad.prior.flat, cal.pps, cal.bal, cal.pps.one, cal.bal.one, cal.opt.one),
         se.cal = cbind(cal.prior.st.se, cal.prior.flat.se, cal.bad.prior.st.se, cal.bad.prior.flat.se, cal.pps.se, cal.bal.se, cal.pps.one.se, cal.bal.one.se, cal.opt.one.se),
         prob = cbind(prior.st.prob, prior.flat.prob, prior.bad.st.prob, prior.bad.flat.prob, pps.w2.prob, bal.w2.prob, pps.w2.prob.one, bal.w2.prob.one, opt.w2.prob.one))
  }, error = function(e) {
    matrix(NA, 9, 9)
  }) 
}
parallel::stopCluster(cl)
saveRDS(r, file = paste0("result", sample(1:10e5,1), ".rds"))
