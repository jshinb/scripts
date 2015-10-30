source('~/Desktop/SKAT_comparison/simfuncs_skat.R')
source('/Users/user/Dropbox/SBB_lab/SKAT_comparison_scripts/2x2_cont_table_LR_score_test_stat_derivations.r')

source('/Users/user/Desktop/scripts_general/sourceDir.r')
sourceDir('/Users/user/Desktop/GAW19_cleaned/scripts/pmlr11/R/')
source("/Users/user/Desktop/GAW19_cleaned/scripts/pmlr11/R/2015-04-30_pmlr_js.R")

load('/Volumes/bull_lab/shin/Desktop/pmlr_retrospective/2015-06-10_simpar_nsim_2000_GEindep_with_beta0_additional_mafs.tab')
load('/Volumes/bull_lab/shin/Desktop/pmlr_retrospective/2015-06-10_simpar_nsim_2000_GEdep_with_beta0_additional_mafs.tab')

source('/Volumes/bull_lab/shin/Desktop/pmlr_retrospective/scripts/2015-06-15_test_sim_functions_simdat_corrGX.r')

logl <- function(theta,x,y){
  y <- y
  x <- as.matrix(x)
  beta <- theta[1:ncol(x)]
  
  # Use the log-likelihood of the Bernouilli distribution, where p is
  # defined as the logistic transformation of a linear combination
  # of predictors, according to logit(p)=(x%*%beta)
  loglik <- sum(-y*log(1 + exp(-(x%*%beta))) - (1-y)*log(1 + exp(x%*%beta)))
  return(-loglik)
}

my.logl <- function(b0,b1,b2,datx,daty){
  y <- daty
  x <- as.matrix(datx)
  x <- cbind(1,x)
  beta <- c(b0,b1,b2)
  
  # Use the log-likelihood of the Bernouilli distribution, where p is
  # defined as the logistic transformation of a linear combination
  # of predictors, according to logit(p)=(x%*%beta)
  loglik <- sum(-y*log(1 + exp(-(x%*%beta))) - (1-y)*log(1 + exp(x%*%beta)))
  return(loglik)
}

my.logl2 <- function(b0,b1,b2,datx,daty){
  y <- daty
  x <- as.matrix(datx)
  x <- cbind(1,x)
  beta <- c(b0,b1,b2)
  eta = exp(x%*%beta)
  myp = eta/(1+eta)
  
  # Use the log-likelihood of the Bernouilli distribution, where p is
  # defined as the logistic transformation of a linear combination
  # of predictors, according to logit(p)=(x%*%beta)
  loglik <- sum( y*log(myp) + (1-y)*log(1-myp) )
  return(loglik)
}


# testing begins here
mybeta1s=seq(from=-20,to=20,length.out=1001)

nreps = 20
plm.score <- rao.pvals <- rep(NA,nreps)
cont.tables <- list()
for(irep in 1:nreps){
mybeta1=0
maf=0.005; b0=-3; b1=mybeta1
dat = simdat(nsim=1000,par=NULL,maf=maf,penmod="add",beta0=b0,beta1=b1,beta2=0,seed=NULL)
dat$z.x1 <- dat$x1>0 #2x2 table (dom)
cont.tables[[irep]] <- table(dat$y,dat$z.x1)
my.glm.null = glm(y~x2, data=dat, family=binomial())
my.glm = glm(y~z.x1+x2, data=dat, family=binomial())
summary(my.glm)
print(anova(my.glm.null, my.glm, test="Rao")$P[2])
rao.pvals[irep] = anova(my.glm.null, my.glm, test="Rao")$P[2]

test.log0 <- rep(NA,length(mybeta1s))
for(i in 1:length(mybeta1s)){
  test.log0[i]=my.logl2(b0=my.glm.null$coef["(Intercept)"],b1=mybeta1s[i],b2=my.glm.null$coef["x2"], 
                        datx=dat[,c(4,3)], daty=dat[,1]) #negative log-likelihood  
}

test.log1 <- rep(NA,length(mybeta1s))
for(i in 1:length(mybeta1s)){
  if(!is.na(my.glm$coef["x2"]))
    test.log1[i]=my.logl(b0=my.glm$coef["(Intercept)"],b1=mybeta1s[i],b2=my.glm$coef["x2"], datx=dat[,c(2,3)], daty=dat[,1]) #negative log-likelihood  
  if(is.na(my.glm$coef["x2"]))
    test.log1[i] <- NA
}

par(mfrow=c(1,1))
plot(mybeta1s,test.log0,type="l",ylab="log-likelihood",
     main=paste("Rao's p-value = ",round(rao.pvals[irep],4),sep=""))
abline(v=0,h=logLik(my.glm.null),col="grey",lty=2)
lines(mybeta1s,test.log1,lty=3,col="red");abline(v=my.glm$coef["z.x1TRUE"],h=logLik(my.glm),col="pink",lty=2)

dat$y.f <- factor(dat$y)
plm.S.gen = pmlr(y.f~z.x1, data=dat, method="score")
plm.LR.gen = pmlr(y.f~z.x1, data=dat, method="likelihood")
plm.LR = pmlr(y.f~z.x1+x2, data=dat, method="likelihood")
plm.LR.null = pmlr(y.f~x2, data=dat, method="likelihood")
plm.S = pmlr(y.f~z.x1+x2, data=dat, method="score")
plm.score[irep] = plm.S$test$pvalue[2,1]

temx <- model.matrix(~x1+x2, dat)
temy <- matrix(dat$y,ncol=1)

pllhd.0 <- pllhd.A <- rep(NA,length(mybeta1s))

for (i in 1:length(mybeta1s)){
  temB.0 <- matrix(c(plm.LR.null$coef[,,1]["(Intercept)"],mybeta1s[i],plm.LR.null$coef[,,1]["x2"]),ncol=1)
  temB.A <- matrix(c(plm.LR$coef[,,1]["(Intercept)"],mybeta1s[i],plm.LR$coef[,,1]["x2"]),ncol=1)
  #temB <- matrix(c(plm.S$coef[,,1][1],mybeta1s[i],plm.S.null$coef[,,1][2]),ncol=1)
  test.0 <- test.LR(x=temx, y=temy, wt=rep(1, nrow(temx)), B=temB.0, penalized=T, h0=1)
  test.A <- test.LR(x=temx, y=temy, wt=rep(1, nrow(temx)), B=temB.A, penalized=T, h0=1)
  pllhd.0[i] <- test.0$la
  pllhd.A[i] <- test.A$la
  print(pllhd.0[i])
  print(pllhd.A[i])
  cat(i,"  ")
}

llkd.res = data.frame(pllhd.0,pllhd.A,test.log0, test.log1)
png(paste("maf_",maf,"_b0_",b0,"_b1_",b1,"_likelihood",irep,".png",sep=""),width=800, height=800, pointsize=18)
par(mfrow=c(2,2)) 
plot(mybeta1s, test.log0, type="l", xlab="beta_g",ylab="log-likelihood under H0");abline(v=0, col="blue3")
plot(mybeta1s, test.log1, type="l", xlab="beta_g",ylab="log-likelihood under HA");#abline(v=c(plm.S$coef[,2,1],0), col=c("red3","blue3"))
plot(mybeta1s, pllhd.0, type="l", xlab="beta_g",ylab="log-likelihood under H0");abline(v=0, col="blue3")
plot(mybeta1s, pllhd.A, type="l", xlab="beta_g",ylab="log-likelihood under HA");abline(v=c(plm.S$coef[,2,1],0), col=c("red3","blue3"))

mtext(paste("Penalized score for G: ",round(plm.S$test$score[2],4),sep=""), side=3, outer=T, line=-1.5)
mtext(paste("Penalized score-test p-val for G: ",round(plm.score[irep],4),sep=""), side=3, outer=T, line=-2.5)
mtext(paste("ML score-test p-val for G: ",round(rao.pvals[irep],4),sep=""), side=3, outer=T, line=-3.5)

mtext("PMLE", side=2, line=4.25)
#lines(mybeta1s, pllhd.A, lty=2);abline(v=c(plm.S$coef[,2,1],0), col=c("red3","blue3"))
dev.off()
}
