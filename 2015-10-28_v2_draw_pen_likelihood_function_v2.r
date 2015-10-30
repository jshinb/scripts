obtain.lrt.stat <- function(glm.null,glm1){
  ret = anova(glm.null,glm1,test="Chi")
  ret = ret$Deviance[2]
  ret
}

##
nreps = 1000
seeds = floor(runif(n=nreps,min=-1000000,max=1000000))
mybeta1s=seq(from=-20,to=20,length.out=1001)
test.log1.mat <- test.log0.mat <- test.plog1.mat <- test.plog0.mat <- matrix(NA,nrow=length(mybeta1s),ncol=nreps)

mybeta1=0
maf=0.005; b0=-3; b1=mybeta1; nsample=2000
plm.betag <- plm.score <- rao.pvals <- rep(NA,nreps)
cont.tables <- list()
glm0.list <- glm1.list <- list()

for(irep in 1:nreps){
  set.seed(seed=seeds[irep])
  dat = simdat(nsim=nsample,par=NULL,maf=maf,penmod="add",beta0=b0,beta1=b1,beta2=0,seed=NULL)
  dat$z.x1 <- dat$x1>0 #2x2 table (dom)
  cont.tables[[irep]] <- table(dat$y,dat$z.x1)
  my.glm.null = glm(y~x2, data=dat, family=binomial())
  my.glm = glm(y~z.x1+x2, data=dat, family=binomial())
  summary(my.glm)
  rao.pvals[irep] = anova(my.glm.null, my.glm, test="Rao")$P[2]
  #cat("irep-",irep,"; ", "Rao Score test P-value = ", round(anova(my.glm.null, my.glm, test="Rao")$P[2],4),"\n",sep="")
  
  dat$y.f <- factor(dat$y)
  plm.LR = pmlr(y.f~z.x1+x2, data=dat, method="likelihood")
  plm.LR.null = pmlr(y.f~x2, data=dat, method="likelihood")
  plm.S = pmlr(y.f~z.x1+x2, data=dat, method="score")
  plm.score[irep] = plm.S$test$pvalue[2,1]
  plm.betag[irep] = plm.S$coef[,,1]["z.x1TRUE"]
  temx <- model.matrix(~z.x1+x2, dat)
  temy <- matrix(dat$y,ncol=1)
  
  # calculating log-likelihoods for the given vector of betags ranging from -20 to 20
  test.log0  <- test.log1 <- rep(NA,length(mybeta1s))
  for(i in 1:length(mybeta1s)){
    # standard log-likelihood
    test.log0[i]=my.logl2(b0=my.glm.null$coef["(Intercept)"],b1=mybeta1s[i],b2=my.glm.null$coef["x2"], 
                          datx=dat[,c(4,3)], daty=dat[,1]) #negative log-likelihood  
    if(!is.na(my.glm$coef["x2"]))
      test.log1[i]=my.logl(b0=my.glm$coef["(Intercept)"],b1=mybeta1s[i],b2=my.glm$coef["x2"], datx=dat[,c(2,3)], daty=dat[,1]) #negative log-likelihood  
    if(is.na(my.glm$coef["x2"]))
      test.log1[i] <- NA
  }
  glm0.list[[irep]] <- my.glm.null
  glm1.list[[irep]] <- my.glm
  test.log0.mat[,irep] <- test.log0
  test.log1.mat[,irep] <- test.log1
  
  # penalized log-likelihood
  pllhd.0 <- pllhd.A <- rep(NA,length(mybeta1s))
  for(i in 1:length(mybeta1s)){
    temB.0 <- matrix(c(plm.LR.null$coef[,,1]["(Intercept)"],mybeta1s[i],plm.LR.null$coef[,,1]["x2"]),ncol=1)
    temB.A <- matrix(c(plm.LR$coef[,,1]["(Intercept)"],mybeta1s[i],plm.LR$coef[,,1]["x2"]),ncol=1)
    #temB <- matrix(c(plm.S$coef[,,1][1],mybeta1s[i],plm.S.null$coef[,,1][2]),ncol=1)
    test.0 <- test.LR(x=temx, y=temy, wt=rep(1, nrow(temx)), B=temB.0, penalized=T, h0=1)
    test.A <- test.LR(x=temx, y=temy, wt=rep(1, nrow(temx)), B=temB.A, penalized=T, h0=1)
    pllhd.0[i] <- test.0$la
    pllhd.A[i] <- test.A$la
    cat(irep,"th irep; ", i,"th beta1, ",sep="")
  }
  test.plog0.mat[,irep] <- pllhd.0
  test.plog1.mat[,irep] <- pllhd.A
  write.table('')
}

lrt.stats = NULL
for(i in 1:irep){
  lrt.stats = c(lrt.stats,obtain.lrt.stat(glm0.list[[i]],glm1.list[[i]]))
}

par(mfrow=c(1,2))
hist(rao.pvals)
lrt.pvals = pchisq(lrt.stats,df=1,lower.tail=F)
plot(1:irep,-log10(lrt.pvals),pch=20,ylim=c(0,7),type="b")
points(1:irep,-log10(rao.pvals),pch=1,col="red",type="b")
abline(v=c(1:irep),lty=2,col="grey")
abline(v=39,lty=1,col="grey")

plot(-log10(rao.pvals),-log10(pchisq(lrt.stats,df=1,lower.tail=F)),pch=20)
abline(0,1)

sum(pchisq(lrt.stats,df=1,lower.tail=F) < 0.05)/length(pchisq(lrt.stats,df=1,lower.tail=F))
sum(rao.pvals < 0.05)/length(rao.pvals)
sum(lrt.pvals < 0.05)/length(rao.pvals)

for(i in which(rao.pvals<0.05))
  print(cont.tables[[i]])

cont.n11.2 <- function(cont.table){
  ret = cont.table[2,2] == 2
  ret
}

cont.n11.3 <- function(cont.table){
  ret = cont.table[2,2] == 3
  ret
}

tem = sapply(cont.tables,cont.n11.2)
print(sum(tem)/length(tem))

tem3 = sapply(cont.tables,cont.n11.3)
print(sum(tem3/length(tem3)))

par(mfrow=c(1,2),ask=T)
for(rep.id in 1:irep){#which(rao.pvals < 0.05)){
  rep.id = 1
  mytitle = paste("irep-",rep.id,"; ", "Rao Score test P-value = ", 
                  round(anova(glm0.list[[rep.id]], glm1.list[[rep.id]], test="Rao")$P[2],4),sep="")
  mytitlep = paste("irep-",rep.id,"; ", "Penalized Score test P-value = ", 
                   round(plm.score[rep.id],4),sep="")
  cat(mytitle)
  print(cont.tables[[rep.id]])
  print(sum(cont.tables[[rep.id]][,2]))
  plot(mybeta1s,test.log0.mat[,rep.id],type="l",ylab="log-likelihood",main=mytitle)
  abline(v=0,h=logLik(glm0.list[[rep.id]]),col="grey",lty=2)
  lines(mybeta1s,test.log1.mat[,rep.id],lty=3,col="red");
  abline(v=glm1.list[[rep.id]]$coef["z.x1TRUE"],h=logLik(glm1.list[[rep.id]]),col="pink",lty=2)
  
  plot(mybeta1s,test.plog0.mat[,rep.id],type="l",ylab="log-likelihood",main=mytitlep)
  abline(v=0,h=logLik(glm0.list[[rep.id]]),col="grey",lty=2)
  lines(mybeta1s,test.plog1.mat[,rep.id],lty=3,col="red");
  abline(v=plm.betag[rep.id],h=logLik(glm1.list[[rep.id]]),col="pink",lty=2)
}

i=1
ns = as.vector(cont.tables[[i]])
for(i in 2:irep) ns = rbind(ns,as.vector(cont.tables[[i]]))
which(tem & rao.pvals>=0.05)
which(rao.pvals < 0.05)
cont.tables[[4]]
cont.tables[[8]]

sum(rao.pvals < 0.01)/irep
sum(rao.pvals < 0.05)/irep
