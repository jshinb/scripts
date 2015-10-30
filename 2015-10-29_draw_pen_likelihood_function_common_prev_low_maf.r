setwd('/Volumes/bull_lab/shin/Desktop/GitHub_jshinb/results')

obtain.lrt.stat <- function(glm.null,glm1){
  ret = anova(glm.null,glm1,test="Chi")
  ret = ret$Deviance[2]
  ret
}

# Load devtools package to read in R files from github
# taken from - http://christophergandrud.blogspot.ca/2012/07/sourcing-code-from-github.html
require(devtools)#to source from urls
pmlr11_R_files = read.table('/Users/user/Desktop/GAW19_cleaned/scripts/pmlr11/pmlr11_R_filenames.out',sep=" ",
                            stringsAsFactors = F)[,1]
github.repos.name = 'https://raw.githubusercontent.com/jshinb/pmlr11/master/'

pmlr11_R_files = paste(github.repos.name,pmlr11_R_files,sep="")
for(i in 1:length(pmlr11_R_files)){
  source_url(pmlr11_R_files[i])
}

##
nreps = 1000
seeds = floor(runif(n=nreps,min=-1000000,max=1000000))
mybeta1s=seq(from=-20,to=20,length.out=1001)
test.log1.mat <- test.log0.mat <- test.plog1.mat <- test.plog0.mat <- matrix(NA,nrow=length(mybeta1s),ncol=nreps)
save(seeds,file='2015-10-29_seeds_common_prev.RData')
save(mybeta1s,file='2015-10-29_mybeta1s_common_prev.RData')

mybeta1=0
maf=0.005; b0=0; b1=mybeta1; nsample=2000
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
  write.table(test.log0.mat,'2015-10-29_draw_likelihood_test.log0_common_prev.mat',
              quote=F,row.names=F,col.names=F)
  write.table(test.log1.mat,'2015-10-29_draw_likelihood_test.log1_common_prev.mat',
              quote=F,row.names=F,col.names=F)
  
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
  write.table(test.plog0.mat,'2015-10-29_draw_likelihood_test.plog0_common_prev.mat',
              quote=F,row.names=F,col.names=F)
  write.table(test.plog0.mat,'2015-10-29_draw_likelihood_test.plog1_common_prev.mat',
              quote=F,row.names=F,col.names=F)

  save(rao.pvals,file='2015-10-29_rao.pvals_common_prev.RData')
  save(plm.score,file='2015-10-29_plm.score_common_prev.RData')
  save(plm.betag,file='2015-10-29_plm.betag_common_prev.RData')
  save(cont.tables,file="2015-10-29_cont.tables_common_prev.RData")
  save(glm1.list,file="2015-10-29_glm1.list_common_prev.RData")
  save(glm0.list,file="2015-10-29_glm0.list_common_prev.RData")
}

