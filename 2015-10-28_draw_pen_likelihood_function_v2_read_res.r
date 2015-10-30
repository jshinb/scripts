setwd('/Users/user/Desktop/GitHub_jshinb/results')

load("2015-10-29_cont.tables.RData")
load("2015-10-29_glm0.list.RData")
load("2015-10-29_glm1.list.RData")
load("2015-10-29_mybeta1s.RData")
load("2015-10-29_plm.betag.RData")
load("2015-10-29_plm.score.RData")
load("2015-10-29_rao.pvals.RData")
load("2015-10-29_seeds.RData")

test.log0.mat = read.table("2015-10-29_draw_likelihood_test.log0.mat", header=F, sep=" ", stringsAsFactors=F)
test.log1.mat = read.table("2015-10-29_draw_likelihood_test.log1.mat", header=F, sep=" ", stringsAsFactors=F)
test.plog0.mat = read.table("2015-10-29_draw_likelihood_test.plog0.mat", header=F, sep=" ", stringsAsFactors=F)
test.plog1.mat = read.table("2015-10-29_draw_likelihood_test.plog1.mat", header=F, sep=" ", stringsAsFactors=F)

par(mfrow=c(1,2),ask=T)
for(rep.id in 1:length(cont.tables)){
  mytitle = paste("irep-",rep.id,"; ", "Rao Score test P-value = ",
                  round(anova(glm0.list[[rep.id]], glm1.list[[rep.id]], test="Rao")$P[2],4),sep="")

  plot(mybeta1s,test.log0.mat[,rep.id],type="l",ylab="log-likelihood",main=mytitle)
  lines(mybeta1s,test.log1.mat[,rep.id],lty=3,col="red");
  abline(v=glm1.list[[rep.id]]$coef["z.x1TRUE"],h=logLik(glm1.list[[rep.id]]),col="pink",lty=2)
  plot(mybeta1s,test.plog0.mat[,rep.id],type="l",ylab="log-likelihood",main=mytitlep)
  abline(v=0,h=logLik(glm0.list[[rep.id]]),col="grey",lty=2)
  lines(mybeta1s,test.plog1.mat[,rep.id],lty=3,col="red");
  abline(v=plm.betag[rep.id],h=logLik(glm1.list[[rep.id]]),col="pink",lty=2)
}
