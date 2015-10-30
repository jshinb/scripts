setwd('/Volumes/bull_lab/shin/Desktop/GitHub_jshinb/results')

load("2015-10-29_cont.tables_common_prev.RData")
load("2015-10-29_glm0.list_common_prev.RData")
load("2015-10-29_glm1.list_common_prev.RData")
load("2015-10-29_mybeta1s_common_prev.RData")
load("2015-10-29_plm.betag_common_prev.RData")
load("2015-10-29_plm.score_common_prev.RData")
load("2015-10-29_rao.pvals_common_prev.RData")
load("2015-10-29_seeds_common_prev.RData")

test.log0.mat = read.table("2015-10-29_draw_likelihood_test.log0_common_prev.mat", header=F, sep=" ", stringsAsFactors=F)
test.log1.mat = read.table("2015-10-29_draw_likelihood_test.log1_common_prev.mat", header=F, sep=" ", stringsAsFactors=F)
test.plog0.mat = read.table("2015-10-29_draw_likelihood_test.plog0_common_prev.mat", header=F, sep=" ", stringsAsFactors=F)
test.plog1.mat = read.table("2015-10-29_draw_likelihood_test.plog1_common_prev.mat", header=F, sep=" ", stringsAsFactors=F)

par(mfrow=c(1,2),ask=T)
for(rep.id in 1:length(cont.tables)){
  mytitle = paste("irep-",rep.id,"; ", "Rao Score test P-value = ",
                  round(anova(glm0.list[[rep.id]], glm1.list[[rep.id]], test="Rao")$P[2],4),sep="")
  mytitlep = paste("irep-",rep.id,"; ", "Penalized Score test P-value = ",
                   round(plm.score[rep.id],4),sep="")
  
  plot(mybeta1s,test.log0.mat[,rep.id],type="l",ylab="log-likelihood",main=mytitle)
  abline(v=0,h=logLik(glm0.list[[rep.id]]),col="grey",lty=2)
  lines(mybeta1s,test.log1.mat[,rep.id],lty=3,col="red");
  abline(v=glm1.list[[rep.id]]$coef["z.x1TRUE"],h=logLik(glm1.list[[rep.id]]),col="pink",lty=2)
  
  plot(mybeta1s,test.plog0.mat[,rep.id],type="l",ylab="log-likelihood",main=mytitlep)
  abline(v=0,h=logLik(glm0.list[[rep.id]]),col="grey",lty=2)
  lines(mybeta1s,test.plog1.mat[,rep.id],lty=3,col="red");
  abline(v=plm.betag[rep.id],h=logLik(glm1.list[[rep.id]]),col="pink",lty=2)
}
