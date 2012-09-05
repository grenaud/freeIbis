data <- read.table("SVMlight_models.index",header=F,as.is=T,sep='\t')
rates <- data[,c(seq(6,51,3),52)]*100
nucleotides <- c(1:4)
nletters <- c("A","C","G","T")

accuracy <- rates[,c(nucleotides+(nucleotides-1)*4,17)]
accuracy[,1:4] <- 100-accuracy[,1:4]
rates <- rates[,c(nucleotides+(nucleotides-1)*4,17)*-1]

pdf("error_profile.pdf",width=18,height=12)
matplot(accuracy,type="l",col=c("green","blue","darkgrey","red","black"),lty=c(rep(1,4),2),lwd=c(rep(1,4),2),main="Control lane prediction accuracy",xlab="Cycles",ylab="Error rate [%]")
legend("topleft",c(nletters,"Average"),fill=c("green","blue","darkgrey","red","black"))

header <- c(unlist(sapply(nletters,FUN=function(first){list(sapply(nletters,FUN=function(second) {if (first != second) {paste(first,"->",second,sep="")}}))})))
matplot(rates,type="l",lty=rep(1,12),col=rainbow(12),main="Control lane missclassification rates by substitution",xlab="Cycles",ylab="Rate per base [%]")
legend("topleft",header,fill=rainbow(12))
dev.off()
