#write.csv (dat, 'dat.csv')
dat <- read.csv ('dat.csv')
rownames (dat) <- dat[,1]
dat <- dat[,-1]
library (gmodels)
library (gtools)
targets <- dat[grep ('targ', rownames(dat)),]
n <- as.character (targets$Target)
levels (targets$Target) <- c(as.character(unique(targets$Target)))
targets$Target <- n
targetsV <- aggregate (targets[,(1:(ncol(targets)-3))],  by = list (targets$Target), mean)
tdat<- mixedsort(as.character(targetsV[,1]))
targetsV <- targetsV[match (tdat, targetsV[,1]),]
getCi <- function (x) {
  x <- as.data.frame (x)
  return (apply (x, 2, ci))
}
targetsCi <- aggregate (targets[,(1:(ncol(targets)-3))],  by = list (targets$Target), getCi)

targetsCi <- targetsCi[match (tdat, targetsCi[,1]),]
dat <- t (targetsCi[,-1])
colnames (dat) <- targetsCi$Group.1

rownames(dat) <- gsub ('1','es', rownames (dat))
rownames(dat)<- gsub ('2','ci.low', rownames (dat))
rownames(dat)<- gsub ('3','ci.up', rownames (dat))
rownames(dat)<- gsub ('4','sd.er', rownames (dat))


cl1 <- dat[grep ('BANK', rownames(dat)),]
barplot (cl1, col=c("chartreuse", "blue4", 'lightblue', 'lightsalmon'))

#Estimate   CI lower   CI upper Std. Error 



t(getCi (targets[which(targets$Target == 'Target: 1'),(1:(ncol(targets)-3))]))


targetsCi <- targetsCi[match (tdat, targetsCi[,1]),]
rownames (targetsCi) <- NULL

#p <- barplot (targetsPlot[,1]) #only first target
targetsPlot <- t(targetsV[,-1])
colnames(targetsPlot)<- (targetsV[,1])
barplot (targetsPlot)

unique (targets$Target)

