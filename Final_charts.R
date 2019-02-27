#write.csv (dat, 'dat.csv')
dat <- read.csv ('dat.csv')
rownames (dat) <- dat[,1]
dat <- dat[,-1]
dat <- dat[grep('target', rownames(dat)),]
dat <- melt (dat)
tryCatch({
  dat <- dat[which (dat$Target == input$targetPlot),]
  var <<- dat
  ggplot(dat, aes(factor(variable), value, colour = variable)) +
    geom_violin(trim = FALSE) + geom_jitter(height = 0, width = 0.1) + geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
    facet_wrap(variable, ncol = 2, scales = "free")
}, warning = function(cond) {}, error = function(cond) {})











library (gmodels)
library (gtools)
targets <- dat[grep ('targ', rownames(dat)),]
n <- as.character (targets$Target)
levels (targets$Target) <- c(as.character(unique(targets$Target)))
targets$Target <- n
targetsV <- aggregate (targets[,(1:(ncol(targets)-3))],  by = list (targets$Target), mean)
tdat<- mixedsort(as.character(targetsV[,1]))
targetsV <- targetsV[match (tdat, targetsV[,1]),]

options(warn=-1)
targetsCi <- aggregate (targets[,(1:(ncol(targets)-3))],  by = list (targets$Target), ci)
options(warn=0)

targetsCi <- targetsCi[match (tdat, targetsCi[,1]),]
rownames(targetsCi) <- seq (1,nrow(targetsCi))
dat <- t (targetsCi[,-1])
colnames (dat) <- targetsCi$Group.1

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

mybarPlot <- dat[grep ('Estimate', rownames(dat)),]
mybarError <- dat[grep ('Error', rownames (dat)),]
ze_barplot = barplot(mybarPlot , beside=T, ylim=c(0,1))
error.bar(ze_barplot,mybarPlot, mybarError)
#First I create a smell function that takes...in entry










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


A=c(rep("sorgho" , 10) , rep("poacee" , 10) )
B=rnorm(20,10,4)
C=rnorm(20,8,3)
D=rnorm(20,5,4)
data=data.frame(A,B,C,D)
colnames(data)=c("specie","cond_A","cond_B","cond_C")


#Let's calculate the average value for each condition and each specie with the *aggregate* function
bilan=aggregate(cbind(cond_A,cond_B,cond_C)~specie , data=data , mean)
rownames(bilan)=bilan[,1]
bilan=as.matrix(bilan[,-1])



#Then I calculate the standard deviation for each specie and condition :
stdev=aggregate(cbind(cond_A,cond_B,cond_C)~specie , data=data , sd)
rownames(stdev)=stdev[,1]
stdev=as.matrix(stdev[,-1]) * 1.96 / 10



#Then it is easy to make a classical barplot :
lim=1.2*max(bilan)
ze_barplot = barplot(bilan , beside=T , legend.text=T , col=c("blue" , "skyblue") , ylim=c(0,lim))

#I becomes a bit more tricky when we want to add the error bar representing the confidence interval.

#First I create a smell function that takes...in entry
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#Then I calculate the standard deviation for each specie and condition :
stdev=aggregate(cbind(cond_A,cond_B,cond_C)~specie , data=data , sd)
rownames(stdev)=stdev[,1]
stdev=as.matrix(stdev[,-1]) * 1.96 / 10

# I am ready to add the error bar on the plot using my "error bar" function !
#   ze_barplot = barplot(bilan , beside=T , legend.text=T,col=c("blue" , "skyblue") , ylim=c(0,lim) , ylab="height")
# error.bar(ze_barplot,bilan, stdev)
