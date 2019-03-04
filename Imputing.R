#Timur Sabitov
#Imputing data
#clean work space
rm(list=ls(all=T)) 
library (zCompositions)

Impute<-function(x){
  
  #Define dimentions of the input data
  D<-dim (x)
  
  #Create matrix that will store detection limits
  mdl <- data.frame(matrix(0,ncol=D[2],nrow=D[1]))
  
  #Get column names and transfer names to the matrix
  Df_names<-names(x)
  colnames(mdl)<-Df_names[1:length(Df_names)]
  
  #Create data frame to store detection limits and values
  Dl<-data.frame(Df_names[1:46])
  colnames(Dl)<-NULL
  colnames(Dl)<-"Elements"
  
  #Define elements names and detection limits corresponding to elements
  Dl$Elements<-c("Al","As","Ba","Be","Bi","C_tot", "C_ing", "C_org",
                 "Ca","Cd","Ce","Co","Cr","Cs", "Cu","Fe","Ga",
                 "Hg","In","K","La","Li","Mg","Mn","Mo","Na","Nb",
                 "Ni","P","Pb","Rb","S","Sb","Sc","Se","Sn","Sr","Te", "Th",
                 "Ti","Tl","U","V","W","Y","Zn")
  Dl$DetLim<-c(0.01,0.06,5,0.01,0.04,0.04,0.2,0,0.01,0.1,0.05,0.1,1,5,0.5,0.01,0.08,0.01,0.02,
               0.01,0.5,1,0.01,5,0.05,0.01,0.1,0.5,50,0.5,0.2,0.01,0.05,0.1,0.2,0.1,0.5,
               0.1,0.2,0.01,0.1,0.1,1,0.1,0.2,1)
  

  #Covariance matrix with corresponding detection limits;
  for (i in seq (1:46)){
    mdl[which(colnames(mdl) == Dl$Elements[i])]<-Dl$DetLim[i]
  }
  
  #Create matrix from original values
  Input_x_mat<-x
  Input_x_mat[] <- lapply(Input_x_mat, function(x) as.numeric(as.character(x)))
  
  #Create vector that will store names of elements to drop from matrix
  drops<-array(0,46)
  
  #Finds the column with number of observations less then 60% of total observations and keeps the name of this column
  for (i in seq (1:dim(Input_x_mat)[2])){
    if ((length(which(is.na(Input_x_mat[,i]))=="True"))/nrow(Input_x_mat)>0.6)
      {drops[i]<-colnames(Input_x_mat)[i]}
  }
  
  #Get all columns except the one that was dropped
  Input_x_mat<-Input_x_mat[ , !(colnames(Input_x_mat) %in% drops)]
  mdl<-mdl[ , !(colnames(mdl) %in% drops)]
  
  #Replace all NA with 0
  Input_x_mat[is.na(Input_x_mat)]<-0
  Input_x_mat[]<-lapply(Input_x_mat,function(x) as.numeric(x))
  Input_x_mat<-as.matrix(Input_x_mat)
  mdl<-as.matrix(mdl)
  
  
  datF <- cbind (colnames(mdl), sample(Dl$DetLim, length(colnames(mdl))))
  datF[,2] <- as.numeric(datF[,2])
  datF <- as.data.frame (datF)
  #match (colnames(mdl),datF[,1])
  datF[,2] <- as.numeric(datF[,2])
  #Input_x_mat[,match (colnames(Input_x_mat),datF[,1])] <- datF[,2][match (colnames(Input_x_mat),datF[,1])]
  
  
  
  ##CHANGED 3/3/2019
  ##
  ##
  ##
  ##
  # for (i in  seq (1, ncol(mdl))){
  #   ind <- which (colnames(mdl)[i] == datF[,1][i])
  #   mdl[,i] <- datF[,2][i]
  # }
  # mdl <- as.matrix(mdl)
  
  #Impute data using z-composition package
  Input_x_lrEM <- lrEM(Input_x_mat,label=0,dl=mdl,ini.cov="multRepl")
  return (Input_x_lrEM)
}

#a function to find and give NA values from the data frame (stackexchange)
which.nonnum <- function(x) {
  badNum <- is.na(suppressWarnings(as.numeric(as.character(x))))
  which(badNum & !is.na(x))
}


#Data frame input and rect.hclust type
DT<-function(input,class,k){
  Clust_tree<-class
  Df<-input
  grp <- cutree(Clust_tree, k = k)
  z<-rownames(Df)
  for (i in seq (1:k)){
    x<-rownames(Df)[grp==i]
    for (j in seq (1:length(x)))
      Df$Class[which(x[j]==z)]<-i
  }
  return (Df)
}


logitTransform <- function(p) { log(p/(1-p)) }#Used for proportion data
squared <- function(p) { p^2 }
cubroot<- function(p) { (p^1/3) }
log_neg<-function(p) { log(p+1)}
one.orig<-function(p) {1/p}
one.sqrt<-function (p) {1/sqrt(p)}
arc<-function(p) { asin(sqrt(p))}#Used for percentage data



TR<-function (x,c1,c2){
  Input<-x
  Temp1<-c("squared","cubroot","log_neg","log","one.orig","sqrt","one.sqrt","arc")
  Temp.Df<-as.data.frame(matrix(nrow = dim(Input)[1],ncol = 8))
  Transf.Df<-as.data.frame(matrix(nrow= dim(Input)[1], ncol = 4))
  colnames(Temp.Df)<-Temp1
  colnames(Transf.Df)<-c("Element","p-val","Method","Result")
  k<-1
  
  for (i in seq (c1,c2)){
    c.loc<-0
    p<-0
    p_value<-array(0,8)
    p<-shapiro.test(Input[,i])
    if (p$p.value<0.05){
      Temp.Df$squared<-squared(Input[,i])
      Temp.Df$cubroot<-cubroot(Input[,i])
      Temp.Df$log_neg<-log_neg(Input[,i])
      Temp.Df$log<-log(Input[,i])
      Temp.Df$one.orig<-one.orig(Input[,i])
      Temp.Df$sqrt<-sqrt(Input[,i])
      Temp.Df$one.sqrt<-1/(sqrt(Input[,i]))
      Temp.Df$arc<-arc(Input[,i])
      Temp.Df[Temp.Df==Inf]<-NaN
      Temp.Df[Temp.Df==-Inf]<-NaN
      
      for (j in seq (1,8)){
        if (length (which(complete.cases(Temp.Df[,j])=="FALSE"))>0){p$p.value<-0}
        else {
          p<-shapiro.test(Temp.Df[,j])
          p_value[j]<-p$p.value
        }
      }
      c.loc<-which(p_value==(max(p_value)))
      Input[,i]<-Temp.Df[,c.loc]
      Transf.Df$Element[k]<-colnames(Input[i])
      Transf.Df$'p-val'[k]<-max(p_value)
      Transf.Df$Method[k]<-colnames(Temp.Df[c.loc])
      
      if (max(p_value)>0.05) {
        Transf.Df$Result[k]<-"Normalized"
        k<-k+1
      }
      else{
        Transf.Df$Result[k]<-"Not normalized"
        k<-k+1
      }
    }
    else {
      Input[,i]<-Input[,i]
      Transf.Df$Element[k]<-colnames(Input[i])
      Transf.Df$'p-val'[k]<-p$p.value
      Transf.Df$Method[k]<-"No transform"
      Transf.Df$Result[k]<-"Normalized"
      k<-k+1
    }
  }
  result<-list()
  result[[1]]<-Transf.Df
  result[[2]]<-Input
  return (result)
}

TR_percent<-function (x,c1,c2,Transf.Df,k){
  for (i in seq (c1,c2)){
    Temp.Df<-as.data.frame(matrix(nrow = dim(x)[1], ncol = 2))
    V1<-arc(x[i])
    Temp.Df$V1<-V1
    #Temp.Df$V2<-V2
    #colnames(Temp.Df)<-NULL
    #names(Temp.Df)<-c("Arc","Logit")
    names(Temp.Df)<-c("Arc","Logit")
    p_value<-array(0,2)
    x1.res<-shapiro.test(Temp.Df$Arc[,1])
    #x2.res<-shapiro.test(Temp.Df$Logit)
    p_value[1]<-x1.res$p.value
    #p_value[2]<-x2.res$p.value
    c.loc<-which(p_value==(max(p_value)))
    x[,i]<-Temp.Df[,c.loc]
    Transf.Df$Element[k]<-"Arc"
    Transf.Df$'p-val'[k]<-max(p_value)
    Transf.Df$Method[k]<-colnames(Temp.Df[c.loc])
    if (max(p_value)>0.05) {
      Transf.Df$Result[k]<-"Normalized"
      k<-k+1
    }
    else{
      Transf.Df$Result[k]<-"Not normalized"
      k<-k+1
    }
  }
  result<-list()
  result[[1]]<-Transf.Df
  result[[2]]<-Input
  return (result)
}

# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


Best.clust.method<-function(dis){
  di<-dis
  Clust_methods<- c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
  r<-0
  for (i in seq (1:length(Clust_methods))){
    tree<-hclust(di, method = Clust_methods[i])
    Clust_cor<-cor(di, cophenetic(tree))
    if (Clust_cor>r){
      r<-Clust_cor
      Clust_tree<-tree
    }
    return (Clust_tree)
  }
}

# @aL3xa's function
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}


