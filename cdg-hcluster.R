genes <- read.csv("~/qiulab/cdg/data/cdgSNPmatrix-reduced.csv", header = F, row.names = 1)
cdg <- read.csv("~/qiulab/cdg/data/cdgTable.csv2", header = T, na.strings = T)
header <- read.csv("~/qiulab/cdg/data/header.csv", sep = ",", header = F)

#Reduce dimension by cluster
d.snp <- dist(genes, method = "manhattan") #manhattan b/c 0/1
hc.snp <- hclust(d.snp)
# plot(hc.snp, las=2, main = "SNP groups")

grp.snps <- lapply(2:10, function(x) {cutree(hc.snp, k = x)}) #make groups

#target
cdg.mean <- tapply(cdg$logcdg, cdg$strains, mean)
# cdg.norm <- (cdg.mean - mean(cdg.mean))/(sd(cdg.mean))
cdg.norm <- scale(cdg.mean)

#Hyperparameters
eta <- 0.1 #learning rate
w <- list() #for 9 lists of weights
b <- runif(9) #bias for each group

#Split the training & prediction(target)
training <- 20
target <- 10

#Empty lists/matrices/vectors
names <- cbind(li.names=c("snp","m.snp","target.snp","weights"),
               num.names=c("acc_train", "acc_target", "group_size","rep"),
               matrices=c("a","y","e","accuracy"))
for (i in 1:4){
  assign(names[i,1], list())
  assign(names[i,2], numeric())
  ifelse(grepl("^ac",names[i,3]) == TRUE,
         assign(names[i,3], matrix(NA, nrow = 1000, ncol = 9)),
         assign(names[i,3], matrix(NA, nrow = training, ncol = 9)))
}


target.y <- list()
target.acc <- numeric()
for(reps in 1:10){
  for(l in 1:length(grp.snps)){
    group <- length(table(grp.snps[[l]]))  
    snp[[l]] <- matrix(0,nrow = 30, ncol = group)
    for(state in 1:group){
      snp[[l]] <- t(genes[sample(nrow(genes[which(grp.snps[[l]]==state),]),group),])
    }
    rownames(snp[[l]]) <- as.character(as.matrix(header[-1]))
    m.snp[[l]] <- snp[[l]][sample(training),] #change training model number
    m.snp[[l]] <- cbind(cdg.norm[row.names(m.snp[[1]])], m.snp[[l]])
    w[[l]] <- runif(group, 1e-3, 1e-2) #random weights for each list
    weights[[l]] <- matrix(0,nrow =1000, ncol = state)
    target.snp[[l]] <- snp[[l]][!rownames(snp[[1]]) %in% rownames(m.snp[[1]]),]
    target.snp[[l]] <- cbind(cdg.norm[row.names(target.snp[[1]])], target.snp[[l]])
    for(epoch in 1:1000){
      a[,l] <- m.snp[[l]][,2:c(group+1)] %*% w[[l]] #activation
      y[,l] <- 1/(1+exp(-a[,l]-b[l])) #output
      e[,l] <- m.snp[[l]][,1] - y[,l] #backpropogation
      w[[l]] <- w[[l]] - eta * -colSums(m.snp[[l]][,2:c(group+1)] * e[,l]) #update weights
      b[l] <- b[l] - eta * sum(-e[,l]) #update bias
      weights[[l]][epoch,] <- w[[l]]
      accuracy[epoch,l] <- cor(y[,l], m.snp[[l]][,1])
      # accuracy[epoch,l] <- sqrt(mean((y[,l]- m.snp[[l]][,1])^2))
    }
    acc_train <- c(acc_train,accuracy[1000,l])
    group_size <- c(group_size, group)
    rep <- c(rep, reps)
  }
  for(i in 1:9){
    target.y[[i]] <- 1/(1+exp(-(target.snp[[i]][,2:c(i+2)] %*% w[[i]])-b[i]))
    # target.acc <- c(target.acc, sqrt(mean((target.y[[i]]- target.snp[[l]][,1])^2)))
    target.acc[i] <- cor(target.y[[i]],target.snp[[i]][,1])
    acc_target <- c(acc_target, target.acc[i])
  }
decision <- data.frame(acc_train,acc_target,group_size, rep,stringsAsFactors = F)
}


library(reshape2)
df <- melt(decision, measure.vars = c("acc_train","acc_target"))
library(ggplot2)
ggplot(df, aes(x=factor(group_size), y=value, color=variable, group=variable)) + geom_line() + facet_wrap(~rep) + ggtitle("CDG Level Prediction w/ Single Neuron")

ggplot(df, aes(x=factor(group_size), y=value, color=variable, group=variable)) + geom_boxplot() + facet_wrap(~rep) + ggtitle("CDG Level Prediction w/ Single Neuron")

xyplot(acc_target ~ group_size | rep, data = decision, type="b", main = "cdg level prediction (hi/lo) by SNL -- training model of 20 strains permuted in each 9 groups in 10 epochs")
