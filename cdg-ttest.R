genes <- read.csv("~/qiulab/cdg/data/genes-pythonver.csv", header = T, row.names = 1)


p.value <- sapply(2:ncol(genes), function(x){
  if((sum(genes[,x]) > 1) & (!sum(genes[,x]) == 29)){
    t.test(genes[,1] ~ genes[,x])$p.value
  }else if(sum(genes[,x]) == 29){
    t.test(genes[,1], mu = 1)$p.value
  }else {
    t.test(genes[,1], mu= genes[which(genes[,x] == 1),][,1])$p.value
  }
})

p.value <- as.data.frame(p.value, row.names=colnames(genes[,2:1079]))

ttest.reduced <- as.matrix(genes[,c(1,which(p.value < 0.00001))])

w <- runif(ncol(ttest.reduced)-1, 1e-3, 1e-2)
b <- runif(1)

eta <- 0.1

weights <- matrix(0,nrow = 1000, ncol = ncol(ttest.reduced)-1) #equivalent to # of SNP position/state
accuracy <- numeric()

training <- 20
test <- 10

train.x <- ttest.reduced[sample(training),]
test.x <- ttest.reduced[!row.names(ttest.reduced) %in% rownames(train.x),]

for(epoch in 1:1000){
  a <- train.x[,2:ncol(train.x)] %*% w #weight each feature/linear comb
  y <- 1/(1+exp(-a-b)) #activation function
  e <- train.x[,1] - y #backpropogation
  w <- w - eta * -colSums(train.x[,2:ncol(train.x)] * e[,1]) #update weights
  b <- b - eta * sum(-e) #update bias
  for(k in 1:ncol(weights)){
    weights[epoch,k] <- w[k]
  }
  #accuracy <- c(accuracy, cor(y, train.x[,1]))
  accuracy <- c(accuracy, mean((train.x[,1] - y)^2))
}

#test
test.y <- 1/(1+exp(-(test.x[,2:ncol(test.x)] %*% w)-b))
# test.acc <- cor(test.y,test.x[,1])
test.acc <- mean((test.x[,1] - test.y)^2)

boxplot(accuracy, type = "p", main = "Accuracy of prediction for train.x(n=20)/model")
plot(test.acc, type = "p", main = "Accuracy of prediction")
plot(NA, xlim = c(1,1000), ylim = c(-10, 10), ylab = "weights range", main = "Weights Variation")
for(i in 1:ncol(ttest.reduced)-1){
  lines(x = 1:1000, y = weights[,i], col=sample(colors(distinct = T),1), lty=1)
}
