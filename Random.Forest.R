all.data=read.csv("all.data.csv", header=T, row.names=1)
all.data$Quantile=as.factor(all.data$Quantile)
all.data$Range.Size=as.factor(all.data$Range.Size)

# traits with quantiles
all.trait.data=all.data[,c(8:10,18,25:31)]
# impute missing values
# imputed values are just the median value for that column
all.trait.data.2=na.roughfix(all.trait.data)

traits.in.pc.data=all.data[,c(18,25:30)]
pc.traits.data=all.data[,c(18:24)]

# traits with Range.size narrow/wide
range.all.trait.data=all.data[,c(8:10,25:32)]
# impute missing values
# imputed values are just the median value for that column
range.all.trait.data.2=na.roughfix(range.all.trait.data)

range.traits.in.pc.data=all.data[,c(25:30,32)]
range.pc.traits.data=all.data[,c(19:24,32)]

library(randomForest)
library(caret)

set.seed(3)
id=sample(2, nrow(all.trait.data.2), prob=c(0.7,0.3), replace=TRUE)
train=all.trait.data.2[id==1,]
test=all.trait.data.2[id==2,]

# Predict Quantile based on all traits

# Tuning randomForest for the optimal mtry parameter
# mtry is the number of variables randomly sampled as candidates at each split
bestmtry=tuneRF(train, train$Quantile, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 3, OOB error=8.86%

forest=randomForest(Quantile~., data=train, importance=TRUE, ntree=1000)

pred.train=predict(forest, train, type="class")
table(pred.train, train$Quantile)

pred.valid=predict(forest, test, type="class")
mean(pred.valid==test$Quantile) 
table(pred.valid, test$Quantile)

RF.all.trait.data.import=importance(forest)
write.csv(RF.all.trait.data.import, file="RF.all.trait.data.import.csv")
# VLA most important, MDA decrease 12.4%, then log.height decrease 9.9%
varImpPlot(forest)

install.packages("e1071")
library(e1071)

confusionMatrix(data=pred.valid, reference=test$Quantile)
# Accuracy = 0.1515
# No Information Rate = 0.3333
# P-value = 0.9948

print(forest)
#OOB = 74.68%
# No.variables tried at each split = 3

# Predict Quantile based on traits in pc included individually

set.seed(3)
id.2=sample(2, nrow(traits.in.pc.data), prob=c(0.7,0.3), replace=TRUE)
train.2=traits.in.pc.data[id.2==1,]
test.2=traits.in.pc.data[id.2==2,]

bestmtry.2=tuneRF(train.2, train.2$Quantile, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 2, OOB = 8.86%

forest.2=randomForest(Quantile~., data=train.2, importance=TRUE, mtry=2)
# OOB = 79.75%
# 2 variables tried at each split

pred.train.2=predict(forest.2, train.2, type="class")
table(pred.train.2, train.2$Quantile)

pred.valid.2=predict(forest.2, test.2, type="class")
mean(pred.valid.2==test.2$Quantile)
#0.2121
table(pred.valid.2, test.2$Quantile)

RF.traits.in.pc.data.import=importance(forest.2)
write.csv(RF.traits.in.pc.data.import, file="RF.traits.in.pc.data.import.csv")
# height most important, MDA decrease 7.26%

varImpPlot(forest.2)

confusionMatrix(data=pred.valid.2, reference=test.2$Quantile)
# Accuracy = 0.2121, 21.21%
# P-value = 0.9566

# Predict Quantile based on 6 pc axes scores

set.seed(3)
id.3=sample(2, nrow(pc.traits.data), prob=c(0.7,0.3), replace=TRUE)
train.3=pc.traits.data[id.3==1,]
test.3=pc.traits.data[id.3==2,]

bestmtry.3=tuneRF(train.3, train.3$Quantile, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 2
# OOB = 12.66%

forest.3=randomForest(Quantile~., data=train.3, importance=TRUE, mtry=2)
# OOB = 77.22%

pred.train.3=predict(forest.3, train.3, type="class")
table(pred.train.3, train.3$Quantile)

pred.valid.3=predict(forest.3, test.3, type="class")
mean(pred.valid.3==test.3$Quantile)
#0.2424
table(pred.valid.3, test.3$Quantile)

RF.pc.traits.data.import=importance(forest.3)
write.csv(RF.pc.traits.data.import, file="RF.pc.traits.data.import.csv")
# Comp1 most important, MDA decrease 2.35%

varImpPlot(forest.3)

confusionMatrix(data=pred.valid.3, reference=test.3$Quantile)
# Accuracy = 0.2424, 24.24%
# P.value = 0.9050

# Predict Wide/Narrow based on all traits

set.seed(3)
id.4=sample(2, nrow(range.all.trait.data.2), prob=c(0.7,0.3), replace=TRUE)
train.4=range.all.trait.data.2[id.4==1,]
test.4=range.all.trait.data.2[id.4==2,]

bestmtry.4=tuneRF(train.4, train.4$Range.Size, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 3
# OOB error = 0%

forest.4=randomForest(Range.Size~., data=train.4, importance=TRUE, mtry=3)
# OOB = 50.63%

pred.train.4=predict(forest.4, train.4, type="class")
table(pred.train.4, train.4$Range.Size)

pred.valid.4=predict(forest.4, test.4, type="class")
mean(pred.valid.4==test.4$Range.Size)
# 0.4848
table(pred.valid.4, test.4$Range.Size)

RF.range.all.trait.data.2.import=importance(forest.4)
write.csv(RF.range.all.trait.data.2.import, file="RF.range.all.trait.data.2.import.csv")
# delta13c most important, MDA decrease 6.87%

varImpPlot(forest.4)

confusionMatrix(data=pred.valid.4, reference=test.4$Range.Size)
# Accuracy = 0.4848

# Predict Wide/Narrow based on traits included in pc, individually

set.seed(3)
id.5=sample(2, nrow(range.traits.in.pc.data), prob=c(0.7,0.3), replace=TRUE)
train.5=range.traits.in.pc.data[id.5==1,]
test.5=range.traits.in.pc.data[id.5==2,]

bestmtry.5=tuneRF(train.5, train.5$Range.Size, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 2
# OOB = 0%

forest.5=randomForest(Range.Size~., data=train.5, importance=TRUE, mtry=2)
# OOB = 48.1%

pred.train.5=predict(forest.5, train.5, type="class")
table(pred.train.5, train.5$Range.Size)

pred.valid.5=predict(forest.5, test.5, type="class")
mean(pred.valid.5==test.5$Range.Size)
#0.4242
table(pred.valid.5, test.5$Range.Size)

RF.range.traits.in.pc.data.import=importance(forest.5)
write.csv(RF.range.traits.in.pc.data.import, file="RF.range.traits.in.pc.data.import.csv")
# log.SLA most important, MDA decrease 6.53%

varImpPlot(forest.5)

confusionMatrix(data=pred.valid.5, reference=test.5$Range.Size)
# Accuracy 0.4242

# Predict Wide/Narrow based on pc axes scores

set.seed(3)
id.6=sample(2, nrow(range.pc.traits.data), prob=c(0.7,0.3), replace=TRUE)
train.6=range.pc.traits.data[id.6==1,]
test.6=range.pc.traits.data[id.6==2,]

bestmtry.6=tuneRF(train.6, train.6$Range.Size, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 2
# OOB = 0%

forest.6=randomForest(Range.Size~., data=train.6, importance=TRUE, mtry=2)
# OOB = 55.7%

pred.train.6=predict(forest.6, train.6, type="class")
table(pred.train.6, train.6$Range.Size)

pred.valid.6=predict(forest.6, test.6, type="class")
mean(pred.valid.6==test.6$Range.Size)
#0.6060
table(pred.valid.6, test.6$Range.Size)

RF.range.pc.traits.data.import=importance(forest.6)
write.csv(RF.range.pc.traits.data.import, file="RF.range.pc.traits.data.import.csv")
# Comp2 most important, MDA decrease 6.31%

varImpPlot(forest.6)

confusionMatrix(data=pred.valid.6, reference=test.6$Range.Size)
# Accuracy 0.6061
