library(randomForest)
library(caret)
library(e1071)

final.data=read.csv("./Data/Original.Data/final.data.csv", header = T)

final.data$gbif.quantile[final.data$gbif.range<=799.4]="Q1"
final.data$gbif.quantile[final.data$gbif.range>799.4 & final.data$gbif.range<=903.8]="Q2"
final.data$gbif.quantile[final.data$gbif.range>903.8 & final.data$gbif.range<=1015.8]="Q3"
final.data$gbif.quantile[final.data$gbif.range>1015.8 & final.data$gbif.range<=1172.8]="Q4"
final.data$gbif.quantile[final.data$gbif.range>1172.8]="Q5"

final.data$gbif.quantile=as.factor(final.data$gbif.quantile)

# all traits (10) with quantiles

shapiro.test(final.data$vla)
shapiro.test(final.data$delta13c)
shapiro.test(final.data$c)

all.trait.data=final.data[,c(13:16,23:28,36,37)] # original traits logged and scaled
all.trait.data.no.log=final.data[,c(13:22,37)] # original 6 traits no logged

# impute missing values
# imputed values are just the median value for that column
all.trait.data.2=na.roughfix(all.trait.data)
all.trait.data.2.no.log=na.roughfix(all.trait.data.no.log)

# 6 traits in pc logged
traits.in.pc.data=final.data[,c(23:28,37)]

# 6 traits in pc not logged
traits.in.pc.data.no.log=final.data[,c(17:22,37)]

# 3 PC axes used in FRic analyses
pc.traits.data=final.data[,c(29:31,37)]


# Predict Quantile based on all traits, original traits logged

set.seed(3)
id=sample(2, nrow(all.trait.data.2), prob=c(0.7,0.3), replace=TRUE)
train=all.trait.data.2[id==1,]
test=all.trait.data.2[id==2,]

# Tuning randomForest for the optimal mtry parameter
# mtry is the number of variables randomly sampled as candidates at each split
bestmtry=tuneRF(train, train$gbif.quantile, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 3, OOB error=13.92%

forest=randomForest(gbif.quantile~., data=train, importance=TRUE, ntree=1000)
# 81.01

pred.train=predict(forest, train, type="class")
table(pred.train, train$gbif.quantile)

pred.valid=predict(forest, test, type="class")
mean(pred.valid==test$gbif.quantile) 
table(pred.valid, test$gbif.quantile)

RF.all.trait.data.import=importance(forest)
#write.csv(RF.all.trait.data.import, file="./Data/gbif.outputs/Random.Forest.gbif/RF.all.trait.data.import.csv")
# log height, MDA decrease 1.4%
varImpPlot(forest)

confusionMatrix(data=pred.valid, reference=test$gbif.quantile)
# Accuracy = 0.1818
# No Information Rate = 0.3636
# P-value = 0.9933

print(forest)
#OOB = 81.01%
# No.variables tried at each split = 3

# Predict Quantile based on all traits, not traits logged

set.seed(3)
id.2=sample(2, nrow(all.trait.data.2.no.log), prob=c(0.7,0.3), replace=TRUE)
train.2=all.trait.data.2.no.log[id.2==1,]
test.2=all.trait.data.2.no.log[id.2==2,]

bestmtry.2=tuneRF(train.2, train.2$gbif.quantile, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 3, OOB = 0.1392%

forest.2=randomForest(gbif.quantile~., data=train.2, importance=TRUE, mtry=2)
# OOB = 82.28%
# 2 variables tried at each split

pred.train.2=predict(forest.2, train.2, type="class")
table(pred.train.2, train.2$gbif.quantile)

pred.valid.2=predict(forest.2, test.2, type="class")
mean(pred.valid.2==test.2$gbif.quantile)
#0.1515
table(pred.valid.2, test.2$gbif.quantile)

RF.all.trait.data.2.no.log.import=importance(forest.2)
#write.csv(RF.all.trait.data.2.no.log.import, file="./Data/gbif.outputs/Random.Forest.gbif/RF.all.trait.data.2.no.log.import.csv")
# height most important, MDA decrease 1.5%

varImpPlot(forest.2)

confusionMatrix(data=pred.valid.2, reference=test.2$gbif.quantile)
# Accuracy = 0.1515, 21.21%
# P-value = 0.9981

# Predict Quantile based on traits in pc included individually

set.seed(3)
id.2=sample(2, nrow(traits.in.pc.data), prob=c(0.7,0.3), replace=TRUE)
train.2=traits.in.pc.data[id.2==1,]
test.2=traits.in.pc.data[id.2==2,]

bestmtry.2=tuneRF(train.2, train.2$gbif.quantile, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 2, OOB = 15.19%

forest.2=randomForest(gbif.quantile~., data=train.2, importance=TRUE, mtry=2)
# OOB = 78.48%
# 2 variables tried at each split

pred.train.2=predict(forest.2, train.2, type="class")
table(pred.train.2, train.2$gbif.quantile)

pred.valid.2=predict(forest.2, test.2, type="class")
mean(pred.valid.2==test.2$gbif.quantile)
#0.2121
table(pred.valid.2, test.2$gbif.quantile)

RF.traits.in.pc.data.import=importance(forest.2)
#write.csv(RF.traits.in.pc.data.import, file="./Data/gbif.outputs/Random.Forest.gbif/RF.traits.in.pc.data.import.csv")
# nitrogen most important, MDA decrease 3.4%

varImpPlot(forest.2)

confusionMatrix(data=pred.valid.2, reference=test.2$gbif.quantile)
# Accuracy = 0.2121, 21.21%
# P-value = 0.98050

# Predict Quantile based on 3 pc axes scores used in FRic analyses

pc.traits.data.2 = pc.traits.data[,c(1:3,7)]

set.seed(3)
id.3=sample(2, nrow(pc.traits.data.2), prob=c(0.7,0.3), replace=TRUE)
train.3=pc.traits.data.2[id.3==1,]
test.3=pc.traits.data.2[id.3==2,]

bestmtry.3=tuneRF(train.3, train.3$gbif.quantile, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 2
# OOB = 1.27%

forest.3=randomForest(gbif.quantile~., data=train.3, importance=TRUE, mtry=2)
# OOB = 73.42%

pred.train.3=predict(forest.3, train.3, type="class")
table(pred.train.3, train.3$gbif.quantile)

pred.valid.3=predict(forest.3, test.3, type="class")
mean(pred.valid.3==test.3$gbif.quantile)
#0.2121
table(pred.valid.3, test.3$gbif.quantile)

RF.pc.traits.data.import=importance(forest.3)
write.csv(RF.pc.traits.data.import, file="./Data/gbif.outputs/Random.Forest.gbif/RF.3pc.traits.data.import.csv")
# Comp1 most important, MDA decrease 3.4%

varImpPlot(forest.3)

confusionMatrix(data=pred.valid.3, reference=test.3$gbif.quantile)
# Accuracy = 0.2121
# P.value = 0.7019


# Predict Wide/Narrow based on all traits

set.seed(3)
id.4=sample(2, nrow(range.all.trait.data.2), prob=c(0.7,0.3), replace=TRUE)
train.4=range.all.trait.data.2[id.4==1,]
test.4=range.all.trait.data.2[id.4==2,]

bestmtry.4=tuneRF(train.4, train.4$gbif.range.size, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 3
# OOB error = 0%

forest.4=randomForest(gbif.range.size~., data=train.4, importance=TRUE, mtry=3)
# OOB = 56.96%

pred.train.4=predict(forest.4, train.4, type="class")
table(pred.train.4, train.4$gbif.range.size)

pred.valid.4=predict(forest.4, test.4, type="class")
mean(pred.valid.4==test.4$gbif.range.size)
# 0.5757
table(pred.valid.4, test.4$gbif.range.size)

RF.range.all.trait.data.2.import=importance(forest.4)
write.csv(RF.range.all.trait.data.2.import, file="./Data/gbif.outputs/Random.Forest.gbif/RF.range.all.trait.data.2.import.csv")
# delta13c most important, MDA decrease 1.76%

varImpPlot(forest.4)

confusionMatrix(data=pred.valid.4, reference=test.4$gbif.range.size)
# Accuracy = 0.5758

# Predict Wide/Narrow based on traits included in pc, individually

set.seed(3)
id.5=sample(2, nrow(range.traits.in.pc.data), prob=c(0.7,0.3), replace=TRUE)
train.5=range.traits.in.pc.data[id.5==1,]
test.5=range.traits.in.pc.data[id.5==2,]

bestmtry.5=tuneRF(train.5, train.5$gbif.range.size, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 2
# OOB = 1.27%

forest.5=randomForest(gbif.range.size~., data=train.5, importance=TRUE, mtry=2)
# OOB = 56.96%

pred.train.5=predict(forest.5, train.5, type="class")
table(pred.train.5, train.5$gbif.range.size)

pred.valid.5=predict(forest.5, test.5, type="class")
mean(pred.valid.5==test.5$gbif.range.size)
#0.5757
table(pred.valid.5, test.5$gbif.range.size)

RF.range.traits.in.pc.data.import=importance(forest.5)
write.csv(RF.range.traits.in.pc.data.import, file="./Data/gbif.outputs/Random.Forest.gbif/RF.range.traits.in.pc.data.import.csv")
# log.height most important, MDA decrease 1.9%

varImpPlot(forest.5)

confusionMatrix(data=pred.valid.5, reference=test.5$gbif.range.size)
# Accuracy 0.5758

# Predict Wide/Narrow based on pc axes scores

set.seed(3)
id.6=sample(2, nrow(range.pc.traits.data), prob=c(0.7,0.3), replace=TRUE)
train.6=range.pc.traits.data[id.6==1,]
test.6=range.pc.traits.data[id.6==2,]

bestmtry.6=tuneRF(train.6, train.6$gbif.range.size, stepFactor=1.2, improve=0.01, trace=T, plot=T)
# mtry = 2
# OOB = 0.01265823%

forest.6=randomForest(gbif.range.size~., data=train.6, importance=TRUE, mtry=2)
# OOB = 54.43%

pred.train.6=predict(forest.6, train.6, type="class")
table(pred.train.6, train.6$gbif.range.size)

pred.valid.6=predict(forest.6, test.6, type="class")
mean(pred.valid.6==test.6$gbif.range.size)
#0.5151
table(pred.valid.6, test.6$gbif.range.size)

RF.range.pc.traits.data.import=importance(forest.6)
write.csv(RF.range.pc.traits.data.import, file="./Data/gbif.outputs/Random.Forest.gbif/RF.range.pc.traits.data.import.csv")
# Comp2 most important, MDA decrease 6.1%

varImpPlot(forest.6)

confusionMatrix(data=pred.valid.6, reference=test.6$gbif.range.size)
# Accuracy 0.5152
