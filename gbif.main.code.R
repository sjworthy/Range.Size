library(vegan)
library(geometry)
devtools::install_version("rangemodelR", "1.0.4")
library(rangemodelR)

# code to redo analysis based on gbif data

final.data=read.csv("./Data/Original.Data/final.data.csv", header = T)

#### Putting Species into Range Size Quantiles ####

### Divide the range sizes from 0 to 1 by .2 or 20%
quantile(final.data$gbif.range, seq(0, 1, by=.2))

### Determine the number of species in each quantile ####

sum(final.data$gbif.range<=799.4)
[1] 23
sum(final.data$gbif.range>799.4 & final.data$gbif.range<=903.8)
[1] 22
sum(final.data$gbif.range>903.8 & final.data$gbif.range<=1015.8)
[1] 22
sum(final.data$gbif.range>1015.8 & final.data$gbif.range<=1172.8)
[1] 22
sum(final.data$gbif.range>1172.8)
[1] 23

### If the value in column range of all.data meets the criteria stated, place the quantile name in column Quantile

final.data$gbif.quantile[final.data$gbif.range<=799.4]="Q1"
final.data$gbif.quantile[final.data$gbif.range>799.4 & final.data$gbif.range<=903.8]="Q2"
final.data$gbif.quantile[final.data$gbif.range>903.8 & final.data$gbif.range<=1015.8]="Q3"
final.data$gbif.quantile[final.data$gbif.range>1015.8 & final.data$gbif.range<=1172.8]="Q4"
final.data$gbif.quantile[final.data$gbif.range>1172.8]="Q5"

### If the value in column range of sonadora.csv meets the criteria stated, place the quantile name in column Quantile
### Get a list of the species that are in each quantile
### Add Q1-Q5 for each of these species to the Quantile column

sonadora=read.csv("./Data/sonadora.csv", header=T)
colnames(sonadora)[7]="SpeciesCode"

# add gbif range to sonoadora

sonadora.2=merge(sonadora, final.data[,c("SpeciesCode","gbif.range","gbif.quantile","gbif.quartile")], by="SpeciesCode")

#### Percentage of abundance in each quantile ####

table(sonadora.2$gbif.quantile)
# Q1
405/8693 = 4.66%
# Q2
709/8693 = 8.16%
# Q3
1562/8693 = 17.99%
# Q4
1086/8693 = 12.49%
# Q5
4931/8693 = 56.72%

#### Making a Community Data Matrix ####

### Make a community data matrix with number of species in each elevation plot
### Make a blank community data matrix with 111 columns, one for each species, and 16 rows, one for each plot in the transect (250m-1000m)

abund.cdm=matrix(data=NA, ncol=112, nrow=16)
colnames(cdm)=all.data$SpeciesCode
rownames(cdm)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
                700, 750, 800, 850, 900, 950, 1000)

pa.cdm=decostand(abund.cdm, "pa")

### Fill in the community data matrix
## Subset the abundance data set by the plot elevation and take that subsetted object and table the SpeciesCode column.
# Enter the values manually because it won't auto place due to missing columns.

elevation=subset(sonadora, sonadora$plotElevation==1000)
output=table(elevation$stemSpeciesCode)
output

# abund.cdm=read.csv("./Data/gbif.outputs/abund.cdm.csv",header=T,row.names = 1)

### Make a community data matrix with just narrow species and one with just wide-ranging species ####
### Species with range < median (970.5) are narrow with species > median range being wide-ranging.

mean(final.data$gbif.range) # 950 m
median(final.data$gbif.range) #970.5

sum(final.data$gbif.range<950)
[1] 54
sum(final.data$gbif.range>950)
[1] 57
sum(final.data$gbif.range == 950)
[1] 1

sum(final.data$gbif.range<970.5)
[1] 56
sum(final.data$gbif.range>970.5)
[1] 56

narrow.species=as.character(final.data$SpeciesCode[final.data$gbif.range < 970.5])
# 56 narrow species 
wide.species=as.character(final.data$SpeciesCode[final.data$gbif.range > 970.5])
# 56 wide-ranging species 

narrow.cdm=abund.cdm[,narrow.species]
wide.cdm=abund.cdm[,wide.species]

#write.csv(narrow.cdm, file="narrow.cdm.gbif.csv")
#write.csv(wide.cdm, file="wide.cdm.gbif.csv")

final.data$gbif.range.size[final.data$SpeciesCode %in% narrow.species]="Narrow"
final.data$gbif.range.size[final.data$SpeciesCode %in% wide.species]="Wide"

### Make a community data matrix where quantiles are the rows/communities ####
### Make a blank community data matrix with 112 columns, one for each species, and 5 rows, one for each quantile.

quant.cdm=matrix(data=NA, nrow=5, ncol=112)
colnames(quant.cdm)=final.data$SpeciesCode
rownames(quant.cdm)=c("Q1", "Q2", "Q3", "Q4", "Q5")

quant.cdm = read.csv("./Data/gbif.outputs/quant.cdm.gbif.csv", row.names = 1)

quant.pa.cdm=decostand(quant.cdm, "pa")

# write.csv(quant.pa.cdm, file="./Data/gbif.outputs/quant.cdm.pa.gbif.csv")

### Fill in the community data matrix
### Subset the abundance data set by the quantile and take that subsetted object and table the SpeciesCode
# Enter the values manually because it won't auto place due to missing columns.

quantile=subset(sonadora.2, sonadora.2$gbif.quantile=="Q1")
output=table(quantile$SpeciesCode)
output

quartile=subset(sonadora.2, sonadora.2$gbif.quartile=="Q4")
output=table(quartile$SpeciesCode)
output

#### Range size correlation with elevation ####
# Hypothesis 1
# correlate range size of each individual with elevation of each individual, not species level

cor.test(sonadora.2$gbif.range, sonadora.2$plotElevation)
# r = 0.3726173, p < 0.0005
lm.test=lm(gbif.range~plotElevation, data=sonadora.2)

plot(sonadora.2$plotElevation,sonadora.2$gbif.range)

# for each elevation need unique species list and combine into one dataframe
# percent of each RSQ in each plot

sp.250 = plot250[!duplicated(plot250$SpeciesCode),]
nrow(sp.250)
table(sp.250$gbif.quantile)/31
# Q5: 0.23%, highest
sp.300 = plot300[!duplicated(plot300$SpeciesCode),]
nrow(sp.300)
table(sp.300$gbif.quantile)/50
# Q4: 0.28%, highest
sp.350 = plot350[!duplicated(plot350$SpeciesCode),]
nrow(sp.350)
table(sp.350$gbif.quantile)/42
# Q1: 0.24%, highest
sp.400 = plot400[!duplicated(plot400$SpeciesCode),]
nrow(sp.400)
table(sp.400$gbif.quantile)/46
# Q4: 0.26%, highest
sp.450 = plot450[!duplicated(plot450$SpeciesCode),]
nrow(sp.450)
table(sp.450$gbif.quantile)/39
# Q3 and Q4: 0.28%, highest
sp.500 = plot500[!duplicated(plot500$SpeciesCode),]
nrow(sp.500)
table(sp.500$gbif.quantile)/44
# Q3 and Q5: 0.25%, highest
sp.550 = plot550[!duplicated(plot550$SpeciesCode),]
nrow(sp.550)
table(sp.550$gbif.quantile)/27
# Q5: 0.33%, highest
sp.600 = plot600[!duplicated(plot600$SpeciesCode),]
nrow(sp.600)
table(sp.600$gbif.quantile)/30
# Q5: 0.30%, highest
sp.650 = plot650[!duplicated(plot650$SpeciesCode),]
nrow(sp.650)
table(sp.650$gbif.quantile)/38
# Q5: 0.32%, highest
sp.700 = plot700[!duplicated(plot700$SpeciesCode),]
nrow(sp.700)
table(sp.700$gbif.quantile)/37
# Q5: 0.35%, highest
sp.750 = plot750[!duplicated(plot750$SpeciesCode),]
nrow(sp.750)
table(sp.750$gbif.quantile)/30
# Q5: 0.40%, highest
sp.800 = plot800[!duplicated(plot800$SpeciesCode),]
nrow(sp.800)
table(sp.800$gbif.quantile)/30
# Q5: 0.43%, highest
sp.850 = plot850[!duplicated(plot850$SpeciesCode),]
nrow(sp.850)
table(sp.850$gbif.quantile)/25
# Q5: 0.44%, highest
sp.900 = plot900[!duplicated(plot900$SpeciesCode),]
nrow(sp.900)
table(sp.900$gbif.quantile)/26
# Q5: 0.46%, highest
sp.950 = plot950[!duplicated(plot950$SpeciesCode),]
nrow(sp.950)
table(sp.950$gbif.quantile)/17
# Q5: 0.53%, highest
sp.1000 = plot1000[!duplicated(plot1000$SpeciesCode),]
nrow(sp.1000)
table(sp.1000$gbif.quantile)/20
# Q5: 0.50%, highest

all.unique.sonadora = rbind(sp.250,sp.300,sp.350,sp.400,sp.450,sp.500,sp.550,sp.600,sp.650,
                            sp.700,sp.750,sp.800,sp.850,sp.900,sp.950,sp.1000)

cor.test(all.unique.sonadora$gbif.range, all.unique.sonadora$plotElevation)
# r = 0.2834519, p < 0.0005

### Make community data matrices for each quantile separately ####
### Get a list of the species in each quantile

Q1=as.character(final.data$SpeciesCode[final.data$gbif.quantile=="Q1"])
Q2=as.character(final.data$SpeciesCode[final.data$gbif.quantile=="Q2"])
Q3=as.character(final.data$SpeciesCode[final.data$gbif.quantile=="Q3"])
Q4=as.character(final.data$SpeciesCode[final.data$gbif.quantile=="Q4"])
Q5=as.character(final.data$SpeciesCode[final.data$gbif.quantile=="Q5"])

### Make empty cdms
q1.cdm=matrix(data=NA, nrow=16, ncol=23)
q2.cdm=matrix(data=NA, nrow=16, ncol=22)
q3.cdm=matrix(data=NA, nrow=16, ncol=22)
q4.cdm=matrix(data=NA, nrow=16, ncol=22)
q5.cdm=matrix(data=NA, nrow=16, ncol=23)

colnames(q1.cdm)=Q1
colnames(q2.cdm)=Q2
colnames(q3.cdm)=Q3
colnames(q4.cdm)=Q4
colnames(q5.cdm)=Q5

rownames(q1.cdm)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
                   700, 750, 800, 850, 900, 950, 1000)
rownames(q2.cdm)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
                   700, 750, 800, 850, 900, 950, 1000)
rownames(q3.cdm)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
                   700, 750, 800, 850, 900, 950, 1000)
rownames(q4.cdm)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
                   700, 750, 800, 850, 900, 950, 1000)
rownames(q5.cdm)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
                   700, 750, 800, 850, 900, 950, 1000)

Q1.subset=sonadora.2[sonadora.2$gbif.quantile=="Q1",]
Q2.subset=sonadora.2[sonadora.2$gbif.quantile=="Q2",]
Q3.subset=sonadora.2[sonadora.2$gbif.quantile=="Q3",]
Q4.subset=sonadora.2[sonadora.2$gbif.quantile=="Q4",]
Q5.subset=sonadora.2[sonadora.2$gbif.quantile=="Q5",]

#write.csv(Q1.subset,file="Q1.subset.csv")
#write.csv(Q2.subset,file="Q2.subset.csv")
#write.csv(Q3.subset,file="Q3.subset.csv")
#write.csv(Q4.subset,file="Q4.subset.csv")
#write.csv(Q5.subset,file="Q5.subset.csv")

Q1.subset=sonadora.2[sonadora.2$gbif.quartile=="Q1",]
Q2.subset=sonadora.2[sonadora.2$gbif.quartile=="Q2",]
Q3.subset=sonadora.2[sonadora.2$gbif.quartile=="Q3",]
Q4.subset=sonadora.2[sonadora.2$gbif.quartile=="Q4",]

q1.cdm = read.csv("./Data/gbif.outputs/q1.cdm.gbif.csv", row.names = 1)
q2.cdm = read.csv("./Data/gbif.outputs/q2.cdm.gbif.csv", row.names = 1)
q3.cdm = read.csv("./Data/gbif.outputs/q3.cdm.gbif.csv", row.names = 1)
q4.cdm = read.csv("./Data/gbif.outputs/q4.cdm.gbif.csv", row.names = 1)
q5.cdm = read.csv("./Data/gbif.outputs/q5.cdm.gbif.csv", row.names = 1)

### Make community data matrix where communities are quantiles of each plot and columns are species in each plot ####

elevation250=subset(sonadora.2, sonadora.2$plotElevation=="250")
elevation300=subset(sonadora.2, sonadora.2$plotElevation=="300")
elevation350=subset(sonadora.2, sonadora.2$plotElevation=="350")
elevation400=subset(sonadora.2, sonadora.2$plotElevation=="400")
elevation450=subset(sonadora.2, sonadora.2$plotElevation=="450")
elevation500=subset(sonadora.2, sonadora.2$plotElevation=="500")
elevation550=subset(sonadora.2, sonadora.2$plotElevation=="550")
elevation600=subset(sonadora.2, sonadora.2$plotElevation=="600")
elevation650=subset(sonadora.2, sonadora.2$plotElevation=="650")
elevation700=subset(sonadora.2, sonadora.2$plotElevation=="700")
elevation750=subset(sonadora.2, sonadora.2$plotElevation=="750")
elevation800=subset(sonadora.2, sonadora.2$plotElevation=="800")
elevation850=subset(sonadora.2, sonadora.2$plotElevation=="850")
elevation900=subset(sonadora.2, sonadora.2$plotElevation=="900")
elevation950=subset(sonadora.2, sonadora.2$plotElevation=="950")
elevation1000=subset(sonadora.2, sonadora.2$plotElevation=="1000")

quantile=subset(elevation1000, elevation1000$gbif.quantile=="Q5")
output=table(quantile$SpeciesCode)
output

quartile=subset(elevation1000, elevation1000$gbif.quartile=="Q4")
output=table(quartile$SpeciesCode)
output

### Make trait data matrix for each quantile ####
### Make character string for species within each quantile

tdm = read.csv("./Data/Trait.Matrices/tdm.csv", header=T, row.names = 1)
tdm.scaled=apply(log(tdm[,c(1:7)]), MARGIN=2, scale)
row.names(tdm.scaled)=row.names(tdm)
colnames(tdm.scaled)=c("log.height.ft", "log.la", "log.sla",
                       "log.wood", "log.n", "log.p", "log.seed.mass")

write.csv(tdm.scaled, file="tdm.scaled.csv")
tdm.scaled=read.csv("tdm.scaled.csv", header = T, row.names = 1)


# Perform a principal components analysis (PCA) on traits
### Performing a PCA will reduce trait redundancy and produce orthogonal axes of function.
# PCA can only be performed for traits with all values 

pca.traits=princomp(tdm.scaled[,c(1:6)])

### Determine the proportion of the total variance explained by each axis of PCA

summary(pca.traits)

# PC1 41.7%
# PC2 20%
# PC3 16%
# PC1-PC3 = 78%

### Determine which traits are most heavily weighted on these axes

print(pca.traits$loadings)

### Axis 1: Negatively associated with Leaf Nitrogen and Leaf Phosphorus
### Axis 2: Positively associated with Maximum Height and Leaf Area
### Axis 3: Positively associated with Wood Specific Gravity and Specific Leaf Area

### Get pca scores for first 3 axes and put them into new matrix
### The first 3 axes explain 78% of the variation in the data.
### Give row names from the original trait matrix

pca.scores=pca.traits$scores[,1:3]
rownames(pca.scores)=rownames(tdm.scaled)
pca.scores.all=pca.traits$scores
rownames(pca.scores.all)=rownames(tdm.scaled)

#write.csv(pca.scores, file="pca.scores.csv")
#write.csv(pca.scores.all, file="pca.scores.all.csv")

tdm.scaled=read.csv("./Data/Trait.Matrices/tdm.scaled.csv", header=T, row.names = 1)

Q1=as.character(final.data$SpeciesCode[final.data$gbif.quantile=="Q1"])
Q2=as.character(final.data$SpeciesCode[final.data$gbif.quantile=="Q2"])
Q3=as.character(final.data$SpeciesCode[final.data$gbif.quantile=="Q3"])
Q4=as.character(final.data$SpeciesCode[final.data$gbif.quantile=="Q4"])
Q5=as.character(final.data$SpeciesCode[final.data$gbif.quantile=="Q5"])

### Subset original trait data matrix by the species in each quantile character string

Q1.traits=tdm.scaled[Q1,]
Q2.traits=tdm.scaled[Q2,]
Q3.traits=tdm.scaled[Q3,]
Q4.traits=tdm.scaled[Q4,]
Q5.traits=tdm.scaled[Q5,]

#write.csv(Q1.traits, file="q1.tdm.gbif.csv")
#write.csv(Q2.traits, file="q2.tdm.gbif.csv")
#write.csv(Q3.traits, file="q3.tdm.gbif.csv")
#write.csv(Q4.traits, file="q4.tdm.gbif.csv")
#write.csv(Q5.traits, file="q5.tdm.gbif.csv")

Q1=as.character(final.data$SpeciesCode[final.data$gbif.quartile=="Q1"])
Q2=as.character(final.data$SpeciesCode[final.data$gbif.quartile=="Q2"])
Q3=as.character(final.data$SpeciesCode[final.data$gbif.quartile=="Q3"])
Q4=as.character(final.data$SpeciesCode[final.data$gbif.quartile=="Q4"])

### Subset original trait data matrix by the species in each quantile character string

Q1.traits=tdm.scaled[Q1,]
Q2.traits=tdm.scaled[Q2,]
Q3.traits=tdm.scaled[Q3,]
Q4.traits=tdm.scaled[Q4,]

#write.csv(Q1.traits, file="q1.tdm.gbif.quartile.csv")
#write.csv(Q2.traits, file="q2.tdm.gbif.quartile.csv")
#write.csv(Q3.traits, file="q3.tdm.gbif.quartile.csv")
#write.csv(Q4.traits, file="q4.tdm.gbif.quartile.csv")

q1.tdm = read.csv("./Data/gbif.outputs/q1.tdm.gbif.csv", row.names = 1)
q2.tdm = read.csv("./Data/gbif.outputs/q2.tdm.gbif.csv", row.names = 1)
q3.tdm = read.csv("./Data/gbif.outputs/q3.tdm.gbif.csv", row.names = 1)
q4.tdm = read.csv("./Data/gbif.outputs/q4.tdm.gbif.csv", row.names = 1)
q5.tdm = read.csv("./Data/gbif.outputs/q5.tdm.gbif.csv", row.names = 1)


#### Calculating Functional Richness ####
### Because dbFD function reduces the data set's dimensionality automatically, it uses a different number of dimensions for each community. To keep number of dimensions constant, must determine convex hull volume manually for each community.

### FRic for each elevation plot for entire data set.
### pca only includes 1st three axes

hull.matrix.all=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
  q.hull=convhulln(pca.scores[names(abund.cdm[i,abund.cdm[i,]>0]),], options="FA")
  q.hull.vol=q.hull$vol
  hull.matrix.all[i,1]=q.hull.vol
}

### FRic for each elevation plot excluding wide ranging species
### pca only includes 1st three axes of narrow ranging species

pca.scores=read.csv("./Data/Trait.Matrices/pca.scores.csv", row.names = 1)

names=as.factor(colnames(narrow.cdm))
narrow.pca=subset(pca.scores, rownames(pca.scores)%in%names)
narrow.pca=narrow.pca[,1:3]

hull.matrix.narrow=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
  q.hull=convhulln(narrow.pca[names(narrow.cdm[i,narrow.cdm[i,]>0]),], options="FA")
  q.hull.vol=q.hull$vol
  hull.matrix.narrow[i,1]=q.hull.vol
}

# 950 don't work because only 2 species

### Get volume when there are only 3 species
### Get the names of species in the community that have an abundance greater than 0
### Subset pc scores for the names in that community
### Get the product of the max values minus the min values 

com.950.names=names(narrow.cdm[15, narrow.cdm[15,]>0])
prod(apply(narrow.pca[com.950.names,], MARGIN=2,max)-
       apply(narrow.pca[com.950.names,], MARGIN=2, min))

# write.csv(hull.matrix.narrow, file="./Data/gbif.outputs/hull.matrix.narrow.pca.gbif.csv")

### FRic for each elevation plot excluding narrow ranging species
### pca only includes 1st three axes of wide ranging species

names.2=as.factor(colnames(wide.cdm))
wide.pca=subset(pca.scores, rownames(pca.scores)%in%names.2)
wide.pca=wide.pca[,1:3]

hull.matrix.wide=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
  q.hull=convhulln(wide.pca[names(wide.cdm[i,wide.cdm[i,]>0]),], options="FA")
  q.hull.vol=q.hull$vol
  hull.matrix.wide[i,1]=q.hull.vol
}

# write.csv(hull.matrix.wide, file="./Data/gbif.outputs/hull.matrix.wide.pca.gbif.csv")

### Abundance of narrow and wide ranging individuals in each elevation ####

rowSums(narrow.cdm)
rowSums(wide.cdm)

narrow.cdm.pa=decostand(narrow.cdm, "pa")
wide.cdm.pa=decostand(wide.cdm, "pa")

rowSums(narrow.cdm.pa)
rowSums(wide.cdm.pa)

#write.csv(narrow.cdm.pa, file="./Data/gbif.outputs/narrow.cdm.pa.gbif.csv")
#write.csv(wide.cdm.pa, file="./Data/gbif.outputs/wide.cdm.pa.gbif.csv")

#### Determine the average range size for each plot ####

plot250=subset(sonadora.2, sonadora.2$plotElevation==250)
plot300=subset(sonadora.2, sonadora.2$plotElevation==300)
plot350=subset(sonadora.2, sonadora.2$plotElevation==350)
plot400=subset(sonadora.2, sonadora.2$plotElevation==400)
plot450=subset(sonadora.2, sonadora.2$plotElevation==450)
plot500=subset(sonadora.2, sonadora.2$plotElevation==500)
plot550=subset(sonadora.2, sonadora.2$plotElevation==550)
plot600=subset(sonadora.2, sonadora.2$plotElevation==600)
plot650=subset(sonadora.2, sonadora.2$plotElevation==650)
plot700=subset(sonadora.2, sonadora.2$plotElevation==700)
plot750=subset(sonadora.2, sonadora.2$plotElevation==750)
plot800=subset(sonadora.2, sonadora.2$plotElevation==800)
plot850=subset(sonadora.2, sonadora.2$plotElevation==850)
plot900=subset(sonadora.2, sonadora.2$plotElevation==900)
plot950=subset(sonadora.2, sonadora.2$plotElevation==950)
plot1000=subset(sonadora.2, sonadora.2$plotElevation==1000)

# second value is where each species only considered once so species-level not individual level
mean(plot250$gbif.range) = 880.956
mean(unique(plot250$gbif.range)) = 954
mean(plot300$gbif.range) = 992.3054
mean(unique(plot300$gbif.range)) = 924.2766
mean(plot350$gbif.range) = 1076.953
mean(unique(plot350$gbif.range)) = 927.2308
mean(plot400$gbif.range) = 1059.351
mean(unique(plot400$gbif.range)) = 956.6364
mean(plot450$gbif.range) = 1059.995
mean(unique(plot450$gbif.range)) = 1022.053
mean(plot500$gbif.range) = 1101.338
mean(unique(plot500$gbif.range)) = 1012
mean(plot550$gbif.range) = 1109.118
mean(unique(plot550$gbif.range)) = 1059.269
mean(plot600$gbif.range) = 1132.853
mean(unique(plot600$gbif.range)) = 1054.069
mean(plot650$gbif.range) = 1076.295
mean(unique(plot650$gbif.range)) = 1040.75
mean(plot700$gbif.range) = 1119.377
mean(unique(plot700$gbif.range)) = 1056.611
mean(plot750$gbif.range) = 1188.019
mean(unique(plot750$gbif.range)) = 1088.567
mean(plot800$gbif.range) = 1172.045
mean(unique(plot800$gbif.range)) = 1076.467
mean(plot850$gbif.range) = 1124.88
mean(unique(plot850$gbif.range)) = 1097.64
mean(plot900$gbif.range) = 1183.781
mean(unique(plot900$gbif.range)) = 1097
mean(plot950$gbif.range) = 1255.433
mean(unique(plot950$gbif.range)) = 1139.312
mean(plot1000$gbif.range) = 1211.722
mean(unique(plot1000$gbif.range)) = 1104.1

# correlation between average range size per plot and elevation
aver.range.plot=c(880.956,992.3054,1076.953,1059.351,1059.995,1101.338,1109.118,1132.853,
                  1076.295,1119.377,1188.019,1172.045,1124.88,1183.781,1255.433,1211.722)

aver.range.plot.unique=c(954,924.2766,927.2308,956.6364,1022.053,1012,1059.269,1054.069,1040.75,
                  1056.611,1088.567,1076.467,1097.64,1097,1139.312,1104.1)


elevation = c(250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000)

cor.test(aver.range.plot, elevation)
# r = 0.43, p = 0.10 # old
# r = 0.8825424, p = 0.00000603 #gbif

cor.test(aver.range.plot.unique, elevation)
# r = 0.938, p < 0.0001

# correlation between range size and elevation

#### percentage of abundance and richness made up of wide and narrow ranging species ####
percent.wide.abund=rowSums(wide.cdm)/rowSums(abund.cdm)
percent.wide.rich=rowSums(wide.cdm.pa)/rowSums(pa.cdm)

percent.narrow.abund=rowSums(narrow.cdm)/rowSums(abund.cdm)
percent.narrow.rich=rowSums(narrow.cdm.pa)/rowSums(pa.cdm)

#write.csv(percent.wide.abund, file="./Data/gbif.outputs/percent.wide.abund.gbif.csv")
#write.csv(percent.wide.rich, file="./Data/gbif.outputs/percent.wide.rich.gbif.csv")
#write.csv(percent.narrow.abund, file="./Data/gbif.outputs/percent.narrow.abund.gbif.csv")
#write.csv(percent.narrow.rich, file="./Data/gbif.outputs/percent.narrow.rich.gbif.csv")

#### correlation test between the percentage of richness and abundance of wide and narrow ranging species and elevation ####
wide.abund=as.numeric(percent.wide.abund)
cor.test(wide.abund, elevation)
# r = 0.48, p = 0.06 # old
# r = 0.87, p = 0.00001165 # gbif
wide.rich=as.numeric(percent.wide.rich)
cor.test(wide.rich, elevation)
# r = -0.12, p = 0.66 # old
# r = 0.90, p = 0.000002209 # gbif

narrow.abund=as.numeric(percent.narrow.abund)
cor.test(narrow.abund, elevation)
# r = -0.48, p = 0.06 # old
# r = -0.87, p = 0.00001165 # gbif

narrow.rich=as.numeric(percent.narrow.rich)
cor.test(narrow.rich, elevation)
# r = 0.12, p = 0.66 # old
# r = -90, p = 0.000002209 # gbif

#### correlation between abundance made of wide ranging species and elevation ####
wide.abundance=as.numeric(rowSums(wide.cdm))
cor.test(elevation, wide.abundance)
# r = 0.61, p = 0.01 # old
# r = 0.77, p = 0.0005 # gbif

plot(wide.abundance~elevation, xlab="Elevation (m)", ylab="Wide-Ranging Species Abundance", pch=19,cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(lm(wide.abundance~elevation))

#### correlation between abundance of narrow ranging species and elevation #####
narrow.abundance=as.numeric(rowSums(narrow.cdm))
cor.test(elevation, narrow.abundance)
# r = 0.21, p = 0.45 # old
# r = -0.55, p = 0.03 # gbif

plot(narrow.abundance~elevation, xlab="Elevation (m)", ylab="Narrow-Ranging Species Abundance", pch=19,cex.lab=1.5, cex.axis=1.5)
abline(lm(narrow.abundance~elevation))

#### Correlation between narrow-ranging species richness and elevation ####
narrow.rich=as.numeric(rowSums(narrow.cdm.pa))
cor.test(elevation, narrow.rich)
# r = -0.72, p = 0.002 # old
# r = -0.86, p = 0.000023 # gbif

plot(narrow.rich~elevation, xlab="Elevation (m)", ylab="Narrow-Ranging Species Richness", pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(lm(narrow.rich~elevation))

#### Correlation between wide-ranging species richness and elevation ####
wide.rich=as.numeric(rowSums(wide.cdm.pa))
cor.test(elevation, wide.rich)
# r = -0.73, p = 0.001 #old
# r = -0.49, p = 0.05 #gbif

plot(wide.rich~elevation, xlab="Elevation (m)", ylab="Wide-Ranging Species Richness", pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(lm(wide.rich~elevation))

#### Correlation between species range size and abundance #####
cor.test(final.data$gbif.range, final.data$abundance)
# r = 0.09, p = 0.37 # old
# r = 0.32, p = 0.0005 # gbif

#### Correlation plot between functional richness and elevation ####
# Hypothesis 2

hull.matrix.pca=read.csv("./Data/gbif.outputs/hull.matrix.pca.gbif.csv", header=T)
hull.matrix.data=hull.matrix.pca[c(1:16),c(2:4)]

cor.test(elevation, hull.matrix.data$Allsp.pca)
# r = -0.66, p = 0.006
cor.test(elevation, hull.matrix.data$Narrowsp.pca)
# r = -0.78, p = 0.0003196
cor.test(elevation, hull.matrix.data$Widesp.pca)
# r = -0.43, p = 0.09978

hull.matrix.data.2 = hull.matrix.pca[c(1:16),c(14:18)]

cor.test(elevation, hull.matrix.data.2$Fric.minus.Q1)
# r = -0.55, p = 0.02695
cor.test(elevation, hull.matrix.data.2$Fric.minus.Q2)
# r = 0.19, p = 0.4663
cor.test(elevation, hull.matrix.data.2$Fric.minus.Q3)
# r = -0.09, p = 0.7278
cor.test(elevation, hull.matrix.data.2$Fric.minus.Q4)
# r = -0.62, p = 0.01034
cor.test(elevation, hull.matrix.data.2$Fric.minus.Q5)
# r = 0.79, p = 0.0002569

#### correlation test between species richness of each quantile and overall species richness ####

q1.cdm = read.csv("./Data/gbif.outputs/q1.cdm.gbif.csv", row.names = 1)
q2.cdm = read.csv("./Data/gbif.outputs/q2.cdm.gbif.csv", row.names = 1)
q3.cdm = read.csv("./Data/gbif.outputs/q3.cdm.gbif.csv",row.names = 1)
q4.cdm = read.csv("./Data/gbif.outputs/q4.cdm.gbif.csv",row.names = 1)
q5.cdm = read.csv("./Data/gbif.outputs/q5.cdm.gbif.csv", row.names = 1)

q1.cdm.pa = decostand(q1.cdm, "pa")
q2.cdm.pa = decostand(q2.cdm, "pa")
q3.cdm.pa = decostand(q3.cdm, "pa")
q4.cdm.pa = decostand(q4.cdm, "pa")
q5.cdm.pa = decostand(q5.cdm, "pa")

## Table.2
cor.test(rowSums(pa.cdm), rowSums(q1.cdm.pa))
# r = 0.47, p = 0.06 # old
# r = 0.81, p = 0.0001274 # gbif
cor.test(rowSums(pa.cdm), rowSums(q2.cdm.pa))
# r = 0.74, p = 0.001 #old
# r = 0.81, p = 0.0001661 # gbif
cor.test(rowSums(pa.cdm), rowSums(q3.cdm.pa))
# r = 0.87, p < 0.001 #old
# r = 0.75, p = 0.0008297 # gbif
cor.test(rowSums(pa.cdm), rowSums(q4.cdm.pa))
# r = 0.90, p < 0.001 #old
# r = 0.84, p = 0.0004687 # gbif
cor.test(rowSums(pa.cdm), rowSums(q5.cdm.pa))
# r = 0.64, p = 0.01 #old
# r = -0.20, p = 0.4497 # gbif

##### Correlation between species richness and range size across elevation ####
# Hypothesis 1

elevation=as.numeric(row.names(q1.cdm.pa))

cor.test(rowSums(q1.cdm.pa), elevation)
# r = -0.84, p = 0.00005338 # gbif
cor.test(rowSums(q2.cdm.pa), elevation)
# r = -0.50, p = 0.04668 # gbif
cor.test(rowSums(q3.cdm.pa), elevation)
# r = -0.52, p = 0.03527 # gbif
cor.test(rowSums(q4.cdm.pa), elevation)
# r = -0.79, p = 0.0002367 # gbif
cor.test(rowSums(q5.cdm.pa), elevation)
# r = 0.60, p = 0.01 # gbif

#### linear model between species richness and average range size for each plot along elevation ####
# compare to linear model between FRic and average range size for each plot along elevation
# Hypothesis 3
species.richness = rowSums(pa.cdm)

rich.model = lm(species.richness ~ aver.range.plot)
# coef = -0.06 p = 0.012006, adjusted R2 = 0.328

Fric.model = lm(hull.matrix.all ~ aver.range.plot)
# coef = -0.06 p = 0.02517, adjusted R2 = 0.2604

#### FRic across elevation for each quantile ####

hull.matrix.minus.q1=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
  q.hull=convhulln(pca.scores[names(q1.cdm.pa[i,q1.cdm.pa[i,]>0]),], options="FA")
  q.hull.vol=q.hull$vol
  hull.matrix.minus.q1[i,1]=q.hull.vol
}

# 550 - 1000 doesn't work

com.550.names=names(q1.cdm.pa[12, q1.cdm.pa[12,]>0])
hull.matrix.minus.q1[16,1]=prod(apply(pca.scores[com.550.names,], MARGIN=2,max)-
       apply(pca.scores[com.550.names,], MARGIN=2, min))

com.750.names=names(q1.cdm.pa[11, q1.cdm.pa[11,]>0])

hull.matrix.minus.q1[11,1]=0.53247270--0.33774486
hull.matrix.minus.q1[13,1]=0
hull.matrix.minus.q1[14,1]=0.53247270--0.33774486
hull.matrix.minus.q1[15,1]=0.12466779--3.15839779
hull.matrix.minus.q1[16,1]=0.12466779--3.15839779

# write.csv(hull.matrix.minus.q1, file="./Data/gbif.outputs/hull.matrix.minus.q1.gbif.csv")

hull.matrix.minus.q2=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
  q.hull=convhulln(pca.scores[names(q2.cdm.pa[i,q2.cdm.pa[i,]>0]),], options="FA")
  q.hull.vol=q.hull$vol
  hull.matrix.minus.q2[i,1]=q.hull.vol
}

q.hull=convhulln(pca.scores[names(q2.cdm.pa[15,q2.cdm.pa[15,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.minus.q2[15,1]=q.hull.vol

com.550.names=names(q2.cdm.pa[16, q2.cdm.pa[16,]>0])
hull.matrix.minus.q2[16,1]=prod(apply(pca.scores[com.550.names,], MARGIN=2,max)-
                                  apply(pca.scores[com.550.names,], MARGIN=2, min))

hull.matrix.minus.q2[15,1]=1.41993368--0.38969741

write.csv(hull.matrix.minus.q2, file="./Data/gbif.outputs/hull.matrix.minus.q2.gbif.csv")

hull.matrix.minus.q3=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
  q.hull=convhulln(pca.scores[names(q3.cdm.pa[i,q3.cdm.pa[i,]>0]),], options="FA")
  q.hull.vol=q.hull$vol
  hull.matrix.minus.q3[i,1]=q.hull.vol
}

q.hull=convhulln(pca.scores[names(q3.cdm.pa[16,q3.cdm.pa[16,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.minus.q3[16,1]=q.hull.vol

com.550.names=names(q3.cdm.pa[15, q3.cdm.pa[15,]>0])
hull.matrix.minus.q3[15,1]=prod(apply(pca.scores[com.550.names,], MARGIN=2,max)-
                                  apply(pca.scores[com.550.names,], MARGIN=2, min))


write.csv(hull.matrix.minus.q3, file="./Data/gbif.outputs/hull.matrix.minus.q3.gbif.csv")

hull.matrix.minus.q4=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
  q.hull=convhulln(pca.scores[names(q4.cdm.pa[i,q4.cdm.pa[i,]>0]),], options="FA")
  q.hull.vol=q.hull$vol
  hull.matrix.minus.q4[i,1]=q.hull.vol
}

q.hull=convhulln(pca.scores[names(q4.cdm.pa[15,q4.cdm.pa[15,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.minus.q4[15,1]=q.hull.vol

com.550.names=names(q4.cdm.pa[16, q4.cdm.pa[16,]>0])
hull.matrix.minus.q4[16,1]=prod(apply(pca.scores[com.550.names,], MARGIN=2,max)-
                                  apply(pca.scores[com.550.names,], MARGIN=2, min))


write.csv(hull.matrix.minus.q4, file="./Data/gbif.outputs/hull.matrix.minus.q4.gbif.csv")

hull.matrix.minus.q5=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
  q.hull=convhulln(pca.scores[names(q5.cdm.pa[i,q5.cdm.pa[i,]>0]),], options="FA")
  q.hull.vol=q.hull$vol
  hull.matrix.minus.q5[i,1]=q.hull.vol
}

write.csv(hull.matrix.minus.q5, file="./Data/gbif.outputs/hull.matrix.minus.q5.gbif.csv")

#### Testing for Mid-Domain Effect ####
# Have to use archived version of package

pa.cdm=read.csv("pa.cdm.csv", header=T, row.names=1)

mde.output=matrix(data=NA, ncol = 999, nrow = 16)

for(i in 1:999){
  mde.output.rand=rangemod1d(pa.cdm, rsize= "observed", reps = 1)
  mde.output[,i]=mde.output.rand[,1]
}

ses.all=(total.rich - apply(mde.output, MARGIN = 1, mean))/apply(mde.output, MARGIN = 1, sd)

p.val.all=apply(cbind(total.rich, mde.output), MARGIN = 1, rank)[1,]/1000

rangemod.data=rangemod1d(pa.cdm, rsize= "observed", reps = 1000, degen = TRUE)

# Presence of MDE was checked by plotting 95% CI of the null model (1000 replicates) against observed speices richness
# Statistical sig. checked using spearman correlation in R

plot(total.rich~elevation, xlab="Elevation (m)", ylab="Species Richness", pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5)
par(new=T)
plot(rangemod.data$out.df$mod.rich~elevation, type="l",ylab="", xlab="", yaxt="n", xaxt="n")
par(new=T)
plot(rangemod.data$out.df$q2.5~elevation, type="l",ylab="", xlab="", yaxt="n", xaxt="n")
par(new=T)
plot(rangemod.data$out.df$q97.5~elevation, type="l",ylab="", xlab="", yaxt="n", xaxt="n")

test.table=cbind(total.rich,rangemod.data$out.df$mod.rich)
test.table=as.data.frame(test.table)
colnames(test.table)[1]="obs.rich"
colnames(test.table)[2]="exp.rich"
row.names(test.table)=elevation

chisq.test(test.table)
# X = 22.294, df = 15, p-value=0.1003

cor.test(test.table$obs.rich, test.table$exp.rich)
# r = 0.25, p = 0.34

test=lm(test.table$obs.rich~test.table$exp.rich)
# slope = 0.39, p = 0.317, adjusted R2 = 0.01

# Observed sp. richness distribution is not different from expected sp. richness under MDE

plot(total.rich,rangemod.data$out.df$mod.rich)
abline(lm(total.rich~rangemod.data$out.df$mod.rich))

test=lm(total.rich~rangemod.data$out.df$mod.rich)















#### Revision Figures ####

library(ggplot2)
library(cowplot)

# Hypothesis 1A
# correlate range size of each species with elevation of elevation

cor.test(all.unique.sonadora$gbif.range, all.unique.sonadora$plotElevation)
# r = 0.2834519, p < 0.0005

test = lm(all.unique.sonadora$gbif.range~all.unique.sonadora$plotElevation)

plot.1 = ggplot(all.unique.sonadora, aes(x = plotElevation, y = gbif.range)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE,color = "#F6955E")+
  theme_classic(base_size = 15)+
  labs(x = "Elevation (m)", y = "Elevational Range Size (m)")

# Rapoport's rule: 

cor.test(aver.range.plot.unique, elevation)
# r = 0.938, p < 0.0001

rap.df = as.data.frame(cbind(aver.range.plot.unique,elevation))

plot.2 = ggplot(rap.df, aes(x = elevation, y = aver.range.plot.unique)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "#F6955E")+
  theme_classic(base_size = 15)+
  labs(x = "Elevation (m)", y = "Average Range Size per Plot (m)")

# Hypothesis 1B

cor.test(rowSums(q1.cdm.pa), elevation)
# r = -0.84, p = 0.00005338 # gbif
cor.test(rowSums(q2.cdm.pa), elevation)
# r = -0.50, p = 0.04668 # gbif
cor.test(rowSums(q3.cdm.pa), elevation)
# r = -0.52, p = 0.03527 # gbif
cor.test(rowSums(q4.cdm.pa), elevation)
# r = -0.79, p = 0.0002367 # gbif
cor.test(rowSums(q5.cdm.pa), elevation)
# r = 0.60, p = 0.01 # gbif

# Hypothesis 2

cor.test(elevation, hull.matrix.data$Allsp.pca)
# r = -0.66, p = 0.006

fric.df = as.data.frame(cbind(hull.matrix.data$Allsp.pca,elevation))

plot.3 = ggplot(fric.df, aes(x = elevation, y = V1)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "#F6955E")+
  theme_classic(base_size = 15)+
  labs(x = "Elevation (m)", y = "Functional Richness")

Figure1 = plot_grid(plot.1,plot.2,plot.3, labels = c("a","b","c"), align = "hv")

ggsave("./Data/gbif.outputs/Figure1.pdf",height=10, width = 12)


cor.test(elevation, hull.matrix.data.2$Fric.minus.Q1)
# r = -0.55, p = 0.02695
cor.test(elevation, hull.matrix.data.2$Fric.minus.Q2)
# r = 0.19, p = 0.4663
cor.test(elevation, hull.matrix.data.2$Fric.minus.Q3)
# r = -0.09, p = 0.7278
cor.test(elevation, hull.matrix.data.2$Fric.minus.Q4)
# r = -0.62, p = 0.01034
cor.test(elevation, hull.matrix.data.2$Fric.minus.Q5)
# r = 0.79, p = 0.0002569

# Hypothesis 3

rich.model = lm(species.richness ~ aver.range.plot.unique)
# coef = -0.119 p = 0.012006, adjusted R2 = 0.6827

Fric.model = lm(hull.matrix.all ~ aver.range.plot.unique)
# coef = -0.11 p = 0.000714, adjusted R2 = 0.54

rich.predict = predict(rich.model)
H3.df = as.data.frame(cbind(rich.predict,aver.range.plot.unique))
fric.predict = predict(Fric.model)
H3.df$fric.predict = fric.predict
H3.df$sp.richness = species.richness
H3.df$fric.obs = hull.matrix.all

rich.plot = ggplot(H3.df, aes(x = aver.range.plot.unique, y = species.richness)) +
  geom_point()+
  geom_line(aes(x = aver.range.plot.unique, y = rich.predict),color = "#F6955E")+
  theme_classic(base_size = 15)+
  scale_x_continuous(limits = c(900,1150))+
  labs(x = "Average Range Size \n per plot across elevation (m)", y = "Species Richness")

fric.plot = ggplot(H3.df, aes(x = aver.range.plot.unique, y = fric.obs)) +
  geom_point()+
  geom_line(aes(x = aver.range.plot.unique, y = fric.predict),color = "#F6955E")+
  theme_classic(base_size = 15)+
  scale_x_continuous(limits = c(900,1150))+
  labs(x = "Average Range Size \n per plot across elevation (m)", y = "Functional Richness")


Figure2 = plot_grid(rich.plot,fric.plot,labels = c("a","b"), align = "hv")

ggsave("./Data/gbif.outputs/Figure2.pdf",height=8, width = 12)

