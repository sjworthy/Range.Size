# Differences in patterns of species and functional richness for narrow and wide-ranging species along an elevation gradient

setwd("~/Documents/Range.size/Range.size/Data")

#########################################################################################
#########################################################################################

all.data=read.csv("all.data.csv", header = T, row.names = 1)

# Species Codes not included: ALOBRY, CASDEC, CYAARB,GESSIN, HELCAR,PHYTOL,THECAC

# Putting Species into Range Size Quantiles

### Divide the range sizes from 0 to 1 by .2 or 20%

quantile(all.data$range, seq(0, 1, by=.2))


### all.data is file final.data.csv
### Determine the number of species in each quantile
### The range given is slightly adjusted from that given by the quantile split of the data so that the data is more evenly distributed into quantiles.

sum(all.data$range<=505)
[1] 23
sum(all.data$range>505 & all.data$range<=749)
[1] 21
sum(all.data$range>749 & all.data$range<=846)
[1] 23
sum(all.data$range>846 & all.data$range<=950) # should be <=949 but = 18 so more uneven distribution
[1] 24
sum(all.data$range>950)
[1] 21

### If the value in column range of all.data meets the criteria stated, place the quantile name in column Quantile

all.data$Quantile[all.data$range<=505]="Q1"
all.data$Quantile[all.data$range>505 & all.data$range<=749]="Q2"
all.data$Quantile[all.data$range>749 & all.data$range<=846]="Q3"
all.data$Quantile[all.data$range>846 & all.data$range<=950]="Q4"
all.data$Quantile[all.data$range>950]="Q5"

### If the value in column range of sonadora.csv meets the criteria stated, place the quantile name in column Quantile
### Get a list of the species that are in each quantile
### Add Q1-Q5 for each of these species to the Quantile column

sonadora=read.csv("sonadora.csv", header=T)

Q1=as.character(all.data$SpeciesCode[all.data$range<=505])
Q2=as.character(all.data$SpeciesCode[all.data$range>505 &
all.data$range<=749])
Q3=as.character(all.data$SpeciesCode[all.data$range>749 &
all.data$range<=846])
Q4=as.character(all.data$SpeciesCode[all.data$range>846 &
all.data$range<=950])
Q5=as.character(all.data$SpeciesCode[all.data$range>950])


sonadora$Quantile[sonadora$stemSpeciesCode %in% Q1]="Q1"
sonadora$Quantile[sonadora$stemSpeciesCode %in% Q2]="Q2"
sonadora$Quantile[sonadora$stemSpeciesCode %in% Q3]="Q3"
sonadora$Quantile[sonadora$stemSpeciesCode %in% Q4]="Q4"
sonadora$Quantile[sonadora$stemSpeciesCode %in% Q5]="Q5"

# Percentage of abundance in each quantile

table(sonadora$Quantile)
# Q1
799/8693 = 9.19%
# Q2
1614/8693 = 18.57%
# Q3
3411/8693 = 39.24%
# Q4
777/8693 = 8.94%
# Q5
2092/8693 = 24.07%

#########################################################################################
#########################################################################################
# Making a Community Data Matrix

### Make a community data matrix with number of species in each elevation plot
### Make a blank community data matrix with 111 columns, one for each species, and 16 rows, one for each plot in the transect (250m-1000m)

library(vegan)

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

abund.cdm=read.csv("abund.cdm.csv",header=T,row.names = 1)

#########################################################################################
### Make a community data matrix with just narrow species and one with just wide-ranging species
### Species with range < mean (758.4732) are narrow with species > mean range being wide-ranging.

narrow.species=as.character(all.data$SpeciesCode[all.data$range< 758.4732])
#50 narrow species 
wide.species=as.character(all.data$SpeciesCode[all.data$range> 758.4732])
# 62 wide-ranging species 

narrow.cdm=abund.cdm[,narrow.species]
wide.cdm=abund.cdm[,wide.species]

write.csv(narrow.cdm, file="narrow.cdm.csv")
write.csv(wide.cdm, file="wide.cdm.csv")


#########################################################################################
### Make a community data matrix where quantiles are the rows/communities
### Make a blank community data matrix with 112 columns, one for each species, and 5 rows, one for each quantile.

quant.cdm=matrix(data=NA, nrow=5, ncol=112)
colnames(q.cdm)=all.data$SpeciesCode
rownames(q.cdm)=c("Q1", "Q2", "Q3", "Q4", "Q5")

quant.pa.cdm=decostand(quant.cdm, "pa")

### Fill in the community data matrix
### Subset the abundance data set by the quantile and take that subsetted object and table the SpeciesCode
# Enter the values manually because it won't auto place due to missing columns.

quantile=subset(sonadora, sonadora$Quantile=="Q5")
output=table(quantile$stemSpeciesCode)
output

quant.cdm=read.csv("quant.cdm.csv",header = T, row.names = 1)


#########################################################################################
### Abundance correlation with range

cor.test(all.data$range, all.data$abundance)
# r = 0.09, p = 0.3714

narrow.subset=all.data[all.data$Range.Size=="Narrow",]
wide.subset=all.data[all.data$Range.Size=="Wide",]

cor.test(narrow.subset$range, narrow.subset$abundance)
# r = 0.23, p = 0.1124
cor.test(wide.subset$range, wide.subset$abundance)
# r = 0.04, p = 0.7308

#########################################################################################
### Make community data matrices for each quantile separately.
### Get a list of the species in each quantile


Q1=as.character(all.data$SpeciesCode[all.data$Quantile=="Q1"])
Q2=as.character(all.data$SpeciesCode[all.data$Quantile=="Q2"])
Q3=as.character(all.data$SpeciesCode[all.data$Quantile=="Q3"])
Q4=as.character(all.data$SpeciesCode[all.data$Quantile=="Q4"])
Q5=as.character(all.data$SpeciesCode[all.data$Quantile=="Q5"])

### Make empty cdms
q1.cdm=matrix(data=NA, nrow=16, ncol=23)
q2.cdm=matrix(data=NA, nrow=16, ncol=21)
q3.cdm=matrix(data=NA, nrow=16, ncol=23)
q4.cdm=matrix(data=NA, nrow=16, ncol=24)
q5.cdm=matrix(data=NA, nrow=16, ncol=21)


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

Q1.subset=sonadora[sonadora$Quantile=="Q1",]
Q2.subset=sonadora[sonadora$Quantile=="Q2",]
Q3.subset=sonadora[sonadora$Quantile=="Q3",]
Q4.subset=sonadora[sonadora$Quantile=="Q4",]
Q5.subset=sonadora[sonadora$Quantile=="Q5",]

write.csv(Q1.subset,file="Q1.subset.csv")
write.csv(Q2.subset,file="Q2.subset.csv")
write.csv(Q3.subset,file="Q3.subset.csv")
write.csv(Q4.subset,file="Q4.subset.csv")
write.csv(Q5.subset,file="Q5.subset.csv")

Q1.df=read.csv("Q1.subset.csv",header=T,row.names=1)
Q2.df=read.csv("Q2.subset.csv",header=T,row.names=1)
Q3.df=read.csv("Q3.subset.csv",header=T,row.names=1)
Q4.df=read.csv("Q4.subset.csv",header=T,row.names=1)
Q5.df=read.csv("Q5.subset.csv",header=T,row.names=1)

### Subset the quantile specific abundance data frames by the elevation and take that subsetted object and table the SpeciesCode
# Enter the values manually because it won't auto place due to missing columns.

elevation=subset(Q1.df, Q1.df$plotElevation==250)
output=table(elevation$stemSpeciesCode)
output

elevation=subset(Q2.df, Q2.df$plotElevation==250)
output=table(elevation$stemSpeciesCode)
output

elevation=subset(Q3.df, Q3.df$plotElevation==1000)
output=table(elevation$stemSpeciesCode)
output

elevation=subset(Q4.df, Q4.df$plotElevation==250)
output=table(elevation$stemSpeciesCode)
output

elevation=subset(Q5.df, Q5.df$plotElevation==250)
output=table(elevation$stemSpeciesCode)
output

q1.cdm=read.csv("q1.cdm.csv", header=T, row.names=1)
q2.cdm=read.csv("q2.cdm.csv", header=T, row.names=1)
q3.cdm=read.csv("q3.cdm.csv", header=T, row.names=1)
q4.cdm=read.csv("q4.cdm.csv", header=T, row.names=1)
q5.cdm=read.csv("q5.cdm.csv", header=T, row.names=1)


#########################################################################################
### Make community data matrix where communities are quantiles of each plot and columns are species in each plot
#Make empty cdms
q250=matrix(data=NA, nrow=5, ncol=31)
q300=matrix(data=NA, nrow=5, ncol=50)
q350=matrix(data=NA, nrow=5, ncol=42)
q400=matrix(data=NA, nrow=5, ncol=46)
q450=matrix(data=NA, nrow=5, ncol=39)
q500=matrix(data=NA, nrow=5, ncol=44)
q550=matrix(data=NA, nrow=5, ncol=27)
q600=matrix(data=NA, nrow=5, ncol=30)
q650=matrix(data=NA, nrow=5, ncol=38)
q700=matrix(data=NA, nrow=5, ncol=37)
q750=matrix(data=NA, nrow=5, ncol=30)
q800=matrix(data=NA, nrow=5, ncol=30)
q850=matrix(data=NA, nrow=5, ncol=25)
q900=matrix(data=NA, nrow=5, ncol=26)
q950=matrix(data=NA, nrow=5, ncol=17)
q1000=matrix(data=NA, nrow=5, ncol=20)

rownames(q250)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q300)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q350)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q400)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q450)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q500)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q550)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q600)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q650)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q700)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q750)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q800)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q850)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q900)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q950)=c("Q1", "Q2", "Q3", "Q4", "Q5")
rownames(q1000)=c("Q1", "Q2", "Q3", "Q4", "Q5")

elevation250=subset(sonadora, sonadora$plotElevation=="250")
elevation300=subset(sonadora, sonadora$plotElevation=="300")
elevation350=subset(sonadora, sonadora$plotElevation=="350")
elevation400=subset(sonadora, sonadora$plotElevation=="400")
elevation450=subset(sonadora, sonadora$plotElevation=="450")
elevation500=subset(sonadora, sonadora$plotElevation=="500")
elevation550=subset(sonadora, sonadora$plotElevation=="550")
elevation600=subset(sonadora, sonadora$plotElevation=="600")
elevation650=subset(sonadora, sonadora$plotElevation=="650")
elevation700=subset(sonadora, sonadora$plotElevation=="700")
elevation750=subset(sonadora, sonadora$plotElevation=="750")
elevation800=subset(sonadora, sonadora$plotElevation=="800")
elevation850=subset(sonadora, sonadora$plotElevation=="850")
elevation900=subset(sonadora, sonadora$plotElevation=="900")
elevation950=subset(sonadora, sonadora$plotElevation=="950")
elevation1000=subset(sonadora, sonadora$plotElevation=="1000")

output1=unique(elevation250$stemSpeciesCode)
output2=unique(elevation300$stemSpeciesCode)
output3=unique(elevation350$stemSpeciesCode)
output4=unique(elevation400$stemSpeciesCode)
output5=unique(elevation450$stemSpeciesCode)
output6=unique(elevation500$stemSpeciesCode)
output7=unique(elevation550$stemSpeciesCode)
output8=unique(elevation600$stemSpeciesCode)
output9=unique(elevation650$stemSpeciesCode)
output10=unique(elevation700$stemSpeciesCode)
output11=unique(elevation750$stemSpeciesCode)
output12=unique(elevation800$stemSpeciesCode)
output13=unique(elevation850$stemSpeciesCode)
output14=unique(elevation900$stemSpeciesCode)
output15=unique(elevation950$stemSpeciesCode)
output16=unique(elevation1000$stemSpeciesCode)

colnames(q250)=output1
colnames(q300)=output2
colnames(q350)=output3
colnames(q400)=output4
colnames(q450)=output5
colnames(q500)=output6
colnames(q550)=output7
colnames(q600)=output8
colnames(q650)=output9
colnames(q700)=output10
colnames(q750)=output11
colnames(q800)=output12
colnames(q850)=output13
colnames(q900)=output14
colnames(q950)=output15
colnames(q1000)=output16

write.csv(elevation250, file="elevation250.subset.csv")
write.csv(elevation300, file="elevation300.subset.csv")
write.csv(elevation350, file="elevation350.subset.csv")
write.csv(elevation400, file="elevation400.subset.csv")
write.csv(elevation450, file="elevation450.subset.csv")
write.csv(elevation500, file="elevation500.subset.csv")
write.csv(elevation550, file="elevation550.subset.csv")
write.csv(elevation600, file="elevation600.subset.csv")
write.csv(elevation650, file="elevation650.subset.csv")
write.csv(elevation700, file="elevation700.subset.csv")
write.csv(elevation750, file="elevation750.subset.csv")
write.csv(elevation800, file="elevation800.subset.csv")
write.csv(elevation850, file="elevation850.subset.csv")
write.csv(elevation900, file="elevation900.subset.csv")
write.csv(elevation950, file="elevation950.subset.csv")
write.csv(elevation1000, file="elevation1000.subset.csv")

### Subset the elevation specific abundance data frames by the quantile and take that subsetted object and table the SpeciesCode
# Enter the values manually because it won't auto place due to missing columns.

quantile=subset(elevation1000, elevation1000$Quantile=="Q3")
output=table(quantile$stemSpeciesCode)
output

q250=read.csv("q250.cdm.csv", header=T, row.names=1)
q300=read.csv("q300.cdm.csv", header=T, row.names=1)
q350=read.csv("q350.cdm.csv", header=T, row.names=1)
q400=read.csv("q400.cdm.csv", header=T, row.names=1)
q450=read.csv("q450.cdm.csv", header=T, row.names=1)
q500=read.csv("q500.cdm.csv", header=T, row.names=1)
q550=read.csv("q550.cdm.csv", header=T, row.names=1)
q600=read.csv("q600.cdm.csv", header=T, row.names=1)
q650=read.csv("q650.cdm.csv", header=T, row.names=1)
q700=read.csv("q700.cdm.csv", header=T, row.names=1)
q750=read.csv("q750.cdm.csv", header=T, row.names=1)
q800=read.csv("q800.cdm.csv", header=T, row.names=1)
q850=read.csv("q850.cdm.csv", header=T, row.names=1)
q900=read.csv("q900.cdm.csv", header=T, row.names=1)
q950=read.csv("q950.cdm.csv", header=T, row.names=1)
q1000=read.csv("q1000.cdm.csv", header=T, row.names=1)

#########################################################################################
### Presence/Absence Community data matrices for each quantile
library(vegan)

q5.cdm.pa=decostand(q5.cdm, "pa")
q4.cdm.pa=decostand(q4.cdm, "pa")
q3.cdm.pa=decostand(q3.cdm, "pa")
q2.cdm.pa=decostand(q2.cdm, "pa")
q1.cdm.pa=decostand(q1.cdm, "pa")

write.csv(q1.cdm.pa, file="q1.cdm.pa.csv")
write.csv(q2.cdm.pa, file="q2.cdm.pa.csv")
write.csv(q3.cdm.pa, file="q3.cdm.pa.csv")
write.csv(q4.cdm.pa, file="q4.cdm.pa.csv")
write.csv(q5.cdm.pa, file="q5.cdm.pa.csv")

#########################################################################################
#########################################################################################
# Making a trait data matrix
### Make a trait data matrix where the number of rows (112) is equal to the number of species and the number of columns (10) is equal to the number of traits.
## Only 6 traits for all species = height.ft, la, sla, wood, n, p

tdm=matrix(data=NA, nrow=112, ncol=10)
colnames(tdm)=c("height.ft", "la", "sla", "wood", "n", "p", "seed.mass", "c", "delta13c", "vla")
rownames(tdm)=all.data$SpeciesCode

### Cells of trait data matrix contain trait values for species

tdm[,1]=all.data$height_ft
tdm[,2]=all.data$la
tdm[,3]=all.data$sla
tdm[,4]=all.data$wood
tdm[,5]=all.data$n
tdm[,6]=all.data$p
tdm[,7]=all.data$seed.mass
tdm[,8]=all.data$c
tdm[,9]=all.data$delta13c
tdm[,10]=all.data$vla

write.csv(tdm, file="tdm.csv")

# Test for normality of traits

shapiro.test(tdm[,10])

# Transform trait data to approximate a normal distribution
# Traits that have normal distribution = c, delta13c, vla

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

write.csv(pca.scores, file="pca.scores.csv")
write.csv(pca.scores.all, file="pca.scores.all.csv")


### Make trait data matrix for each quantile
### Make character string for species within each quantile

tdm.scaled=read.csv("tdm.scaled.csv", header=T, row.names = 1)

Q1=as.character(all.data$SpeciesCode[all.data$range<=505])
Q2=as.character(all.data$SpeciesCode[all.data$range>505 & 
                                       all.data$range<=749])
Q3=as.character(all.data$SpeciesCode[all.data$range>749 & 
                                       all.data$range<=846])
Q4=as.character(all.data$SpeciesCode[all.data$range>846 & 
                                       all.data$range<=950])
Q5=as.character(all.data$SpeciesCode[all.data$range>950])

### Subset original trait data matrix by the species in each quantile character string

Q1.traits=tdm.scaled[Q1,]
Q2.traits=tdm.scaled[Q2,]
Q3.traits=tdm.scaled[Q3,]
Q4.traits=tdm.scaled[Q4,]
Q5.traits=tdm.scaled[Q5,]

write.csv(Q1.traits, file="q1.tdm.csv")
write.csv(Q2.traits, file="q2.tdm.csv")
write.csv(Q3.traits, file="q3.tdm.csv")
write.csv(Q4.traits, file="q4.tdm.csv")
write.csv(Q5.traits, file="q5.tdm.csv")

# Calculating Functional Richness
### Because dbFD function reduces the data set's dimensionality automatically, it uses a different number of dimensions for each community. To keep number of dimensions constant, must determine convex hull volume manually for each community.
install.packages("geometry")
library(geometry)

### FRic for each elevation plot for entire data set.
### pca only includes 1st three axes.

hull.matrix.all=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
q.hull=convhulln(pca.scores[names(abund.cdm[i,abund.cdm[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.all[i,1]=q.hull.vol
}

write.csv(hull.matrix.all, file="hull.matrix.pca.csv")

# FRic based on p/a cdm = same values as abund.cdm

hull.matrix.pa.all=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
  q.hull=convhulln(pca.scores[names(pa.cdm[i,pa.cdm[i,]>0]),], options="FA")
  q.hull.vol=q.hull$vol
  hull.matrix.pa.all[i,1]=q.hull.vol
}

### FRic for each elevation plot excluding wide ranging species
### pca only includes 1st three axes of narrow ranging species

names=as.factor(colnames(narrow.cdm))
narrow.pca=subset(pca.scores, rownames(pca.scores)%in%names)
narrow.pca=narrow.pca[,1:3]

hull.matrix.narrow=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
q.hull=convhulln(narrow.pca[names(narrow.cdm[i,narrow.cdm[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.narrow[i,1]=q.hull.vol
}

write.csv(hull.matrix.narrow, file="hull.matrix.narrow.pca.csv")

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

write.csv(hull.matrix.wide, file="hull.matrix.wide.pca.csv")

### Abundance of narrow and wide ranging individuals in each elevation

rowSums(narrow.cdm)
rowSums(wide.cdm)

narrow.cdm.pa=decostand(narrow.cdm, "pa")
wide.cdm.pa=decostand(wide.cdm, "pa")

rowSums(narrow.cdm.pa)
rowSums(wide.cdm.pa)

write.csv(narrow.cdm.pa, file="narrow.cdm.pa.csv")
write.csv(wide.cdm.pa, file="wide.cdm.pa.csv")

### Determine the average range size for each plot

sonadora=read.csv("sonadora.csv", header=T, row.names = 1)

plot250=subset(sonadora, sonadora$plotElevation==250)
plot300=subset(sonadora, sonadora$plotElevation==300)
plot350=subset(sonadora, sonadora$plotElevation==350)
plot400=subset(sonadora, sonadora$plotElevation==400)
plot450=subset(sonadora, sonadora$plotElevation==450)
plot500=subset(sonadora, sonadora$plotElevation==500)
plot550=subset(sonadora, sonadora$plotElevation==550)
plot600=subset(sonadora, sonadora$plotElevation==600)
plot650=subset(sonadora, sonadora$plotElevation==650)
plot700=subset(sonadora, sonadora$plotElevation==700)
plot750=subset(sonadora, sonadora$plotElevation==750)
plot800=subset(sonadora, sonadora$plotElevation==800)
plot850=subset(sonadora, sonadora$plotElevation==850)
plot900=subset(sonadora, sonadora$plotElevation==900)
plot950=subset(sonadora, sonadora$plotElevation==950)
plot1000=subset(sonadora, sonadora$plotElevation==1000)

mean(plot250$Range) = 683.9197
mean(plot300$Range) = 780.9163
mean(plot350$Range) = 860.7509
mean(plot400$Range) = 784.538
mean(plot450$Range) = 812.1144
mean(plot500$Range) = 822.6497
mean(plot550$Range) = 740.1176
mean(plot600$Range) = 791.0245
mean(plot650$Range) = 728.8867
mean(plot700$Range) = 805.9375
mean(plot750$Range) = 853.6667
mean(plot800$Range) = 850.832
mean(plot850$Range) = 725.1356
mean(plot900$Range) = 823.2753
mean(plot950$Range) = 950.8442
mean(plot1000$Range) = 823.9358

# correlation between average range size per plot and elevation
aver.range.plot=c(683.9197,780.9163,860.7509,784.538,812.1144,822.6497,740.1176,791.0245,
                  728.8867,805.9375,853.6667,850.832,725.1356,823.2753,950.8442,823.9358)
cor.test(aver.range.plot, elevation)
# r = 0.43, p = 0.10

#########################################################################################
#########################################################################################
### Shifting Ranges
# shift ranges in 50 m increments

og.shift.range=all.data[,c(1, 5,6)]
shift.range.50=og.shift.range
shift.range.50[,2:3]=shift.range.50[,2:3]+50
shift.range.100=shift.range.50
shift.range.100[,2:3]=shift.range.100[,2:3]+50
shift.range.150=shift.range.100
shift.range.150[,2:3]=shift.range.150[,2:3]+50
shift.range.200=shift.range.150
shift.range.200[,2:3]=shift.range.200[,2:3]+50
shift.range.250=shift.range.200
shift.range.250[,2:3]=shift.range.250[,2:3]+50
shift.range.300=shift.range.250
shift.range.300[,2:3]=shift.range.300[,2:3]+50
shift.range.350=shift.range.300
shift.range.350[,2:3]=shift.range.350[,2:3]+50
shift.range.400=shift.range.350
shift.range.400[,2:3]=shift.range.400[,2:3]+50
shift.range.450=shift.range.400
shift.range.450[,2:3]=shift.range.450[,2:3]+50
shift.range.500=shift.range.450
shift.range.500[,2:3]=shift.range.500[,2:3]+50
shift.range.550=shift.range.500
shift.range.550[,2:3]=shift.range.550[,2:3]+50

write.csv(shift.range.1, file="shift.range.50.csv")

### Need cdm which changes every increase in 50 m

all.data=read.csv("all.data.csv", header=T, row.names = 1)

cdm.shift.range=matrix(data=NA, ncol=112, nrow=16)
colnames(cdm.shift.range)=all.data$SpeciesCode
rownames(cdm.shift.range)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
700, 750, 800, 850, 900, 950, 1000)

### Fill in the community data matrix
### Make a list object containing the elevation of each plot

listy=c(250, 300, 350, 400, 450, 500, 550, 600, 650, 700,750,
800, 850, 900, 950, 1000)

### For loop that loops through each element in listy, each elevation. Then, subset the abundance data set by the plot elevation and take that subsetted object and table the SpeciesCode column and take that object and place those values into row i of cdm object.

for(j in 1:112){
for(i in 1:length(listy)){
cdm.shift.range[i,j]=ifelse(og.shift.range[j,2] <= listy[i] & og.shift.range[j,3] >= listy[i], 1, 0) 
}
}

write.csv(cdm.shift.range, file="cdm.shift.range.550.csv")

### Calculate FRic

cdm.shift.range.og=read.csv("cdm.shift.range.og.csv", header=T, row.names = 1)
cdm.shift.range.50=read.csv("cdm.shift.range.50.csv", header=T, row.names = 1)
cdm.shift.range.100=read.csv("cdm.shift.range.100.csv", header=T, row.names = 1)
cdm.shift.range.150=read.csv("cdm.shift.range.150.csv", header=T, row.names = 1)
cdm.shift.range.200=read.csv("cdm.shift.range.200.csv", header=T, row.names = 1)
cdm.shift.range.250=read.csv("cdm.shift.range.250.csv", header=T, row.names = 1)
cdm.shift.range.300=read.csv("cdm.shift.range.300.csv", header=T, row.names = 1)
cdm.shift.range.350=read.csv("cdm.shift.range.350.csv", header=T, row.names = 1)
cdm.shift.range.400=read.csv("cdm.shift.range.400.csv", header=T, row.names = 1)
cdm.shift.range.450=read.csv("cdm.shift.range.450.csv", header=T, row.names = 1)
cdm.shift.range.500=read.csv("cdm.shift.range.500.csv", header=T, row.names = 1)
cdm.shift.range.550=read.csv("cdm.shift.range.550.csv", header=T, row.names = 1)


hull.matrix=matrix(data=NA, nrow=16,ncol=1)

for(i in 1:16){
q.hull=convhulln(pca.scores[names(cdm.shift.range.250[i,cdm.shift.range.250[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix[i,1]=q.hull.vol
}

hull.matrix=as.data.frame(hull.matrix)
hull.matrix[,2]=rowSums(cdm.shift.range.250)
write.csv(hull.matrix, file="hull.matrix.shift.range.250.csv")

hull.matrix=matrix(data=NA, nrow=16,ncol=1)
for(i in 2:16){
q.hull=convhulln(pca.scores[names(cdm.shift.range.300[i,cdm.shift.range.300[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix[i,1]=q.hull.vol
}

hull.matrix=as.data.frame(hull.matrix)
hull.matrix[,2]=rowSums(cdm.shift.range.300)
write.csv(hull.matrix, file="hull.matrix.shift.range.300.csv")

hull.matrix=matrix(data=NA, nrow=16,ncol=1)
for(i in 3:16){
q.hull=convhulln(pca.scores[names(cdm.shift.range.350[i,cdm.shift.range.350[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix[i,1]=q.hull.vol
}

hull.matrix=as.data.frame(hull.matrix)
hull.matrix[,2]=rowSums(cdm.shift.range.350)
write.csv(hull.matrix, file="hull.matrix.shift.range.350.csv")

hull.matrix=matrix(data=NA, nrow=16,ncol=1)
for(i in 4:16){
q.hull=convhulln(pca.scores[names(cdm.shift.range.400[i,cdm.shift.range.400[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix[i,1]=q.hull.vol
}
hull.matrix=as.data.frame(hull.matrix)
hull.matrix[,2]=rowSums(cdm.shift.range.400)
write.csv(hull.matrix, file="hull.matrix.shift.range.400.csv")

hull.matrix=matrix(data=NA, nrow=16,ncol=1)
for(i in 5:16){
q.hull=convhulln(pca.scores[names(cdm.shift.range.450[i,cdm.shift.range.450[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix[i,1]=q.hull.vol
}

hull.matrix=as.data.frame(hull.matrix)
hull.matrix[,2]=rowSums(cdm.shift.range.450)
write.csv(hull.matrix, file="hull.matrix.shift.range.450.csv")

hull.matrix=matrix(data=NA, nrow=16,ncol=1)
for(i in 6:16){
q.hull=convhulln(pca.scores[names(cdm.shift.range.500[i,cdm.shift.range.500[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix[i,1]=q.hull.vol
}

hull.matrix=as.data.frame(hull.matrix)
hull.matrix[,2]=rowSums(cdm.shift.range.500)
write.csv(hull.matrix, file="hull.matrix.shift.range.500.csv")

hull.matrix=matrix(data=NA, nrow=16,ncol=1)
for(i in 7:16){
q.hull=convhulln(pca.scores[names(cdm.shift.range.550[i,cdm.shift.range.550[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix[i,1]=q.hull.vol
}

hull.matrix=as.data.frame(hull.matrix)
hull.matrix[,2]=rowSums(cdm.shift.range.550)
write.csv(hull.matrix, file="hull.matrix.shift.range.550.csv")

plot(1, type="b", xlab="Species Richness", ylab="FRic", xlim=rev(range(c(0,105))), ylim=c(0,85))

x=seq(0,110,by=10)
y=seq(0,110, by=10)
plot(x,y, xlim=(rev(range(x))))
plot(all.hull$R.250, all.hull$E.250, xlim=(rev(range(all.hull$R.250))), type="b")

### Fill in the community data matrix
### Make a list object containing the elevation of each plot

og.time=all.data[,c(1, 5,6)]
time.2030=og.time
time.2030[,2:3]=time.2030[,2:3]+250
time.2050=og.time
time.2050[,2:3]=og.time[,2:3]+333
time.2070=og.time
time.2070[,2:3]=og.time[,2:3]+417
time.2090=og.time
time.2090[,2:3]=og.time[,2:3]+500
time.2100=og.time
time.2100[,2:3]=og.time[,2:3]+667

### Fill in the community data matrix
### Make a list object containing the elevation of each plot

listy=c(250, 300, 350, 400, 450, 500, 550, 600, 650, 700,750,
        800, 850, 900, 950, 1000)

### For loop that loops through each element in listy, each elevation. Then, subset the abundance data set by the plot elevation and take that subsetted object and table the SpeciesCode column and take that object and place those values into row i of cdm object.

cdm.time.og=matrix(data=NA, ncol=112, nrow=16)
colnames(cdm.time.og)=all.data$SpeciesCode
rownames(cdm.time.og)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
                                 700, 750, 800, 850, 900, 950, 1000)


for(j in 1:112){
  for(i in 1:length(listy)){
    cdm.time.og[i,j]=ifelse(og.time[j,2] <= listy[i] & og.time[j,3] >= listy[i], 1, 0)
}
}

write.csv(cdm.time.og, file="cdm.time.og.csv")

### For loop that loops through each element in listy, each elevation. Then, subset the abundance data set by the plot elevation and take that subsetted object and table the SpeciesCode column and take that object and place those values into row i of cdm object.

cdm.time.2030=matrix(data=NA, ncol=112, nrow=16)
colnames(cdm.time.2030)=all.data$SpeciesCode
rownames(cdm.time.2030)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
                        700, 750, 800, 850, 900, 950, 1000)

for(j in 1:112){
for(i in 1:length(listy)){
cdm.time.2030[i,j]=ifelse(time.2030[j,2] <= listy[i] & time.2030[j,3] >= listy[i], 1, 0) 
}
}

write.csv(cdm.time.2030, file="cdm.time.2030.csv")


cdm.time.2050=matrix(data=NA, ncol=112, nrow=16)
colnames(cdm.time.2050)=all.data$SpeciesCode
rownames(cdm.time.2050)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
700, 750, 800, 850, 900, 950, 1000)


for(j in 1:112){
for(i in 1:length(listy)){
cdm.time.2050[i,j]=ifelse(time.2050[j,2] <= listy[i] & time.2050[j,3] >= listy[i], 1, 0) 
}
}

write.csv(cdm.time.2050, file="cdm.time.2050.csv")

cdm.time.2070=matrix(data=NA, ncol=112, nrow=16)
colnames(cdm.time.2070)=all.data$SpeciesCode
rownames(cdm.time.2070)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
700, 750, 800, 850, 900, 950, 1000)


for(j in 1:112){
for(i in 1:length(listy)){
cdm.time.2070[i,j]=ifelse(time.2070[j,2] <= listy[i] & time.2070[j,3] >= listy[i], 1, 0) 
}
}

write.csv(cdm.time.2070, file="cdm.time.2070.csv")

cdm.time.2090=matrix(data=NA, ncol=112, nrow=16)
colnames(cdm.time.2090)=all.data$SpeciesCode
rownames(cdm.time.2090)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
700, 750, 800, 850, 900, 950, 1000)


for(j in 1:112){
for(i in 1:length(listy)){
cdm.time.2090[i,j]=ifelse(time.2090[j,2] <= listy[i] & time.2090[j,3] >= listy[i], 1, 0) 
}
}

write.csv(cdm.time.2090, file="cdm.time.2090.csv")

cdm.time.2100=matrix(data=NA, ncol=112, nrow=16)
colnames(cdm.time.2100)=all.data$SpeciesCode
rownames(cdm.time.2100)=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
700, 750, 800, 850, 900, 950, 1000)


for(j in 1:112){
for(i in 1:length(listy)){
cdm.time.2100[i,j]=ifelse(time.2100[j,2] <= listy[i] & time.2100[j,3] >= listy[i], 1, 0) 
}
}

write.csv(cdm.time.2100, file="cdm.time.2100.csv")

hull.matrix.2030=matrix(data=NA, nrow=16,ncol=1)

for(i in 1:16){
q.hull.2030=convhulln(pca.scores[names(cdm.time.2030[i,cdm.time.2030[i,]>0]),], options="FA")
q.hull.vol=q.hull.2030$vol
hull.matrix.2030[i,1]=q.hull.vol
}

hull.matrix.2030=as.data.frame(hull.matrix.2030)
hull.matrix.2030[,2]=rowSums(cdm.time.2030)
write.csv(hull.matrix.2030, file="hull.matrix.2030.csv")


hull.matrix.2050=matrix(data=NA, nrow=16,ncol=1)

for(i in 3:16){
q.hull.2050=convhulln(pca.scores[names(cdm.time.2050[i,cdm.time.2050[i,]>0]),], options="FA")
q.hull.vol=q.hull.2050$vol
hull.matrix.2050[i,1]=q.hull.vol
}

hull.matrix.2050=as.data.frame(hull.matrix.2050)
hull.matrix.2050[,2]=rowSums(cdm.time.2050)
write.csv(hull.matrix.2050, file="hull.matrix.2050.csv")

hull.matrix.2070=matrix(data=NA, nrow=16,ncol=1)

for(i in 5:16){
q.hull.2070=convhulln(pca.scores[names(cdm.time.2070[i,cdm.time.2070[i,]>0]),], options="FA")
q.hull.vol=q.hull.2070$vol
hull.matrix.2070[i,1]=q.hull.vol
}

hull.matrix.2070=as.data.frame(hull.matrix.2070)
hull.matrix.2070[,2]=rowSums(cdm.time.2070)
write.csv(hull.matrix.2070, file="hull.matrix.2070.csv")

hull.matrix.2090=matrix(data=NA, nrow=16,ncol=1)

for(i in 6:16){
q.hull.2090=convhulln(pca.scores[names(cdm.time.2090[i,cdm.time.2090[i,]>0]),], options="FA")
q.hull.vol=q.hull.2090$vol
hull.matrix.2090[i,1]=q.hull.vol
}

hull.matrix.2090=as.data.frame(hull.matrix.2090)
hull.matrix.2090[,2]=rowSums(cdm.time.2090)
write.csv(hull.matrix.2090, file="hull.matrix.2090.csv")


hull.matrix.2100=matrix(data=NA, nrow=16,ncol=1)

for(i in 10:16){
q.hull.2100=convhulln(pca.scores[names(cdm.time.2100[i,cdm.time.2100[i,]>0]),], options="FA")
q.hull.vol=q.hull.2100$vol
hull.matrix.2100[i,1]=q.hull.vol
}

hull.matrix.2100=as.data.frame(hull.matrix.2100)
hull.matrix.2100[,2]=rowSums(cdm.time.2100)
write.csv(hull.matrix.2100, file="hull.matrix.2100.csv")


## Make cdm excluding species in each quantile one at a time

q1.sp=subset(all.data, all.data$Quantile == "Q1")
listy=q1.sp$SpeciesCode
q1.pa.cdm=pa.cdm[,-which(names(pa.cdm)%in% listy)]
hull.matrix.minus.q1=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
q.hull=convhulln(pca.scores[names(q1.pa.cdm[i,q1.pa.cdm[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.minus.q1[i,1]=q.hull.vol
}

write.csv(hull.matrix.minus.q1, file="hull.matrix.minus.q1.csv")


q2.sp=subset(all.data, all.data$Quantile == "Q2")
listy=q2.sp$SpeciesCode
q2.pa.cdm=pa.cdm[,-which(names(pa.cdm)%in% listy)]
hull.matrix.minus.q2=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
q.hull=convhulln(pca.scores[names(q2.pa.cdm[i,q2.pa.cdm[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.minus.q2[i,1]=q.hull.vol
}

write.csv(hull.matrix.minus.q2, file="hull.matrix.minus.q2.csv")

q3.sp=subset(all.data, all.data$Quantile == "Q3")
listy=q3.sp$SpeciesCode
q3.pa.cdm=pa.cdm[,-which(names(pa.cdm)%in% listy)]
hull.matrix.minus.q3=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
q.hull=convhulln(pca.scores[names(q3.pa.cdm[i,q3.pa.cdm[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.minus.q3[i,1]=q.hull.vol
}

write.csv(hull.matrix.minus.q3, file="hull.matrix.minus.q3.csv")


q4.sp=subset(all.data, all.data$Quantile == "Q4")
listy=q4.sp$SpeciesCode
q4.pa.cdm=pa.cdm[,-which(names(pa.cdm)%in% listy)]
hull.matrix.minus.q4=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
q.hull=convhulln(pca.scores[names(q4.pa.cdm[i,q4.pa.cdm[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.minus.q4[i,1]=q.hull.vol
}

write.csv(hull.matrix.minus.q4, file="hull.matrix.minus.q4.csv")

q5.sp=subset(all.data, all.data$Quantile == "Q5")
listy=q5.sp$SpeciesCode
q5.pa.cdm=pa.cdm[,-which(names(pa.cdm)%in% listy)]
hull.matrix.minus.q5=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
q.hull=convhulln(pca.scores[names(q5.pa.cdm[i,q5.pa.cdm[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.minus.q5[i,1]=q.hull.vol
}

write.csv(hull.matrix.minus.q5, file="hull.matrix.minus.q5.csv")

q1.2.sp=subset(all.data, subset=Quantile %in% c("Q1", "Q2"))
listy=q1.2.sp$SpeciesCode
q1.2.pa.cdm=pa.cdm[,-which(names(pa.cdm)%in% listy)]
hull.matrix.minus.q1.2=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
q.hull=convhulln(pca.scores[names(q1.2.pa.cdm[i,q1.2.pa.cdm[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.minus.q1.2[i,1]=q.hull.vol
}

write.csv(hull.matrix.minus.q1.2, file="hull.matrix.minus.q1.q2.csv")


q4.5.sp=subset(all.data, subset=Quantile %in% c("Q4", "Q5"))
listy=q4.5.sp$SpeciesCode
q4.5.pa.cdm=pa.cdm[,-which(names(pa.cdm)%in% listy)]
hull.matrix.minus.q4.5=matrix(data=NA, nrow=16,ncol=1)
for(i in 1:16){
q.hull=convhulln(pca.scores[names(q4.5.pa.cdm[i,q4.5.pa.cdm[i,]>0]),], options="FA")
q.hull.vol=q.hull$vol
hull.matrix.minus.q4.5[i,1]=q.hull.vol
}

write.csv(hull.matrix.minus.q4.5, file="hull.matrix.minus.q4.q5.csv")

# correlation test between species richness of each quantile and overall species richness
## Table.2
cor.test(rowSums(pa.cdm), rowSums(q1.cdm.pa))
# r = 0.47, p = 0.06
cor.test(rowSums(pa.cdm), rowSums(q2.cdm.pa))
# r = 0.74, p = 0.001
cor.test(rowSums(pa.cdm), rowSums(q3.cdm.pa))
# r = 0.87, p < 0.001
cor.test(rowSums(pa.cdm), rowSums(q4.cdm.pa))
# r = 0.90, p < 0.001
cor.test(rowSums(pa.cdm), rowSums(q5.cdm.pa))
# r = 0.64, p = 0.01


## Figure.1

elevation=c(250, 300, 350, 400, 450, 500, 550, 600, 650,
                                  700, 750, 800, 850, 900, 950, 1000)

# Correlation between species richness and elevation
pa.cdm=read.csv("pa.cdm.csv", header=T, row.names=1)
total.rich=as.numeric(rowSums(pa.cdm))
cor.test(elevation, total.rich)
# r = -0.78, p < 0.001

plot(total.rich~elevation, xlab="Elevation (m)", ylab="Species Richness", pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(lm(total.rich~elevation))

# Correlation between species abundance and elevation
abund.cdm=read.csv("abund.cdm.csv", header=T, row.names=1)
total.abund=as.numeric(rowSums(abund.cdm))
cor.test(elevation, total.abund)
# r = 0.59, r = 0.02

plot(total.abund~elevation, xlab="Elevation (m)", ylab="Species Abundance", pch=19,cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(lm(total.abund~elevation))

# percentage of abundance and richness made up of wide and narrow ranging species
percent.wide.abund=rowSums(wide.cdm)/rowSums(abund.cdm)
percent.wide.rich=rowSums(wide.cdm.pa)/rowSums(pa.cdm)

percent.narrow.abund=rowSums(narrow.cdm)/rowSums(abund.cdm)
percent.narrow.rich=rowSums(narrow.cdm.pa)/rowSums(pa.cdm)

write.csv(percent.wide.abund, file="percent.wide.abund.csv")
write.csv(percent.wide.rich, file="percent.wide.rich.csv")
write.csv(percent.narrow.abund, file="percent.narrow.abund.csv")
write.csv(percent.narrow.rich, file="percent.narrow.rich.csv")

# correlation test between the percentage of richness and abundance of wide and narrow ranging species and elevation
wide.abund=as.numeric(percent.wide.abund)
cor.test(wide.abund, elevation)
# r = 0.48, p = 0.06
wide.rich=as.numeric(percent.wide.rich)
cor.test(wide.rich, elevation)
# r = -0.12, p = 0.66

narrow.abund=as.numeric(percent.narrow.abund)
cor.test(narrow.abund, elevation)
# r = -0.48, p = 0.06
narrow.rich=as.numeric(percent.narrow.rich)
cor.test(narrow.rich, elevation)
# r = 0.12, p = 0.66


# correlation between abundance made of wide ranging species and elevation
wide.abundance=as.numeric(rowSums(wide.cdm))
cor.test(elevation, wide.abundance)
# r = 0.61, p = 0.01

plot(wide.abundance~elevation, xlab="Elevation (m)", ylab="Wide-Ranging Species Abundance", pch=19,cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(lm(wide.abundance~elevation))

# correlation between abundance of narrow ranging species and elevation
narrow.abundance=as.numeric(rowSums(narrow.cdm))
cor.test(elevation, narrow.abundance)
# r = 0.21, p = 0.45

plot(narrow.abundance~elevation, xlab="Elevation (m)", ylab="Narrow-Ranging Species Abundance", pch=19,cex.lab=1.5, cex.axis=1.5)
abline(lm(narrow.abundance~elevation))

# Correlation between narrow-ranging species richness and elevation
narrow.rich=as.numeric(rowSums(narrow.cdm.pa))
cor.test(elevation, narrow.rich)
# r = -0.72, p = 0.002

plot(narrow.rich~elevation, xlab="Elevation (m)", ylab="Narrow-Ranging Species Richness", pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(lm(narrow.rich~elevation))

# Correlation between wide-ranging species richness and elevation
wide.rich=as.numeric(rowSums(wide.cdm.pa))
cor.test(elevation, wide.rich)
# r = -0.73, p = 0.001

plot(wide.rich~elevation, xlab="Elevation (m)", ylab="Wide-Ranging Species Richness", pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(lm(wide.rich~elevation))

# Correlation between species range size and abundance
cor.test(all.data$range, all.data$abundance)
# r = 0.09, p = 0.37


# Plots of correlations
par(mfrow=c(2,2))
plot(total.rich~elevation, xlab="Elevation (m)", ylab="Species Richness", pch=19)
abline(lm(total.rich~elevation))
plot(total.abund~elevation, xlab="Elevation (m)", ylab="Species Abundance", pch=19)
abline(lm(total.abund~elevation))
plot(narrow.abundance~elevation, xlab="Elevation (m)", ylab="Narrow-Ranging Species Abundance", pch=19)
abline(lm(narrow.abundance~elevation))
plot(wide.abundance~elevation, xlab="Elevation (m)", ylab="Wide-Ranging Species Abundance", pch=19)
abline(lm(wide.abundance~elevation))

### Plot of number of species per quantile

all.data$Quantile=as.factor(all.data$Quantile)
plot(all.data$Quantile, ylab="Number of Species", xlab="Quintiles",
main="Number of Species per Quintile")

## Correlation plot between functional richness and elevation
# Figure 3

hull.matrix.pca=read.csv("hull.matrix.pca.csv", header=T)
hull.matrix.data=hull.matrix.pca[c(1:16),c(1:4)]

cor.test(elevation, hull.matrix.data$Allsp.pca)
# r = -0.66, p = 0.006
cor.test(elevation, hull.matrix.data$Narrowsp.pca)
# r = -0.63, p = 0.01
cor.test(elevation, hull.matrix.data$Widesp.pca)
# r = -0.50, p = 0.05


plot(hull.matrix.data$Allsp.pca~elevation, xlab="Elevation (m)", ylab="Functional Richness", pch=19,cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(lm(hull.matrix.data$Allsp.pca~elevation))
plot(hull.matrix.data$Narrowsp.pca~elevation, xlab="Elevation (m)", ylab="Narrow-Ranging Functional Richness", pch=19,cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(lm(hull.matrix.data$Narrowsp.pca~elevation))
plot(hull.matrix.data$Widesp.pca~elevation, xlab="Elevation (m)", ylab="Wide-Ranging Functional Richness", pch=19,cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(lm(hull.matrix.data$Widesp.pca~elevation))

par(mfrow=c(3,1))
plot(fric.all$All~elevation, xlab="Elevation (m)", ylab="Functional Richness", pch=19)
abline(lm(fric.all$All~elevation))
plot(fric.all$Narrow~elevation, xlab="Elevation (m)", ylab="Narrow-Ranging Functional Richness", pch=19)
abline(lm(fric.all$Narrow~elevation))
plot(fric.all$Wide~elevation, xlab="Elevation (m)", ylab="Wide-Ranging Functional Richness", pch=19)
abline(lm(fric.all$Wide~elevation))



# Plot of functional volume of wide ranging vs. narrow ranging species in each plot
# Figure 4
row.names(hull.matrix.data)=hull.matrix.data[,1]
hull.matrix.data.2=hull.matrix.data[,c(3:4)]
hull.matrix.data.3=t(hull.matrix.data.2)

barplot(as.matrix(hull.matrix.data.3),beside=TRUE, legend=c("Narrow", "Wide"),
args.legend=list(x="topright"), xlab="Elevation (m)",
ylab="Functional Volume", las=2)

# Plot of species and abundance percent loss
## Figure 5A
# Species is considered lost if the lowest elevation it occurs is > 1000 m above sea level which is the top of the mountain.

sonadora=read.csv("sonadora.csv", header=T, row.names=1)

table(sonadora$Lost)
table(sonadora$stemSpeciesCode, sonadora$Lost)

# Percent richness lost over time
# 2030 = 1/112 = 0.89%
# 2050 = 2/112 = 1.79%
# 2070 = 5/112 = 4.46%
# 2090 = 7/112 = 6.25%
# 2100 = 30/112 = 26.79%

# Percent abundance lost over time
# 2030 = 52/8693 = 0.60%
# 2050 = 130/8693 = 1.50%
# 2070 = 164/8693 = 1.89%
# 2090 = 292/8693 = 3.36%
# 2100 = 2891/8693 = 33.26%

per.loss.plot=read.csv("per.loss.plot.csv", header=T, row.names=1)

plot(per.loss.plot$per.rich.loss~row.names(per.loss.plot), pch=19, col="black", ylim=c(0,35), type="b", cex=2, xaxt="n", ylab="Percent Lost", xlab="Year")
par(new=TRUE)
plot(per.loss.plot$per.abund.loss~row.names(per.loss.plot), pch=19, col="red",ylim=c(0,35),type="b", cex=2, xaxt="n", ylab="Percent Lost", xlab="Year")
axis(side=1, at=rownames(per.loss.plot), labels=c("2030", "2050", "2070", "2090", 
"2100"))
legend("topleft", legend=c("Species Lost", "Abundance Lost"), col=c("black", "red"), lty=1, pch=19)

# Plot of percent of quantiles at each elevation
# Figure 2

percent.q.table=read.csv("percent.q.table.csv", header=T, row.names=1)

q250.pa=decostand(q250, "pa")
rowSums(q250.pa)
q300.pa=decostand(q300, "pa")
rowSums(q300.pa)
q350.pa=decostand(q350, "pa")
rowSums(q350.pa)
q400.pa=decostand(q400, "pa")
rowSums(q400.pa)
q450.pa=decostand(q450, "pa")
rowSums(q450.pa)
q500.pa=decostand(q500, "pa")
rowSums(q500.pa)
q550.pa=decostand(q550, "pa")
rowSums(q550.pa)
q600.pa=decostand(q600, "pa")
rowSums(q600.pa)
q650.pa=decostand(q650, "pa")
rowSums(q650.pa)
q700.pa=decostand(q700, "pa")
rowSums(q700.pa)
q750.pa=decostand(q750, "pa")
rowSums(q750.pa)
q800.pa=decostand(q800, "pa")
rowSums(q800.pa)
q850.pa=decostand(q850, "pa")
rowSums(q850.pa)
q900.pa=decostand(q900, "pa")
rowSums(q900.pa)
q950.pa=decostand(q950, "pa")
rowSums(q950.pa)
q1000.pa=decostand(q1000, "pa")
rowSums(q1000.pa)

library("RColorBrewer")

percent.q.table.matrix=as.matrix(percent.q.table)
colnames(percent.q.table.matrix)=elevation

barplot(percent.q.table.matrix, beside=TRUE, ylab="Percent of Individuals in each Quantile", xlab="Elevation (m)", las=2, ylim=c(0,50), legend=c("Quantile.1", "Quantile.2","Quantile.3", "Quantile.4", "Quantile.5"), args.legend=list(x="topleft"), col=brewer.pal(n=5, name="RdBu"))

# Plot of functional richness changes at range shift time points
## Figure 5B

hull.time=read.csv("hull.time.csv", header=T, row.names=1)
hull.time.2=t(hull.time)

plot(hull.time.2[1,]~colnames(hull.time.2), type="b", xlab="Elevation (m)", ylab="Functional Richness", ylim=c(20,90), xaxt="n", pch=19)
xticks=seq(250, 1000, by=50)
axis(side=1, at=xticks, labels=elevation, las=2)
par(new=T)
plot(hull.time.2[3,]~colnames(hull.time.2), type="b", xlab="", ylab="", axes=FALSE, col="red", xaxt="n", pch=19)
par(new=T)
plot(hull.time.2[5,]~colnames(hull.time.2), type="b", xlab="", ylab="", axes=FALSE, col="blue", xaxt="n", pch=19)
par(new=T)
plot(hull.time.2[7,]~colnames(hull.time.2), type="b", xlab="", ylab="", axes=FALSE, col="yellow", xaxt="n", pch=19)
par(new=T)
plot(hull.time.2[9,]~colnames(hull.time.2), type="b", xlab="", ylab="", axes=FALSE, col="gray", xaxt="n", pch=19)
par(new=T)
plot(hull.time.2[11,]~colnames(hull.time.2), type="b", xlab="", ylab="", axes=FALSE, col="green", xaxt="n", pch=19)

legend("bottomright", legend=c("Present", "Year 2030", "Year 2050", "Year 2070", "Year 2090", "Year 2100"), col=c("black", "red", "blue", "yellow", "gray", "green"), pch=19, cex=0.8)

# Testing for Mid-Domain Effect
# Have to use archived version of package

devtools::install_version("rangemodelR", "1.0.4")
library(rangemodelR)

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
# slope = 0.37, p = 0.344, adjusted R2 = -0.002

# Observed sp. richness distribution is not different from expected sp. richness under MDE

plot(total.rich,rangemod.data$out.df$mod.rich)
abline(lm(total.rich~rangemod.data$out.df$mod.rich))

test=lm(total.rich~rangemod.data$out.df$mod.rich)