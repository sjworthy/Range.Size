# estimating range sizes from gbif

install.packages("rgbif")
library(rgbif)
library(stringr)
library(tidyr)
library(dplyr)
install.packages("usethis")
usethis::edit_r_environ()
library(geonames)
library(elevatr)
library(rgdal)
library(progress)

all.data=read.csv("./Data/all.data.csv")

all.data$binomial=str_c(all.data$Genus," ",all.data$Species)

# getting taxonomic key from gbif
taxon.key=name_backbone_checklist(all.data$binomial)

higherrank.issue=name_backbone("Tetragastris balsamifera", verbose = TRUE)
higherrank.issue$verbatim_index = c(0,17,0,0)
# remove Tetragastris balsamifera usage key

taxon.key.2=taxon.key[c(1:16,18:112),]

# add new taxon key info for Tetragastris balsamifera

taxon.key.3=rbind(taxon.key.2, higherrank.issue[2,])

higherrank.issue.2=name_backbone("Miconia tetrandra", verbose=TRUE)

# get taxon key
taxon.key.all=taxon.key.3$usageKey
taxon.key.all.2=as.numeric(as.vector(taxon.key.all))

occ_download(pred("hasCoordinate", TRUE),pred("occurrenceStatus","PRESENT"),
                  pred_in("taxonKey",taxon.key.all.2),
                  pred_or(pred_lt("coordinateUncertaintyInMeters",10000),pred_isnull("coordinateUncertaintyInMeters")),
                  format = "SIMPLE_CSV")

all.gbif=occ_download_get('0163342-230224095556074', path = "/Users/samanthaworthy/Documents/Range.size/Range.size/Data") %>%
  occ_download_import()

# citation: https://www.gbif.org/citation-guidelines
# DOI: 10.15468/dl.gbe8zm
# Citation: GBIF Occurrence Download https://doi.org/10.15468/dl.gbe8zm Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-04-12

tetra=occ_search(pred("taxonKey", 8660243), format = "SIMPLE_CSV")

all.t.balsamifera=occ_download_get('0163788-230224095556074', path = "/Users/samanthaworthy/Documents/Range.size/Range.size/Data") %>%
  occ_download_import()

all.t.balsamifera.2 = subset(all.t.balsamifera, all.t.balsamifera$stateProvince %in% "Puerto Rico" | all.t.balsamifera$countryCode %in% "PR")
# no lat long for this species

# getting elevation based on lat/long from gbif

# some have 0,0 for coordinates so need to remove

all.gbif.2=subset(all.gbif, all.gbif$decimalLatitude != 0)


# subset for Puerto Rico

all.gbif.3 = subset(all.gbif.2, all.gbif.2$stateProvince %in% "Puerto Rico" | all.gbif.2$countryCode %in% "PR")

options(geonamesUsername="sjworthy")
options(geonamesHost="api.geonames.org")

test=elevation(all.gbif.3,  username = options(geonamesUsername="sjworthy"))

# try to get elevation a different way

test.2=test[,c(23,22)]
new.elev=get_elev_point(test.2, prj="EPSG:4326", src = "aws")
test$elevatr.elev=new.elev$elevation

# geonames seem most accurate, but need to change negative values

neg.elev=subset(test, test$elevation_geonames < 0 & test$elevatr.elev < 0)

# subset when both elevation_geonames and elevatr.elev are negative

test.3 = subset(test, test$elevation_geonames > 0 & test$elevatr.elev > 0)
test.3$row.number = 1:5639

test.3=test.3[-5307,]

# test.3=test[c(1:41,43:154,156:372,374:778,780:990),]

species.elevation = test.3 %>%
  group_by(species) %>%
  summarise(ele.min=min(elevation, na.rm=TRUE),
            ele.max=max(elevation, na.rm=TRUE),
            geo.ele.min=min(elevation_geonames, na.rm=TRUE),
            geo.ele.max=max(elevation_geonames, na.rm=TRUE),
            elevar.ele.min=min(elevatr.elev, na.rm=TRUE),
            elevar.ele.max=max(elevatr.elev, na.rm = TRUE))

# write.csv(species.elevation, file="gbif.elevation.2.csv")






