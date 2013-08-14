#restart the metdata analysis using Git in the new directory
#redo first_to_merge_class_issues, 
#this time set characterAsFactor=FALSE when import

library(BiocInstaller)
library(XML)
library(XMLSchema)
library(SSOAP)
library(rols)
library(reshape2)
library(plyr)


setwd("~/EarthMicrobiomeProject/R/EMPmetadataReport")

def.par <- par(no.readonly = TRUE) # save default, for resetting...

#set directory to location of all mapping files
setwd("~/EarthMicrobiomeProject/QIIME_metadata_download/mapping_files")

#list mapping files
map.file.names<-list.files()
#create list of data frames for each mapping file
all.maps<-list()

#import all mapping files, name as individual objects and place in list
#note characterAsFactors=FALSE to reduce class issues
#though this only removes issues of character -> factor and vice versa issues
#does not resolve charcter-> integer, numeric, date etc...
for (i in 1:length(map.file.names)){
	all.maps[[i]]<-assign(paste("map", i, sep=""), read.delim(map.file.names[i], quote="", stringsAsFactors=FALSE))
}

#reset working directory 
setwd("~/EarthMicrobiomeProject/R/EMPmetadataReport")

#check merge
emp.map<-Reduce(function(x, y) merge(x, y, all=TRUE), all.maps)
#well no more warnings, so character/ factor resolved


#pick through each study and find which common columns have different classes
map.class.tmp<-list()
map.class<-list()

for(i in 1:length(all.maps)){
	for(j in 1:length(all.maps)){
		all.maps.cnames<-colnames(all.maps[[i]]
															[which(unlist(lapply(all.maps[i], colnames)) %in%  
																		 	unlist(lapply(all.maps[j], colnames)))])
		dat1.class<-lapply(all.maps[[i]][all.maps.cnames], class)
		dat2.class<-lapply(all.maps[[j]][all.maps.cnames], class)
		dat1.class<-unlist(dat1.class)
		dat2.class<-unlist(dat2.class)
		map.class.tmp[[j]]<-cbind(dat1.class[which(!dat1.class==dat2.class)], 
													 dat2.class[which(!dat2.class==dat1.class)])
	}
	map.class[[i]]<-map.class.tmp
}
map.class
#haha, and this works