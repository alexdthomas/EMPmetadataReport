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

#can still find columns with same names and different classes?
emp.col.names<-unique(unlist(lapply(all.maps, colnames)))

#create subset to test code
few.maps<-list(all.maps[[1]], all.maps[[2]], all.maps[[3]])
lapply(few.maps, head)
emp.col.names<-unique(unlist(lapply(few.maps, colnames)))

head(few.maps[[1]])
class(few.maps[[1]])
identical(class(few.maps[[1]][,1]), class(few.maps[[2]][,1]))

dim(
	few.maps[[1]][,which(unlist(lapply(few.maps[1], colnames)) %in%  unlist(lapply(few.maps[2], colnames)))]
)

dim(few.maps[[1]])

which(unlist(lapply(few.maps[1], colnames)) %in%  
				unlist(lapply(few.maps[2], colnames)))

few.maps.cnames<-colnames(few.maps[[1]]
													[,which(unlist(lapply(few.maps[1], colnames)) %in%  
																		unlist(lapply(few.maps[2], colnames)))])

summary(lapply(few.maps[[3]][,c(few.maps.cnames)], class) %in%
					lapply(few.maps[[4]][,c(few.maps.cnames)], class))

#test loop methods
dat1<-data.frame(col1=as.integer(c(1,2,3)), col3=as.integer(c(4,5,6)))
dat2<-data.frame(col1=as.numeric(c(1,2,3)), col3=as.integer(c(4,5,6)))
dat3<-data.frame(col4=as.integer(c(1,2,3)), col3=as.character(c(4,5,6)))

dat<-list(dat1, dat2, dat3)

dat.names<-colnames(dat[[1]]
										[which(unlist(lapply(dat[1], colnames)) %in%  
													 	unlist(lapply(dat[2], colnames)))])

dat1.class<-lapply(dat[[1]][dat.names], class)
dat2.class<-lapply(dat[[2]][dat.names], class)

dat.class1<-which(!dat1.class %in% dat2.class)
dat.class2<-which(!dat2.class %in% dat1.class)

dat1.class==dat2.class
isTRUE(all.equal(dat1.class, dat2.class))
duplicated(dat1.class, dat2.class)
duplicated(dat1.class[2],dat2.class[2])
dat1.class[2] 
dat2.class[2]
all.equal(dat1.class[2],dat2.class[2])
intersect(dat1.class, dat2.class)
dat1.class[intersect(names(dat1.class), names(dat2.class))] == dat2.class[intersect(names(dat1.class), names(dat2.class))]

submaps1<-which(lapply(dat1.class, function(x) isTRUE(all.equal(x, dat2.class))))

#found this example here
#http://stackoverflow.com/questions/12958335/compare-two-character-vectors-matching-names
#modified as list, realizes issue was using lists
x <- list("a", "b", "c", "d", "e")
names(x) <- c("foo", "bar", "baz", "qux", "grault")

y <- list("c", "a", "d", "b")
names(y) <- c("bar", "foo", "qux", "corge")
x
x[intersect(names(x), names(y))] == y[intersect(names(x), names(y))]

#just unlist
dat1.class<-unlist(dat1.class)
dat2.class<-unlist(dat2.class)

#and this works
dat1.class==dat2.class

#can get fancy
cbind(dat1.class[which(!dat1.class==dat2.class)], dat2.class[which(!dat2.class==dat1.class)])

#this only works to compare one dataset to the next one
map.class1<-list()
map.class2<-list()

for(i in 1:(length(all.maps)-1)){
	all.maps.cnames<-colnames(all.maps[[i]]
														[,which(unlist(lapply(all.maps[i], colnames)) %in%  
																			unlist(lapply(all.maps[-i], colnames)))])
	
	dat1.class<-lapply(all.maps[[i]][,c(all.maps.cnames)], class)
	dat2.class<-lapply(all.maps[[-i]][,c(all.maps.cnames)], class)
	submaps1<-which(!dat1.class %in% dat2.class) 
	submaps2<-which(!dat2.class %in% dat1.class)
	map.class1[[i]]<-colnames(all.maps[[i]][submaps1])
	map.class2[[i]]<-colnames(all.maps[[i+1]][submaps2])
	
}

map.class1
map.class2

#should go through each study and compare to all othe studies
all.maps.cnames<-colnames(dat[[1]]
													[which(unlist(lapply(dat[1], colnames)) %in%  
																 	unlist(lapply(dat[2], colnames)))])
dat1.class<-lapply(dat[[1]][all.maps.cnames], class)
dat2.class<-lapply(dat[[2]][all.maps.cnames], class)

submaps1<-which(lapply(dat1.class, function(x) isTRUE(all.equal(x, dat2.class)))==FALSE)
submaps2<-which(lapply(dat2.class, function(x) isTRUE(all.equal(x, dat1.class)))==FALSE)


submaps1<-which(!dat1.class %in% dat2.class) 
submaps2<-which(!dat2.class %in% dat1.class)
map.class1[[i]]<-colnames(dat[[1]][submaps1])
map.class2[[i]]<-colnames(dat[[2]][submaps2])

all.maps.cnames<-colnames(all.maps[[1]]
													[which(unlist(lapply(all.maps[1], colnames)) %in%  
																 	unlist(lapply(all.maps[45], colnames)))])
dat1.class<-lapply(all.maps[[1]][all.maps.cnames], class)
dat2.class<-lapply(all.maps[[45]][all.maps.cnames], class)

isTRUE(all.equal(dat1.class, dat2.class))
lapply(dat1.class, function(x) isTRUE(all.equal(x, dat2.class)))

which(sapply(dat1.class, function(x) isTRUE(all.equal(x, dat2.class)))) 

submaps1<-which(lapply(dat1.class, function(x) isTRUE(all.equal(x, dat1.class)))==FALSE)
submaps2<-which(lapply(dat2.class, function(x) isTRUE(all.equal(x, dat1.class)))==FALSE)

submaps1<-which(!dat1.class %in% dat2.class) 
submaps2<-which(!dat2.class %in% dat1.class)
map.class1[[i]]<-colnames(dat[[1]][submaps1])
map.class2[[i]]<-colnames(dat[[45]][submaps2])

#ok, start with first dataset
map.class1<-list()

for(i in 2:length(all.maps)){
	all.maps.cnames<-colnames(all.maps[[1]]
														[which(unlist(lapply(all.maps[1], colnames)) %in%  
																	 	unlist(lapply(all.maps[i], colnames)))])
	dat1.class<-lapply(all.maps[[1]][all.maps.cnames], class)
	dat2.class<-lapply(all.maps[[i]][all.maps.cnames], class)
	dat1.class<-unlist(dat1.class)
	dat2.class<-unlist(dat2.class)
	map.class1[[i]]<-cbind(dat1.class[which(!dat1.class==dat2.class)], 
												 dat2.class[which(!dat2.class==dat1.class)])
	
}
map.class1
#haha, and this works