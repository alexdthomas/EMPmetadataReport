#restart the metdata analysis using Git in the new directory

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
for (i in 1:length(map.file.names)){
	all.maps[[i]]<-assign(paste("map", i, sep=""), read.delim(map.file.names[i], quote=""))
}

#reset working directory 
setwd("~/EarthMicrobiomeProject/R/EMPmetadataReport")

#explore the list of all studies before reducing to one dataframe
dim(all.maps[[1]])
head(all.maps[[1]])
sum(is.na(all.maps[[1]]))
boxplot(is.na(all.maps[[1]]))
sort(apply(all.maps[[1]], 2, function(x){length(which(is.na(x)))}))
all.maps[[1]][, c("KEY_SEQ", "REGION", "AGE_IN_YEARS")]

#so can make a list of how many values are NA out of total
unlist(lapply(all.maps, function(x) paste(length(which(is.na(x))), 
																					nrow(x)*ncol(x), sep="/")))

#mean NA?
sort(unlist(lapply(all.maps, function(x) mean(is.na(x)))))
#one of them is 25% NA?
which(unlist(lapply(all.maps, function(x) mean(is.na(x))))>0.24)

#explore that one
head(all.maps[[41]])
sort(apply(all.maps[[41]], 2, function(x){length(which(is.na(x)))}))

#can get study titles
lapply(all.maps, function(x) unique(x[,"TITLE"]))
#great, some studies gave each sample a unique  title... annoying
which(lapply(all.maps, function(x) length(unique(x[,"TITLE"])))>1)
#oh, only 2

#copy to edit
all.maps2<-all.maps

#should probably fix this for later reporting...
sort(all.maps2[[9]][,"TITLE"])
#so this is an Excel 'draw down error', 
#started with "Intertidal microbes 16s for 2009 and 2010"
#have to juggle a little as TITLE is factor
all.maps2[[9]][,"TITLE"]<-as.character(all.maps2[[9]][,"TITLE"])
all.maps2[[9]][,"TITLE"]<-"Intertidal microbes 16s for 2009 and 2010"
all.maps2[[9]][,"TITLE"]<-as.factor(all.maps2[[9]][,"TITLE"])

#the other one
sort(all.maps2[[16]][,"TITLE"])
#same issue
all.maps2[[16]][,"TITLE"]<-as.character(all.maps2[[16]][,"TITLE"])
all.maps2[[16]][,"TITLE"]<-"EPOCA_Svalbard2018"
all.maps2[[16]][,"TITLE"]<-as.factor(all.maps2[[16]][,"TITLE"])

#and titles
unlist(lapply(all.maps2, function(x) unique(x[,"TITLE"])))
#ah good

#make a little table...
study_count_NA_table<-
	data.frame(TITLE=unlist(lapply(all.maps2, function(x) unique(x[,"TITLE"]))),
						 No_SAMPLES=unlist(lapply(all.maps2, nrow)),
						 No_COLUMNS=unlist(lapply(all.maps2, ncol)),
						 No_NA=unlist(lapply(all.maps2, function(x) length(which(is.na(x))))),
						 No_VALUES=unlist(lapply(all.maps2, function(x) nrow(x)*ncol(x))),
						 PERCENT_NA=(round(unlist(lapply(all.maps2, function(x) mean(is.na(x)))),2)*100)
	)
#sort by TITLE
study_count_NA_table<-study_count_NA_table[order(as.character(study_count_NA_table$TITLE)),]

#export to subfolder of working directory, remove row names
write.csv(study_count_NA_table, 
					file.path(paste(getwd(), 
													"outputs/study_count_NA_table.csv", sep="/")), 
					row.names=FALSE)


#now a little more depth...
#check out column names
study_column_table<-sort(table(unlist(lapply(all.maps2, colnames))), decreasing=TRUE)

#clean up for exporting
study_column_table<-as.data.frame(study_column_table)
#study names are currently row names, add as column
study_column_table$COLUMN_NAMES<-rownames(study_column_table)
#rename count column
colnames(study_column_table)[1]<-"No_STUDIES"
#replace row names
rownames(study_column_table)<-as.character(1:nrow(study_column_table))
#reorder
study_column_table<-study_column_table[,c(2,1)]
#subset to fit on MSword page, place side by side
study_column_table<-cbind(study_column_table[1:40,], study_column_table[41:80,])

#add a little note, how many columns not shown...
study_column_table<-rbind(study_column_table, 
													c(paste("*", 
																	paste(as.character(length(unique(unlist(lapply(all.maps2, colnames))))-nrow(study_column_table)*2), 
																				"columns not shown", sep=" "), 
																	sep=""), "", "", ""))

#this is handy...
write.csv(study_column_table, 
					file.path(paste(getwd(), 
													"outputs/study_column_table.csv", sep="/")), 
					row.names=FALSE)

#explore most common column values
com.col<-sort(table(unlist(lapply(all.maps2, colnames))), decreasing=TRUE)
#subset columns that are in every study
com.col<-com.col[which(com.col==49)]
com.col<-names(com.col)
lapply(all.maps2, function(x) head(x[ ,which(colnames(x) %in% com.col)], 2))
#interesting, but kind of a pain...

#kinda want to screen for outliers, easiest with everything in one dataframe...

#merge all edited mapping files into one dataframe
#emp.map<-Reduce(function(x, y) merge(x, y, all=TRUE), all.maps2)

#finally realized what the warning messages I was getting with
#above code were for!
#data frames with same column name, but different class were
#not merging correctly and values were being converted to NA
#here http://stackoverflow.com/questions/1632772/appending-rows-to-a-dataframe-the-factor-problem

#so need to find which columns in which data frames have different classes...
#uh...