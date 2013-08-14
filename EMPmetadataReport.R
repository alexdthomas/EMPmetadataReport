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

#explore data
#make a list of how many values are NA out of total
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


#fix titles
sort(all.maps[[9]][,"TITLE"])
#so this is an Excel 'draw down error', 
#started with "Intertidal microbes 16s for 2009 and 2010"
all.maps[[9]][,"TITLE"]<-"Intertidal microbes 16s for 2009 and 2010"

#the other one
sort(all.maps[[16]][,"TITLE"])
#same issue
all.maps[[16]][,"TITLE"]<-"EPOCA_Svalbard2018"

#and titles
unlist(lapply(all.maps, function(x) unique(x[,"TITLE"])))
#ah good

#make a little table...
study_count_NA_table<-
	data.frame(TITLE=unlist(lapply(all.maps, function(x) unique(x[,"TITLE"]))),
						 No_SAMPLES=unlist(lapply(all.maps, nrow)),
						 No_COLUMNS=unlist(lapply(all.maps, ncol)),
						 No_NA=unlist(lapply(all.maps, function(x) length(which(is.na(x))))),
						 No_VALUES=unlist(lapply(all.maps, function(x) nrow(x)*ncol(x))),
						 PERCENT_NA=(round(unlist(lapply(all.maps, function(x) mean(is.na(x)))),2)*100)
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
study_column_table<-sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE)

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
																	paste(as.character(length(unique(unlist(lapply(all.maps, colnames))))-nrow(study_column_table)*2), 
																				"columns not shown", sep=" "), 
																	sep=""), "", "", ""))

#this is handy...
write.csv(study_column_table, 
					file.path(paste(getwd(), 
													"outputs/study_column_table.csv", sep="/")), 
					row.names=FALSE)

#explore most common column values
com.col<-sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE)
#subset columns that are in every study
com.col<-com.col[which(com.col==49)]
com.col<-names(com.col)
lapply(all.maps, function(x) head(x[ ,which(colnames(x) %in% com.col)], 2))

#pick through each study and find which common columns have different classes
map.class.tmp<-list()
map.class<-list()

for(i in 1:length(all.maps)){
	for(j in i:length(all.maps)){
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
map.class[[1]]
#volumous, but works

#now can pick through and see what the differences are...
map.class[[1]][6]
head(all.maps[[1]][c("RUN_DATE", "DEPTH", "EXPERIMENT_DESIGN_DESCRIPTION")])
head(all.maps[[6]][c("RUN_DATE", "DEPTH", "EXPERIMENT_DESIGN_DESCRIPTION")])

#see if can get all the colnames with problems
class(map.class[[1]][6])
rownames(map.class[[1]][[6]])

tmp.rows<-list()
col.class.dif<-list()

for(i in 1:length(map.class)){
	tmp.map<-map.class[[i]]
	for(j in 1:length(map.class[[i]])){
		tmp.rows[[i]]<-list(rownames(tmp.map[[j]]))
	}
	col.class.dif[[i]]<-list(tmp.rows[[i]])
}
unlist(col.class.dif)
unique(unlist(col.class.dif))
sort(table(unlist(col.class.dif)))
#so it is mostly depth and a few others

#look at just these columns across studes
unique(unlist(col.class.dif))
col.class.dif<-unique(unlist(col.class.dif))

#too volumous to look at all, select one at a time
check.find<-lapply(all.maps, function(x) {
	isTRUE(col.class.dif[2] %in% colnames(x))})

check.val<-lapply(all.maps, function(x) {head(x[colnames(x) %in% 
			col.class.dif[2]], 1)})
#compare classes
check.class<-lapply(all.maps, function(x) {lapply(x[colnames(x) %in% 
					col.class.dif[2]], class)})

length(unlist(check.find)) 
length(unlist(check.val))
length(unlist(check.class)) 

cbind(check.find, c(check.val))

names(check.val)<-paste("map", 1:49, sep="")
check.val[which(lapply(check.val, length)==0)]<-NA
unlist(check.val)

study_class_diff_table<-data.frame(TITLE=unlist(lapply(all.maps, function(x) unique(x[,"TITLE"]))))
study_class_diff_tableM<-cbind(study_class_diff_table, 
															 paste(col.class.dif[1], ".ex", sep="")=unlist(check.val))

#do this in a loop
#look at just these columns across studes
l1 <- vector('list', 11)

for(i in 1:length(col.class.dif)){
	
	#get example of values in column for each data frame
	check.val<-lapply(all.maps, function(x) {head(x[colnames(x) %in% 
							unique(unlist(col.class.dif))[i]], 1)})
	#give elements names
	#names(check.val)<-paste("map", 1:49, sep="")
	#replace empty data frames with NA
	check.val[which(lapply(check.val, length)==0)]<-NA
	
	#build table
	l1[[i]]<-unlist(check.val)
	names(l1[i])<-paste(col.class.dif[i], ".ex", sep="")
													 
}


l2 <- vector('list', 11)
for(i in 1:length(col.class.dif)){
	
	#get class of each column for each data frame
	check.class<-lapply(all.maps, function(x) {lapply(x[colnames(x) %in% 
								col.class.dif[i]], class)})
	#convert empty values to NA
	check.class[which(lapply(check.class, length)==0)]<-NA
	
	l2[[i]]<-unlist(check.class)
	names(l2[i])<-paste(col.class.dif[i], ".class", sep="")
	
}

#create table
l1<-data.frame(l1)
colnames(l1)<-paste(col.class.dif, ".ex", sep="")

l2<-data.frame(l2)
colnames(l2)<-paste(col.class.dif, ".class", sep="")

study_class_diff_table<-data.frame(TITLE=unlist(lapply(all.maps, function(x) unique(x[,"TITLE"]))))
study_class_diff_table<-cbind(study_class_diff_table, l1, l2)
	
#sort by column name 
study_class_diff_table<-study_class_diff_table[, c("TITLE", sort(colnames(study_class_diff_table)[2:23]))]
#export
write.csv(study_class_diff_table, 
					file.path(paste(getwd(), 
													"outputs/study_class_diff_table.csv", sep="/")), 
					row.names=FALSE)

