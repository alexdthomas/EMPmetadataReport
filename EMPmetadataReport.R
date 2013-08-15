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

#search colnames
colnames(emp.map)[grep('contact', colnames(emp.map), ignore.case=TRUE)]

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
col.class.dif<-unique(unlist(col.class.dif))

l1 <- vector('list', 11)

for(i in 1:length(col.class.dif)){
	
	#get example of values in column for each data frame
	check.val<-lapply(all.maps, function(x) {head(x[colnames(x) %in% 
							unique(unlist(col.class.dif))[i]], 1)})
	#give elements names
	#names(check.val)<-paste("map", 1:49, sep="")
	#replace empty data frames with NA
	check.val[which(lapply(check.val, length)==0)]<-"No Field"
	
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
	check.class[which(lapply(check.class, length)==0)]<-"No Field"
	
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

#nice thought, not much better...
#study_class_diff_table.melt<-melt(data.frame(TITLE=study_class_diff_table[,"TITLE"], l2), id.vars="TITLE")

#study_class_diff_table.cast<-dcast(study_class_diff_table.melt, variable~ TITLE, value.var="value")
#head(study_class_diff_table.cast)

#investigate all of the most common fields using this script...
col.comm<-c(study_column_table[1:40,1], study_column_table[,3])
col.comm<-col.comm[1:43]
#add contact fields (mostly missing...?)
col.comm<-c(col.comm, "PRINCIPAL_INVESTIGATOR_CONTACT", 
						"LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")

l1 <- vector('list', 43)

for(i in 1:length(col.comm)){
	
	#get example of values in column for each data frame
	check.val<-lapply(all.maps, function(x) {head(x[colnames(x) %in% 
										col.comm[i]], 1)})
	#give elements names
	#names(check.val)<-paste("map", 1:49, sep="")
	#replace empty data frames with NA
	check.val[which(lapply(check.val, length)==0)]<-"No Field"
	
	#build table
	l1[[i]]<-unlist(check.val)
	names(l1[i])<-paste(col.comm[i], ".ex", sep="")
	
}


l2 <- vector('list', 43)
for(i in 1:length(col.comm)){
	
	#get class of each column for each data frame
	check.class<-lapply(all.maps, function(x) {lapply(x[colnames(x) %in% 
																												col.comm[i]], class)})
	#convert empty values to NA
	check.class[which(lapply(check.class, length)==0)]<-"No Field"
	
	l2[[i]]<-unlist(check.class)
	names(l2[i])<-paste(col.comm[i], ".class", sep="")
	
}

#create table
l1<-data.frame(l1)
colnames(l1)<-paste(col.comm, ".ex", sep="")

l2<-data.frame(l2)
colnames(l2)<-paste(col.comm, ".class", sep="")

study_class_comm_col_table<-data.frame(TITLE=unlist(lapply(all.maps, function(x) unique(x[,"TITLE"]))))
study_class_comm_col_table<-cbind(study_class_comm_col_table, l1, l2)

#sort by column name 
study_class_comm_col_table<-study_class_comm_col_table[, c("TITLE", sort(colnames(study_class_comm_col_table)[2:87]))]
#export
write.csv(study_class_comm_col_table, 
					file.path(paste(getwd(), 
													"outputs/study_class_comm_col_table.csv", sep="/")), 
					row.names=FALSE)

#can create list of studies and key data to export as report?
test<-list(list(TITLE="TITLE1", CONTACT=c("Contact", "No Field", NA), 
								RUN_DATE=2007, PRIMERS=NA), 
					 list(TITLE="TITLE1", CONTACT=c("Contact", "No Field", NA), 
					 		 RUN_DATE="8/24/2007", PRIMERS="FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT"))

test<-lapply(all.maps, function(x) list(TITLE=unique(x[,"TITLE"]),			
	CONTACTS=paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "),
	RUN_DATE=ifelse(TRUE %in% (colnames(x) %in% "RUN_DATE"), paste(range(x[,"RUN_DATE"]), collapse="-"), "No Field"),
	RUN_DATE.class=ifelse(TRUE %in% (colnames(x) %in% "RUN_DATE"), class(x[,"RUN_DATE"]), "No Field"),
	DEPTH=ifelse(TRUE %in% (colnames(x) %in% "DEPTH"), paste(range(x[,"DEPTH"]), collapse="-"), "No Field"),
	DEPTH.class=ifelse(TRUE %in% (colnames(x) %in% "DEPTH"), class(x[,"DEPTH"]), "No Field"),
	BREAK=""																	
))


test<-lapply(all.maps, function(x) list(			
			CONTACTS=paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "),
			RUN_DATE=ifelse(TRUE %in% (colnames(x) %in% "RUN_DATE"), paste(range(x[,"RUN_DATE"]), collapse="-"), "No Field"),
			RUN_DATE.class=ifelse(TRUE %in% (colnames(x) %in% "RUN_DATE"), class(x[,"RUN_DATE"]), "No Field"),
			DEPTH=ifelse(TRUE %in% (colnames(x) %in% "DEPTH"), paste(range(x[,"DEPTH"]), collapse="-"), "No Field"),
			DEPTH.class=ifelse(TRUE %in% (colnames(x) %in% "DEPTH"), class(x[,"DEPTH"]), "No Field"),
			BREAK=""																	
))

names(test)<-unlist(lapply(all.maps, function(x) unique(x[,"TITLE"])))

ifelse(TRUE %in% (colnames(all.maps[[1]]) %in% "RUN_DATE"), paste(range(all.maps[[1]][,"RUN_DATE"]), collapse="-"), "No Field")
summary(all.maps[[1]][,"RUN_DATE"])
all.maps[[1]][,"RUN_DATE"]
table(all.maps[[1]][,"RUN_DATE"])
cat(range(all.maps[[1]][,"RUN_DATE"]))
head(unique(all.maps[[2]][,"TITLE"]),1)

#export
wb = loadWorkbook("EMP_metadata_issues.xlsx", create = TRUE)
# Create a new sheet
createSheet(wb, name = "mysheet")
# cumulative length (rows) of matrices
# +2 = 1 for list names, 1 for header row
cumlen = cumsum(c(1, head(sapply(test, length), n = -1) + 2))
#cumlen = cumsum(c(1, rep(2,48) + 2))

# Write data rows (implicitly vectorized!)
writeWorksheet(wb, data = test, sheet = "mysheet", startRow = cumlen + 1, header = TRUE, rownames=lapply(test, rownames))
# Write list names
writeWorksheet(wb, data = as.list(names(test)), sheet = "mysheet", startRow = cumlen, header = FALSE)
saveWorkbook(wb)

#some other attempts
write.table(unlist(test))
class(unlist(test))
test.names<-unlist(test)
test.names<-unlist(test.names$names)
test<-cbind(test.names, unlist(test))
class(test.bind)

#test<-data.frame(test)
test.bind<-cbind(unlist(test$names), unlist(test))
test.bind<-cbind(rownames(test.bind), test.bind)
colnames(test.bind)<-c("FIELD", "VALUE")
rownames(test.bind)<-seq(1:343)

write.table(test, 
					file.path(paste(getwd(), 
													"outputs/EMP_metadata_issues.csv", sep="/")), 
					row.names=FALSE, sep=",")

test.data.frame<-do.call("rbind", test)
write.csv(test.data.frame, 
					file.path(paste(getwd(), 
													"outputs/EMP_metadata_issues.csv", sep="/")), 
					row.names=TRUE)

test<-lapply(all.maps, function(x) data.frame(Value=c(
	paste(unique(x[,"TITLE"]), collapse=""),
	paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "),
	ifelse(TRUE %in% (colnames(x) %in% "RUN_DATE"), paste(range(x[,"RUN_DATE"]), collapse="-"), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "RUN_DATE"), class(x[,"RUN_DATE"]), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "DEPTH"), paste(range(x[,"DEPTH"]), collapse="-"), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "DEPTH"), class(x[,"DEPTH"]), "No Field"),
	paste("")),																						
																							row.names=c("TITLE","CONTACTS", "RUN_DATE","RUN_DATE.class", "DEPTH","DEPTH.class", "break")  																	
))


test.data.frame<-do.call("rbind", test)
class(test.data.frame)
test.data.frame$Field<-rep(c("TITLE", "CONTACTS", "RUN_DATE","RUN_DATE.class", "DEPTH","DEPTH.class", ""), 49)
test.data.frame<-test.data.frame[,c("Field", "Value")]

write.csv(test.data.frame, 
					file.path(paste(getwd(), 
													"outputs/EMP_metadata_issues.csv", sep="/")), 
					row.names=FALSE)