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

#name studies in list
names(all.maps)<-unlist(lapply(all.maps, function(x) unique(x[,"TITLE"])))

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
study_column_table<-cbind(study_column_table[1:38,], study_column_table[39:76,])

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

l1 <- vector('list', length(col.class.dif))

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


l2 <- vector('list', length(col.class.dif))
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

#nice though, not much better...
#study_class_diff_table.melt<-melt(data.frame(TITLE=study_class_diff_table[,"TITLE"], l2), id.vars="TITLE")

#study_class_diff_table.cast<-dcast(study_class_diff_table.melt, variable~ TITLE, value.var="value")
#head(study_class_diff_table.cast)

#investigate all of the most common fields using this script...
#select fields present in greater than half of the studies
col.comm<-names(which(sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE)>length(all.maps)/2))

#add contact fields (mostly missing...?)
col.comm<-c(col.comm, "PRINCIPAL_INVESTIGATOR_CONTACT", 
						"LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")

l1 <- vector('list', length(col.comm))

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


l2 <- vector('list', length(col.comm))
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

EMP_metadata_issues<-lapply(all.maps, function(x) data.frame(Value=c(
	paste(unique(x[,"TITLE"]), collapse=" "),
	paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "),
	ifelse(TRUE %in% (colnames(x) %in% "COLLECTION_DATE"), 
				 paste(paste(range(x[,"COLLECTION_DATE"], finite=TRUE), collapse=" to "),
				 			paste(length(which(is.na(x[,"COLLECTION_DATE"]))), length(x[,"COLLECTION_DATE"]), sep=" NA of "), sep="; "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "COLLECTION_DATE"), class(x[,"COLLECTION_DATE"]), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "SAMP_SIZE"), 
				 paste(paste(range(x[,"SAMP_SIZE"], finite=TRUE), collapse=" to "),
						paste(length(which(is.na(x[,"SAMP_SIZE"]))), length(x[,"SAMP_SIZE"]), sep=" NA of "), sep="; "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "DEPTH"),  
				 paste(paste(range(x[,"DEPTH"], finite=TRUE), collapse=" to "),
				 			paste(length(which(is.na(x[,"DEPTH"]))), length(x[,"DEPTH"]), sep=" NA of "), sep="; "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "DEPTH"), class(x[,"DEPTH"]), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "RUN_DATE"), 
				 paste(paste(range(x[,"RUN_DATE"], finite=TRUE), collapse=" to "),
				 			paste(length(which(is.na(x[,"RUN_DATE"]))), length(x[,"RUN_DATE"]), sep=" NA of "), sep="; "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "RUN_DATE"), class(x[,"RUN_DATE"]), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "SEQUENCING_METH"), 
				 paste(unique(x[which(!is.na(x[,"SEQUENCING_METH"])),"SEQUENCING_METH"]), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "LIBRARY_CONSTRUCTION_PROTOCOL"), 
				 paste(unique(x[which(!is.na(x[,"LIBRARY_CONSTRUCTION_PROTOCOL"])),"LIBRARY_CONSTRUCTION_PROTOCOL"]), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "PCR_PRIMERS"), 
				 paste(unique(x[which(!is.na(x[,"PCR_PRIMERS"])),"PCR_PRIMERS"]), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "PLATFORM"), 
				 paste(unique(x[which(!is.na(x[,"PLATFORM"])),"PLATFORM"]), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "RUN_CENTER"), 
				 paste(unique(x[which(!is.na(x[,"RUN_CENTER"])),"RUN_CENTER"]), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "SAMPLE_CENTER"), 
				 paste(unique(x[which(!is.na(x[,"SAMPLE_CENTER"])),"SAMPLE_CENTER"]), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "SAMPLE_LOCATION"), 
				 paste(unique(x[which(!is.na(x[,"SAMPLE_LOCATION"])),"SAMPLE_LOCATION"]), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "TARGET_GENE"),
				 paste(unique(x[which(!is.na(x[,"TARGET_GENE"])),"TARGET_GENE"]), collapse=" "), "No Field"),
	paste("")),	
	#FIELD=c("TITLE","CONTACTS", "RUN_DATE","RUN_DATE.class", "DEPTH","DEPTH.class", "break"),  																	
	row.names=c("TITLE", 
							"CONTACTS", 
							"COLLECTION_DATE_range", 
							"COLLECTION_DATE_class", 
							"SAMP_SIZE_range",
							"DEPTH_range",
							"DEPTH.class", 
							"RUN_DATE_range",
							"RUN_DATE.class",
							"SEQUENCING_METH", 
							"LIBRARY_CONSTRUCTION", 
					  	"PCR_PRIMERS", 
							"PLATFORM", 
							"RUN_CENTER",  
							"SAMPLE_CENTER", 
							"SAMPLE_LOCATION", 
							"TARGET_GENE", 
							"")  																	
))

names(EMP_metadata_issues)<-unlist(lapply(all.maps, function(x) unique(x[,"TITLE"])))

#test
ifelse(TRUE %in% (colnames(all.maps[[1]]) %in% "RUN_CENTER"), paste(head(all.maps[[1]][,"RUN_CENTER"]), 1), "No Field")
ifelse(TRUE %in% (colnames(all.maps[[1]]) %in% "RUN_CENTER"), paste(range(all.maps[[1]][,"RUN_CENTER"]), collaps= ""), "No Field")

ifelse(TRUE %in% (colnames(all.maps[[1]]) %in% "RUN_CENTER"), class(all.maps[[1]][,"RUN_CENTER"]), "No Field")

#lapply(all.maps, function(x) paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "))
#all.maps[["EPOCA_Svalbard2018"]][,colnames(all.maps[["EPOCA_Svalbard2018"]]) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]

#convert to data frame and export
EMP_metadata_issues.df<-do.call("rbind", EMP_metadata_issues)
EMP_metadata_issues.df$Field<-rep(c("TITLE", 
																		"CONTACTS", 
																		"COLLECTION_DATE_range", 
																		"COLLECTION_DATE_class", 
																		"SAMP_SIZE_range",
																		"DEPTH_range",
																		"DEPTH.class", 
																		"RUN_DATE_range",
																		"RUN_DATE.class",
																		"SEQUENCING_METH", 
																		"LIBRARY_CONSTRUCTION", 
																		"PCR_PRIMERS", 
																		"PLATFORM", 
																		"RUN_CENTER",  
																		"SAMPLE_CENTER", 
																		"SAMPLE_LOCATION", 
																		"TARGET_GENE", 
																		""), 49)

EMP_metadata_issues.df<-EMP_metadata_issues.df[,c("Field", "Value")]

write.csv(EMP_metadata_issues.df, 
					file.path(paste(getwd(), 
													"outputs/EMP_metadata_issues.csv", sep="/")), 
					row.names=FALSE)

#this seemed useful, but in the end not, may want later though
#http://stackoverflow.com/questions/13006909/export-a-list-of-matrices-nicely-to-the-same-worksheet-in-excel

#wb = loadWorkbook("matrix.xlsx", create = TRUE)
# Create a new sheet
#createSheet(wb, name = "mysheet")
# cumulative length (rows) of matrices
# +2 = 1 for list names, 1 for header row
#cumlen = cumsum(c(1, head(sapply(EMP_metadata_issues, nrow), n = -1) + 2))
# Write data rows (implicitly vectorized!)
#writeWorksheet(wb, data = EMP_metadata_issues, sheet = "mysheet", startRow = cumlen + 1, header = FALSE, rownames=lapply(EMP_metadata_issues, rownames))
# Write list names
#writeWorksheet(wb, data = as.list(names(EMP_metadata_issues)), sheet = "mysheet", startRow = cumlen, header = FALSE)
#saveWorkbook(wb)

#check out some really bad ones...
#tibetan_plateau_salt_lake_sediment
head(all.maps[["tibetan_plateau_salt_lake_sediment"]])
head(all.maps[[22]])
colnames(all.maps[["tibetan_plateau_salt_lake_sediment"]])
colnames(all.maps[[1]])
all.maps[["tibetan_plateau_salt_lake_sediment"]]["COLLECTION_DATE"]

#Kilauea geothermal soils and biofilms
head(all.maps[["Kilauea geothermal soils and biofilms"]])
colnames(all.maps[["Kilauea geothermal soils and biofilms"]])

dim(all.maps[["Environmental metagenomic interrogation of Thar desert microbial communities"]])
all.maps[["Environmental metagenomic interrogation of Thar desert microbial communities"]]

colnames(all.maps[["Fermilab_spatial_study"]])
head(all.maps[["Fermilab_spatial_study"]])

colnames(all.maps[["Great Lake Microbiome"]])
head(all.maps[["Great Lake Microbiome"]])

write.csv(head(all.maps[["tibetan_plateau_salt_lake_sediment"]]), 
					file.path(paste(getwd(), 
													"outputs/tibetan_plateau_salt_lake_sediment_head.csv", sep="/")), 
					row.names=FALSE)

write.csv(head(all.maps[["Kilauea geothermal soils and biofilms"]]), 
					file.path(paste(getwd(), 
													"outputs/Kilauea_geothermal_soils_and_biofilms_head.csv", sep="/")), 
					row.names=FALSE)

write.csv(head(all.maps[["Fermilab_spatial_study"]]), 
					file.path(paste(getwd(), 
													"outputs/Fermilab_spatial_study_head.csv", sep="/")), 
					row.names=FALSE)

write.csv(head(all.maps[["Great Lake Microbiome"]]), 
					file.path(paste(getwd(), 
													"outputs/Great_Lake_Microbiome_head.csv", sep="/")), 
					row.names=FALSE)

#try to evaluate all other columns across studies
col.n.comm<-names(which(sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE)<=length(all.maps)/2))

#as all spaces between words are '_' can split all column names
#by '_' and find common words for fields might indicate a pattern to data
col.n.comm.sp<-strsplit(col.n.comm, "_")
head(col.n.comm.sp)
col.n.comm.sp<-unlist(col.n.comm.sp)
col.n.comm.sp<-as.matrix(sort(table(col.n.comm.sp), decreasing=TRUE))
#this is useful...
#frequencies...
table(col.n.comm.sp)
#start with most frequent
col.n.comm.sp[col.n.comm.sp>8, ]
length(col.n.comm.sp[col.n.comm.sp>8, ])

col.n.comm[grep(paste(names(col.n.comm.sp[col.n.comm.sp>8, ]), collapse="|"), col.n.comm, ignore.case=TRUE)]
#99 columns contain the 9 most frequent words

#search colnames
col.n.comm[grep('TOT', col.n.comm, fixed=TRUE)]
#well here's some chemistry
col.chem<-col.n.comm[grep('TOT', col.n.comm, fixed=TRUE)]
col.chem<-col.chem[which(!col.chem %in% c("TOTAL_PHYTOPLANKTON_COUNT", "TOTAL_BACTERIA_COUNT", "DISTANCETOTRANSECTA_CM"))]
#look at metadata with chemistry
chem.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.chem), head(x[,"TITLE"], 1), "No Field"))
chem.maps<-names(chem.maps[which(!chem.maps %in% "No Field")])
lapply(all.maps[chem.maps], function(x) head(x[colnames(x) %in% col.chem], 2))
#not many columns shared across studies, see if can focus
length(chem.maps) #14 studies have any of these columns
as.matrix(table(unlist(lapply(all.maps[chem.maps], function(x) colnames(x[colnames(x) %in% col.chem])))))
#10 have TOT_NITRO, only 7 have TOT_N_METH and only 1 has TOT_NITRO_UNITS
lapply(all.maps[chem.maps], function(x) head(x[colnames(x) %in% c("TOT_NITRO", "TOT_N_METH", "TOT_NITRO_UNITS")], 10))
#well, that instersting, values don't mean much without units and methods...

#export table for EMP chemistry fields
EMP_chem_fq<-as.matrix(table(unlist(lapply(all.maps[chem.maps], function(x) colnames(x[colnames(x) %in% col.chem])))))

EMP_chem_table<-lapply(all.maps[chem.maps], function(x) cbind(
	head(x[,"TITLE"], 1),
	ifelse(TRUE %in% (colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")), 
														 paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "), 
														 "No Contacts")
))

EMP_chem_table<-do.call("rbind", EMP_chem_table)

for(i in 1:length(col.chem)){
	EMP_chem_table<- cbind(EMP_chem_table, unlist(lapply(all.maps[chem.maps], function(x) ifelse(TRUE %in% (colnames(x) %in% col.chem[i]), paste(summary(x[,col.chem[i]]), collapse=" "), "No Field")), use.names=FALSE))
	#colnames(chem.maps.df)[i]<-col.chem[i]
}
colnames(EMP_chem_table)<-c("TITLE", "CONTACTS", col.chem)
EMP_chem_table<-as.data.frame(EMP_chem_table)
#this is actually pretty hard to read, maybe a better way

#well this is the easy way to do what I did above... 
EMP_chem_table<-lapply(all.maps[chem.maps], function(x) colnames(x[colnames(x) %in% col.chem]))
EMP_chem_table<-melt(EMP_chem_table)
EMP_chem_table<-dcast(EMP_chem_table, L1~value)
EMP_chem_table[is.na(EMP_chem_table)]<-"No Field"
#without contacts... but can add
#can't change what is in field though...
EMP_chem_table$CONTACTS<-do.call("rbind", lapply(all.maps[chem.maps], function(x) cbind(
	ifelse(TRUE %in% (colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")), 
				 paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "), 
				 "No Contacts")
)))
#fix title and order
colnames(EMP_chem_table)[1]<-"TITLE"
EMP_chem_table<-EMP_chem_table[, c("TITLE", "CONTACTS", col.chem)]

colnames(EMP_chem_table)


#just nitrogen?
col.n.comm[grep('NITRO', col.n.comm, fixed=TRUE)]

#look at metadata with nitrogen, make special list
col.nitro<-c("TOT_NITRO", "TOT_NITRO_PERCENT", "TOT_NITRO_UNIT", "TOT_NITRO_UNITS", 
						 "TOT_N_METH")

nitro.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.nitro), head(x[,"TITLE"], 1), "No Field"))
nitro.maps<-names(nitro.maps[which(!nitro.maps %in% "No Field")])
lapply(all.maps[nitro.maps], function(x) head(x[colnames(x) %in% col.nitro], 2))
#not many columns shared across studies, see if can focus
length(nitro.maps) #11 studies have some nitrogen
as.matrix(table(unlist(lapply(all.maps[nitro.maps], function(x) colnames(x[colnames(x) %in% col.nitro])))))
#10 have TOT_NITRO, only 7 have TOT_N_METH and only 1 has TOT_NITRO_UNITS

#export table for EMP nitrogen fields
EMP_nitro_fq<-as.matrix(table(unlist(lapply(all.maps[nitro.maps], function(x) colnames(x[colnames(x) %in% col.nitro])))))

EMP_nitro_table<-lapply(all.maps[nitro.maps], function(x) cbind(
	head(x[,"TITLE"], 1),
	ifelse(TRUE %in% (colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")), 
				 paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "), 
				 "No Contacts")
))

EMP_nitro_table<-do.call("rbind", EMP_nitro_table)
EMP_nitro_table<-as.data.frame(EMP_nitro_table)
colnames(EMP_nitro_table)<-c("TITLE", "CONTACTS")

#add NITRO fields
EMP_nitro_table$TOT_NITRO<-unlist(lapply(all.maps[nitro.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "TOT_NITRO"), 
				 paste(paste(range(x[,"TOT_NITRO"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"TOT_NITRO"]))), length(x[,"TOT_NITRO"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_nitro_table$TOT_NITRO_PERCENT<-unlist(lapply(all.maps[nitro.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "TOT_NITRO_PERCENT"),  
				 paste(paste(range(x[,"TOT_NITRO_PERCENT"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"TOT_NITRO_PERCENT"]))), length(x[,"TOT_NITRO_PERCENT"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_nitro_table$TOT_NITRO_UNIT<-unlist(lapply(all.maps[nitro.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "TOT_NITRO_UNIT"), 
				 paste(unique(x[which(!is.na(x[,"TOT_NITRO_UNIT"])),"TOT_NITRO_UNIT"]), collapse=" "), "No Field")), use.names=FALSE)


EMP_nitro_table$TOT_NITRO_UNITS<-unlist(lapply(all.maps[nitro.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "TOT_NITRO_UNITS"), 
				 paste(unique(x[which(!is.na(x[,"TOT_NITRO_UNITS"])),"TOT_NITRO_UNITS"]), collapse=" "), "No Field")), use.names=FALSE)

EMP_nitro_table$TOT_N_METH<-unlist(lapply(all.maps[nitro.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "TOT_N_METH"), 
				 paste(unique(x[which(!is.na(x[,"TOT_N_METH"])),"TOT_N_METH"]), collapse=" "), "No Field")), use.names=FALSE)


#just nitrogen
write.csv(EMP_nitro_table, 
					file.path(paste(getwd(), "outputs/EMP_nitro_table.csv", sep="/")), 
					row.names=FALSE)

#other chemistry besides nitrogen and carbon
col.chem<-colnames(EMP_chem_table)[which(!colnames(EMP_chem_table) %in% colnames(EMP_nitro_table))]
col.chem<-col.chem[which(!col.chem %in% col.chem[grep("CARB", col.chem)])]
col.chem<-col.chem[which(!col.chem %in% c("TOT_DEPTH_WATER_COL", "TOT_ORG_C_METH"))]
#look at metadata with chemistry
chem.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.chem), head(x[,"TITLE"], 1), "No Field"))
chem.maps<-names(chem.maps[which(!chem.maps %in% "No Field")])
lapply(all.maps[chem.maps], function(x) head(x[colnames(x) %in% col.chem], 2))
#not many columns shared across studies, see if can focus
length(chem.maps) #3 studies have any of these columns
as.matrix(table(unlist(lapply(all.maps[chem.maps], function(x) colnames(x[colnames(x) %in% col.chem])))))
#not a lot of repeats
lapply(all.maps[chem.maps], function(x) head(x[colnames(x) %in% col.chem], 10))
#well, that instersting, values don't mean much without units and methods...

#all other chemistry fields from Environmental metagenomic interrogation of Thar desert microbial communities
#note TOT_MASS and TOT_MASS_UNIT from `Green Iguana hindgut microbiome`

#search colnames
col.n.comm[grep('DATE', col.n.comm, fixed=TRUE)]
#yup, dates
col.date<-col.n.comm[grep('DATE', col.n.comm, fixed=TRUE)]
date.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.date), head(x[,"TITLE"], 1), "No Field"))
date.maps<-names(date.maps[which(!date.maps %in% "No Field")])
length(date.maps)
#well, all 49 have some date field (knew that)
sort(table(unlist(lapply(all.maps[date.maps], function(x) colnames(x[colnames(x) %in% col.date])))), decreasing=TRUE)
#well, basically only COLLECTION_DATE and RUN_DATE, most of the rest look like Juan's

#search colnames
col.n.comm[grep('METH', col.n.comm, fixed=TRUE)]
#pretty diverse, already looked at SEQUENCING_METH and chemistry *METH
#some soil ones are ineresting
#WATER_CONTENT_SOIL_METH, SOIL_TYPE_METH, TEXTURE_METH
sort(table(unlist(lapply(all.maps, function(x) colnames(x[colnames(x) %in% c("WATER_CONTENT_SOIL_METH", "SOIL_TYPE_METH", "TEXTURE_METH")])))), decreasing=TRUE)
lapply(all.maps, function(x) head(x[colnames(x) %in% c("WATER_CONTENT_SOIL_METH", "SOIL_TYPE_METH", "TEXTURE_METH")]))
#well it's something...

#search colnames
col.n.comm[grep('SOIL', col.n.comm, fixed=TRUE)]
#I like soil...
col.soil<-col.n.comm[grep('SOIL', col.n.comm, fixed=TRUE)]
soil.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.soil), head(x[,"TITLE"], 1), "No Field"))
soil.maps<-names(soil.maps[which(!soil.maps %in% "No Field")])
length(soil.maps)
#hmm, 9 studies have any fields with soil...
as.matrix(sort(table(unlist(lapply(all.maps[soil.maps], function(x) colnames(x[colnames(x) %in% col.soil])))), decreasing=TRUE))
#and most of them only occur once 
#again most common is WATER_CONTENT_SOIL (7), but only 4 have a method...
lapply(all.maps[soil.maps], function(x) head(x[colnames(x) %in% col.soil], 2))

EMP_soil_table<-lapply(all.maps[soil.maps], function(x) cbind(
	head(x[,"TITLE"], 1),
	ifelse(TRUE %in% (colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")), 
				 paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "), 
				 "No Contacts")
))

EMP_soil_table<-do.call("rbind", EMP_soil_table)
EMP_soil_table<-as.data.frame(EMP_soil_table)
colnames(EMP_soil_table)<-c("TITLE", "CONTACTS")
#add soil fields
EMP_soil_table$WATER_CONTENT_SOIL<-unlist(lapply(all.maps[soil.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "WATER_CONTENT_SOIL"),  
				 paste(paste(range(x[,"WATER_CONTENT_SOIL"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"WATER_CONTENT_SOIL"]))), length(x[,"WATER_CONTENT_SOIL"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_soil_table$WATER_CONTENT_SOIL_METH<-unlist(lapply(all.maps[soil.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "WATER_CONTENT_SOIL_METH"), 
				 paste(unique(x[which(!is.na(x[,"WATER_CONTENT_SOIL_METH"])),"WATER_CONTENT_SOIL_METH"]), collapse=" "), "No Field")), use.names=FALSE)

EMP_soil_table$SOIL_TYPE<-unlist(lapply(all.maps[soil.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "SOIL_TYPE"), 
				 paste(unique(x[which(!is.na(x[,"SOIL_TYPE"])),"SOIL_TYPE"]), collapse=" "), "No Field")), use.names=FALSE)

EMP_soil_table$SOIL_TYPE_METH<-unlist(lapply(all.maps[soil.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "SOIL_TYPE_METH"), 
				 paste(unique(x[which(!is.na(x[,"SOIL_TYPE_METH"])),"SOIL_TYPE_METH"]), collapse=" "), "No Field")), use.names=FALSE)

EMP_soil_table$WATER_CONTENT_SOIL_UNIT<-unlist(lapply(all.maps[soil.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "WATER_CONTENT_SOIL_UNIT"), 
				 paste(unique(x[which(!is.na(x[,"WATER_CONTENT_SOIL_UNIT"])),"WATER_CONTENT_SOIL_UNIT"]), collapse=" "), "No Field")), use.names=FALSE)

EMP_soil_table$SOIL_TEMP<-unlist(lapply(all.maps[soil.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "SOIL_TEMP"),  
				 paste(paste(range(x[,"SOIL_TEMP"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"SOIL_TEMP"]))), length(x[,"SOIL_TEMP"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

write.csv(EMP_soil_table, 
					file.path(paste(getwd(), "outputs/EMP_soil_table.csv", sep="/")), 
					row.names=FALSE)

#search colnames
col.n.comm[grep('CARB', col.n.comm, fixed=TRUE)]
#hmm, maybe TOT_ORG_CAB and ORG_CARB are duplicates
col.carb<-col.n.comm[grep('CARB', col.n.comm, fixed=TRUE)]
#add method
col.carb<-c(col.carb, "TOT_ORG_C_METH")
carb.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.carb), head(x[,"TITLE"], 1), "No Field"))
carb.maps<-names(carb.maps[which(!carb.maps %in% "No Field")])
length(carb.maps)
#15 studies have fields related to carbon
as.matrix(sort(table(unlist(lapply(all.maps[carb.maps], function(x) colnames(x[colnames(x) %in% col.carb])))), decreasing=TRUE))
#the most common is TOT_ORG_CARB (8), however only 1 TOT_ORG_CARB_UNITS and no methods...
lapply(all.maps[carb.maps], function(x) head(x[colnames(x) %in% col.carb], 2))

#export carbon table
EMP_carb_table<-lapply(all.maps[carb.maps], function(x) cbind(
	head(x[,"TITLE"], 1),
	ifelse(TRUE %in% (colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")), 
				 paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "), 
				 "No Contacts")
))

EMP_carb_table<-do.call("rbind", EMP_carb_table)
EMP_carb_table<-as.data.frame(EMP_carb_table)
colnames(EMP_carb_table)<-c("TITLE", "CONTACTS")
#add carb fields
EMP_carb_table$TOT_CARB<-unlist(lapply(all.maps[carb.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "TOT_CARB"),  
				 paste(paste(range(x[,"TOT_CARB"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"TOT_CARB"]))), length(x[,"TOT_CARB"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_carb_table$TOT_ORG_CARB<-unlist(lapply(all.maps[carb.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "TOT_ORG_CARB"),  
				 paste(paste(range(x[,"TOT_ORG_CARB"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"TOT_ORG_CARB"]))), length(x[,"TOT_ORG_CARB"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_carb_table$ORG_CARB<-unlist(lapply(all.maps[carb.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "ORG_CARB"),  
				 paste(paste(range(x[,"ORG_CARB"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"ORG_CARB"]))), length(x[,"ORG_CARB"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_carb_table$TOTAL_INORG_CARB<-unlist(lapply(all.maps[carb.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "TOTAL_INORG_CARB"),  
				 paste(paste(range(x[,"TOTAL_INORG_CARB"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"TOTAL_INORG_CARB"]))), length(x[,"TOTAL_INORG_CARB"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_carb_table$NITRO_ORG_CARB_UNIT<-unlist(lapply(all.maps[carb.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "NITRO_ORG_CARB_UNIT"), 
				 paste(unique(x[which(!is.na(x[,"NITRO_ORG_CARB_UNIT"])),"NITRO_ORG_CARB_UNIT"]), collapse=" "), "No Field")), use.names=FALSE)

EMP_carb_table$TOT_ORG_CARB_UNITS<-unlist(lapply(all.maps[carb.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "TOT_ORG_CARB_UNITS"), 
				 paste(unique(x[which(!is.na(x[,"TOT_ORG_CARB_UNITS"])),"TOT_ORG_CARB_UNITS"]), collapse=" "), "No Field")), use.names=FALSE)

EMP_carb_table$CARB_NITRO_RATIO<-unlist(lapply(all.maps[carb.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "CARB_NITRO_RATIO"),  
				 paste(paste(range(x[,"CARB_NITRO_RATIO"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"CARB_NITRO_RATIO"]))), length(x[,"CARB_NITRO_RATIO"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_carb_table$TOT_ORG_C_METH<-unlist(lapply(all.maps[carb.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "TOT_ORG_C_METH"), 
				 paste(unique(x[which(!is.na(x[,"TOT_ORG_C_METH"])),"TOT_ORG_C_METH"]), collapse=" "), "No Field")), use.names=FALSE)

write.csv(EMP_carb_table, 
					file.path(paste(getwd(), "outputs/EMP_carb_table.csv", sep="/")), 
					row.names=FALSE)

#search colnames
col.n.comm[grep('COUNT', col.n.comm, fixed=FALSE)]
#note "COUNTY" is less frequent, not "COUNTRY"
#includes many microbial counts, maybe interesting, but not right now...
col.count<-col.n.comm[grep('COUNT', col.n.comm, fixed=TRUE)]
count.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.count), head(x[,"TITLE"], 1), "No Field"))
count.maps<-names(count.maps[which(!count.maps %in% "No Field")])
length(count.maps)
#15 studies have fields related to counton
as.matrix(sort(table(unlist(lapply(all.maps[count.maps], function(x) colnames(x[colnames(x) %in% col.count])))), decreasing=TRUE))
#the most common is TOT_ORG_count (8), however only 1 TOT_ORG_count_UNITS and no methods...
lapply(all.maps[count.maps], function(x) head(x[colnames(x) %in% col.count], 2))

#search colnames
col.n.comm[grep('TEMP', col.n.comm, fixed=FALSE)]
#hey, here are some potential duplicates
col.temp<-col.n.comm[grep('TEMP', col.n.comm, fixed=FALSE)]
temp.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.temp), head(x[,"TITLE"], 1), "No Field"))
temp.maps<-names(temp.maps[which(!temp.maps %in% "No Field")])
length(temp.maps)
#20 studies have temp data
as.matrix(sort(table(unlist(lapply(all.maps[temp.maps], function(x) colnames(x[colnames(x) %in% col.temp])))), decreasing=TRUE))
#what is differenve between "TEMP" and all other "*TEMP*"
#ANNUAL_SEASON_TEMP vs. ANNUAL_TEMP vs. SEASON_TEMP
#SOIL_TEMP vs. MEAN_SOIL_TEMP_DAY
lapply(all.maps[temp.maps], function(x) head(x[colnames(x) %in% col.temp], 2))

EMP_temp_table<-lapply(all.maps[temp.maps], function(x) cbind(
	head(x[,"TITLE"], 1),
	ifelse(TRUE %in% (colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")), 
				 paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "), 
				 "No Contacts")
))

EMP_temp_table<-do.call("rbind", EMP_temp_table)
EMP_temp_table<-as.data.frame(EMP_temp_table)

######### shoot, if any NA's this returns range NA to NA...
#for(i in 1:length(col.temp)){
#	EMP_temp_table<- cbind(EMP_temp_table, unlist(lapply(all.maps[temp.maps], function(x) ifelse(TRUE %in% (colnames(x) %in% col.temp[i]), paste(range(x[,col.temp[i]]), collapse=" to "), "No Field")), use.names=FALSE))
#	#colnames(chem.maps.df)[i]<-col.chem[i]
#}

EMP_temp_table$TEMP<-unlist(lapply(all.maps[temp.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "TEMP"),  
				 paste(paste(range(x[,"TEMP"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"TEMP"]))), length(x[,"TEMP"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_temp_table$ANNUAL_SEASON_TEMP<-unlist(lapply(all.maps[temp.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "ANNUAL_SEASON_TEMP"),  
				 paste(paste(range(x[,"ANNUAL_SEASON_TEMP"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"ANNUAL_SEASON_TEMP"]))), length(x[,"ANNUAL_SEASON_TEMP"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_temp_table$AIR_TEMP<-unlist(lapply(all.maps[temp.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "AIR_TEMP"),  
				 paste(paste(range(x[,"AIR_TEMP"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"AIR_TEMP"]))), length(x[,"AIR_TEMP"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_temp_table$SOIL_TEMP<-unlist(lapply(all.maps[temp.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "SOIL_TEMP"),  
				 paste(paste(range(x[,"SOIL_TEMP"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"SOIL_TEMP"]))), length(x[,"SOIL_TEMP"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)


EMP_temp_table<-as.data.frame(EMP_temp_table)
colnames(EMP_temp_table)<-c("TITLE", "CONTACTS", "TEMP", "ANNUAL_SEASON_TEMP", "AIR_TEMP", "SOIL_TEMP")

#add temp fields

write.csv(EMP_temp_table, 
					file.path(paste(getwd(), "outputs/EMP_temp_table.csv", sep="/")), 
					row.names=FALSE)

#search colnames
col.n.comm[grep('SEASON', col.n.comm, fixed=FALSE)]
#some duplicates from TEMP, and other climate data...
col.season<-col.n.comm[grep('SEASON', col.n.comm, fixed=FALSE)]
season.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.season), head(x[,"TITLE"], 1), "No Field"))
season.maps<-names(season.maps[which(!season.maps %in% "No Field")])
length(season.maps)
#8 studies have season data
as.matrix(sort(table(unlist(lapply(all.maps[season.maps], function(x) colnames(x[colnames(x) %in% col.season])))), decreasing=TRUE))
lapply(all.maps[season.maps], function(x) head(x[colnames(x) %in% col.season], 2))
#wonder what the methods were to find these seasonal values, but interesting

#create table for export
EMP_season_table<-lapply(all.maps[season.maps], function(x) cbind(
	head(x[,"TITLE"], 1),
	ifelse(TRUE %in% (colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")), 
				 paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "), 
				 "No Contacts")
))

EMP_season_table<-do.call("rbind", EMP_season_table)

for(i in 1:length(col.season)){
	EMP_season_table<- cbind(EMP_season_table, unlist(lapply(all.maps[season.maps], function(x) ifelse(TRUE %in% (colnames(x) %in% col.season[i]), paste(range(x[,col.season[i]]), collapse=" to "), "No Field")), use.names=FALSE))
	#colnames(chem.maps.df)[i]<-col.chem[i]
}
EMP_season_table<-as.data.frame(EMP_season_table)
colnames(EMP_season_table)<-c("TITLE", "CONTACTS", col.season)

#add SEASON_ENVIRONMENT fields
EMP_season_table$SEASON_ENVIRONMENT<-unlist(lapply(all.maps[season.maps], function(x) ifelse(TRUE %in% (colnames(x) %in% "SEASON_ENVIRONMENT"), paste(head(x[,"SEASON_ENVIRONMENT"], 1), collapse=" "), "No Field")), use.names=FALSE)
	
write.csv(EMP_season_table, 
					file.path(paste(getwd(), "outputs/EMP_season_table.csv", sep="/")), 
					row.names=FALSE)

#search colnames
col.n.comm[grep('UNIT', col.n.comm, fixed=FALSE)]
#some duplicates from TEMP, and other climate data...
col.unit<-col.n.comm[grep('UNIT', col.n.comm, fixed=FALSE)]
unit.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.unit), head(x[,"TITLE"], 1), "No Field"))
unit.maps<-names(unit.maps[which(!unit.maps %in% "No Field")])
length(unit.maps)
#8 studies have unit data
as.matrix(sort(table(unlist(lapply(all.maps[unit.maps], function(x) colnames(x[colnames(x) %in% col.unit])))), decreasing=TRUE))
lapply(all.maps[unit.maps], function(x) head(x[colnames(x) %in% col.unit], 2))
#can at least find sudies that have units...

#manually adding fields to the 
#Summary of potential metadata issues across studies table

#inspect to make sure problems exist 
head(all.maps[["Friedman_alaska_peat_soils"]])
summary(is.na(all.maps[["Friedman_alaska_peat_soils"]][["TEMP"]]))

all.maps[["Friedman_alaska_peat_soils"]][["TOT_N_METH"]]
range(all.maps[["Bergen Ocean Acidification Mesocosms"]][["ORG_NITRO"]], finite=TRUE)
summary(is.na(all.maps[["Bergen Ocean Acidification Mesocosms"]][["ORG_NITRO"]]))