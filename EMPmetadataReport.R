#restart the metdata analysis using an additional 15 studies
#also very cleaned up script as of 8/27/2013
#look out for issues with additional studies and/or 
#different order in all.maps list

# library(BiocInstaller)
# library(XML)
# library(XMLSchema)
# library(SSOAP)
# library(rols)
# library(reshape2)
# library(plyr)
library(rgeos)
library(rgdal)

#set working directory
setwd("~/EarthMicrobiomeProject/R/EMPmetadataReport")

#save defualt plotting options
def.par <- par(no.readonly = TRUE) # save default, for resetting...

#list mapping files
map.file.names<-list.files(path.expand("~/EarthMicrobiomeProject/QIIME_metadata_download"), pattern="_map.txt", recursive=TRUE)

#create list of data frames for each mapping file
all.maps<-list()

#import all mapping files, 
#place each one as a seperate data frame in the list 'all.maps'
#the arguement stringsAsFactors=False 
#removes issues of character -> factor and vice versa issues
#does not resolve integer, numeric, date -> charcter etc...
all.maps<-lapply(map.file.names, function(x) read.delim(file.path(paste(path.expand("~/EarthMicrobiomeProject/QIIME_metadata_download"), x, sep="/")), quote="", stringsAsFactors=FALSE))

# this is a slicker way to import all the mapping files than old line
# but changes the order of the 'all.maps' list
# some lines are (unfortunately) hard coded to this order… so look out!
	
#check merge
emp.map<-Reduce(function(x, y) merge(x, y, all=TRUE), all.maps)
#, so character/ factor resolved
#there are not 17072 obs x 539 fields...

#explore NA
#mean NA?
sort(unlist(lapply(all.maps, function(x) mean(is.na(x)))))
#some are >25%?
which(unlist(lapply(all.maps, function(x) mean(is.na(x))))>0.24)

#take a look at those ones...
head(all.maps[[7]])
dim(all.maps[[7]])
#which columns have the most NAs?
sort(apply(all.maps[[7]], 2, function(x){length(which(is.na(x)))}))

head(all.maps[[28]])
dim(all.maps[[28]])
sort(apply(all.maps[[28]], 2, function(x){length(which(is.na(x)))}))

#can get study titles
lapply(all.maps, function(x) unique(x[,"TITLE"]))
#some studies have more than one unique  title, ID those
which(lapply(all.maps, function(x) length(unique(x[,"TITLE"])))>1)
#3

#fix titles
sort(all.maps[[40]][,"TITLE"])
#so this is an Excel 'draw down error', 
#looks like started with "EPOCA_Svalbard2018" 
all.maps[[40]][,"TITLE"]<-"EPOCA_Svalbard2018"

#next title issue
sort(all.maps[[54]][,"TITLE"])
#also this is an Excel 'draw down error', 
#started with "Intertidal microbes 16s for 2009 and 2010"
all.maps[[54]][,"TITLE"]<-"Intertidal microbes 16s for 2009 and 2010"

#next title issue
sort(all.maps[[63]][,"TITLE"])
#this is more complicated
table(all.maps[[63]][,"TITLE"])
unique(all.maps[[63]][,"TITLE"])
#this mapping file appears to contain 11 different studies  
#useful to have only one title for labeling purposes, 
#can add column of id another column to use...
sort(colnames(all.maps[[63]]))
#called "Thomas_sponge_communities" on http://microbio.me/qiime
head(all.maps[[63]][grep("Thomas_sponge_communities", all.maps[[63]])])
#ah, just switch TITLE and ExPERIMENT_TITLE
which(colnames(all.maps[[63]]) %in% c("TITLE", "EXPERIMENT_TITLE"))
colnames(all.maps[[63]])[c(8,32)]<-c("EXPERIMENT_TITLE", "TITLE")
#check
unique(all.maps[[63]][,"TITLE"])
unique(all.maps[[63]][,"EXPERIMENT_TITLE"])
#ok, that worked

#check titles again
unlist(lapply(all.maps, function(x) unique(x[,"TITLE"])))
#now have one unique title per study

#name studies in list by unique title
names(all.maps)<-unlist(lapply(all.maps, function(x) unique(x[,"TITLE"])))
#note that titles are different from the name of the QIIME file
qiime.file.names<-list.dirs(path.expand("~/EarthMicrobiomeProject/QIIME_metadata_download"))
qiime.file.names<-qiime.file.names[-1]
qiime.file.names<-lapply(qiime.file.names, function(x)
	substr(x, nchar("C:/Users/asus4/Documents/EarthMicrobiomeProject/QIIME_metadata_download//"), nchar(x)))								
qiime.file.names<-unlist(qiime.file.names)

#Earth Microbiome Project (EMP) metadata summary
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

#get list of all unique columns
all.col<-unique(unlist(lapply(all.maps, colnames)))

#search colnames 
all.col[grep('elev', all.col, ignore.case=TRUE)]

#create a sorted table of column name frequency 
#(how many studies have this metadata field)
study_column_table<-sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE)
#convert from array to data frame for ease of use and exporting as table
study_column_table<-data.frame(study_column_table)
#column names are currently row names, add as column
study_column_table$COLUMN_NAMES<-rownames(study_column_table)
#rename count column
colnames(study_column_table)[1]<-"No_STUDIES"
#replace row names with # (easier to read)
rownames(study_column_table)<-as.character(1:nrow(study_column_table))
#reorder
study_column_table<-study_column_table[,c(2,1)]
#subset to fit on MSword page, place side by side
study_column_table<-cbind(study_column_table[1:38,], study_column_table[39:76,])

#add a little note, how many columns not shown...
#don't have to hard code this, it is based on total number of columns and 
#number of rows in table with blanks added to make the table readily presentable
study_column_table<-rbind(study_column_table, 
													c(paste("*", 
																	paste(as.character(length(unique(unlist(lapply(all.maps, colnames))))-nrow(study_column_table)*2), 
																				"columns not shown", sep=" "), 
																	sep=""), "", "", ""))

#export table
write.csv(study_column_table, 
					file.path(paste(getwd(), 
													"outputs/study_column_table.csv", sep="/")), 
					row.names=FALSE)

#how many fields are only in a single study?
summary(data.frame(sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE))<2)

#NA values in core items of the MIMARKS checklists Environment 
missing.crit<-lapply(1:length(all.maps), function(x) all.maps[[x]][which(is.na(all.maps[[x]][["LATITUDE"]]) |
																																				 is.na(all.maps[[x]][["LONGITUDE"]]) |
																																				 is.na(all.maps[[x]][["DEPTH"]]) |
																																				 is.na(all.maps[[x]][["ALTITUDE"]]) |	
																																				 is.na(all.maps[[x]][["COLLECTION_DATE"]]) |
																																				 is.na(all.maps[[x]][["ENV_FEATURE"]]) |
																																				 is.na(all.maps[[x]][["ENV_BIOME"]]) |
																																				 is.na(all.maps[[x]][["ENV_MATTER"]])), 
	c("X.SampleID", "TITLE", "LATITUDE", "LONGITUDE", "DEPTH", "ALTITUDE", "COLLECTION_DATE", 'ENV_FEATURE', 'ENV_BIOME', "ENV_MATTER")])

#convert to data frame
missing.crit<-do.call("rbind", missing.crit)
#aggregate by study title and count no. of NA's in all critical fields
missing.crit.tab<-aggregate(missing.crit, by=list(missing.crit$TITLE), FUN=function(x) length(which(is.na(x))))
#check
missing.crit[missing.crit$TITLE == "Bird Egg Shells from Spain", "LATITUDE"]

#clean
colnames(missing.crit.tab)
missing.crit.tab<-missing.crit.tab[,-which(colnames(missing.crit.tab) %in% c("X.SampleID", "TITLE"))]
colnames(missing.crit.tab)[1]<-"TITLE"

#add ELEVATION because not all studies have it, need to account for that
#gotta make a temp list (subset studies with NA)
tmp<-all.maps[missing.crit.tab$TITLE]
#add the ELEVATION field
missing.crit.tab$ELEVATION<-unlist(lapply(1:length(tmp), function(x) ifelse(isTRUE("ELEVATION" %in% colnames(tmp[[x]])), length(which(is.na(tmp[[x]][["ELEVATION"]]))), "No Field")))
#add number of samples in that study
missing.crit.tab$No.Samples<-unlist(lapply(tmp, nrow))

#which study does not have ELEVATION?
unlist(lapply(all.maps, function(x) which(!"ELEVATION"%in% colnames(x))))
nrow(all.maps[["Comparison of groundwater samples from karst sinkholes (cenotes) from the Yucatan Peninsula, Mexico"]])
#manually add this to table

#export
write.csv(missing.crit.tab, 
					file.path(paste(getwd(), 
													"outputs/missing_crit.csv", sep="/")), 
					row.names=FALSE)

###################
#Geographic Data
#check lat long coordinates against COUNTRY field
library(rgdal)
#import natural earth administrative boundaries (country borders map)
borders<-readOGR(dsn="C:/Users/asus4/Documents/GIS/Data/Natural_Earth", 
								 layer="ne_50m_admin_0_countries")

#subset samples with lat/long (according to missing.crit should be all but 31)
#use the emp.map data frame (not good for investigating data, but should be ok for plotting coordinates)
emp.gis<-subset(emp.map, complete.cases(emp.map[,c("LONGITUDE", "LATITUDE")]))
nrow(emp.map)-nrow(emp.gis)

#check class, make sure lat/long usable
class(emp.gis[, "LONGITUDE"])
class(emp.gis[, "LATITUDE"])

#remove "GAZ:" from COUNTRY
emp.gis$COUNTRY<-substr(emp.gis$COUNTRY, 5, nchar(as.character(emp.gis$COUNTRY)))

#convert to SpatialPointsDataFrame
coordinates(emp.gis) = c("LONGITUDE", "LATITUDE")

#need to check and set spatial projection
proj4string(borders)
proj4string(emp.gis)
proj4string(emp.gis)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#overlay EMP spatial points on country map
emp.border<-over(emp.gis, borders, returnList=FALSE)

#check overlay
colnames(data.frame(emp.border))
dim(data.frame(emp.border))

summary(emp.border$admin)
summary(emp.gis$COUNTRY == emp.border$admin)
#and 1745 FALSE and 2196 NA's

#so the NA's are points not in polygons (off land)
#the FALSE are points whose country does not match the name of the polygon

#try to match NA to nearest country
emp.gis.na<-data.frame(emp.gis[which(is.na(emp.border$admin)), ]) 
coordinates(emp.gis.na)= c("LONGITUDE", "LATITUDE")
proj4string(emp.gis.na)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#this was helpful https://stat.ethz.ch/pipermail/r-sig-geo/2011-June/012008.html

#rgeos gives warnings that object is not projected
#this would be an issue, except I am just interetsted in the clostest
#polygon, the measure does not need to be accurate 
is.projected(emp.gis.na)
is.projected(borders)
#not projected

#get the distance of each point to polygons
na.dist<- gDistance(emp.gis.na, borders, byid=TRUE)
class(na.dist)
dim(na.dist)
dim(borders)
#so there is a measure from each point to every polygon
attributes(na.dist)
#so polygons named 0:240, same as rownames of borders
#points labeled 1:2196, same as rownames for data.frame(emp.gis.na)

#find which polygon is closest to points
na.dist.min<- apply(na.dist, 2, function(x) which(x==min(x)))
#convert to regular data frame
na.dist.min<-data.frame(na.dist.min)

dim(na.dist.min)
dimnames(na.dist.min)
#now rows are points and the value in only column is rowname of country polygon
#match the row name of polygon to the min distance to point

#convert to data.frame without factors
borders.df <- data.frame(lapply(data.frame(borders), as.character), stringsAsFactors=FALSE)

#test
borders.df[rownames(borders.df) %in% as.character(na.dist.min[1,]), "admin"]
borders.df["77", "admin"]
borders.df[, "admin"]

emp.gis.na.admin<-lapply(1:nrow(na.dist.min), function(x) borders.df[rownames(borders.df) %in% as.character(na.dist.min[x,]), "admin"])
head(emp.gis.na.admin)
#looks good

#collapse list
emp.gis.na.admin<-do.call("rbind", emp.gis.na.admin)
head(emp.gis.na.admin)
dim(emp.gis.na.admin)
#so each row should correspond to a point that was NA and the Country is nearest polygon

#and compare
emp.gis.na.df <- data.frame(lapply(data.frame(emp.gis.na), as.character), stringsAsFactors=FALSE)

summary(as.character(emp.gis.na.admin) == as.character(emp.gis.na$COUNTRY))
#so 1872 COUNTRY fields correclty match the nearest country polygon
#324 do not
unique(cbind(as.character(emp.gis.na.admin), as.character(emp.gis.na$COUNTRY)))

#which studies are these 324 samples in?
unique(emp.gis.na.df[which(!as.character(emp.gis.na$COUNTRY) == as.character(emp.gis.na.admin)), "TITLE"])
#3 studies have lat/long coordinates that do not place them nearest 
#to the country border in the COUNTRY field
unique(cbind(emp.gis.na.df[which(!as.character(emp.gis.na$COUNTRY) == as.character(emp.gis.na.admin)), "COUNTRY"],
						 emp.gis.na.admin[which(!as.character(emp.gis.na$COUNTRY) == as.character(emp.gis.na.admin)), ]))
#ah, and find what the mismatches are, only 13 and somewhat reasonable

#now investigate points where lat/long does not place point in
#country specified in COUNTRY field
emp.gis.f<-data.frame(emp.gis)
emp.gis.f$ne.admin<-emp.border$admin
#subset the rows where COUNTRY does not = admin
emp.gis.f<-emp.gis.f[which(!as.character(emp.gis.f$COUNTRY) == as.character(emp.gis.f$ne.admin)), ]
#correct number according to original summary
dim(emp.gis.f)
unique(emp.gis.f[, "TITLE"])
unique(emp.gis.f[, c("COUNTRY", "ne.admin")])
#some reasonable Scotland vs. United Kindtom, 
#Republic of South Africa vs. South Africa
#some unreasonable United States of America vs. China
#Unisted States of America vs. Mexico
unique(emp.gis.f[, c("TITLE", "COUNTRY", "ne.admin")])
#United States of America vs. China (the Oregon Alder_Fir study...)

#subset just problem ones
emp.gis.f<-emp.gis.f[c(which(emp.gis.f$COUNTRY=="United States of America" & emp.gis.f$ne.admin =="Russia"),
						which(emp.gis.f$COUNTRY=="United States of America" & emp.gis.f$ne.admin =="China"),
						which(emp.gis.f$COUNTRY=="United States of America" & emp.gis.f$ne.admin =="Mexico")), ]
dim(emp.gis.f)
#inspect
emp.gis.f[emp.gis.f$TITLE=="Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample", ]


#export to GIS to investigate
emp.gis.wrong<-data.frame(emp.gis[which(!emp.gis$COUNTRY == emp.border$admin | is.na(emp.border$admin)), ]) 
emp.gis.wrong<-emp.gis.wrong[,c("X.SampleID", "TITLE", "LATITUDE", "LONGITUDE", "COUNTRY","Description")]
emp.gis.wrong<-unique(emp.gis.wrong)
dim(emp.gis.wrong)
class(emp.gis.wrong)

#add the matching natural earth admin column
emp.gis.wrong$ne_admin<-emp.border[which(!emp.gis$COUNTRY == emp.border$admin | is.na(emp.border$admin)), "admin"]
coordinates(emp.gis.wrong)<-c("LONGITUDE", "LATITUDE")

#export these points to 
writeOGR(emp.gis.wrong, getwd(), "emp_gis_wrong", driver="ESRI Shapefile", overwrite_layer=TRUE)

#export just studies where NA points COUNTRY does not match nearest 
#country border
emp.gis.na.bad<-emp.gis.na.df[which(!as.character(emp.gis.na$COUNTRY) == as.character(emp.gis.na.admin)), c("X.SampleID", "TITLE", "LATITUDE", "LONGITUDE", "COUNTRY","Description")]
emp.gis.na.bad$ne.admin<-emp.gis.na.admin[which(!as.character(emp.gis.na$COUNTRY) == as.character(emp.gis.na.admin)), ]
dim(emp.gis.na.bad)
#somewhere lat/long got converted to numeric, fix now but 
#should fix where that happened later
emp.gis.na.bad$LATITUDE<-as.numeric(emp.gis.na.bad$LATITUDE)
emp.gis.na.bad$LONGITUDE<-as.numeric(emp.gis.na.bad$LONGITUDE)
#convert to sp class
coordinates(emp.gis.na.bad)<-c("LONGITUDE", "LATITUDE")
#export these points to 
writeOGR(emp.gis.na.bad, getwd(), "emp_gis_NA_bad", driver="ESRI Shapefile", overwrite_layer=TRUE)

#export just studies where lat/long matches country NOT in COUNTRY field
emp.gis.f<-emp.gis.f[, c("X.SampleID", "TITLE", "LATITUDE", "LONGITUDE", "COUNTRY","Description", "ne.admin")]
coordinates(emp.gis.f)<-c("LONGITUDE", "LATITUDE")
writeOGR(emp.gis.f, getwd(), "emp_gis_f", driver="ESRI Shapefile", overwrite_layer=TRUE)

#make summary table of GIS issues
#first NA for matching COUNTRY
emp.gis.na.bad.df<-data.frame(emp.gis.na.bad)
head(emp.gis.na.bad.tab)
colnames(emp.gis.na.bad.tab)
unique(emp.gis.na.bad.tab[,c("TITLE", "COUNTRY", "ne.admin", "Description")])
emp.gis.na.bad.tab<-unique(emp.gis.na.bad.tab[,c("TITLE", "COUNTRY", "ne.admin", "Description")])

#add no. of samples with one of these unique combinations?
emp.gis.na.bad.tab$No.Samples<-unlist(lapply(1:nrow(emp.gis.na.bad.tab), function(x) nrow(emp.gis.na.bad.df[which(emp.gis.na.bad.df$TITLE == emp.gis.na.bad.tab[x, "TITLE"] &
																																																										emp.gis.na.bad.df$COUNTRY == emp.gis.na.bad.tab[x, "COUNTRY"] &
																																																										emp.gis.na.bad.df$ne.admin == emp.gis.na.bad.tab[x, "ne.admin"] &
																																																										emp.gis.na.bad.df$Description == emp.gis.na.bad.tab[x, "Description"]), ])))
#total no. of samples in that study
emp.gis.na.bad.tab$Tot_No._Study_Samples<-unlist(lapply(1:nrow(emp.gis.na.bad.tab), function(x) nrow(emp.map[which(emp.map$TITLE == emp.gis.na.bad.tab[x, "TITLE"]), ])))
#write table
write.csv(emp.gis.na.bad.tab, 
					file.path(paste(getwd(), 
													"outputs/emp_gis_na_bad_tab.csv", sep="/")), 
					row.names=FALSE)

#table of mismatch between COUNTRY and country of lat/long coordinate 
emp.gis.f.df<-data.frame(emp.gis.f)
colnames(emp.gis.f.df)
emp.gis.f.tab<-unique(emp.gis.f.df[, c("TITLE", "COUNTRY", "Description", "ne.admin")])
#No. of samples with mismatch
emp.gis.f.tab$No.Samples<-unlist(lapply(1:nrow(emp.gis.f.tab), function(x) nrow(emp.gis.f.df[which(emp.gis.f.df$TITLE == emp.gis.f.tab[x, "TITLE"] &
																																																					emp.gis.f.df$COUNTRY == emp.gis.f.tab[x, "COUNTRY"] &
																																																					emp.gis.f.df$ne.admin == emp.gis.f.tab[x, "ne.admin"] &
																																																					emp.gis.f.df$Description == emp.gis.f.tab[x, "Description"]), ])))
#No. of samples in study
emp.gis.f.tab$Tot_No._Study_Samples<-unlist(lapply(1:nrow(emp.gis.f.tab), function(x) nrow(emp.map[which(emp.map$TITLE == emp.gis.f.tab[x, "TITLE"]), ])))
#write table
write.csv(emp.gis.f.tab, 
					file.path(paste(getwd(), 
													"outputs/emp_gis_f_tab.csv", sep="/")), 
					row.names=FALSE)

##############################################
#go into more detail
#pick through each study and find which common columns have different classes
#had written a really ugly nested loop to do this and it was not that helpful
#below is much cleaner and more helpful

#investigate all of the most common fields
#select fields present in greater than half of the studies
col.comm<-names(which(sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE)>length(all.maps)/2))

#add contact fields (mostly missing...?)
col.comm<-c(col.comm, "PRINCIPAL_INVESTIGATOR_CONTACT", 
						"LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")


#apply an ifelse statement for each column in 'col.comm' to each 
#mapping file in 'all.maps'
check.class<-lapply(all.maps, function(i) 
	unlist(lapply(col.comm, function(x) ifelse(x %in% colnames(i), class(i[ ,x]), "No Field"))))
#the first lapply is to iterate through each mapping file in the list of data frames
#the second lapply iterates through the list of common fields (columns)
#so each mapping file is checked to see if it contains each common field 
#by ifelse(x %in% colnames(i))
#if that particular field is in that mapping file, provide the class (data type)
#if it is not 'No Field' should be stated
#as lapply produces a list nesting it produces a list of lists 
#so unlist the nested list before saving
#this allows for an easier conversion to a flat table (row x col) instead of a list of lists

#now convert the list to a table
check.class<-data.frame(check.class)
#change the row names
rownames(check.class)<-paste(col.comm, ".class", sep="")
#transpose the table so it is study x field (instead of field x study)
check.class<-t(check.class)
#what types of classes are present?
unique(unlist(apply(check.class, 2, unique)))

#use a very similar function to produce the first value in the column
#or an example of what the field loooks like
#or state no field
check.val<-lapply(all.maps, function(i) 
	unlist(lapply(col.comm, function(x) ifelse(x %in% colnames(i), head(i[ ,x], 1), "No Field"))))
check.val<-data.frame(check.val)
rownames(check.val)<-paste(col.comm, ".ex", sep="")
check.val<-t(check.val)

#now make a new table of just study names
study_class_comm_col_table<-data.frame(TITLE=unlist(lapply(all.maps, function(x) unique(x[,"TITLE"]))))
#bind the class and examples tables to this table
study_class_comm_col_table<-cbind(study_class_comm_col_table, check.class, check.val)
dim(study_class_comm_col_table)
colnames(study_class_comm_col_table)
#reorder the columns to place classes and examples of each column side by side
study_class_comm_col_table<-study_class_comm_col_table[, c("TITLE", sort(colnames(study_class_comm_col_table)[2:93]))]

#export as .csv to easily view and format
write.csv(study_class_comm_col_table, 
					file.path(paste(getwd(), 
													"outputs/study_class_comm_col_table.csv", sep="/")), 
					row.names=FALSE)

#now that I have picked through that table and identified particularly 
#problematic fields can collect just those into a table
#can create list of studies and key data to export as report?

#create two lists, one for numeric and one for string data
meta.num<-c("COLLECTION_DATE", "SAMP_SIZE", "DEPTH", "RUN_DATE")
meta.ex<-c("SEQUENCING_METH", "LIBRARY_CONSTRUCTION_PROTOCOL", "PCR_PRIMERS", 
					 "PLATFORM", "RUN_CENTER", "SAMPLE_CENTER", "SAMPLE_LOCATION", "TARGET_GENE")

#create new table of study title, contacts and field ranges, 
#classes and values of interest
ifelse(isTRUE(TRUE %in% (colnames(all.maps[["Polluted Polar Coastal Sediments"]]) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT"))),
			 paste(unique(all.maps[["Polluted Polar Coastal Sediments"]][,colnames(all.maps[["Polluted Polar Coastal Sediments"]]) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "),
			 "No Contacts")

EMP.meta.ls<-lapply(all.maps, function(x) data.frame(
	Field=c("TITLE", 
					"CONTACTS", 
					"COLLECTION_DATE_range",
					"SAMP_SIZE_range",
					"DEPTH_range",
					"RUN_DATE_range",
					"COLLECTION_DATE_class", 
					"DEPTH.class", 
					"RUN_DATE.class",
					"SEQUENCING_METH", 
					"LIBRARY_CONSTRUCTION", 
					"PCR_PRIMERS", 
					"PLATFORM", 
					"RUN_CENTER",  
					"SAMPLE_CENTER", 
					"SAMPLE_LOCATION", 
					"TARGET_GENE",
					""),
	Value=c(
	#study title first
	paste(unique(x[,"TITLE"]), collapse=" "),
	#any values that may indicate contact information
	ifelse(isTRUE(TRUE %in% (colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT"))),
					paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "),
				 "No Contacts"),
	#paste the range and #NA out of total for all numeric columns 
	#(listed in meta.num)
	unlist(lapply(meta.num, function(i) 
	unlist(ifelse(isTRUE(i %in% colnames(x)), paste(paste(range(x[,i], finite=TRUE), collapse=" to "),
																									paste(length(which(is.na(x[,i]))), length(x[,i]), sep=" NA of "), sep="; "), "No Field")
	))),
	#paste class of certain fields
	ifelse(TRUE %in% (colnames(x) %in% "COLLECTION_DATE"), class(x[,"COLLECTION_DATE"]), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "DEPTH"), class(x[,"DEPTH"]), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "RUN_DATE"), class(x[,"RUN_DATE"]), "No Field"),
	#paste first non NA value for listed fields 
	unlist(lapply(meta.ex, function(i)
		unlist(ifelse(isTRUE(i %in% colnames(x)), 
									paste(unique(x[which(!is.na(x[,i])),i]), collapse=" "), "No Field")))),
	#add a break between studies
	paste("")),
	stringsAsFactors=FALSE
))	

#maybe could figure out how to reorder rows of each data frame in the list
#so far promising but trick...
EMP.meta.ls[["Polluted Polar Coastal Sediments"]][["Field"]]["TITLE", "CONTACTS", sort(EMP.meta.ls[["Polluted Polar Coastal Sediments"]][["Field"]][
	which(!EMP.meta.ls[["Polluted Polar Coastal Sediments"]][["Field"]] %in% c("TITLE", "CONTACTS", ""))
	]), ""]

test<-c("TITLE", "CONTACTS", sort(EMP.meta.ls[["Polluted Polar Coastal Sediments"]][["Field"]][
	which(!EMP.meta.ls[["Polluted Polar Coastal Sediments"]][["Field"]] %in% c("TITLE", "CONTACTS", ""))
	]), "")

sort(EMP.meta.ls[["Polluted Polar Coastal Sediments"]][["Field"]][
	which(!EMP.meta.ls[["Polluted Polar Coastal Sediments"]][["Field"]] %in% c("TITLE", "CONTACTS", ""))
	])


#first look at just the new studies
#got this list from my notebook
new<-c("Caporaso_illumina_time_series",
					 "CaporasoIlluminaPNAS2011_5prime", 
					 "Dominguez_Puerto_Rico_Cohort",
					 "Ezenwa_Cape_Buffalo",
					 "Gasser_MWC_catchment_microbes",
					 "Grossart_German_lake_water_sediment",
					 "Haig_WaterPurif_temp_spat",
					 "Hultman_Geochemical_Landscapes_permafrost",
					 "Jansson_Alaskan_fire_chrono_Tanana",
					 "Jurelivicius_Antarctic_cleanup",
					 "Metcalf_SanDiegoZoo_folivorus_primate",
					 "Moore_Yucatan_cenotes",
					 "Rees_VulcanoIsland_seawaterMedSeA",
					 "Spirito_Monensin_Cow_Hindgut",
					 "Thomas_sponge_communities")
new<-which(qiime.file.names %in% new)

EMP.meta.ls.new<-EMP.meta.ls[new]
length(EMP.meta.ls.new)
#convert to flatten to data frame 
EMP.meta.df.new<-do.call("rbind", EMP.meta.ls.new)
#replace row names
rownames(EMP.meta.df.new)<-seq(1:nrow(EMP.meta.df.new))

#export whole data frame to visually inspect 
write.csv(EMP.meta.df.new, 
					file.path(paste(getwd(), 
													"outputs/EMP_metadata_issues_new.csv", sep="/")), 
					row.names=FALSE)

#new and intersting issues to investigate
all.maps[["Geochemical landscapes"]]["COLLECTION_DATE"]
sort(all.maps[["Geochemical landscapes"]][["COLLECTION_DATE"]])
head(all.maps[["Geochemical landscapes"]], 1)

#now work with whole data frame of all studies
#add study title as list title
names(EMP.meta.ls)<-unlist(lapply(all.maps, function(x) unique(x[,"TITLE"])))

#convert to flatten to data frame 
EMP.meta.df<-do.call("rbind", EMP.meta.ls)

#replace row names
rownames(EMP.meta.df)<-seq(1:nrow(EMP.meta.df))

#18 values per study
which(EMP.meta.df$Value=="Polluted Polar Coastal Sediments")
#this is the example of good(-ish) metadata 
EMP.meta.df[595:612, ]

#previously edited manually to find issues, think can subset now 
#should probably export later after adding additional fields (chemistry, etc...)
#but test now

#subset titles, contacts, No Field, NA to NA and integer dates

#this method scrambles the order of everything
EMP.meta.df.sub<-EMP.meta.df[c(which(EMP.meta.df$Field==c("TITLE", "CONTACTS") |
																		 	EMP.meta.df$Value== "No Field"),
															 grep("NA to NA", EMP.meta.df$Value),
															 grep(glob2rx("20* to 20*"), EMP.meta.df$Value)), ]

dim(EMP.meta.df.sub)

#experiment with indexing dataframes in lists
EMP.meta.ls[[1]][which(EMP.meta.ls[[1]][["Field"]] %in% "TITLE"|
											 	EMP.meta.ls[[1]][["Value"]] %in% "No Field"),]
#note format list[[dataframe]][row#, ]
#this makes sense '[[' for single element, '[' for multiple elements
#and ',' to note rows, not columns

#check if additional fields work...
EMP.meta.ls[[1]][c(which(EMP.meta.ls[[1]][["Field"]] %in% c("TITLE", "CONTACTS")|
											 	EMP.meta.ls[[1]][["Value"]] %in% "No Field"),
									 grep("NA to NA", EMP.meta.ls[[1]][["Value"]]),
									 grep(glob2rx("20* to 20*"), EMP.meta.ls[[1]][["Value"]]),
									 nrow(EMP.meta.ls[[1]])
									 ) ,]

#ok, so for the number of studies in the EMP.meta.ls
#subset the TITLE, CONTACTS, any fields are 'No Field', 
#any values are NA to NA
#any dates that are just year
#and the final row (the break between studies)
EMP.meta.ls.sub<-lapply(1:length(EMP.meta.ls), function(x) EMP.meta.ls[[x]][c(which(EMP.meta.ls[[x]][["Field"]] %in% c("TITLE", "CONTACTS") |
																														EMP.meta.ls[[x]][["Value"]] %in% "No Field"),
																										grep("NA to NA", EMP.meta.ls[[x]][["Value"]]),
																										grep(glob2rx("20* to 20*"), EMP.meta.ls[[x]][["Value"]]),
																										nrow(EMP.meta.ls[[x]])													
									) ,] )

#convert to flatten to data frame 
EMP.meta.df.sub<-do.call("rbind", EMP.meta.ls.sub)

#replace row names
rownames(EMP.meta.df.sub)<-seq(1:nrow(EMP.meta.df.sub))

#write the table
write.csv(EMP.meta.df.sub, 
					file.path(paste(getwd(), 
													"outputs/EMP_metadata_issues.csv", sep="/")), 
					row.names=FALSE)
#remember to copy and paste the 
#Geochemical landscapes COLLECTION_DATE_range (starts spring 1951...)

#this seemed useful, but in the end not, may want later though
#http://stackoverflow.com/questions/13006909/export-a-list-of-matrices-nicely-to-the-same-worksheet-in-excel

#previously exported the 4 studies missing most/all sequencing data as csv
#in fact do not have any fields that appear to contain this data
#those studies are "tibetan_plateau_salt_lake_sediment", 
# "Kilauea geothermal soils and biofilms", "Fermilab_spatial_study"
# "Great Lake Microbiome"

#############################################try to evaluate all columns in more than half of the studies
#next
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
col.n.comm[grep(paste(c("NITRO", "_N_"), collapse="|"), col.n.comm, ignore.case=TRUE)]

#look at metadata with nitrogen, make special list
col.nitro<-c("TOT_NITRO", "TOT_NITRO_PERCENT", "TOT_NITRO_UNIT", "TOT_NITRO_UNITS", 
						 "TOT_N_METH")

col.nitro<-c(col.n.comm[grep(paste(c("NITRO", "_N_"), collapse="|"), col.n.comm, ignore.case=TRUE)])

nitro.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.nitro), head(x[,"TITLE"], 1), "No Field"))
nitro.maps<-names(nitro.maps[which(!nitro.maps %in% "No Field")])
lapply(all.maps[nitro.maps], function(x) head(x[colnames(x) %in% col.nitro], 2))
#not many columns shared across studies, see if can focus
length(nitro.maps) #18 studies have some nitrogen
as.matrix(table(unlist(lapply(all.maps[nitro.maps], function(x) colnames(x[colnames(x) %in% col.nitro])))))
#11 have TOT_NITRO, only 7 have TOT_N_METH and only 1 has TOT_NITRO_UNITS/TOT_NITRO_UNIT

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

EMP_temp_table$SAMP_STORE_TEMP<-unlist(lapply(all.maps[temp.maps], function(x) 
	ifelse(TRUE %in% (colnames(x) %in% "SAMP_STORE_TEMP"),  
				 paste(paste(range(x[,"SAMP_STORE_TEMP"], finite=TRUE), collapse=" to "), 
				 			paste(length(which(is.na(x[,"SAMP_STORE_TEMP"]))), length(x[,"SAMP_STORE_TEMP"]), sep=" NA of "), sep="; "), "No Field")), use.names=FALSE)

EMP_temp_table<-as.data.frame(EMP_temp_table)
colnames(EMP_temp_table)<-c("TITLE", "CONTACTS", "TEMP", "ANNUAL_SEASON_TEMP", "AIR_TEMP", "SOIL_TEMP", "SAMP_STORE_TEMP")

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
	
# write.csv(EMP_season_table, 
# 					file.path(paste(getwd(), "outputs/EMP_season_table.csv", sep="/")), 
# 					row.names=FALSE)

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

range(all.maps[["Great Lake Microbiome"]][["NITROGEN_SATURATION"]], finite=TRUE)
summary(all.maps[["Great Lake Microbiome"]][["NITROGEN_SATURATION"]])
length(all.maps[["Great Lake Microbiome"]][["NITROGEN_SATURATION"]])

#other column fragments
col.n.comm.sp
col.n.comm[grep('NAME', col.n.comm, fixed=TRUE)]
col.n.comm[grep('BODY', col.n.comm, fixed=TRUE)]
col.n.comm[grep('HOST', col.n.comm, fixed=TRUE)]
col.n.comm[grep('ORG', col.n.comm, fixed=TRUE)]
col.n.comm[grep('DAY', col.n.comm, fixed=TRUE)]
col.n.comm[grep('DENSITY', col.n.comm, fixed=TRUE)]
col.n.comm[grep('DISS', col.n.comm, fixed=TRUE)]
col.n.comm[grep('GROWTH', col.n.comm, fixed=TRUE)]
col.n.comm[grep('ID', col.n.comm, fixed=TRUE)]
col.n.comm[grep('LOG', col.n.comm, fixed=TRUE)]
col.n.comm[grep('NUMBER', col.n.comm, fixed=TRUE)]
col.n.comm[grep('PRECPT', col.n.comm, fixed=TRUE)]
col.n.comm[grep('SAMP', col.n.comm, fixed=TRUE)]
col.n.comm[grep('SAMPLE', col.n.comm, fixed=TRUE)]
col.n.comm[grep('SITE', col.n.comm, fixed=TRUE)]
#none of these look very necessary to investigate
#what about just the study specific fields
col.n.comm
col.only.1<-sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE)
col.only.1<-col.only.1[which(col.only.1<2)]
col.only.1<-names(col.only.1)
col.only.1.sp<-strsplit(names(col.only.1), "_")
col.only.1.sp<-as.matrix(table(unlist(col.only.1.sp)))
col.only.1.sp<-col.only.1.sp[order(col.only.1.sp[,1], decreasing=TRUE), ]
as.matrix(col.only.1.sp)
col.only.1[grep('TOT', col.only.1, fixed=TRUE)]
col.only.1[grep(names(col.only.1.sp[1]), col.only.1, fixed=TRUE)]

col.only.1.tab<-lapply(1:length(col.only.1.sp), function(x) col.only.1[grep(names(col.only.1.sp[x]), col.only.1, fixed=TRUE)])
names(col.only.1.tab)<-names(col.only.1.sp)
col.only.1.tab[1:50]

#######################
#in regards to the "Have you seen these OTUs" e-mail
#realized I can check the OTU tables from EMP that I downloaded
#when I got all 62 studies

library(phyloseq)
#emp.otu<-import_biom("C:/Users/asus4/Documents/EarthMicrobiomeProject/QIIME_metadata_download/8_17_2013_all_qiime_proccessed_studies_n62/meta_analysis_tmpiuSyWGySPR1izselizZR_otu_table.biom", 
#										 parseFunction=parse_taxonomy_greengenes)
#don't do this, crashed all running programs

#bring in the mapping file for all 62 studies
emp.map62<-read.delim("C:/Users/asus4/Documents/EarthMicrobiomeProject/QIIME_metadata_download/8_17_2013_all_qiime_proccessed_studies_n62/meta_analysis_tmpiuSyWGySPR1izselizZR_map.txt", quote="", stringsAsFactors=FALSE)

#get all titles
unique(emp.map62[, "TITLE"])

#fix known issues
emp.map62[grep("Intertidal microbes 16s for 2009 and", emp.map62$TITLE), "TITLE"]<-"Intertidal microbes 16s for 2009 and 2010"
emp.map62[grep("EPOCA_Svalbard", emp.map62$TITLE), "TITLE"]<-"EPOCA_Svalbard2018"
#now there are only 61 unique titles...
#which means two QIIME processed studies have the same title...

#are study IDs unique?
colnames(emp.map62)[grep("STUDY", colnames(emp.map62))]
length(unique(emp.map62[,"STUDY_ID"]))
#and there are 63 unique study IDs...

unique(emp.map62[,c("TITLE", "STUDY_ID")])
dim(unique(emp.map62[,c("TITLE", "STUDY_ID")]))
#and there are 64 unique combinations of TITLe and STUDY_ID
#ah, two studies have NA for TITLE, hence repeat...

head(emp.map62[91, ], 2)[!is.na(head(emp.map62[91, ], 2))]
colnames(head(emp.map62[91, ], 2))[!is.na(head(emp.map62[91, ], 2))]
#this could work as TITLE EMP_Thomas_sponges_16S_L007
which(emp.map62$RUN_PREFIX=="EMP_Thomas_sponges_16S_L007")
test<-emp.map62[which(emp.map62$RUN_PREFIX=="EMP_Thomas_sponges_16S_L007"), ]
head(test[,which(sapply(test, function(x) !any(is.na(x))))],2)

emp.map62[which(emp.map62$RUN_PREFIX=="EMP_Thomas_sponges_16S_L007"), "TITLE"]<-"EMP_Thomas_sponges_16S_L007"


head(emp.map62[92, ], 2)[!is.na(head(emp.map62[92, ], 2))]
colnames(head(emp.map62[92, ], 2))[!is.na(head(emp.map62[92, ], 2))]
head(emp.map62[92, ], 2)

unique(emp.map62[, "LONGITUDE"])
#oh, this is wrong... hope it is only issue with combining all studies
#in one file... then is a QIIME issue
#can go through list and find ones I don't have yet...

#######################
#new idea for "Have you seen these OTUs" e-mail...

#list biom files
biom.file.names<-list.files(path.expand("~/EarthMicrobiomeProject/QIIME_metadata_download"), pattern=".biom", recursive=TRUE)


library(phyloseq)

#test importing .biom table and searching for 'have you seen these OTUs'
test<-import_biom(paste(path.expand("~/EarthMicrobiomeProject/QIIME_metadata_download"), biom.file.names[1], sep="/"), parseFunction=parse_taxonomy_greengenes)
#importing biom files takes a really long time...
test
otu_table(test)
taxa_names(test)
isTRUE(c(124309, 141145) %in% taxa_names(test))
rm(test)

isTRUE(c(124309, 141145) %in% taxa_names(import_biom(paste(path.expand("~/EarthMicrobiomeProject/QIIME_metadata_download"), biom.file.names[2], sep="/"), parseFunction=parse_taxonomy_greengenes)))

#create list of data frames for each mapping file
have.seen.OTUs<-list()

#import all mapping files, 
# have.seen.OTUs<-lapply(biom.file.names, function(x) 
# 	isTRUE(c(124309, 141145) %in% taxa_names(
# 		import_biom(file.path(paste(path.expand("~/EarthMicrobiomeProject/QIIME_metadata_download"), x, sep="/")), 
# 								parseFunction=parse_taxonomy_greengenes))))

#this also crashed my computer