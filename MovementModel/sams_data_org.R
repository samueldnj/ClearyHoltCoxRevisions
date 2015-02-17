# This code takes simulated data from a population with known release and recovery locations over a series of years and organizes
# it into a sum of the number of recovered fish released in each area given the number of years it was at liberty. The
# end result is a list where each component is the recovery matrix for that many years at liberty (e.g: recovery_subsets[[2]] will
#give a recovery matrux for all fish at libety two years where the rows are the release area, the columns are the recovery areas


source("Sam_Michelle_simulatorV2.r") 	#source in the function to simulate the data to be re-organized
sablefish_example<-Srecoveries(Nyears=15,Nsims=1, Release.matrix=Release.matrix) # Run the function to make data
as.data.frame(sablefish_example$recoveries)->recovery_sim 	# The function produces all sorts of info in a list, so just take the recoveries which are needed for this model

##### Start of Data Organization############
recovery_subsets<-list() ## Make a storage space for sorting the data by the liberty years.
Nareas<-length(AreaNames) # make an index for the number of areas
libyears<-5 #decide the number of liberty years you want to make subsets for
for (i in 1:libyears)  # for each liberty year the loop will run
	{
	releases_and_recoveries<-matrix(nrow=Nareas, ncol=Nareas) #define storage space for the liberty year recovery matrix
	rownames(releases_and_recoveries)<-AreaNames  # add in the rownames 
	colnames(releases_and_recoveries)<-AreaNames #add in the column names
	temp<-subset(recovery_sim,liberty.year==i) # subset the data by the ith liberty year 
	      for (j in 1:length(AreaNames)) # for each release area
		  {
			for (k in 1:length(AreaNames)) #and for each recovery area
			{
			 releases_and_recoveries[j,k]<-nrow(subset(temp,release_area==AreaNames[j] & recovery_area==AreaNames[k])) #determines the number of individuals released in area j and recovered in area k who were at liberty i years
	       }
	     } 
		recovery_subsets[[i]]<-releases_and_recoveries #save the matrix as the ith component of all recovery subsets
	}

recovery_subsets # print out the recovery subsets

# list recovery_subsets