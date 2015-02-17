#Sablefish Recovery Matrix Simulator
# This simulator uses a pre-determined releases matrix to generate a simulated recovery matrix. .
# The default code assumes: mortality is constant over time, migration rates are constant over time, and the proportion of migrants and residents is constant over time
#############   SIMULATION SET-UP###################################
############# Input data  ##########################################
A<-c(rep(100, 10))           #number of tags released in A each year (10 years)
B<-c(rep(100, 10))           # of tags released in B each year
C<-c(rep(100, 10))           #of tags released in c each year
D<-c(rep(100,10))
E<-c(rep(0,10))
F<-c(rep(0,10))
G<-c(rep(0,10))
H<-c(rep(0,10))
I<-c(rep(0,10))

releases<-(cbind(A,B,C,D,E,F,G,H,I))

AreaNames<-c("A","B","C","D","E","F","G","H","I")                                        #vector of #tags released each year in each area
Release.matrix<-matrix(data=releases, nrow=10, byrow=FALSE, dimnames=list(c(seq(1:10)),AreaNames)) # places release data into a matrix
#####################################################################################################################################################
#Nsims=25
Nyears=15
set.seed(123086)
Srecoveries<-function (Nyears=15,Nsims=25, Release.matrix=Release.matrix)                                                                          #numbr of siulations to run
 {  ### Set up Parameters
      AreaNames<-colnames(Release.matrix)      #define strata
      NAreas<-length(AreaNames)
      Ntags<- sum(colSums(Release.matrix))                                             # Total number of tags that will be introduced throughout entire time series
      nareas=9

              Migrant.Rates<-matrix(data=rep(1/nareas,nareas*nareas),nareas,nareas)
              for (i in (1:nareas))
              {
                Migrant.Rates[i,i]=0 
               }

      #Migrant.Rates<-matrix(data=c(.1,.45,.45,.45,.1,.45,.45,.45,.1), ncol=9,nrow=9,byrow=TRUE)  #Prob a migrant moves from i to j (PAA,PAB,PAC,PBA,PBB,PBC,PCA,PCB,PCC), assumes constant over time
      MigrantProbs<-matrix(data=c(rep(1,3*10)),ncol=nareas,nrow=10)                 #probability something is a migrant in each area in each year    
      p.Rec<-.20                                                                  # Probability of observing tags, currently written to be constant between places and time
      mortality.rates<-rep(.1,nareas)
      release.year<-vector(length=Ntags)
      Recovery.year<-vector(length=Ntags) 
                                                                              #vector of release year for each tag
                                                                              #vector of recovery year for each tag
      #######Storage needed for function to run###############
      all.sim.tracking<-list()                                                        #used to store positions and fates of all tags at all time steps
      all.sim.recoveries<-list()                                                      #used to store recovery matricies
      resident.types<-vector(length=Ntags)                                            #identifies the M/R type for each tagged individual
      sightings<-paste("Y",AreaNames,sep="")                                          #identifies observed individs
      total.tags.per.year<-colSums(t(Release.matrix))                                #index to enter in new tags each year of the simulation
      cumsum(total.tags.per.year)-total.tags.per.year[1]+1->release.matrix.index  #provides an index to later assign which year the initial release of individuals came from
      Recovery.Matrix.New <-matrix(nrow = Ntags, ncol =Nyears+4, byrow=FALSE,
      dimnames=list(c(seq(1:Ntags)),c("tag.num",seq(1:Nyears),"release_year","release_area","recovery_area")))                     #Recovery.Matrix= matrix to store recoveries from each simulation
       Recovery.Matrix.New[,1]=seq(1:Ntags)                                                     #assigns tag number and places it in the first column of the recovery matrix
         for (k in 1:nrow(Release.matrix))                                                      #loop that will circle through all years of recovery data
           {
             temp<-vector(length=total.tags.per.year[k])                           #vector to store area names until they are assigned to the recovery matrix
             counter<-rep(1,length=length(AreaNames)+1)                            #index to place area names within the temp vector during each cycle of the loop
             for (j in c(1:length(AreaNames)))                                 # Loop that generates the area names to put in the recovery matrix
               {
                 c(rep(AreaNames[j],Release.matrix[k,j])) ->for.temp            #repetes the area name for each release in this area during year k
                 # browser()
                 if(length(for.temp)==0) counter[j+1] =counter[j+1] 
                 if (length(for.temp)!=0)
                 {
                 length(for.temp)+counter[j]->counter[j+1] 
                                      #updates the counter
                 for.temp->temp[counter[j]:c(counter[j]+length(for.temp)-1)]    #puts the area names in the correct place with temp, which can then be assigned to the recovery matrix
                  #browser()
                 }
                 }

             Recovery.Matrix.New[release.matrix.index[k]:(release.matrix.index[k]+        # puts the area names into the correct portion of the recovery matrix
             length(temp)-1),k+1]=temp
              
             Recovery.Matrix.New[release.matrix.index[k]:(release.matrix.index[k]+        # puts the area names into the correct portion of the recovery matrix
             length(temp)-1),Nyears+2]=as.numeric(k)
               Recovery.Matrix.New[release.matrix.index[k]:(release.matrix.index[k]+        # puts the area names into the correct portion of the recovery matrix
             length(temp)-1),Nyears+3]=temp
              }
for (p in 1:Nsims)                                          
 {  Recovery.Matrix.New->Recovery.Matrix
     print(p)
      #remove(Recovery.Matrix) 
 
################################## Start Simulation ####################################3
   ### Assign each indidividuals migratory status
#for (p in 1:Nsims)
  #{
     for (l in 1:length(AreaNames))                           # For each area
         {                                                    

            for (i in 1:nrow(Recovery.Matrix))               # go through each tag
     {
        for (j in 1:(nrow(Release.matrix)) )                # and for each year with releases
       {
          ifelse (Recovery.Matrix[i,j+1]==AreaNames[l], sample(x=c("M","R"), size=1, prob=c(MigrantProbs[j,l],(1-MigrantProbs[j,l])))->resident.types[i],resident.types[i]->resident.types[i] )
       }} }                                                                   #assign the tagged individual as a migrant or a resident
Tracking.Matrix<-cbind(Recovery.Matrix,resident.types)                          
Tracking.Matrix


 for (k in 2:Nyears)
    {
      if (k<=nrow(Release.matrix)) Tracking.Matrix[1:(release.matrix.index[k]-1),]->old.releases         #Make sure the previous years releases are added into the popultaion
      if (k>nrow(Release.matrix)) Tracking.Matrix[1:Ntags,]->old.releases                              #if there are no more new releases, simply start where you left off
           for (i in 1:nrow(old.releases))                                                       ###Survival Loop
            {  if(old.releases[i,k]=="0")   old.releases[i,k+1]<-0                               #if an animal was dead last year, make it dead this year
             for (l in 1:length(AreaNames))
               {
                 if(old.releases[i,k]==AreaNames[l]) old.releases[i,k+1]<-sample(c(0,1),size=1,prob=c(mortality.rates[l],(1-mortality.rates[l]))) ##determine if the animal survived
                 if(old.releases[i,k]==sightings[l])  old.releases[i,k+1]<-0                         ###if an animal was seen last year, mak sure it is dead this year
               }
            }

         for (i in 1:nrow(old.releases))
             {
             for (l in 1:length(AreaNames))
                { if(old.releases[i,ncol(old.releases)]=="R" && old.releases[i,k+1]=="1") old.releases[i,k+1]<-old.releases[i,k]        ##if this is a living resident, put it in the same spot
                  if (old.releases[i,k]==AreaNames[l] && old.releases[i,k+1]=="1") old.releases[i,k+1]<-sample((AreaNames), size=1, prob=Migrant.Rates[l,])
          }
    }
 ###observe.Indviduals
   for (i in 1:nrow(old.releases))
     {
        for (l in 1:length(AreaNames))
          {
        if (old.releases[i,k+1]==AreaNames[l]) old.releases[i,k+1]<-sample(c(sightings[l],AreaNames[l]), size=1, prob=c(p.Rec,1-p.Rec))     #put either a sighting flag or same area flag into the old releases based on p recovery

          }
        }
       old.releases[,k+1] ->Tracking.Matrix[1:nrow(old.releases),k+1]
    }
as.data.frame(Tracking.Matrix)->all.sim.tracking[[p]]
 #Tracking.Matrix   to Recovery Matrix

 for (k in 1:(ncol(Tracking.Matrix)-1))
     {
        for (i in 1:nrow(Tracking.Matrix))
          {
            for( l in 1:length(AreaNames))
             {
               ifelse( Tracking.Matrix[i,k]==sightings[l], Recovery.Matrix[i,k]<-AreaNames[l], Recovery.Matrix[i,k]<-Recovery.Matrix[i,k]) # if it was sighted put in the area name, otherwise leave it alone
               ifelse( Tracking.Matrix[i,k]==sightings[l], Recovery.year[i]<-k-1, Recovery.year[i]<-Recovery.year[i]) #if it was sighted mark down the recapture year, otherwise leave it alone
               ifelse( Tracking.Matrix[i,k]==sightings[l], Recovery.Matrix[i,Nyears+4]<-AreaNames[l], Recovery.Matrix[i,Nyears+4]<-Recovery.Matrix[i,Nyears+4]) #if it was recovered add the recovery area name to the right column for later subsets that pool all years
          }
        }
      }
 Recovery.Matrix[is.na(Recovery.Matrix)] <- 0
 as.data.frame(Recovery.Matrix)->test
test$recovery.year<-Recovery.year
test$liberty.year<-as.numeric(as.character(test$recovery.year))-as.numeric(as.character(test$release_year)) # calculate the liberty years for later subsetting and data concatination. 
test[test$liberty.year<=0,]$liberty.year<-0
test->all.sim.recoveries[[p]]


}
return(list(recoveries=all.sim.recoveries,tracking=all.sim.tracking,Recovery.year=Recovery.year)) }

Srecoveries(Nyears=15,Nsims=1, Release.matrix=Release.matrix)->sablefish_example   #numbr of simulations to run