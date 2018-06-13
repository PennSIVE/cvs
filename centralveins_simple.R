# Central Vein Detection algorithm - Simple
# Code written by Jordan Dworkin (jdwor@upenn.edu)

#Variable descriptions:
# **probmap** is an image of class nifti, containing the probability that each voxel
# is a lesion voxel.
#---
# **binmap** is a nifti mask in which voxels are classified as either lesion voxels
# or not lesion voxels. Note that mask should be in the same space as the probmap volume.
#---
# **veinmap** is an image of class nifti, containing the vesselness score for each voxel.
#---
# **tissue** is an image of class nifti, where candidate voxels (e.g., deep white matter) are labeled 1,
# and non-candidate voxels are labeled 0.
#---
# **parallel** is a logical value that indicates whether the user's computer
# is Linux or Unix (i.e. macOS), and should run the code in parallel.
#---
# **cores** (if parallel = TRUE) is an integer value that indicates how many cores
# the function should be run on.
#---

# Function output:
# A list containing candidate.lesions (a nifti file with labeled lesions evaluated for CVS),
# cvs.probmap (a nifti file in which candidate lesions are labeled with their CVS probability), and
# cvs.biomarker (a numeric value representing the average CVS probability of a subject's lesions).

# Required packages:
library(neurobase) # for reading and writing nifti files
library(ANTsRCore) # for cluster labeling
library(extrantsr) # for transitioning between ants and nifti
library(parallel) # for working in parallel
library(pbmcapply) # for working in parallel


source("helperfunctions.R") # load necessary helper functions

# Simple CVS detection function:
centralveins_simple=function(probmap,binmap,veinmap,tissue,parallel=F,cores=2){

  ###############################################################
  ####### Split confluent lesions for individual analysis #######
  ###############################################################
  les=lesioncenters(probmap,binmap,parallel=parallel,cores=cores)$lesioncenters
  
  ############################################
  ####### Remove non-candidate lesions #######
  ############################################
  nctis=(tissue==0)
  lables=ants2oro(labelClusters(oro2ants(les>0),minClusterSize=27))
  for(j in 1:max(lables)){
    if(sum(nctis[lables==j])>0){
      lables[lables==j]<-0
    }
  }
  les=lables>0
  
  ###################################################
  ####### Obtain distance-to-the-boundary map #######
  ###################################################
  dtb=dtboundary(les)
  
  ######################################################################
  ####### Perform permutation procedure to get CVS probabilities #######
  ######################################################################
  lables=ants2oro(labelClusters(oro2ants(les),minClusterSize=27))
  probles=lables
  avprob=NULL
  maxles=max(as.vector(lables))
  for(j in 1:maxles){
    # get true coherence for lesion j
    frangsub=veinmap[lables==j]
    centsub=dtb[lables==j]
    coords=which(lables==j,arr.ind=T)
    prod=frangsub*centsub
    score=sum(prod)
    # get 1000 null coherence values for lesion j
    if(parallel==T){
      nullscores=as.vector(unlist(mclapply(1:1000,getnulldist,
                                           centsub,coords,frangsub,
                                           mc.cores=cores)))
    }else{
      nullscores=as.vector(unlist(lapply(1:1000,getnulldist,
                                         centsub,coords,frangsub)))
    }
    
    # get CVS probability for lesion j
    lesprob=sum(nullscores<score)/length(nullscores)
    avprob=c(avprob,lesprob)
    probles[lables==j]<-lesprob
    
    print(paste0("Done with lesion ",j," of ",maxles))
  }
  
  return(list(candidate.lesions=lables>0,cvs.probmap=probles,cvs.biomarker=mean(avprob)))
}
