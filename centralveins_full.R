# Central Vein Detection algorithm
# Code written by Jordan Dworkin (jdwor@upenn.edu)

#Variable descriptions:
#---
# **epi** is a T2*-EPI volume of class nifti.
#---
# **t1** is a T1-weighted volume of class nifti.
#---
# **flair** is a T2-FLAIR volume of class nifti.
#---
# **probmap** is an image of class nifti, containing the probability that each voxel
# is a lesion voxel. If a probability map is not included, the MIMoSA model will be applied (Valcarcel et al., 2017).
#---
# **binmap** (if probmap!=NULL) is a nifti mask in which voxels are classified as either lesion voxels
# or not lesion voxels. Note that mask should be in the same space as the probmap volume.
#---
# **parallel** is a logical value that indicates whether the user's computer
# is Linux or Unix (i.e. macOS), and should run the code in parallel.
#---
# **cores** (if parallel = TRUE) is an integer value that indicates how many cores
# the function should be run on.
#---
# **skullstripped** is a logical value reflecting whether or not the images have already been skull-stripped.
#---
# **biascorrected** is a logical value reflecting whether or not the images have already been bias-corrected.
#---
# **c3d** is a logical value reflecting whether or not the Convert3D imaging toolbox is installed.
#---

# Function output:
# A list containing candidate.lesions (a nifti file with labeled lesions evaluated for CVS),
# cvs.probmap (a nifti file in which candidate lesions are labeled with their CVS probability), and
# cvs.biomarker (a numeric value representing the average CVS probability of a subject's lesions).

# Required packages:
library(neurobase) # for reading and writing nifti files
library(extrantsr) # for bias correction, skull-stripping, and registration
library(mimosa) # for performing lesion segmentation
library(ANTsRCore) # for cluster labeling
library(fslr) # for smoothing and tissue class segmentation
library(parallel) # for working in parallel
library(pbmcapply) # for working in parallel

source("helperfunctions.R") # load necessary helper functions

# CVS detection function:
centralveins=function(epi,t1,flair,probmap=NULL,binmap=NULL,parallel=F,
                      cores=2,skullstripped=F,biascorrected=F,c3d=F){

  #######################################
  ####### Perform bias correction #######
  #######################################
  if(biascorrected==F){
    epi=bias_correct(epi,correction="N4",reorient=F)
    t1=bias_correct(t1,correction="N4",reorient=F)
    flair=bias_correct(flair,correction="N4",reorient=F)
  }
  
  ####################################
  ####### Register flair to T1 #######
  ####################################
  flair=registration(filename=flair,template.file=t1,typeofTransform="Rigid",
                     remove.warp=FALSE,outprefix="fun")
  flair=flair$outfile
  
  #######################################
  ####### Perform skull stripping #######
  #######################################
  if(skullstripped==F){
    t1=fslbet_robust(t1,correct=F)
    epi=fslbet_robust(epi,correct=F)
    flair[t1==0]<-0
  }
  mask=(t1!=0)
  
  ###########################################
  ####### Perform lesion segmentation #######
  ###########################################
  if(is.null(probmap) & is.null(binmap)){
    mimosa_data = mimosa_data(
      brain_mask = mask,
      FLAIR = flair,
      T1 = t1,
      normalize = 'Z',
      cores = cores,
      verbose = TRUE)
    
    mimosa_df = mimosa_data$mimosa_dataframe
    mimosa_cm = mimosa_data$top_voxels
    rm(mimosa_data)
    
    predictions = predict(mimosa::mimosa_model_No_PD_T2,
                          newdata = mimosa_df,type = 'response')
    probmap = niftiarr(mask, 0)
    probmap[mimosa_cm == 1] = predictions
    probmap = fslsmooth(probmap,sigma = 1.25,mask = mask,
                        retimg = TRUE,smooth_mask = TRUE)
  }
  
  ###############################################################
  ####### Split confluent lesions for individual analysis #######
  ###############################################################
  les=lesioncenters(probmap,probmap>0.3,parallel=parallel,cores=cores)$lesioncenters
  
  ###############################
  ####### Obtain vein map #######
  ###############################
  if(c3d==T){
    frangi=frangifilter(image=epi,mask=epi!=0,parallel=parallel,cores=cores,c3d=c3d)
    frangi[frangi<0]<-0
    
    # Register vein map to T1
    regs=labelreg(epi,frangi,t1)
    epi=regs$image_reg
    frangi=regs$label_reg
  }else{
    epi=registration(filename=epi,template.file=t1,typeofTransform="Rigid",
                       remove.warp=FALSE,outprefix="fun")
    epi=epi$outfile
    
    frangi=frangifilter(image=epi,mask=probmap>0.3,parallel=parallel,cores=cores,c3d=c3d)
    frangi[frangi<0]<-0
  }
  
  ######################################################################
  ####### Find and expand CSF for periventricular lesion removal #######
  ######################################################################
  csf=fast(t1,opts='--nobias')
  csf[csf!=1]<-0
  csf=ants2oro(labelClusters(oro2ants(csf),minClusterSize=300))
  csf[csf>0]<-1
  csf=(csf!=1)
  csf=fslerode(csf, kopts = paste("-kernel boxv",3), verbose = TRUE)
  csf=(csf==0)
  
  ##############################################
  ####### Remove periventricular lesions #######
  ##############################################
  lables=ants2oro(labelClusters(oro2ants(les>0),minClusterSize=27))
  for(j in 1:max(lables)){
    if(sum(csf[lables==j])>0){
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
    frangsub=frangi[lables==j]
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
