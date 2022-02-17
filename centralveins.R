# Central Vein Detection algorithm
# Code written by Jordan Dworkin (jdwor@upenn.edu) and edited by Russell Shinohara, Timothy Robert-Fitzgerald, and Melissa Martin

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
# library(pbmcapply) # for working in parallel
library(WhiteStripe)

source("helperfunctions.R") # load necessary helper functions

centralveins_preproc<-function(epi,t1,flair,skullstripped=F,biascorrected=F){

  # print(".Random.seed = ")
  # print(.Random.seed)
  
  #######################################
  ####### Perform bias correction #######
  #######################################
  if(biascorrected==F){
    epi_n4=bias_correct(epi,correction="N4",reorient=F)
    t1_n4=bias_correct(t1,correction="N4",reorient=F)
    flair_n4=bias_correct(flair,correction="N4",reorient=F)
    writenii(epi_n4, "epi_n4")
    writenii(t1_n4, "t1_n4")
    writenii(flair_n4, "flair_n4")
  } else {
    epi_n4<-epi; t1_n4<-t1; flair_n4<-flair
  }
  
  #######################################
  ####### Perform skull stripping #######
  #######################################
  if(skullstripped==F){
    t1_n4_brain=fslbet_robust(t1,correct=F)
    epi_n4_brain=fslbet_robust(epi,correct=F)
    writenii(epi_n4_brain, "epi_brain")
    writenii(t1_n4_brain, "t1_brain")
  } else { # assume the T1 is skull-stripped
    t1_n4_brain<-t1_n4; epi_n4_brain<-epi_n4
    #epi_n4_brain=fslbet_robust(epi_n4,correct=F)
    #writenii(epi_n4_brain, "epi_brain")
  }
  mask<-1*(t1_n4_brain!=0)

  ############################################
  ###### Register T1 to FLAIR space ##########
  ############################################
  t1_to_flair = registration(filename = t1_n4,
                template.file = flair_n4,
                typeofTransform = "Rigid", remove.warp = FALSE,
                outprefix="t1_reg_to_flair") ### rigid

  t1_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(flair_n4), moving = oro2ants(t1_n4_brain),
            transformlist = t1_to_flair$fwdtransforms, interpolator = "welchWindowedSinc"))
  brainmask_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(flair_n4), moving = oro2ants(mask),
            transformlist = t1_to_flair$fwdtransforms, interpolator = "nearestNeighbor"))
  writenii(t1_reg, "t1_n4_reg_flair")
  writenii(brainmask_reg, "brainmask_reg_flair")
  flair_n4_brain<-flair_n4; flair_n4_brain[brainmask_reg==0]<-0
  writenii(flair_n4_brain, "flair_brain")
  
  ######################################################################
  ####### Find and expand CSF for periventricular lesion removal #######
  ######################################################################
  csf=fast(t1_reg,opts='--nobias')
  csf[csf!=1]<-0
  csf=ants2oro(labelClusters(oro2ants(csf),minClusterSize=300))
  csf[csf>0]<-1
  csf=(csf!=1)
  csf=fslerode(csf, kopts = paste("-kernel boxv",2), verbose = TRUE)
  csf=(csf==0)
  writenii(csf, "csf")

  ###################################################
  ###### Register T1 in FLAIR space to EPI ##########
  ###################################################
  flair_to_epi = registration(filename = flair_n4,
                template.file = epi_n4,
                typeofTransform = "Rigid", remove.warp = FALSE,
                outprefix="flair_to_epi") ### rigid
  save(flair_to_epi,file="flair_to_epi")
  brainmask_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_n4), moving = oro2ants(mask),
            transformlist = c(flair_to_epi$fwdtransforms,t1_to_flair$fwdtransforms), interpolator = "nearestNeighbor"))
  writenii(brainmask_reg_epi, 'brainmask_reg_epi')

  ###############################
  ####### Obtain vein map #######
  ###############################
  # frangi=frangifilter(image=epi_n4,mask=brainmask_reg_epi)
  # frangi[frangi<0]<-0
  # writenii(frangi, 'frangi')

  frangi=frangifilternoc3d(image=epi_n4,mask=brainmask_reg_epi)
  frangi[frangi<0]<-0
  writenii(frangi, 'frangi')
}


centralveins_seg<-function(epi_n4_brain,t1_reg,flair_n4_brain,brainmask_reg,csf,flair_to_epi_filename,mimosa_model = mimosa::mimosa_model_No_PD_T2,probmap=NULL,parallel=F,cores=2,probthresh=0.2){

  ###########################################
  ####### Perform lesion segmentation #######
  ###########################################

### JORDAN'S CODE

  # if(is.null(probmap)){

  #   mimosa_data = mimosa_data(
  #     brain_mask = brainmask_reg,
  #     FLAIR = flair_n4_brain,
  #     T1 = t1_reg,
  #     normalize = 'Z',
  #     cores = cores,
  #     verbose = TRUE)
    
  #   mimosa_df = mimosa_data$mimosa_dataframe
  #   mimosa_cm = mimosa_data$top_voxels
  #   rm(mimosa_data)
    
  #   predictions = predict(mimosa::mimosa_model_No_PD_T2,
  #                         newdata = mimosa_df,type = 'response')
  #   probmap = niftiarr(brainmask_reg, 0)
  #   probmap[mimosa_cm == 1] = predictions
  #   probmap = fslsmooth(probmap,sigma = 1.25,mask = brainmask_reg,
  #                       retimg = TRUE,smooth_mask = TRUE)
  #   writenii(probmap, 'mimosa_prob')
  #   writenii(probmap>probthresh, 'mimosa_mask_30')


### MELISSA'S CODE

  if (is.null(probmap)) {
 
# WhiteStripe normalize data
    ind = whitestripe(t1_reg, "T1")
    t1_n4_reg_brain_ws = whitestripe_norm(t1_reg, ind$whitestripe.ind)
    ind = whitestripe(flair_n4_brain, "T2")
    flair_n4_brain_ws = whitestripe_norm(flair_n4_brain, ind$whitestripe.ind)
    
# Prepare data for MIMoSA
    mimosa = mimosa_data(brain_mask=brainmask_reg, FLAIR=flair_n4_brain_ws, T1=t1_n4_reg_brain_ws, gold_standard=NULL, normalize="no", cores = cores, verbose = TRUE)
    mimosa_df = mimosa$mimosa_dataframe
    cand_voxels = mimosa$top_voxels
    tissue_mask = mimosa$tissue_mask

  # Apply model to test image
    # load("/mimosa_model.RData")
    predictions_WS = predict(mimosa_model, mimosa_df, type="response")
    predictions_nifti_WS = niftiarr(cand_voxels, 0)
    predictions_nifti_WS[cand_voxels==1] = predictions_WS
    probmap = fslsmooth(predictions_nifti_WS, sigma = 1.25, mask=tissue_mask, retimg=TRUE, smooth_mask=TRUE) # probability map
    
    writenii(probmap, 'mimosa_prob')
    writenii(probmap>probthresh, 'mimosa_mask')
  }

  ###############################################################
  ####### Split confluent lesions for individual analysis #######
  ###############################################################
  les=lesioncenters(probmap,probmap>probthresh,parallel=parallel,cores=cores)$lesioncenters
    
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
  
  #################################################
  ###### Map masks in FLAIR space to EPI ##########
  #################################################

  load(flair_to_epi_filename)

  t1_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_n4_brain), moving = oro2ants(t1_reg),
            transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
  flair_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_n4_brain), moving = oro2ants(flair_n4_brain),
            transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
  mimosa_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_n4_brain), moving = oro2ants(probmap),
            transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
  mimosa_mask_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_n4_brain), moving = oro2ants(probmap>probthresh),
            transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
  les_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_n4_brain), moving = oro2ants(les),
            transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))

  writenii(t1_reg_epi, 't1_reg_epi')
  writenii(flair_reg, 'flair_reg')
  writenii(mimosa_reg, 'mimosa_reg')
  writenii(mimosa_mask_reg, 'mimosa_mask_reg')
  writenii(les_reg, 'les_reg')
  
  ###################################################
  ####### Obtain distance-to-the-boundary map #######
  ###################################################
  dtb=dtboundary(les_reg)
  writenii(dtb,"dtb")

}

# CVS detection function:
centralveins<-function(les_reg,frangi,dtb,parallel=F,cores=2){

  # print(".Random.seed = ")
  # print(.Random.seed)
  
  ######################################################################
  ####### Perform permutation procedure to get CVS probabilities #######
  ######################################################################
  lables=ants2oro(labelClusters(oro2ants(les_reg),minClusterSize=27))
  writenii(lables,"lables.nii.gz")
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
    if(parallel==T) {
      nullscores=as.vector(unlist(mclapply(1:1000,getnulldist,
                                          centsub,coords,frangsub,
                                          mc.cores=cores)))
    } else {
      nullscores=as.vector(unlist(lapply(1:1000,getnulldist,centsub,coords,frangsub)))
    }
    
    # get CVS probability for lesion j
    lesprob=sum(nullscores<score)/length(nullscores)
    avprob=c(avprob,lesprob)
    probles[lables==j]<-lesprob

    print(paste0("Done with lesion ",j," of ",maxles))
  }
  
  writenii(probles,"cvs_probmap.nii.gz")
  saveRDS(avprob,"cvs_avprob")
  return(list(candidate.lesions=lables>0,cvs.probmap=probles,cvs.biomarker=mean(avprob),numles=maxles))
}
