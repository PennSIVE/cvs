source("helperfunctions.R")
source("centralveins_presegmented.R")
library(neurobase)
library(ANTsRCore)
library(extrantsr)
library(parallel)
library(pbmcapply)

### Read in lesion probability map
probmap=readnii('mimosa.nii') # see following github repo for lesion segmentation software
# https://github.com/avalcarcel9/mimosa

### Read in vein map (using Frangi filter, or similar vesselness filter)
veinmap=readnii("frangi.nii") # see following link for Frangi filter software
# https://sourceforge.net/p/c3d/git/ci/master/tree/doc/c3d.md#-hessobj-hessian-objectness-hessian-objectness-filter

### Read in map of ventricles to remove periventricular lesions from consideration
ventricles=readnii("ventricles.nii") # can use FSL fast ('fast' function in fslr) for tissue seg
ventricles=1-ventricles
vent_expanded=fslerode(ventricles, kopts = "-kernel boxv 3", verbose = TRUE)
vent_expanded=1-vent_expanded

### Locate distinct lesions within probability map for individual analysis
# set probability threshold for lesions
thresh=0.30

# run lesion center detection code
lesion_centers=lesioncenters(probmap=probmap,binmap=probmap>thresh,
                             c3d=T,parallel=T,cores=4)
lesionmap=lesion_centers$lesioncenters

# remove periventricular lesions
for(i in 1:max(lesionmap)){
  if(sum(vent_expanded[lesionmap==i])>0){
    lesionmap[lesionmap==i]<-0
  }
}
lesionmap[lesionmap!=0]<-1

### Run probabilistic central vein detection
results=centralveins_presegmented(lesionmap,veinmap,parallel=T,cores=4)

# can examine/save results$candidate.lesions to see lesions that were evaluated
# can examine/save results$cvs.probmap to see CVS probabilities for each candidate lesion
# can examine results$cvs.biomarker to see the subject-level average of CVS probabilities





