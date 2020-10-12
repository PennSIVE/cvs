source("centralveins_full.R")
library(neurobase)
model_path = Sys.getenv("MIMOSA_MODEL")
if (model_path == "") {
    mimosa_model = mimosa::mimosa_model_No_PD_T2
} else {
    load(model_path)
}
setwd("/out")


centralveins_preproc(readnii("/epi.nii.gz"), readnii("/t1.nii.gz"), readnii("/flair.nii.gz"), skullstripped = F, biascorrected = F)

epi_n4_brain = readnii('epi_brain.nii.gz')
t1_reg = readnii('t1_n4_reg_flair.nii.gz')
flair_n4_brain = readnii('flair_brain.nii.gz')
brainmask_reg = readnii('brainmask_reg_flair.nii.gz')
csf = readnii('csf.nii.gz')
flair_to_epi_filename = 'flair_to_epi'
centralveins_seg(epi_n4_brain, t1_reg, flair_n4_brain, brainmask_reg, csf, flair_to_epi_filename, mimosa_model = mimosa_model, probmap = NULL, parallel = T, cores = as.numeric(Sys.getenv("NSLOTS")), probthresh = 0.2)

result = centralveins(readnii("les_reg.nii.gz"), readnii("frangi.nii.gz"), readnii("dtb.nii.gz"), parallel = T, cores = as.numeric(Sys.getenv("NSLOTS")))

saveRDS(result, file = "centralveins")
