#!/usr/bin/env Rscript
set.seed(123, kind = "L'Ecuyer-CMRG")
source("centralveins.R")
library(neurobase)
library(argparser)

# Create a parser
p <- arg_parser("Run CVS")
# Add command line arguments
p <- add_argument(p, "--outdir", help = "Output directory", default = "/tmp")
p <- add_argument(p, "--indir", help = "Input directory", default = "/data")
p <- add_argument(p, "--flair", help = "FLAIR image", default = "flair.nii.gz")
p <- add_argument(p, "--t1", help = "T1 image", default = "t1.nii.gz")
p <- add_argument(p, "--epi", help = "EPI/T2star image", default = "epi.nii.gz")
p <- add_argument(p, "--thresh", help = "Threshold to binarize probability map", default = "0.2")
p <- add_argument(p, "--skullstrip", help = "Whether to skull strip inputs", flag = TRUE)
p <- add_argument(p, "--n4", help="Whether to N4 correct input", flag = TRUE)
# Parse the command line arguments
argv <- parse_args(p)


model_path = Sys.getenv("MIMOSA_MODEL")
if (model_path == "") {
    mimosa_model = mimosa::mimosa_model_No_PD_T2
} else {
    load(model_path)
}

setwd(argv$indir)
epi = readnii(argv$epi)
t1 = readnii(argv$t1)
flair = readnii(argv$flair)
setwd(argv$outdir)
centralveins_preproc(epi, t1, flair, skullstripped = argv$skullstrip, biascorrected = argv$n4)


epi_n4_brain = readnii('epi_brain.nii.gz')
t1_reg = readnii('t1_n4_reg_flair.nii.gz')
flair_n4_brain = readnii('flair_brain.nii.gz')
brainmask_reg = readnii('brainmask_reg_flair.nii.gz')
csf = readnii('csf.nii.gz')
flair_to_epi_filename = 'flair_to_epi'
centralveins_seg(epi_n4_brain, t1_reg, flair_n4_brain, brainmask_reg, csf, flair_to_epi_filename, mimosa_model = mimosa_model, probmap = NULL, parallel = T, cores = as.numeric(Sys.getenv("CORES")), probthresh = argv$thresh)

result = centralveins(readnii("les_reg.nii.gz"), readnii("frangi.nii.gz"), readnii("dtb.nii.gz"), parallel = T, cores = as.numeric(Sys.getenv("CORES")))

saveRDS(result, file = "centralveins")
