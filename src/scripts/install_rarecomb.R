#!/R

# # unload RareComb (only needed if RareComb is currently loaded)
# unloadNamespace("RareComb")
# # remove old package
# remove.packages('RareComb')
# install local package
install.packages('/data5/UK_Biobank/analysis_for_oligogenic_paper/software/RareComb', '/data5/deepro/miniconda3/envs/rarecomb/lib/R/library', repos = NULL, type = "source")
# or install from repo
# install.packages('RareComb', repos='http://cran.us.r-project.org')
# load package
library('RareComb')
