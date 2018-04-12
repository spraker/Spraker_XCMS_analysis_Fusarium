##Joe Spraker 3/26/17 - contact: jspraker@email.arizona.edu
##XCMS/MetaXCMS Analyses of Positive mode LCMS datasets from Fusarium response to Ralstonia solanacearum conditioned media 
##Prior to analyses all directories need to be set up set up for pairwise comparisons between groups as described in  Patti GJ, Tautenhahn R, Siuzdak G. 2012. Meta-analysis of untargeted metabolomic data from multiple profiling experiments. Nat Protoc 7:508–16.

library(xcms)
library(multtest)

##Set your working directory to the one contaning data files
setwd("~path/to/directory/media_High_Low")
#optional: reduce print output for checking datasets to first 50 (saves space)
options(max.print=50)
#Peak identification and filtration - one row per dataset output (m/z:number of peaks)
## Add argument nSlaves=x to specify number of cores to be used.
xset<-xcmsSet(method="centWave", mzCenterFun="wMeanApex3", noise=1000, ppm=5, peakwidth=c(5,30), prefilter=c(3,50000), nSlaves=4)
#Perform retention time correction across samples
xset1 <- retcor(xset, method="obiwarp", plottype = c("deviation"))
#Match peaks across samples
xset2 <- group(xset1, bw = 10, minfrac = 0.5, mzwid = 0.025)
#Fill in missing peaks and calculate statistics
xset3 <- fillPeaks(xset2) 
#Generate feature tables and EICs
##change class lavels based on pairwise comparison
dr <- diffreport(xset3, filebase ="mediaH_vs_mediaL", eicmax=100)
write.csv(dr, file = "mediaH_vs_mediaL.csv")

setwd("~path/to/directory/High_GMI_media")
#optional: reduce print output for checking datasets to first 50 (saves space)
options(max.print=50)
#Peak identification and filtration - one row per dataset output (m/z:number of peaks)
## Add argument nSlaves=x to specify number of cores to be used.
xset<-xcmsSet(method="centWave", mzCenterFun="wMeanApex3", noise=1000, ppm=5, peakwidth=c(5,30), prefilter=c(3,50000), nSlaves=4)
#Perform retention time correction across samples
xset1 <- retcor(xset, method="obiwarp", plottype = c("deviation"))
#Match peaks across samples
xset2 <- group(xset1, bw = 10, minfrac = 0.5, mzwid = 0.025)
#Fill in missing peaks and calculate statistics
xset3 <- fillPeaks(xset2) 
#Generate feature tables and EICs
##change class lavels based on pairwise comparison
dr <- diffreport(xset3, filebase ="HighN_GMI_vs_media", eicmax=100)
write.csv(dr, file = "HighN_GMI_vs_media.csv")

setwd("~path/to/directory/High_rmy_media")
#optional: reduce print output for checking datasets to first 50 (saves space)
options(max.print=50)
#Peak identification and filtration - one row per dataset output (m/z:number of peaks)
## Add argument nSlaves=x to specify number of cores to be used.
xset<-xcmsSet(method="centWave", mzCenterFun="wMeanApex3", noise=1000, ppm=5, peakwidth=c(5,30), prefilter=c(3,50000), nSlaves=4)
#Perform retention time correction across samples
xset1 <- retcor(xset, method="obiwarp", plottype = c("deviation"))
#Match peaks across samples
xset2 <- group(xset1, bw = 10, minfrac = 0.5, mzwid = 0.025)
#Fill in missing peaks and calculate statistics
xset3 <- fillPeaks(xset2) 
#Generate feature tables and EICs
##change class lavels based on pairwise comparison
dr <- diffreport(xset3, filebase ="HighN_rmy_vs_media", eicmax=100)
write.csv(dr, file = "HighN_rmy_vs_media.csv")

