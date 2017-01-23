source("supplychaindefs.r")
# A. Process example 
S<-supplychainsimple("default")
E<-example("aimable",S)
cat("\nExample AIMABLE\n")
source("resourceefficiency_analysis.r")
source("resourceuseefficiency_report.r")

E<-example("manureexport",S)
cat("\nExample Manure not recycled but exported\n")
source("resourceefficiency_analysis.r")
source("resourceuseefficiency_report.r")

E<-example("nofoodres",S)
cat("\nExample No Food residues recycling - only food product\n")
source("resourceefficiency_analysis.r")
source("resourceuseefficiency_report.r")

E<-example("norecycling",S)
cat("\nExample No recycling - only food product and manure export\n")
source("resourceefficiency_analysis.r")
source("resourceuseefficiency_report.r")

E<-example(c("largemanurerecycling"),S)
cat("\nExample Virtually all livestock coproduct recycled\n")
source("resourceefficiency_analysis.r")
source("resourceuseefficiency_report.r")
