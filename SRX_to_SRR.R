# Author: Julien Devilliers

#### Library ----
library(rentrez)

#### Function SRX_to_SRR ----
SRX_to_SRR <- function(SRX){
  a=entrez_search(db='sra', term = SRX)
  b=entrez_summary(db="sra", id=a$ids)
  c=strsplit(b$runs, 'acc=')
  c=c[[1]][-1]
  for (j in 1:length(c)){
    c[j]=unlist(strsplit(c[j], '\"'))[2]
  }
  expID=SRX
  runID=paste(as.character(unlist(c)), collapse=", ")
  data.frame(expID, runID)
}


#### Run the function ----

## Input: 
# Accession numbers (SRA ID)
## Output: 
# data.frame SRA ID and run IDs

SRX_list=c("SRX3417227", "ERX2262916", "ERX2262915", "ERX9462204", "SRX8919952", "SRX1211589", "SRX1211611", "SRX1457528", "ERX149935", "ERX043006", "SRX3540532", "ERX095890", "ERX102235", "SRX7283575", "SRX9116083", "SRX1527573", "ERX149934", "ERX047884", "ERX087587", "ERX303999", "ERX304000", "SRX8230501", "SRX8230498", "SRX15182548", "SRX15182545", "SRX15182532", "SRX15182546", "SRX15182541", "SRX15182530", "SRX15182542", "SRX15182544", "SRX15182534", "SRX15182538", "SRX15182540", "SRX15182535", "SRX15182537", "SRX476023", "SRX478896", "SRX5185244", "SRX5228411", "SRX5228406", "SRX5228385", "SRX5228388", "SRX5228380", "SRX5228409", "SRX5228384", "SRX5228407", "SRX5228417", "SRX5228390", "SRX5228427", "SRX5228396", "SRX5228394", "SRX5228418", "SRX5228420", "SRX5228437", "SRX5228410", "SRX5228432", "SRX5228424", "SRX5228382", "SRX5228377", "SRX5228403", "SRX5228441", "SRX5228378", "SRX5228405")
list_SRR=do.call(rbind, lapply(SRX_list, SRX_to_SRR))


  
  




