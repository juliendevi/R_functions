library(readr)

# Put a file from the annotation output folder
folder=dirname(file.choose())

#Function itself
get_information <- function(file){
  #Take ID OG
  ID = as.character(strsplit(file, '[/.]')[[1]][1])
  #Import EggNOG annotation table
  annot_table = read_delim(file, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 5, show_col_types = F)
  #Number of sequences OG
  numberseq = nrow(annot_table)-3
  if(numberseq>0){
    #Put the names of the columns
    colnames(annot_table)=c('query',	'seed_ortholog',	'evalue',	'score',	'eggNOG_OGs',	'max_annot_lvl',	'COG_category',	'Description',	'Preferred_name',	'GOs',	'EC',	'KEGG_ko',	'KEGG_Pathway',	'KEGG_Module',	'KEGG_Reaction',	'KEGG_rclass',	'BRITE',	'KEGG_TC',	'CAZy',	'BiGG_Reaction',	'PFAMs')
    #Get all unique lists
    Description = unique(annot_table$Description)
    Prefname = unique(annot_table$Preferred_name)
    PFAMs=unique(annot_table$PFAMs)
    #Make tables
    if(names(which(table(unlist(Prefname))==max(table(unlist(Prefname)))))[1]!='-'){
      pref_name=paste0(
        names(which(table(unlist(Prefname))==max(table(unlist(Prefname)))))[1], ', ',
        names(which(table(unlist(Prefname))==max(table(unlist(Prefname)))))[2])
    }
    else{
      pref_name=paste0(
        names(which(table(unlist(Prefname))==max(table(unlist(Prefname)))))[1], ', ',
        names(which(table(unlist(Prefname))==max(table(unlist(Prefname)))))[2], ', ',
        names(which(table(unlist(Prefname))==max(table(unlist(Prefname)))))[3], ', ',
        names(which(table(unlist(Prefname))==max(table(unlist(Prefname)))))[4]
      )
    }
    if (names(which(table(unlist(Description))==max(table(unlist(Description)))))[1]!='-'){
      descrip=names(which(table(unlist(Description))==max(table(unlist(Description)))))[1]
    }
    else{
      descrip=names(which(table(unlist(Description))==max(table(unlist(Description)))))[2]
    }
    if(names(which(table(unlist(PFAMs))==max(table(unlist(PFAMs)))))[1]!='-'){
      pfams=names(which(table(unlist(PFAMs))==max(table(unlist(PFAMs)))))[1]
    }
    else{
      pfams=names(which(table(unlist(PFAMs))==max(table(unlist(PFAMs)))))[2]
    }
    #List GO terms of interest
    global=c('GO:0050907', 'GO:0043695', 'GO:0009582', 'GO:0048512')
    chemicalstim=c('GO:0050969', 'GO:0050968', 'GO:0050911', 'GO:0050912')
    taste=c('GO:0001580', 'GO:0001583', 'GO:0001581', 'GO:0001582', 'GO:0046535')
    abiotic=c('GO:0009590', 'GO:0098513', 'GO:0009583', 'GO:0050982', 'GO:0016048', 'GO:0050981')
    #Extract GO characters
    GO_list_unique = unique(unlist(strsplit(annot_table$GOs, '[,]')))
    if (length(intersect(GO_list_unique, global))>0){
      GO_list=unlist(strsplit(annot_table$GOs, '[,]'))
      if (global[1] %in% intersect(GO_list_unique, global)){
        if(length(intersect(GO_list_unique, chemicalstim))>0){
          if(chemicalstim[1] %in% intersect(GO_list_unique, chemicalstim)){
            magnet=length(which(GO_list==chemicalstim[1]))
          }
          else{
            magnet=0
          }
          if(chemicalstim[2] %in% intersect(GO_list_unique, chemicalstim)){
            pain=length(which(GO_list==chemicalstim[2]))
          }
          else{
            pain=0
          }
          if(chemicalstim[3] %in% intersect(GO_list_unique, chemicalstim)){
            smell=length(which(GO_list==chemicalstim[3]))
          }
          else{
            smell=0
          }
          if(length(intersect(GO_list_unique, taste))>0){
            if(taste[1] %in% intersect(GO_list_unique, taste)){
              bitter=length(which(GO_list==taste[1]))
            }
            else{
              bitter=0
            }
            if(taste[2] %in% intersect(GO_list_unique, taste)){
              salty=length(which(GO_list==taste[2]))
            }
            else{
              salty=0
            }
            if(taste[3] %in% intersect(GO_list_unique, taste)){
              sour=length(which(GO_list==taste[3]))
            }
            else{
              sour=0
            }
            if(taste[4] %in% intersect(GO_list_unique, taste)){
              sweet=length(which(GO_list==taste[4]))
            }
            else{
              sweet=0
            }
            if(taste[5] %in% intersect(GO_list_unique, taste)){
              umami=length(which(GO_list==taste[5]))
            }
            else{
              umami=0
            }
          }
          else{
            bitter=salty=sour=sweet=umami=0
          }
        }
        else{
          magnet=pain=smell=bitter=salty=sour=sweet=umami=0
        }
      }
      else{
        magnet=pain=smell=bitter=salty=sour=sweet=umami=0
      }
      if(global[2] %in% intersect(GO_list_unique, global)){
        phero=length(which(GO_list==global[2]))
      }
      else{
        phero=0
      }
      if (global[3] %in% intersect(GO_list_unique, global)){
        if(length(intersect(GO_list_unique, abiotic))>0){
          if(abiotic[1] %in% intersect(GO_list_unique, abiotic)){
            gravity=length(which(GO_list==abiotic[1]))
          }
          else{
            gravity=0
          }
          if(abiotic[2] %in% intersect(GO_list_unique, abiotic)){
            humidity=length(which(GO_list==abiotic[2]))
          }
          else{
            humidity=0
          }
          if(abiotic[3] %in% intersect(GO_list_unique, abiotic)){
            light=length(which(GO_list==abiotic[3]))
          }
          else{
            light=0
          }
          if(abiotic[4] %in% intersect(GO_list_unique, abiotic)){
            mecha=length(which(GO_list==abiotic[4]))
          }
          else{
            mecha=0
          }
          if(abiotic[5] %in% intersect(GO_list_unique, abiotic)){
            temp=length(which(GO_list==abiotic[5]))
          }
          else{
            temp=0
          }
          if(abiotic[6] %in% intersect(GO_list_unique, abiotic)){
            elec=length(which(GO_list==abiotic[6]))
          }
          else{
            elec=0
          }
        }
        else{
          gravity=humidity=light=mecha=temp=elec=0
        }
      }
      else{
        gravity=humidity=light=mecha=temp=elec=0
      }
      if(global[4] %in% intersect(GO_list_unique, global)){
        circa=length(which(GO_list==global[4]))
      }
      else{
        circa=0
      }
    }
    else{
      magnet=pain=smell=bitter=salty=sour=sweet=umami=phero=gravity=humidity=light=mecha=temp=elec=circa=0
    }
    #Fill table
    data.frame(ID, numberseq, pref_name, descrip, pfams, magnet, pain, smell, bitter, salty, sour, sweet, umami, phero, gravity, humidity, light, mecha, temp, elec, circa)
  }
  else{
    pref_name=descrip=pfams=magnet=pain=smell=bitter=salty=sour=sweet=umami=phero=gravity=humidity=light=mecha=temp=elec=circa='###'
    numberseq='#!ERROR!#'
    data.frame(ID, numberseq, pref_name, descrip, pfams, magnet, pain, smell, bitter, salty, sour, sweet, umami, phero, gravity, humidity, light, mecha, temp, elec, circa)
  }
}


#Find all paths
paths=list.files(folder, pattern = "*.fa.eggnog.emapper.annotations")

#Make the table
try=do.call(rbind, lapply(paths, get_information))

library(googlesheets4)
#add it to googledrive
write_sheet(try, 
            ss = as_sheets_id('https://docs.google.com/spreadsheets/d/1gHc8D87coBrwIB56qF8bsNBllC0tqX6di-a60OBWrmo/edit#gid=1400871183'), 
            sheet='9august22')

