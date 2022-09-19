# Author: Julien Devilliers

library(ggplot2)

# Function ----
busco_graph <- function(table, tall=6, type='piechart', leg_pos='right',cols=NULL) {
  if(type=='piechart' | type=='barplot'){
    l = length(1:nrow(table))
    long=4*l
    identif = rep(NA,long)
    group = rep(c('S', 'D', 'F', 'M'),l)
    value = rep(NA, long)
    for (i in 1:nrow(table)){
      a=4*i-3
      b=4*i
      c=a+1
      d=a+2
      identif[a:b] = table[i, 1]
      value[a]=table[i,2]
      value[c]=table[i,3]
      value[d]=table[i,4]
      value[b]=table[i,5]
    }
    busco_table = data.frame(cbind(identif, value, group)) # new table for graph
    busco_table$val=as.numeric(busco_table$value) # set values as numeric (vs. characters)
    if(type=='piechart'){
      fig=ggplot(busco_table, aes(x="", y=val, group=group)) +
        geom_bar(aes(fill = factor(group, levels=c('M', 'F', 'D', 'S'))), width = 1, stat = "identity") +
        ylab(' ') + 
        xlab(' ')+
        coord_polar("y", start=0) + facet_wrap(~ identif, ncol=cols) +
        labs(fill='busco')+
        scale_fill_manual(values = c('brown3', 'yellow2', 'dodgerblue4', 'deepskyblue'), labels=c('Missing', 'Fragmentated', 'Duplicated', 'Single')) +
        theme_void()+
        theme(strip.text = element_text(size=tall), axis.text = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank(), legend.position = leg_pos)
        
    }
    if(type=='barplot'){
      fig = ggplot(busco_table, aes(y=identif, x=val)) + 
        geom_bar(aes(fill = factor(group, levels=c('M', 'F', 'D', 'S'))),stat = "identity",position = 'stack') + 
        ylab(' ') + 
        xlab('Percentage %') + 
        labs(fill='busco') +
        scale_fill_manual(values = c('brown3', 'yellow2', 'dodgerblue4', 'deepskyblue'), labels=c('Missing', 'Fragmentated', 'Duplicated', 'Single')) + 
        theme(axis.text.y = element_text(size=tall,hjust=1), legend.position = leg_pos)
    }
    return(fig) #output = graph
  }
  else{
    warning('type argument non-valid: put \'piechart\' or \'barplot\'')
  }
}

# Input : 
### table = data.frame with 5 columns: ID, S, D, F, M
# Arguments: 
### size = fontsize axis ID 
### type = 'barplot' or 'piechart'
### leg_pos = position legend following style ggplot2
### cols = number of columns for piecharts

busco_graph(busco_table, tall=6, type='piechart', leg_pos='right',cols=NULL)







