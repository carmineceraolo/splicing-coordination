library(nparcomp)
library(dplyr)
?mctp

read.exontable<-function(sample.type=""){
  
  data_dir<-"/project/hfa_work/ceraolo/theco/tables/exon_counts/"
  
  event_types=c("IR","ES","5AS","3AS")
  
  
  df.l<-list()
  
  if (sample.type!="") {
    sample.type<-paste0("_",sample.type)
  }
  
  for (i in 1:4) {
    t<-event_types[i]
    fname<-paste0(data_dir,"exon_count_",t,sample.type,".csv")
    df.l[[i]]<-read.csv(fname)
    df.l[[i]]$type<-t
  }
  
  names(df.l)<-event_types
  
  df<-bind_rows(df.l)
  df$type<-factor(df$type)
  
  return(df)
}

### TESTING FOR ALL THE SAMPLES TOGETHER

ensemble<-read.exontable()

colnames(ensemble)

ensemble.test<-mctp(n_ex_d~type,data=ensemble)
summary(ensemble.test)
plot(ensemble.test)

### TESTING FOR CELL LINES ONLY

cell.line<-read.exontable(sample.type = "cell_line")

cell.line.test<-mctp(n_ex_d~type,data=cell.line)

summary(cell.line.test)
plot(cell.line.test)


### TESTING FOR TISSUES ONLY

tissue<-read.exontable(sample.type = "tissue")

tissue.test<-mctp(n_ex_d~type,data=tissue)

summary(tissue.test)
plot(tissue.test)


### TESTING FOR PRIMARY CELLS ONLY

primary.cells<-read.exontable(sample.type = "primary_cell")

primary.cells.test<-mctp(n_ex_d~type,data=primary.cells)

summary(primary.cells.test)
plot(primary.cells.test)

### TESTING FOR IN VITRO DIFFERENTIATED CELLS ONLY

in.vitro<-read.exontable(sample.type = "in_vitro_differentiated_cells")

in.vitro.test<-mctp(n_ex_d~type,data=in.vitro)

summary(in.vitro.test)
plot(in.vitro.test)
