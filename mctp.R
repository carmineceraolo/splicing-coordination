library(nparcomp)
library(dplyr)
library(splitstackshape)
library(MASS)
library(ggplot2)
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
  df<-expandRows(df,"coverage")
  return(df)
}

filter.df<-function(df){
  df[df$n_ex>=3,]
}

odds<-function(p1,p2){
  (p1/(1-p1))/(p2/(1-p2))
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


# TEST AFTER FILTERING out transcript with less than 3 exons

ensemble.filtered<-filter.df(ensemble)
ensemble.filtered.test<-mctp(n_ex_d~type,data=ensemble.filtered)
summary(ensemble.filtered.test)
plot(ensemble.filtered.test)

# test log 2 nu/nd

ensemble.filtered$ud<-(ensemble.filtered$n_ex_u+1)/(ensemble.filtered$n_ex_d+1)

ud.test<-mctp(ud~type,data=ensemble.filtered)
summary(ud.test)
plot(ud.test)

ud.test$Data.Info

log2(odds(0.5662993,0.4707851)) # IR,ES
log2(odds(0.5662993,0.4408684)) # IR,5AS
log2(odds(0.5662993,0.5220473)) # IR,3AS
     
ggplot(data=ensemble.filtered,aes(x=log2(ud),fill=type))+
  geom_density(alpha=0.6,bw=0.16)+
  facet_grid(type~.)

ggplot(data=ensemble.filtered,aes(x=log2(ud),fill=type))+
  geom_density(alpha=0.6,bw=0.16)

1-0.4653252

# the effect is not so big

norm.pos.test<-mctp(norm_pos~type,data=ensemble.filtered)

summary(norm.pos.test)
plot(norm.pos.test)

log2(((0.5621490)/(1-0.5621490))/((0.4770922)/(1-0.4770922)))

fractions(round(odds(0.5621490,0.4770922),2))

log2(odds(0.5621490,0.4770922))

odds(0.56)
