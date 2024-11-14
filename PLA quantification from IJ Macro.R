


list.files<-list.files(pattern="1d.+.csv")
list.data<-list()
for (i in 1:length(list.files)){
  list.data[[i]]<-read.csv(list.files[i],header=TRUE,sep=",",row.names = 1)
}
names(list.data)<-list.files

list.data_30d[[3]][[1]]
results_vector_1d <- character()
k=0
for (k in 1:length(list.data)){
  results_vector_1d <- append(results_vector_1d, list.data[[k]][[1]])
}


low_flow_events <- results_vector_1d[c(grep("blue", results_vector_1d, invert = T))]
low_flow_events_gsub <- gsub(".*:","",low_flow_events)
write.table(table(low_flow_events_gsub), "1d_ch1_counts_nuclear-events.csv")


list.files_30d<-list.files(pattern="30d.+.csv")
list.data_30d<-list()
for (j in 1:length(list.files_30d)){
  list.data_30d[[j]]<-read.csv(list.files_30d[j],header=TRUE,sep=",",row.names = 1)
}
names(list.data_30d)<-list.files_30d

results_vector_30d <- character()
l=0
for (l in 1:length(list.data_30d)){
  results_vector_30d <- append(results_vector_30d, list.data_30d[[l]][[1]])
}


high_flow_events <- results_vector_30d[c(grep("blue", results_vector_30d, invert = T))]
high_flow_events_gsub <- gsub(".*:","",high_flow_events)
write.table(table(high_flow_events_gsub), "30d_ch1_counts_nuclear-events.csv")

