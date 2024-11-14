cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

#Working Directory
setwd("Z:/data/")

#Total OS

list.files_total_OS<-list.files(pattern="Total.+.OS", recursive = T)
list.filenames_total_OS<-list()
for (i in 1:length(list.files_total_OS)){
  list.filenames_total_OS[[i]]<-read.table(list.files_total_OS[i],header=TRUE,sep="\t",row.names = 1)
}
names(list.filenames_total_OS)<-list.files_total_OS

results_vector_OS_total <- character()
k=0
for (k in 1:length(list.files_total_OS)){
  results_vector_OS_total <- append(results_vector_OS_total, nrow(list.filenames_total_OS[[k]]))
}
results_OS_total <- as.data.frame(t(rbind(list.files_total_OS, results_vector_OS_total)))

write.table(results_OS_total, "PLA_OS_total-events.csv")

#Nuclear OS

list.files_nuclear_OS<-list.files(pattern="Nuclear.+.OS", recursive = T)
list.filenames_nuclear_OS<-list()
for (i in 1:length(list.files_nuclear_OS)){
  list.filenames_nuclear_OS[[i]]<-read.table(list.files_nuclear_OS[i],header=TRUE,sep="\t",row.names = 1)
}
names(list.filenames_nuclear_OS)<-list.files_nuclear_OS

results_vector_OS_nuclear <- character()
k=0
for (k in 1:length(list.files_nuclear_OS)){
  results_vector_OS_nuclear <- append(results_vector_OS_nuclear, nrow(list.filenames_nuclear_OS[[k]]))
}
results_OS_nuclear <- as.data.frame(t(rbind(list.files_nuclear_OS, results_vector_OS_nuclear)))

write.table(results_OS_nuclear, "PLA_OS_nuclear-events.csv")


#Total PS

list.files_Total_PS<-list.files(pattern="Total.+.PS", recursive = T)
list.filenames_Total_PS<-list()
for (i in 1:length(list.files_Total_PS)){
  list.filenames_Total_PS[[i]]<-read.table(list.files_Total_PS[i],header=TRUE,sep="\t",row.names = 1)
}
names(list.filenames_Total_PS)<-list.files_Total_PS

results_vector_PS_Total <- character()
k=0
for (k in 1:length(list.files_Total_PS)){
  results_vector_PS_Total <- append(results_vector_PS_Total, nrow(list.filenames_Total_PS[[k]]))
}
results_PS_Total <- as.data.frame(t(rbind(list.files_Total_PS, results_vector_PS_Total)))

write.table(results_PS_Total, "PLA_PS_Total-events.csv")

#Nuclear PS

list.files_nuclear_PS<-list.files(pattern="Nuclear.+.PS", recursive = T)
list.filenames_nuclear_PS<-list()
for (i in 1:length(list.files_nuclear_PS)){
  list.filenames_nuclear_PS[[i]]<-read.table(list.files_nuclear_PS[i],header=TRUE,sep="\t",row.names = 1)
}
names(list.filenames_nuclear_PS)<-list.files_nuclear_PS

results_vector_PS_nuclear <- character()
k=0
for (k in 1:length(list.files_nuclear_PS)){
  results_vector_PS_nuclear <- append(results_vector_PS_nuclear, nrow(list.filenames_nuclear_PS[[k]]))
}
results_PS_nuclear <- as.data.frame(t(rbind(list.files_nuclear_PS, results_vector_PS_nuclear)))

write.table(results_PS_nuclear, "PLA_PS_nuclear-events.csv")

