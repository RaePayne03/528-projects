data<-read.csv('normalized_expression.csv')
filter1<-read.csv('Project1_4_1.csv')
nrow(data)
ncol(data)
patients<-colnames(data)



step1=data.frame(data)
step2=data.frame(filter1)

step1[,3:135]
vari<-c()
var2<-c()
for(i in 1:nrow(step1)){
  vari<-append(vari,var(as.numeric(step1[i,2:135])))
}

for(i in 1:nrow(step2)){
  var2<-append(var2,var(as.numeric(step2[i,3:136])))
}


med<-median(vari)
med2<-median(var2)


chist<-133*(var2/med2)
chist1<-133*(var2/med2)^2
chim<-133*((var2-med2)^2)/med2
chim1<-133*((vari-mean(vari))^2)/mean(vari)

qchisq(0.01,133,lower.tail=FALSE)



ind<-c()
for(i in 1:length(chist)){
  
  if(chist[i]>qchisq(0.01,133,lower.tail=FALSE)){
    ind<-append(ind,i)
  }
  
}

ind.1<-c()
for(i in 1:length(chist1)){
  
  if(chist1[i]>qchisq(0.01,133,lower.tail=FALSE)){
    ind.1<-append(ind.1,i)
  }
  
}

ind.2<-c()
for(i in 1:length(chim)){
  
  if(chim[i]>qchisq(0.01,133,lower.tail=FALSE)){
    ind.2<-append(ind.2,i)
  }
  
}

ind.3<-c()
for(i in 1:length(chim1)){
  
  if(chim1[i]>qchisq(0.01,133,lower.tail=FALSE)){
    ind.3<-append(ind.3,i)
  }
  
}



covar<-c()
for(i in 1:nrow(step2)){
  covar<-append(covar,var(as.numeric(step2[i,3:136]))/mean((as.numeric(step2[i,3:136]))))
}

ind3<-c()
for(i in 1:length(covar)){
  
  if(covar[i]>0.186){
    ind3<-append(ind3,i)
  }
  
}

fi<-intersect(ind3,ind)
fi.1<-intersect(ind3,ind.1)
fi.2<-intersect(ind3,ind.2)
fi.3<-intersect(ind3,ind.3)



length(intersect(ind3,ind))

length(intersect(ind3,ind.1))

length(intersect(ind3,ind.2))

length(intersect(ind3,ind.3))

filtered<-(step2[intersect(ind3,ind.1),])
dist_mat <- dist(t(filtered[,3:136]), method = 'euclidean')
agnes(x = t(filtered[, 3:136]), method = "average")$ac 
agnes(x = t(filtered[, 3:136]), method = "complete")$ac 
agnes(x = t(filtered[, 3:136]), method = "single")$ac
agnes(x = t(filtered[, 3:136]), method = "ward")$ac 
agnes(x = t(filtered[, 3:136]), method = "ward")
hclust_avg <- hclust(dist_mat, method = 'ward.D')
par(mar=c(1, 1, 1, 1))
plot(hclust_avg)
cut<-cutree(hclust_avg,k=2)
sum(cut==1)
sum(cut==2)
patients<-colnames(filtered)
patients<-patients[3:length(patients)]
cluster1<-c()
cluster2<-c()
for(i in 1:length(patients)){
  if(cut[i]==1){
    cluster1<-append(cluster1,patients[i]) 
  }
  if(cut[i]==2){
    cluster2<-append(cluster2,patients[i]) 
  }
}
cluster.1<-filtered[,c("Unnamed..0" ,cluster1)]
cluster.2<-filtered[,c("Unnamed..0" ,cluster2)]
welchn<-c()
welchstat<-c()
welchp<-c()
welchn<-cluster.1[,1]
for(i in 1:nrow(cluster.1)){
  wel<-t.test(cluster.1[i,2:ncol(cluster.1)],cluster.2[i,2:ncol(cluster.2)])
  welchstat<-append(welchstat,as.numeric(wel['statistic']))
  welchp<-append(welchp,wel['p.value'])
} 
welchpa<-p.adjust(welchp)
welch<-data.frame(matrix(ncol=0,nrow=length(welchn)))
welch$name<-welchn
welch$t.statistic<-welchstat
welch$p.value<-welchp
welch$p.adjusted<-welchpa

welch<-as.data.frame(welch)
weg <- as.data.frame(lapply(welch, unlist))

differ<-c()
for(i in 1:nrow(welch)){
  if(welch[i,'p.adjusted']<0.05){
    differ<-append(differ,i)
  }
}
differential_expression_results<-data.frame(weg[differ,])

write.csv(differential_expression_results,'differential_expression_results.csv')


write.csv(weg,'welch.csv')
write.csv(step2[ind.1,],'output4_5.csv')
write.csv(filtered,'output4_4.csv')


filtered<-(step2[intersect(ind.1,ind.1),])
dist_mat <- dist(t(filtered[,3:136]), method = 'euclidean')
agnes(x = t(filtered[, 3:136]), method = "average")$ac 
agnes(x = t(filtered[, 3:136]), method = "complete")$ac 
agnes(x = t(filtered[, 3:136]), method = "single")$ac
agnes(x = t(filtered[, 3:136]), method = "ward")$ac 
agnes(x = t(filtered[, 3:136]), method = "ward")
hclust_avg <- hclust(dist_mat, method = 'ward.D')
par(mar=c(1, 1, 1, 1))
cut<-cutree(hclust_avg,k=2)
sum(cut==1)
sum(cut==2)
patients<-colnames(filtered)
patients<-patients[3:length(patients)]
cluster1<-c()
cluster2<-c()
for(i in 1:length(patients)){
  if(cut[i]==1){
    cluster1<-append(cluster1,patients[i]) 
  }
  if(cut[i]==2){
    cluster2<-append(cluster2,patients[i]) 
  }
}
cluster.1<-filtered[,c("Unnamed..0" ,cluster1)]
cluster.2<-filtered[,c("Unnamed..0" ,cluster2)]
welchn<-c()
welchstat<-c()
welchp<-c()
welchn<-cluster.1[,1]
for(i in 1:nrow(cluster.1)){
  wel<-t.test(cluster.1[i,2:ncol(cluster.1)],cluster.2[i,2:ncol(cluster.2)])
  welchstat<-append(welchstat,as.numeric(wel['statistic']))
  welchp<-append(welchp,wel['p.value'])
} 
welchpa<-p.adjust(welchp)
welch<-data.frame(matrix(ncol=0,nrow=length(welchn)))
welch$name<-welchn
welch$t.statistic<-welchstat
welch$p.value<-welchp
welch$p.adjusted<-welchpa

welch<-as.data.frame(welch)
weg <- as.data.frame(lapply(welch, unlist))

differ<-c()
for(i in 1:nrow(welch)){
  if( welch[i,'p.adjusted']<0.05){
    differ<-append(differ,i)
  }
}
differential_expression_results_5_6<-data.frame(weg[differ,])
write.csv(differential_expression_results_5_6,'differential_expression_results_5_6.csv')
write.csv(weg,'welch_5_6.csv')

cancer<-read.csv('proj_metadata.csv')
hot<-(as.matrix(filtered[,3:ncol(filtered)]))
colon<-(filtered[,3:ncol(filtered)])
colm<-(cancer[,'cit.coloncancermolecularsubtype'])
colors<-c()
for(i in 1:length(colm)){
  if(colm[i]=='C3'){colors<-append(colors,'red')}
  if(colm[i]!='C3'){colors<-append(colors,'blue')}
}

heatmap.2(hot,dendrogram = 'none',ColSideColors = colors)
# package.install('superheat')
library(superheat)
superheat(hot,
          # scale the matrix columns
          scale = TRUE,
          # add row dendrogram
          n.clusters.cols = 2,
          n.clusters.rows = 2)
