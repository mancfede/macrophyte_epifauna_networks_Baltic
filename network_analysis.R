rm(list=ls())

library(igraph)
library(vegan)
library(veganUtils)
library(DescTools)
library(ggplot2)
library(reshape2)
library(viridis)
library(bipartite)

setwd("~/Desktop/PROJECTS/Macrophyte_animal_networks/Tvärminne_sampling")

AllData<-read.csv("./DATA/INPUT/all_data_ab.csv",header=T)
head(AllData)
names(AllData)

##Summary

nrow(AllData)   #430 macrophyte samples
unique(AllData$macrophyte)   #15 macrophytes species
sum(AllData[,9:38])   #17329 epifauna individuals
ncol(AllData[,9:38])   #30 epifauna taxa

#n. unique associations in study area before filtering
AllMacro<-unique(AllData$macrophyte) 
AllInv<-colnames(AllData[,9:38])
n_m<-length(AllMacro)
n_i<-length(AllInv)

subData<-AllData[,c(5,9:38)]
names(subData)
m_UF<-aggregate(.~macrophyte,data=subData,FUN='mean') 
m_UF<-m_UF[,2:ncol(m_UF)]
m_UF<-1*(m_UF>0)  
sum(m_UF)     #227 total observed associations in study area
  
#create ID for sites (to merge subsites later)
names(AllData)
sites<-unique(AllData$site)
sites

SiteID<-c()
for (i in AllData$site){
  site<-c()
  if (i==sites[1]){site<-"SP"}
  if (i==sites[2]){site<-"HK"}
  if (i==sites[3]){site<-"AB"}
  SiteID<-c(SiteID,site)}

cbind(AllData$site,SiteID)
NetID2<-paste(AllData$date,SiteID,sep='_')
AllData<-cbind(NetID2,AllData)
head(AllData)
names(AllData)

AllData$NetID2
unique(table(AllData$macrophyte,AllData$NetID2))
sum(table(AllData$macrophyte,AllData$NetID2))
AllInv


#######################################

####Null model

#######################################

#function reshuffles matrix keeping row totals and col totals fixed  

null_alg<-function(mat_num){
  row_tots<-rowSums(mat_num)
  col_tots<-colSums(mat_num)
  R<-nrow(mat_num)
  C<-ncol(mat_num)
  
  l<-length(mat_num)
  null_mat<-matrix(rep(0,l),nrow=R)
  tot<-sum(mat_num)   #total epifauna abundance per site/season
  while(sum(null_mat)<tot){
    null_row_tot<-rowSums(null_mat)
    null_col_tot<-colSums(null_mat)
    r_val<-(1:R)[which(!null_row_tot==row_tots)]  #select rows where totals have not reach observed row totals
    c_val<-(1:C)[which(!null_col_tot==col_tots)]  #same
    
    if (length(r_val)>1){
      rand_r<-sample(r_val,1)  #if more than one selected row pick only one at random
    } 
    else {
      rand_r<-r_val[1]   #else pick the only selected row
    }
    
    if (length(c_val)>1){   
      rand_c<-sample(c_val,1)   #same for columns
    } 
    else {
      rand_c<-c_val[1]
    }
    
    null_mat[rand_r,rand_c]<-null_mat[rand_r,rand_c]+1   #increment selected cell by one
    #    if (any(rowSums(null_mat)>row_tots) | any(colSums(null_mat)>col_tots)){break}
    #    print(sum(null_mat))
  }
  return(null_mat)
}


names(AllData)
epi<-AllData[,c(1,6,10:39)]  #net ID (site), macrophytes, epifauna
names(epi)

Nets<-unique(epi$NetID2)

AllMats<-list() 
for (net in Nets){   #for every site/season (network)
  mat<-epi[which(epi$NetID==net),2:ncol(epi)]   #take subset occurrences corresponding to specific site/season
  tot_mac<-aggregate(.~macrophyte,data=mat,FUN='mean')  #mean epifauna abundance on all macrophyte species
  sig_mat<-tot_mac[,2:ncol(tot_mac)]*0
  rownames(sig_mat)<-tot_mac[,1]
  
  mat_num<-as.matrix(mat[,2:ncol(mat)])   #take only epifauna 
  mat_num[which(is.na(mat_num))]<-0
  
  for (rep in 1:1000){
    null_mat<-null_alg(mat_num)  #randomize site/season matrix with null algorithm (considering replicates separately)
    null_mat<-cbind(mat[,1],null_mat)   #add macrophytes
    null_mat<-data.frame(null_mat)
    null_mat[,2:ncol(null_mat)]<-sapply(null_mat[,2:ncol(null_mat)],as.numeric)
    names(null_mat)<-names(mat)
    tot_mac_null<-aggregate(.~macrophyte,data=null_mat,FUN='mean')  #mean epifauna abundance on all macrophyte species (null matrix)
    sig_mat<-sig_mat+1*(tot_mac_null[,2:ncol(tot_mac_null)]<tot_mac[,2:ncol(tot_mac)])
    #print(rep)
  }
  sig_mat<-sig_mat/rep
  sig_mat<-1*(sig_mat>0.95)
  
  delc<-which(colSums(sig_mat)==0)
  if (length(delc)>0){sig_mat<-sig_mat[,-delc]}  #remove empty cols (if any)
  delr<-which(rowSums(sig_mat)==0)
  if (length(delr)>0){sig_mat<-sig_mat[-delr,]}  #remove empty rows (if any)
  dim(sig_mat)
  
  AllMats[[net]]<-sig_mat
  
  write.csv(sig_mat,paste0("./DATA/OUTPUT/matrix_nm2_",net,".csv"))  
  
  print(net)
}

AllMats

AllMacro<-unique(epi$macrophyte)
AllInv<-colnames(epi)[3:ncol(epi)]
n_m<-length(AllMacro)
n_i<-length(AllInv)

###Assemble meta-network

m_U<-matrix(rep(0,n_m*n_i),nrow=n_m)
rownames(m_U)<-AllMacro  #macrophytes
colnames(m_U)<-AllInv   #invertebrates
dim(m_U)

for (net in Nets){
  m<-as.data.frame(AllMats[[net]])
  for (ma in rownames(m)) {
    for (i in colnames(m)) {
      if (m[ma,i]>0) {m_U[ma,i]<-m_U[ma,i]+1}  #association present in metanetwork if observed at least once in local networks
    }
  }
}

m_U<-(m_U>0)*1

#check empty rows or columns
which(colSums(m_U)==0)   #5 epifauna taxa (mostly rare) have no significant associations
length(which(colSums(m_U)==0))/ncol(m_U)  #17% of total n. epifauna taxa 
which(rowSums(sig_mat)==0)  #no macrophyte without associations
dim(m_U)
delc<-which(colSums(m_U)==0)
if (length(delc)>0){m_U<-m_U[,-delc]}  #remove empty cols (if any)
dim(m_U)   #15 macrophyte species, 26 epifauna taxa with statistically significant associations

write.csv(m_U,"./DATA/OUTPUT/metanetwork_nm2.csv")

m_U<-read.csv("./DATA/OUTPUT/metanetwork_nm2.csv",header=T,row.names=1)

AllMacro<-rownames(m_U)
AllInv<-colnames(m_U)
n_m<-length(AllMacro)
n_i<-length(AllInv)
sum(m_U)  #total n. significant associations


#######################################

####Assemble local networks

#######################################

#common network layout
#all invertebrate and macrophyte species in the same position to visualize differences in network structure

lay_x0<-seq(5,n_i-5,length=length(AllMacro))   #space out macrophyte nodes
macro_xy<-list()
for (i in 1:n_m){
  macro_xy[[AllMacro[i]]]<-c(lay_x0[i],0)
}

lay_x1<-1:n_i
inv_xy<-list()
for (i in 1:n_i){
  inv_xy[[AllInv[i]]]<-c(lay_x1[i],1)
}

#function to make specific layout for each network
make_lay<-function(g){
  vn<-length(V(g))
  lay<-c()  #specific layout for each network
  for (i in 1:vn){
    if (V(g)$type[i]==TRUE){
      lay<-rbind(lay,macro_xy[[V(g)$name[i]]])
    }
    else{
      lay<-rbind(lay,inv_xy[[V(g)$name[i]]])
    }
  }  
  return(lay)
}

Nets<-unique(AllData$NetID2)

net_sum<-c()
AllNets<-list()  #store local network matrix, igraph object and layout in list
for (net in Nets){
  
  #m<-AllMats[[net]]
  
  #alternatively import matrices from csv
  m<-as.matrix(read.csv(paste0("./DATA/OUTPUT/matrix_nm2_",net,".csv"),header=T,row.names=1))

  net_sum<-rbind(net_sum,c(net,nrow(m),ncol(m),sum(m)/length(m)))
  
  #2) ASSEMBLE NETWORKS
  #adjust eastethics to save it in pdf format
  g<-graph_from_incidence_matrix(t(m),directed=T,mode=c("in"))  #Note: transpose matrix
  V(g)$type<-FALSE
  V(g)[which(V(g)$name %in% AllMacro)]$type<-TRUE
  V(g)$label.cex<-0.5  #to visualize taxa name
  V(g)$label.color<-"grey"
  #V(g)$label.dist<-2
  #V(g)$label.degree<-pi/2
  #V(g)$label<-""
  V(g)$color<-"darkblue"
  V(g)[V(g)$type==TRUE]$color<-"forestgreen"
  V(g)$size<-10
  V(g)$frame.width<-1
  V(g)$frame.color<-"white"
  E(g)$arrow.size<-0.5
  E(g)$width<-1.2
  
  #3) CREATE CUSTOM LAYOUT FOR NETWORK
  lay<-make_lay(g)
  
  #pdf(paste0("./FIGURES/morpho_network_",net,"_nm2.pdf"),width=10,height=8)
  plot(g,layout=lay,main=net)
  #dev.off()
  
  AllNets[[net]]<-list(m,g,lay)    #store matrices, networks and layouts in list (for every season and site)
  
  print(net)
  
} 

##Summary: n. nodes, connectance
net_sum
net_sum<-data.frame(net_sum)
colnames(net_sum)<-c("net_ID","n_macro","n_inv","connectance")
summary(net_sum)
round(sapply(net_sum,FUN="sd"),2)
cbind(Nets,net_sum)


#######################################

####Nestedness vs. modularity

#######################################

#null model algorithm (keep row and col totals fixed)
FF<-function(m){   
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:(5*l_hp)){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot,l_tot,replace=FALSE,prob=NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
    
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}

NetProp<-c()
for (net in Nets){
  m<-AllNets[[net]][[1]]
  g<-AllNets[[net]][[2]]
  
  nodf0<-nestednodf(m)$statistic[3]   #obs. nestedness
  
  comm<-computeModules(m)
  mod0<-comm@likelihood    #obs. modularity
  
  #compare to nestedness and modularity null networks
  nodf_nm<-c()
  mod_nm<-c()
  for (i in 1:1000){
    nm<-FF(m)   #null matrix
    nodf_nm<-c(nodf_nm,nestednodf(nm)$statistic[3])   #null nestedness
    ncomm<-computeModules(nm)    
    mod_nm<-c(mod_nm,ncomm@likelihood)   #null modularity
  }
  
  nodf_z<-(nodf0-mean(nodf_nm))/sd(nodf_nm)
  nodf_p<-sum(nodf_nm>=nodf0)/1000 
  
  mod_z<-(mod0-mean(mod_nm))/sd(mod_nm) 
  mod_p<-sum(mod_nm>=mod0)/1000
  
  NetProp<-rbind(NetProp,c(net,nodf0,nodf_z,nodf_p,mod0,mod_z,mod_p))
  
  #plot network with different node colors based on module membership
  #memb<-module2constraints(comm)
  #df<-data.frame(c(rownames(m),colnames(m)),memb) 
  #colnames(df)<-c("vertices","modules")
  
  #gc<-g   #make copy of network
  #memb2<-c()
  #for (v in V(g)$name){
  #  memb2<-c(memb2,df[df$vertices==v,2])
  #}
  #V(gc)$color<-viridis(max(memb2))[memb2]
  
  #pdf(paste0("./FIGURES/",net,"_nm2_mod_bipartite.pdf"),width=10,height=8)
  #plot(gc,layout=layout_with_sugiyama(gc),main=paste(net,"mod_z=",round(mod_z,2)))
  #dev.off()
  
  print(net)
}

NetProp
NetProp<-data.frame(NetProp)
str(NetProp)
names(NetProp)<-c("net_ID","nodf0","nodf_z","nodf_p","mod0","mod_z","mod_p")
NetProp[,2:7]<-sapply(NetProp[,2:7],as.numeric)
NetProp[,2:7]<-sapply(NetProp[,2:7],function(x) round(x,3))
NetProp


#######################################

####Robustness

#######################################

#run co-extinction simulations and calculate AUC

#random macrophyte loss
random_coex<-function(g,mat){
  inv<-colnames(mat)
  macro<-rownames(mat)
  coex<-c()
  for(i in 1:100){   #100 simulations with random removal
    g_r<-g  
    rem_step<-0
    coex<-rbind(coex,c(rem_step,1))
    for (sp in sample(macro)){    #remove macrophytes at random 
      to_del<-which(V(g_r)$name==sp)  
      g_r<-delete_vertices(g_r,to_del)  
      to_del<-which(igraph::degree(g_r,mode=c("in"))[V(g_r)$name %in% inv]==0)
      g_r<-delete_vertices(g_r,to_del)  
      inv_surv<-sum(V(g_r)$name %in% inv)/length(inv)
      rem_step<-rem_step+1
      coex<-rbind(coex,c(rem_step,inv_surv))
    }
  }
  mean_coex<-aggregate(coex[,2]~coex[,1],FUN="mean")   
  y<-mean_coex[,2]
  x<-mean_coex[,1]/length(macro)
  #plot(x,y,type='l')
  rob<-AUC(x,y)  #trapezoid (default)
  return (rob)
}

rob<-random_coex(g,m)

#to simulate random associations patterns between epifauna and macrophytes, reshuffle
#links in networks with null model EE, run co-extinction simulations and compute robustness

#randomize both rows and columns
EE<-function(mat){ 
  r<-dim(mat)[1]
  c<-dim(mat)[2]
  n_mat<-matrix(sample(mat),r,c)   #randomize matrix (considered as a vector)
  while (min(c(rowSums(n_mat),colSums(n_mat)))==0){  #process continues until a matrix with no null 
    #rows and columns is created
    n_mat<-matrix(sample(mat),r,c)}
  return(n_mat)}

m
EE(m)

#Apply to all networks
NetRob<-c()
for (net in Nets){
  m<-AllNets[[net]][[1]]
  g<-AllNets[[net]][[2]]
  inv<-colnames(m)
  macro<-rownames(m)
  rob0<-random_coex(g,m)   #observed robustness
  
  rob_nm<-c()
  for (i in 1:100){  #100 null models
    nm<-EE(m)  #null matrix (all links reshuffled with no constraints)
    rownames(nm)<-macro
    colnames(nm)<-inv
    ng<-graph_from_incidence_matrix(t(nm),directed=T,mode=c("in"))  #null graph
    rob_null<-random_coex(ng,nm)   
    rob_nm<-c(rob_nm,rob_null)
  }
  mean_rob_nm<-mean(rob_nm)
  rob_z<-(rob0-mean_rob_nm)/sd(rob_nm)
  rob_p<-sum(rob_nm>=rob0)/100
  
  NetRob<-rbind(NetRob,c(net,rob0,mean_rob_nm,rob_z,rob_p))
  print(net)
}

NetRob<-data.frame(NetRob)
names(NetRob)<-c("net_ID","rob0","mean_nrob","rob_z","rob_p")
NetRob[,2:5]<-sapply(NetRob[,2:5],as.numeric)
NetRob[,2:5]<-sapply(NetRob[,2:5],function(x) round(x,3))
NetRob

#merge datasets network properties
NetProp
NetRob
NetProp<-merge(NetProp,NetRob,by="net_ID",sort=F)

#write.csv(NetProp,"./DATA/OUTPUT/network_properties.csv",row.names=F)


###BARPLOTS

#order NetProp by site for plotting purposes
NetProp$net_ID
NetProp_v2<-NetProp[c(1,4,7,10,2,5,8,12,3,6,9,11),]
NetProp_v2
str(NetProp_v2)

#plot nestedness (z-score)
pdf("./FIGURES/nodf_z_nm2.pdf",height=6,width=10)
range(NetProp_v2$nodf_z)
barplot(NetProp_v2$nodf_z,ylim=c(-5,2))
abline(h=2,col='red',lwd=2)  #red lines indicate significant pattern
abline(h=-2,col='red',lwd=2)  
dev.off()
##Networks with Fucus (ABI, SP) appear less nested than expected by chance

#plot modularity (z-score)   
pdf("./FIGURES/mod_z_nm2.pdf",height=6,width=10)
range(NetProp_v2$mod_z)
barplot(NetProp_v2$mod_z,ylim=c(-2,3))
abline(h=2,col='red',lwd=2)  #red lines indicate significant pattern
abline(h=-2,col='red',lwd=2)  
dev.off()

#plot robustness (z-score)
pdf("./FIGURES/rob_z_nm2.pdf",height=6,width=10)
range(NetProp_v2$rob_z)
barplot(NetProp_v2$rob_z,ylim=c(-2,2))
abline(h=2,col='red',lwd=2)  #red lines indicate significant pattern
abline(h=-2,col='red',lwd=2)  
dev.off()



#######################################

####Link overlap

#######################################

##calculate % shared links 1) out of all the links and 2) out of all the links
#between species that occur in both networks

###TEMPORAL COMPARISON
#Interannual vs. interseasonal comparison

sites<-c("SP","HK","AB")

for (s in sites){
  print(s)
  
  ##Interannual comparison (autumn 2021-autumn 2022)
  print("Interannual")
  mO<-AllNets[[paste0("Oct_2021_",s)]][[1]]
  mS<-AllNets[[paste0("Sept_2022_",s)]][[1]]
  
  #matrix intersection
  l_mO<-melt(mO)
  l_mO<-subset(l_mO,value>0)[,1:2]
  
  l_mS<-melt(mS)
  l_mS<-subset(l_mS,value>0)[,1:2]
  
  sh_l<-merge(l_mO,l_mS,by=c("Var1","Var2"))  #shared links
  print(round(length(which(sh_l$Var1=="Fucus_vesiculosus"))/nrow(sh_l),2))
  print(sh_l)
  
  gO<-AllNets[[paste0("Oct_2021_",s)]][[2]]
  gS<-AllNets[[paste0("Sept_2022_",s)]][[2]]
  
  g_U<-union(gO,gS)    #all links and vertices at the site
  
  #plot network with all links and nodes (union), shared links in black (from intersection)
  V(g_U)$type<-FALSE
  V(g_U)[which(V(g_U)$name %in% AllMacro)]$type<-TRUE
  V(g_U)$label.cex<-0.5  #to visualize taxa name
  V(g_U)$label.color<-"grey"  #to visualize taxa name
  V(g_U)$color<-"darkblue"
  V(g_U)[V(g_U)$type==TRUE]$color<-"forestgreen"
  V(g_U)$size<-10
  V(g_U)$frame.width<-1
  V(g_U)$frame.color<-"white"
  E(g_U)$arrow.size<-0.5
  E(g_U)$width<-1.2
  df_gU<-as_data_frame(g_U)
  sh<-which(df_gU[,1] %in% sh_l[,1] & df_gU[,2] %in% sh_l[,2])  
  
  sh_links<-round((nrow(sh_l)/length(E(g_U)))*100,3)
  
  E(g_U)$color<-"lightgray"
  E(g_U)[sh]$color<-"black"
  
  lay<-make_lay(g_U)   
  
  pdf(paste0("./FIGURES/sh_links_nm2_annual_",s,".pdf"),width=10,height=8)
  plot(g_U,layout=lay,main=paste("% shared links:",sh_links))
  dev.off()
  
  ###
  
  ##Interseasonal
  print("Interseasonal")
  mJ<-AllNets[[paste0("Jul_2022_",s)]][[1]]
  mA<-AllNets[[paste0("Aug_2022_",s)]][[1]]
  mS<-AllNets[[paste0("Sept_2022_",s)]][[1]]
  
  #matrix intersection
  l_mJ<-melt(mJ)
  l_mJ<-subset(l_mJ,value>0)[,1:2]
  
  l_mA<-melt(mA)
  l_mA<-subset(l_mA,value>0)[,1:2]
  
  l_mS<-melt(mS)
  l_mS<-subset(l_mS,value>0)[,1:2]
  
  sh_l1<-merge(l_mJ,l_mA,by=c("Var1","Var2"))  #shared links
  sh_l<-merge(sh_l1,l_mS,by=c("Var1","Var2"))  #shared links
  print(round(length(which(sh_l$Var1=="Fucus_vesiculosus"))/nrow(sh_l),2))
  print(sh_l)
  
  gJ<-AllNets[[paste0("Jul_2022_",s)]][[2]]
  gA<-AllNets[[paste0("Aug_2022_",s)]][[2]]
  gS<-AllNets[[paste0("Sept_2022_",s)]][[2]]
  
  g_U1<-union(gJ,gA)    #all links and vertices in the site
  g_U<-union(g_U1,gS)
  
  sh_links<-round((nrow(sh_l)/length(E(g_U)))*100,3)
  
  #plot network with all links and nodes (union), shared links in black (from intersection)
  V(g_U)$type<-FALSE
  V(g_U)[which(V(g_U)$name %in% AllMacro)]$type<-TRUE
  V(g_U)$label.cex<-0.5  #to visualize taxa name
  V(g_U)$label.color<-"grey"  #to visualize taxa name
  V(g_U)$color<-"darkblue"
  V(g_U)[V(g_U)$type==TRUE]$color<-"forestgreen"
  V(g_U)$size<-10
  V(g_U)$frame.width<-1
  V(g_U)$frame.color<-"white"
  E(g_U)$arrow.size<-0.5
  E(g_U)$width<-1.2
  df_gU<-as_data_frame(g_U)
  sh<-which(df_gU[,1] %in% sh_l[,1] & df_gU[,2] %in% sh_l[,2])  
  E(g_U)$color<-"lightgray"
  E(g_U)[sh]$color<-"black"
  
  lay<-make_lay(g_U)     
  
  pdf(paste0("./FIGURES/sh_links_nm2_seas_",s,".pdf"),width=10,height=8)
  plot(g_U,layout=lay,main=paste("% shared links:",sh_links))
  dev.off()
}


###SPATIAL COMPARISON
#% shared associations between sites for each seasons
#plot links in different color if shared between AB-HK, AB-SP, HK-SP 

Nets
seasons<-c("Oct_2021","Jul_2022","Aug_2022","Sept_2022")

all_sh<-c()
for (se in seasons){
  
  mSP<-AllNets[[paste0(se,"_SP")]][[1]]
  mHK<-AllNets[[paste0(se,"_HK")]][[1]]
  mAB<-AllNets[[paste0(se,"_AB")]][[1]]
  
  l_mSP<-melt(mSP)
  l_mSP<-subset(l_mSP,value>0)[,1:2]
  
  l_mHK<-melt(mHK)
  l_mHK<-subset(l_mHK,value>0)[,1:2]
  
  l_mAB<-melt(mAB)
  l_mAB<-subset(l_mAB,value>0)[,1:2]
  
  #intersection
  sh_l_SPHK<-merge(l_mSP,l_mHK,by=c("Var1","Var2"))  #shared links SP-HK
  sh_l_ABSP<-merge(l_mAB,l_mSP,by=c("Var1","Var2"))  #shared links AB-SP
  sh_l_ABHK<-merge(l_mAB,l_mHK,by=c("Var1","Var2"))  #shared links AB-HK
  sh_l<-merge(sh_l_SPHK,l_mAB,by=c("Var1","Var2"))  #shared links all networks
  
  print(se)
  print(sh_l_SPHK)
  print(sh_l_ABSP)
  print(sh_l_ABHK)
  print(sh_l)
  
  gSP<-AllNets[[paste0(se,"_SP")]][[2]]
  gHK<-AllNets[[paste0(se,"_HK")]][[2]]
  gAB<-AllNets[[paste0(se,"_AB")]][[2]]
  
  #union
  g_U_SPHK<-union(gSP,gHK)   #all links SP-HK
  g_U_ABSP<-union(gAB,gSP)   #all links AB-SP
  g_U_ABHK<-union(gAB,gHK)   #all links AB-HK
  g_U<-union(g_U_SPHK,gAB)   #all links per season
  print(length(E(g_U)))
  
  sh_links_SPHK<-round((nrow(sh_l_SPHK)/length(E(g_U_SPHK)))*100,2)
  sh_links_ABSP<-round((nrow(sh_l_ABSP)/length(E(g_U_ABSP)))*100,2)
  sh_links_ABHK<-round((nrow(sh_l_ABHK)/length(E(g_U_ABHK)))*100,2)
  sh_links<-round((nrow(sh_l)/length(E(g_U)))*100,2)
  
  #plot network with all links and nodes (union), shared links in different colors (from intersection)
  V(g_U)$type<-FALSE
  V(g_U)[which(V(g_U)$name %in% AllMacro)]$type<-TRUE
  V(g_U)$label.cex<-0.5  #to visualize taxa name
  V(g_U)$label.color<-"grey"  #to visualize taxa name
  V(g_U)$color<-"darkblue"
  V(g_U)[V(g_U)$type==TRUE]$color<-"forestgreen"
  V(g_U)$size<-10
  V(g_U)$frame.width<-1
  V(g_U)$frame.color<-"white"
  E(g_U)$arrow.size<-0.5
  E(g_U)$width<-1.2
  df_gU<-as_data_frame(g_U)
  sh_SPHK<-which(df_gU[,1] %in% sh_l_SPHK[,1] & df_gU[,2] %in% sh_l_SPHK[,2]) 
  sh_ABSP<-which(df_gU[,1] %in% sh_l_ABSP[,1] & df_gU[,2] %in% sh_l_ABSP[,2]) 
  sh_ABHK<-which(df_gU[,1] %in% sh_l_ABHK[,1] & df_gU[,2] %in% sh_l_ABHK[,2]) 
  sh<-which(df_gU[,1] %in% sh_l[,1] & df_gU[,2] %in% sh_l[,2])
  E(g_U)$color<-"lightgray"
  E(g_U)[sh_SPHK]$color<-"#8C7AA9"
  E(g_U)[sh_ABSP]$color<-"#E58F65"
  E(g_U)[sh_ABHK]$color<-"#F9E784"
  E(g_U)[sh]$color<-"black"
  
  lay<-make_lay(g_U)     
  
  all_sh<-rbind(all_sh,c(se,sh_links_SPHK,sh_links_ABSP,sh_links_ABHK,sh_links))
  
  #pdf(paste0("./FIGURES/sh_links_nm2_spat_pair_comp_",se,".pdf"),width=10,height=8)
  plot(g_U,layout=lay)
  #dev.off()
  
  print(se)
}

all_sh<-data.frame(all_sh)
all_sh_v2<-all_sh[,2:5]
rownames(all_sh_v2)<-all_sh[,1]
colnames(all_sh_v2)<-c("SP-HK","AB-SP","AB-HK","all_sites")
all_sh_v2<-sapply(all_sh_v2,as.numeric)
t(all_sh_v2)
round(rowMeans(t(all_sh_v2)),2)
sd(t(all_sh_v2)[2,])

