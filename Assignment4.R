#  << Population Genetics  >>
#  The files hapmap3 r2 b36 chr15 CEU.haps and hapmap3 r2 b36 chr20 CEU.haps contain 
#  Single-Nucleotide-Polymorphisms (SNPs) covering 500-kb regions of chromosomes 15 and 20, 
#  respectively, which I have selected from phased haplotypes of 117 individuals with northwest
#  European ancestry that were sampled as part of Phase 3 of the HapMap project
# (http://www.hapmap.org). The files hapmap3 r2 b36 chr15 PopX.haps and hapmap3 
# r2 b36 chr20 PopX.haps contain data at the same SNPs from 174 individuals with unknown ancestry.
# Within each file, rows represent SNPs in this region, and each column is a distinct haplotype, 
# with every two consecutive columns representing the DNA from a single diploid individual.
# The two possible allele types at each SNP are coded as {0,1}. 

library(ggplot2)
ceu2 <- read.table("/Users/richardtsai/Documents/DataScience/Classes/Genetics/HW/Assignment 4/hapmap3_r2_b36_chr15_CEU.haps")
ceu  <- read.table("/Users/richardtsai/Documents/DataScience/Classes/Genetics/HW/Assignment 4/hapmap3_r2_b36_chr20_CEU.haps")
Popx <- read.table("/Users/richardtsai/Documents/DataScience/Classes/Genetics/HW/Assignment 4/hapmap3_r2_b36_chr20_PopX.haps")
Popx2 <- read.table("/Users/richardtsai/Documents/DataScience/Classes/Genetics/HW/Assignment 4/hapmap3_r2_b36_chr15_PopX.haps") 

################# question 1-a finding mutation rate ###################
# For each of chromosomes 15 and 20, infer the mutation rate for each of populations CEU and PopX.

mutation_rate <- function(n,snp){
   a = 0
   for ( i in 1:(n-1))  {
      a = a + 1/i
   }
  print(snp/a) 
}

# chr 20 (ceu) theta =a 73.96  n = 234, snp = 446
# chr 15 (ceu2) theta = 37.15  n = 234 snp =224
# chr 20 (Popx)  theta = 69.38 n= 348  snp = 446
# chr 15 (Popx2) theta = 34.85  n= 348 snp = 224

################## Question1-b  ####################
# Separately for each of chromosomes 15 and 20,
# and for each population CEU and PopX, plot the distribution of p-values when 
# esting each SNP for Hardy- Weinburg-Equilibrium. 

hwe.calc=function(x)
{
  geno.X=apply(cbind(x[seq(1,length(x),2)],x[seq(2,length(x),2)]),1,sum)
  p_X=length(x[x==0])/length(x)
  expected.vals=c(p_X^2,2*p_X*(1-p_X),(1-p_X)^2)*length(geno.X)
  observed.vals=c(length(geno.X[geno.X==0]),length(geno.X[geno.X==1]),length(geno.X[geno.X==2]))
  pearson.chi.sq=sum(((observed.vals-expected.vals)^2)/expected.vals)
  p.val=1-pchisq(pearson.chi.sq,1)
  #return(pearson.chi.sq)
  return(p.val)
}

num.snps = dim(ceu)[2]
hwe.snps.ceu = rep( NA, num.snps)
for (i in 1: num.snps ) hwe.snps.ceu[i] = hwe.calc(ceu[,i])
par(mfrow=c(2,2)) ## plots 2 figures per row, and 2 per column
hist(hwe.snps.ceu,main="20 chr CEU",xlab="HWE scores across SNPs")

num.snps2 = dim(ceu2)[2]
hwe.snps.ceu2 = rep( NA, num.snps2)
for (i in 1: num.snps2 ) hwe.snps.ceu2[i] = hwe.calc(ceu2[,i])
hist(hwe.snps.ceu2,main="15 chr CEU",xlab="HWE scores across SNPs")

num.snps.popx = dim(Popx)[2]
hwe.snps.popx = rep( NA, num.snps.popx)
for (i in 1: num.snps.popx ) hwe.snps.popx[i] = hwe.calc(Popx[,i])
hist(hwe.snps.popx,main="20 chr Popx",xlab="HWE scores across SNPs")

num.snps.popx2 = dim(Popx2)[2]
hwe.snps.popx2 = rep( NA, num.snps.popx2)
for (i in 1: num.snps.popx2 ) hwe.snps.popx2[i] = hwe.calc(Popx2[,i])
hist(hwe.snps.popx2,main="15 chr Popx",xlab="HWE scores across SNPs")

all_CEU = cbind(ceu,ceu2)
all_POPX = cbind(Popx,Popx2)

num.snps.all.popx = dim(all_POPX)[2]
hwe.snps.all.popx = rep( NA, num.snps.all.popx)
for (i in 1: num.snps.all.popx) hwe.snps.all.popx[i] = hwe.calc(all_POPX[,i])
#hist(hwe.snps.all.popx,main="All Popx",xlab="HWE scores across SNPs")

num.snps.all.ceu = dim(all_CEU)[2]
hwe.snps.all.ceu = rep( NA, num.snps.all.ceu)
for (i in 1: num.snps.all.ceu) hwe.snps.all.ceu[i] = hwe.calc(all_CEU[,i])
#hist(hwe.snps.all.ceu,main="All CEU",xlab="HWE scores across SNPs")

par(mfrow=c(1,1)) ## goes back to one figure total
boxplot(hwe.snps.all.ceu, hwe.snps.all.popx,
        ylab="HWE p-values",names=c("CEU","POPX"))

################## Question 1-c  #####################
# Separately within each chromosome and each population, explore linkage disequilibrium (LD) 
# among all pairs of SNPs using both r2 and |D′|. Also, for each population, calculate r2 and 
# |D′| between all pairs of SNPs on different chromosomes (i.e. one SNP on chromosome 15,
# the other on chromosome 20)


d.prime.calc=function(x,y)
{
  D.00=length(x[x==0 & y==0])/length(x)-(length(x[x==0])/
                                           length(x))*(length(y[y==0])/length(y))
  D.minus=min((length(x[x==1])/length(x))*(length(y[y==1])/length(y)),
              (length(x[x==0])/length(x))*(length(y[y==0])/length(y)))
  D.plus=min((length(x[x==1])/length(x))*(length(y[y==0])/length(y)),
             (length(x[x==0])/length(x))*(length(y[y==1])/length(y)))
  if (D.00>=0) D.prime=D.00/D.plus
  if (D.00<0) D.prime=D.00/D.minus
  return(abs(D.prime))
}

##################### D prime for Europe 20 chr  ############################
D.prime.ceu=matrix(NA,nrow=dim(ceu)[2],ncol=dim(ceu)[2])
for (i in 1:(dim(ceu)[2]-1)){
  for(j in (i+1):dim(ceu)[2]){
     D.prime.ceu[i,j]=d.prime.calc(ceu[,i],ceu[,j])
  }  
}
#scatterplot_ for D-prime Europe 20 chr   
sca_distance_ceu20=c()
sca_d_prime_ceu20=c()
for (i in 1:(dim(ceu)[2]-1)){
  for(j in (i+1):dim(ceu)[2]){
    sca_distance_ceu20 = c( sca_distance_ceu20 , j-i ) 
    sca_d_prime_ceu20 =  c( sca_d_prime_ceu20 , D.prime.ceu[i,j])
  }  
}
par(mfrow=c(1,1)) 
ceu20_dprime_data<- data.frame(sca_distance_ceu20,sca_d_prime_ceu20)
ggplot(ceu20_dprime_data, aes(sca_distance_ceu20,sca_d_prime_ceu20)) + geom_point() + geom_smooth()+xlab("Distance between SNP(SNP number)")+ ylab("D prime") + ggtitle(" Northwest European chromosome 20 ")

### Correlation for Europe 20 chr 
correlation.ceu=matrix(NA,nrow=dim(ceu)[2],ncol=dim(ceu)[2])
for (i in 1:(dim(ceu)[2]-1)){
  for(j in (i+1):dim(ceu)[2]){
    correlation.ceu[i,j]= cor(ceu[,i],ceu[,j])^2
  }  
}
#scatterplot_ for r2 Europe 20 chr   
sca_distance_ceu20=c()
sca_r2_ceu20=c()
for (i in 1:(dim(ceu)[2]-1)){
  for(j in (i+1):dim(ceu)[2]){
    sca_distance_ceu20 = c( sca_distance_ceu20 , j-i ) 
    sca_r2_ceu20 =  c( sca_r2_ceu20 , correlation.ceu[i,j])
  }  
}
par(mfrow=c(1,1)) 
ceu20_r2_data<- data.frame(sca_distance_ceu20,sca_r2_ceu20)
ggplot(ceu20_r2_data, aes(sca_distance_ceu20,sca_r2_ceu20)) + geom_point() + geom_smooth()+xlab("Distance between SNP(SNP number)")+ ylab("R^2") + ggtitle(" Northwest European chromosome 20 squared correlation coefficient ")


####################### D prime for Europe 15 chr #############################
D.prime.ceu15=matrix(NA,nrow=dim(ceu2)[2],ncol=dim(ceu2)[2])
for (i in 1:(dim(ceu2)[2]-1)){
  for(j in (i+1):dim(ceu2)[2]){
    D.prime.ceu15[i,j]=d.prime.calc(ceu2[,i],ceu2[,j])
  }  
}
#scatterplot_ for D-prime Europe 15 chr   
distance_ceu15=c()
d_prime_ceu15=c()
for (i in 1:(dim(ceu2)[2]-1)){
  for(j in (i+1):dim(ceu2)[2]){
    distance_ceu15 = c( distance_ceu15 , j-i ) 
    d_prime_ceu15 =  c( d_prime_ceu15 , D.prime.ceu15[i,j])
  }  
}
par(mfrow=c(1,1)) 
Data<- data.frame(distance_ceu15,d_prime_ceu15)
ggplot(Data, aes(distance_ceu15,d_prime_ceu15)) + geom_point() + geom_smooth()+xlab("Distance between SNP(SNP number)")+ ylab("D-prime") + ggtitle(" Northwest European chromosome 15 D_prime / distance ")

### Correlation for Europe 15 chr 
cor.ceu15=matrix(NA,nrow=dim(ceu2)[2],ncol=dim(ceu2)[2])
for (i in 1:(dim(ceu2)[2]-1)){
  for(j in (i+1):dim(ceu2)[2]){
    cor.ceu15[i,j]=cor(ceu2[,i],ceu2[,j])^2
  }  
}

#scatterplot_ for cor Europe 15 chr   
distance_ceu15=c()
cor_ceu15=c()
for (i in 1:(dim(ceu2)[2]-1)){
  for(j in (i+1):dim(ceu2)[2]){
    distance_ceu15 = c( distance_ceu15 , j-i ) 
    cor_ceu15 =  c( cor_ceu15 , cor.ceu15[i,j])
  }  
}
par(mfrow=c(1,1)) 
Data<- data.frame(distance_ceu15,cor_ceu15)
ggplot(Data, aes(distance_ceu15,cor_ceu15)) + geom_point() + geom_smooth()+ geom_smooth()+xlab("Distance between SNP(SNP number)")+ ylab(" r^2 ") + ggtitle(" Northwest European chromosome 15 correlation / distance ")

################################ D prime for popX 20 chr  ####################################
D.prime.popx20=matrix(NA,nrow=dim(Popx)[2],ncol=dim(Popx)[2])
for (i in 1:(dim(Popx)[2]-1)){
  for(j in (i+1):dim(Popx)[2]){
    D.prime.popx20[i,j]=d.prime.calc(Popx[,i],Popx[,j])
  }  
}

#scatterplot_ for dprime popx 20 chr   
distance_pop20=c()
dprime_pop20=c()
for (i in 1:(dim(Popx)[2]-1)){
  for(j in (i+1):dim(Popx)[2]){
    distance_pop20 = c( distance_pop20 , j-i ) 
    dprime_pop20 =  c( dprime_pop20 , D.prime.popx20[i,j])
  }  
}
par(mfrow=c(1,1)) 
Data<- data.frame(distance_pop20,dprime_pop20 )
ggplot(Data, aes(distance_pop20,dprime_pop20 )) + geom_point() + geom_smooth()+xlab("Distance between SNP(SNP number)")+ ylab(" d prime ") + ggtitle(" unknown population chromosome 20 d prime / distance ")

### correlation for popX 20 chr 
cor.popx20=matrix(NA,nrow=dim(Popx)[2],ncol=dim(Popx)[2])
for (i in 1:(dim(Popx)[2]-1)){
  for(j in (i+1):dim(Popx)[2]){
    cor.popx20[i,j]=cor(Popx[,i],Popx[,j])^2
  }  
}

#scatterplot_ for cor popx 20 chr   
distance_pop20=c()
cor_pop20=c()
for (i in 1:(dim(Popx)[2]-1)){
  for(j in (i+1):dim(Popx)[2]){
    distance_pop20 = c( distance_pop20 , j-i ) 
    cor_pop20 =  c( cor_pop20 , cor.popx20[i,j])
  }  
}
par(mfrow=c(1,1)) 
Data<- data.frame(distance_pop20,cor_pop20)
ggplot(Data, aes(distance_pop20,cor_pop20)) + geom_point() + geom_smooth()+xlab("Distance between SNP(SNP number)")+ ylab(" r^2 ") + ggtitle(" unknown population chromosome 20 correlation / distance ")

################ D prime for popX 15 chr 
D.prime.popx15=matrix(NA,nrow=dim(Popx2)[2],ncol=dim(Popx2)[2])
for (i in 1:(dim(Popx2)[2]-1)){
  for(j in (i+1):dim(Popx2)[2]){
    D.prime.popx15[i,j]=d.prime.calc(Popx2[,i],Popx2[,j])
  }  
}

#scatterplot_ for dprime popx 15 chr   
distance_pop15=c()
cor_pop15=c()
for (i in 1:(dim(Popx2)[2]-1)){
  for(j in (i+1):dim(Popx2)[2]){
    distance_pop15 = c( distance_pop15 , j-i ) 
    cor_pop15 =  c( cor_pop15 , D.prime.popx15[i,j])
  }  
}
par(mfrow=c(1,1)) 
Data<- data.frame(distance_pop15,cor_pop15)
ggplot(Data, aes(distance_pop15,cor_pop15)) + geom_point() + geom_smooth()+xlab("Distance between SNP(SNP number)")+ ylab(" d prime ") + ggtitle(" unknown population chromosome 15 d prime / distance ")


### correlation for popX 15 chr 
cor.popx15=matrix(NA,nrow=dim(Popx2)[2],ncol=dim(Popx2)[2])
for (i in 1:(dim(Popx2)[2]-1)){
  for(j in (i+1):dim(Popx2)[2]){
    cor.popx15[i,j]=cor(Popx2[,i],Popx2[,j])^2
  }  
}

#scatterplot_ for r popx 15 chr   
distance_pop15=c()
cor_pop15=c()
for (i in 1:(dim(Popx2)[2]-1)){
  for(j in (i+1):dim(Popx2)[2]){
    distance_pop15 = c( distance_pop15 , j-i ) 
    cor_pop15 =  c( cor_pop15 , cor.popx15[i,j])
  }  
}
par(mfrow=c(1,1)) 
Data<- data.frame(distance_pop15,cor_pop15)
ggplot(Data, aes(distance_pop15,cor_pop15)) + geom_point() + geom_smooth()+xlab("Distance between SNP(SNP number)")+ ylab(" r^2 ") + ggtitle(" unknown population chromosome 15 r^2 / distance ")

### D prime for popX 15chr vs popX 20chr 
D.prime.popx15_20= matrix(NA,nrow=dim(Popx)[2],ncol=dim(Popx2)[2])
for (i in 1:dim(Popx)[2]){
  for(j in 1:dim(Popx2)[2]){
    D.prime.popx15_20[i,j]= d.prime.calc(Popx[,i],Popx2[,j])
  }  
}

pop15_20 <- c( )
for(i in 1:dim(Popx)[2] ){
    for( j in 1: dim(Popx2)[2]) {
      pop15_20 <-  c(pop15_20 , D.prime.popx15_20[i,j])
   }
}
ggplot() + aes(pop15_20)+ geom_histogram(binwidth=0.1, colour="black", fill="black")+xlab("d prime")+ ylab(" frequency ") + ggtitle(" D prime for unknown population 15chr vs 20chr  ")

### Correlation for popX 15chr vs popX 20chr 
cor.popx15_20= matrix(NA,nrow=dim(Popx)[2],ncol=dim(Popx2)[2])
for (i in 1:dim(Popx)[2]){
  for(j in 1:dim(Popx2)[2]){
    cor.popx15_20[i,j]= cor(Popx[,i],Popx2[,j])^2
  }  
}

pop15_20_cor <- c( )
for(i in 1:dim(Popx)[2] ){
  for( j in 1: dim(Popx2)[2]) {
    pop15_20_cor <-  c(pop15_20_cor , cor.popx15_20[i,j])
  }
}
ggplot() + aes(pop15_20_cor)+ geom_histogram(binwidth=0.001, colour="black", fill="black")+xlab("r^2")+ ylab(" frequency ") + ggtitle(" r^2 for unknown population 15chr vs 20chr  ")

### D prime for Ceu 15chr vs popX 20chr 
D.prime.ceu15_20= matrix(NA,nrow=dim(ceu)[2],ncol=dim(ceu2)[2])
for (i in 1:dim(ceu)[2]){
  for(j in 1:dim(ceu2)[2]){
    D.prime.ceu15_20[i,j]= d.prime.calc(ceu[,i],ceu2[,j])
  }  
}

ceu15_20_dprime <- c( )
for(i in 1:dim(ceu)[2] ){
  for( j in 1: dim(ceu2)[2]) {
    ceu15_20_dprime <-  c(ceu15_20_dprime ,  D.prime.ceu15_20[i,j])
  }
}
ggplot() + aes(ceu15_20_dprime)+ geom_histogram(binwidth=0.01, colour="black", fill="black")+xlab("d prime")+ ylab(" frequency ") + ggtitle(" dprime for European 15chr vs 20chr  ")

### Correlation for Ceu 15chr vs popX 20chr 
cor.ceu15_20= matrix(NA,nrow=dim(ceu)[2],ncol=dim(ceu2)[2])
for (i in 1:dim(ceu)[2]){
  for(j in 1:dim(ceu2)[2]){
    cor.ceu15_20[i,j]= cor(ceu[,i],ceu2[,j])^2
  }  
}

ceu15_20_cor <- c( )
for(i in 1:dim(ceu)[2] ){
  for( j in 1: dim(ceu2)[2]) {
    ceu15_20_cor <-  c(ceu15_20_cor , cor.ceu15_20[i,j])
  }
}
ggplot() + aes(ceu15_20_cor)+ geom_histogram(binwidth=0.001, colour="black", fill="black")+xlab("r^2")+ ylab(" frequency ") + ggtitle(" r^2 for European 15chr vs 20chr  ")

################################### Question 2 #########################################
#2. The region of chromosome 15 considered in question one contains the gene SLC24A5, 
# which is believed to have been affected by selection (for light skin pigmentation) in 
# some groups such as CEU. The file ldsel.R includes an R function ldsel that simulates 
# a population undergoing random mating with mutation and selection at two biallelic loci. 
# Notice that while the wf.R models a haploid population, to study selection we need to 
# consider diploid individuals and so nall here refers to the total number of diploid individuals. 
# The default is nall = 2000 individuals, which corresponds to 4000 haplotypes.
# The user-input value init is a vector that gives the relative frequency of the four haplotypes
# {AB,Ab,aB,ab} at generation 0; the default is {0.5,0.01,0.01,0.48}. One of these four
# haplotypes is undergoing selection, with the parameter s controlling the increase in fitness
# for the selected haplotype. (You can fix s=0.01, for example, to figure out which haplotype 
#   is selected.) The output is a vector with 6 elements that gives 
# final results at the end of ngen generations (default ngen is 500 generations), with the first 
# four elements giving the haplotype proportions for {AB,Ab,aB,ab} and the last two elements giving |D′| and r2 between the two loci. As long as to.plot in the ldsel function is set to “yes”, there will also be two plots. One gives the haplotype frequency trajectories over time, as well as the frequency trajectories of the allele frequencies of A and B. The other plot gives the values of linkage disequilibrium measures r2 and |D′| over this same time frame. Note that unlike with wf.R, you can only simulate one population at a time with ldsel, so that the plot can illustrate the frequencies of the different haplotypes. (Hint: You should first set to.plot equal to “no” before copy n’ pasting ldsel into R to run, anytime you want to simulate lots of populations in e.g. a for loop.) 
# Assume that selection has acted on the SLC24A5 gene for 100 generations up to today. Using the
# function ldsel, perform simulations that mimic data for CEU from question one to answer the 
# questions below. For simplicity, always use the default mu (i.e. mu=0) and ngen=100. Also,
# also use the default init unless otherwise noted (i.e. in part (d) below). As ldsel simulates
# data from 2 biallelic loci, it may make sense to compare your results to pairs of SNPs from 
# the CEU data. For simplicity, use pairs of SNPs from the CEU data that are physically close
# to each other, e.g. with no more than 10 other SNPs between them. 
# Be sure to justify your answers, again writing a 2-2.5 page report including figures
# (again 2-2.5 pages in total for ALL of question two). High precision is not needed for 
# these questions, but you need to say roughly how accurate are your answers (which you can explore
# by repetition of e.g. simulations). Explaining the reasoning behind the results you see and 
# choices you make is more important than simply reporting answers. 

ldsel = function(ngen=100,nall=2000,init=c(0.5,0.01,0.01,0.48),rho=0,s=0,mu=0){
  hp = matrix(0,ngen,4)
  p = rep(0,4)
  hp[1,] = init/sum(init)
  ldstat = matrix(0,ngen,2)
  p1 = hp[1,1]+hp[1,2]		# allele prop at loc 1
  p2 = hp[1,1]+hp[1,3]		# allele prop at loc 2
  for(i in 2:ngen){
    # haplotype proportions in next generation after recombination:
    pp = hp[i-1,]*(1-rho)+rho*c(p1*p2,p1*(1-p2),(1-p1)*p2,(1-p1)*(1-p2))
    # now let's have some mutation and selection:
    p[1] = sum(pp*c((1-mu)^2,mu*(1-mu),mu*(1-mu),mu^2))
    p[2] = sum(pp*c(mu*(1-mu),(1-mu)^2,mu^2,mu*(1-mu)))*(1+s)
    p[3] = sum(pp*c(mu*(1-mu),mu^2,(1-mu)^2,mu*(1-mu)))
    p[4] = sum(pp*c(mu^2,mu*(1-mu),mu*(1-mu),(1-mu)^2))
    # sample haplotypes in next generation and record counts
    tmp = sample(1:4,nall,repl=T,prob=p)
    for (j in 1:4) hp[i,j]=sum(tmp==j)/nall
    p1 = hp[i,1]+hp[i,2]		# allele prop at loc 1
    p2 = hp[i,1]+hp[i,3]		# allele prop at loc 2
    # compute D'_00
    D00 = hp[i,1]-p1*p2
    if(D00>0) ldstat[i,1] = D00/min(p1*(1-p2),p2*(1-p1))
    if(D00<=0) ldstat[i,1] = -D00/min(p1*p2,(1-p2)*(1-p1))
    # compute r^2_00
    ldstat[i,2] = D00^2/(p1*p2*(1-p1)*(1-p2))
  }
  pA=apply(hp[,1:2],1,sum)
  pB=apply(hp[,c(1,3)],1,sum)
  to.plot='yes'
  if (to.plot=='yes')
  {
    lwd.val=2
    par(mfrow=c(2,1),mar=c(3,2,2,1))
    matplot(hp[,1:4],type="l",xlim=c(0,ngen*1.2),ylim=c(0,1),lty=1,lwd=lwd.val,xlab="",ylab="haplotype proportion",main=paste("LD: rho=",rho,", s=",s,", mu=",mu,sep=''))
    lines(1:dim(hp)[1],pA,lty=1,lwd=2*lwd.val,col=5)
    lines(1:dim(hp)[1],pB,lty=1,lwd=2*lwd.val,col=6)
    legend(ngen,1,leg=c("AB","Ab","aB","ab","A","B"),lty=1,col=1:6,lwd=2,cex=0.5)
    ylim.plot=c(0,max(c(ldstat)))
    if (sum(is.na(ylim.plot))!=0) ylim.plot=c(0,1)
    matplot(ldstat[-1,],type="l",xlim=c(0,ngen*1.2),ylim=ylim.plot,lty=1,lwd=lwd.val,col=3:4,xlab="",ylab="")
    legend(ngen,max(ylim.plot),leg=c("|D'|","r^2"),lty=1,col=3:4,lwd=2,cex=0.5)
  }
  return(c(hp[ngen,],ldstat[ngen,]))
}

#r2=0 , d_prime=1 recombination is not high,,,, allele frequency few.... AB   r2=0.026 