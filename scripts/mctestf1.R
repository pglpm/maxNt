## this script computes the probability for the conditional frequencies
## frequency(symptom | snp)
## where each symptom can assume values listed in 'symptomvariants'
## and each snp can assume values (alleles) listed in 'snpvariants'
## such values can be tuples
#library('ggplot2')
#library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
library('doParallel')
library('foreach')
library('LaplacesDemon')
#library('dplyr')

#pre <- FALSE
pre <- TRUE

filename <- paste0(if(pre){'pre_'},'x4testmc_p1_')

cores <- 30
thinning <- 10
mciterations <- 100
mcstatus <- 10
mcdiscard <- 10
nn <- 10000
## f <- choose(n,(0:n)+1)
## f <- f/sum(f)
f <- unlist(read.csv('stensola_activity_freqs_bins417641.csv',header=FALSE,sep=','))
ind <- which(f>0)
n <- length(f)-1
L <- 10
T <- 417641
smooth <- nn*10
#prior <- (1-(0:nn)+nn)/((nn+1)*(nn+2))
poste <- unlist(read.csv('stensola_6mom_prior2_N10000.csv',header=FALSE,sep=','))
#dirich <- nn*(1/(nn+1)+unlist(read.csv('stensola_6mom_prior2_N1000.csv',header=FALSE,sep=',')))

nchains <- cores

#if(cores>1){registerDoMC(cores=cores)}


##condp <- function(n,a,nn,aa){exp(lchoose(n,a)+lchoose(nn-n,aa-a)-lchoose(nn,aa))}
## G <- t(sapply(0:n, function(a){
##     exp(lchoose(n,a)+lchoose(nn-n,(0:nn)-a)-lchoose(nn,(0:nn)))
## }))

# this is slower
#condtt <- sapply(0:nn, function(aa){
#    exp(lchoose(n,(0:n))+lchoose(nn-n,aa-(0:n))-lchoose(nn,aa))
#})

if(!pre){
## setup parallel processing
## `%dox%` <- `%do%`
 if(cores>1){
## `%dox%` <- `%dopar%`
 cl <- makeCluster(cores, outfile="")
     registerDoParallel(cl)
 }
print('calculating G...')
G <- foreach(a=(ind-1), .combine='cbind'#, .export=c('nn','n')
               ) %:% 
       foreach(aa=0:nn, .combine='c'#, .export=c('nn','n')
               ) %dopar%
        {exp(lchoose(n,a)+lchoose(nn-n,aa-a)-lchoose(nn,aa))}

print('...done')
if(cores>1){
    stopCluster(cl)
}
}


### Here is the set-up for the Markov-chain Monte Carlo sampling of the parameters
### Uses package LaplacesDemon
if(pre){

    PGF <- function(data){nn <- data$nn
    ff <- 0.9/(nn+1)+0.1*rdirichlet(1,rep(1,nn+1))
    (1-cumsum(ff[-(nn+1)]))/c(1,1-cumsum(ff[-c(nn,nn+1)]))
}

mydata <- list(y=1, PGF=PGF,
               parm.names=paste0('x',1:nn),
               mon.names=c(''),
               nn=nn,
               L=L,
               #f=f[ind],
               #G=G,
               #T=T,
               #poste=poste,
               smooth=smooth
               )

logprob <- function(parm,data){
    parm <- interval(parm,0,1)
    F <- c(1,cumprod(parm))*c(1-parm,1)
    nn <- data$nn
    LP <- #data$T*sum(data$f*log(crossprod(data$G, F)),na.rm=TRUE) 
        -data$L*sum(F*log(F),na.rm=TRUE) +
        sum((nn-(1:(nn-1)))*log(parm[1:(nn-1)])) -
        data$smooth*sum(diff(F,differences=4)^2)
    list(LP=LP, Dev=-2*LP, Monitor=1, yhat=1, parm=parm)
}
} else {
        PGF <- function(data){nn <- data$nn
    ff <- 0.9*poste+0.1*rdirichlet(1,rep(1,nn+1))
    (1-cumsum(ff[-(nn+1)]))/c(1,1-cumsum(ff[-c(nn,nn+1)]))
}

mydata <- list(y=1, PGF=PGF,
               parm.names=paste0('x',1:nn),
               mon.names=c(''),
               nn=nn,
               L=L,
               f=f[ind],
               G=G,
               T=T,
               #poste=poste,
               smooth=smooth
               )

logprob <- function(parm,data){
    parm <- interval(parm,0,1)
    F <- c(1,cumprod(parm))*c(1-parm,1)
    nn <- data$nn
    LP <- data$T*sum(data$f*log(crossprod(data$G, F)),na.rm=TRUE) -
        data$L*sum(F*log(F),na.rm=TRUE) +
        sum((nn-(1:(nn-1)))*log(parm[1:(nn-1)])) -
        data$smooth*sum(diff(F,differences=4)^2)
    list(LP=LP, Dev=-2*LP, Monitor=1, yhat=1, parm=parm)
}
}
        ##Initial.Values <- c(0,0) #log(mf/sum(mf)) #GIV(prob, mydata, n=1000, PGF=T)

print('running Monte Carlo...')
if(cores>1){
 cl <- makeCluster(cores, outfile="")
# registerDoParallel(cl)

    Initial.Values <- foreach(i=1:nchains,.combine='rbind')%do%{GIV(logprob, mydata, n=1000, PGF=TRUE)}
    
        usamplesp <- LaplacesDemon.hpc(logprob, mydata, Initial.Values,
                        Covar=NULL,
                        Thinning=thinning,
                        Iterations=mciterations, Status=mcstatus,
                        Chains=nchains,CPUs=nchains,LogFile=paste0(filename,'_LDlog'), #Packages=c('Matrix'),#Type="MPI",
                        ##Algorithm="RDMH"#, Specs=list(B=list(1:d,d1:d2,d3:dnp))
                        ##Algorithm="Slice", Specs=list(B=NULL, Bounds=c(0,1), m=Inf, Type="Continuous", w=0.001)
                        Algorithm="AFSS", Specs=list(A=mcdiscard, B=NULL, m=100, n=0, w=1)
                        )
        usamples <- Combine(usamplesp,mydata)

        stopCluster(cl)

} else{
        Initial.Values <- GIV(logprob, mydata, n=1000, PGF=TRUE)
        usamples <- LaplacesDemon(logprob, mydata, Initial.Values,
                        Covar=NULL,
                        Thinning=thinning,
                        Iterations=mciterations, Status=mciterations,
                        #Chains=nchains,CPUs=nchains,LogFile=paste0(filename,'_LDlog'), #Packages=c('Matrix'),#Type="MPI",
                        ##Algorithm="RDMH"#, Specs=list(B=list(1:d,d1:d2,d3:dnp))
                        ##Algorithm="Slice", Specs=list(B=NULL, Bounds=c(0,1), m=Inf, Type="Continuous", w=0.001)
                        Algorithm="AFSS", Specs=list(A=mcdiscard, B=NULL, m=100, n=0, w=1)
                        )
                                     }

print('...finished')
samples2 <- usamples$Posterior2
samples1 <- usamples$Posterior1
samplesf <- t(apply(samples1,1,function(v){c(1,cumprod(v))*c(1-v,1)}))
    
print('writing samples to file...')

#write.csv(samples1,'testmcsamples.csv',row.names=FALSE)

    write.table(samples1,file=paste0(filename,'samplesx_N',nn,'_L',L,'_s',dim(samples1)[1],'_a',smooth,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')
    write.table(samplesf,file=paste0(filename,'samplesf_N',nn,'_L',L,'_s',dim(samples1)[1],'_a',smooth,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

print(usamples$Rec.Thinning)
print(Consort(usamples))
