library('doParallel')

cores <- 30
nn <- 10000
n <- 2000

condp <- function(n,a,nn,aa){exp(lchoose(n,a)+lchoose(nn-n,aa-a)-lchoose(nn,aa))}

## setup parallel processing
`%dox%` <- `%do%`
if(cores>1){
`%dox%` <- `%dopar%`
cl <- makeCluster(cores, outfile="")
registerDoParallel(cl)
}

data1a <- unlist(c(read.csv(paste0('stensoladist2mom',10000,'.csv'),header=F)))
data1b <- unlist(c(read.csv(paste0('stensoladist4mom',10000,'.csv'),header=F)))
data2a <- unlist(c(read.csv(paste0('riehledist2mom',10000,'.csv'),header=F)))
data2b <- unlist(c(read.csv(paste0('riehledist4mom',10000,'.csv'),header=F)))


result <- foreach(a=0:n, .combine='cbind') %do% {message(a)
    table <- foreach(aa=0:nn, .combine='c', .export=c('condp','nn','n','a','data1a','data1b','data2a','data2b')) %dox% {condp(n,a,nn,aa)}
    c(sum(table*data1a),
      sum(table*data1b),
      sum(table*data2a),
      sum(table*data2b))
}

write.csv(result[1,],paste0('stensola2mom_margto',n,'.csv'),row.names=F)
write.csv(result[2,],paste0('stensola4mom_margto',n,'.csv'),row.names=F)
write.csv(result[3,],paste0('riehle2mom_margto',n,'.csv'),row.names=F)
write.csv(result[4,],paste0('riehle4mom_margto',n,'.csv'),row.names=F)

if(cores>1){stopCluster(cl)}
