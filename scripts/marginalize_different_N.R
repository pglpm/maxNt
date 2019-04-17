library('doParallel')

cores <- 1

condp <- function(n,a,nn,aa){exp(lchoose(n,a)+lchoose(nn-n,aa-a)-lchoose(nn,aa))}

## setup parallel processing
`%dox%` <- `%do%`
if(cores>1){
`%dox%` <- `%dopar%`
cl <- makeCluster(cores, outfile="")
registerDoParallel(cl)
}

for(nn in c(1000)){message(nn)

    data1a <- unlist(c(read.csv(paste0('stensoladist2mom',nn,'.csv'),header=F)))
    data1b <- unlist(c(read.csv(paste0('stensoladist4mom',nn,'.csv'),header=F)))
    n <- 65

    result <- foreach(a=0:n, .combine='cbind') %do% {
    table <- foreach(aa=0:nn, .combine='c', .export=c('condp','nn','n','a','data1a','data1b')) %dox% {condp(n,a,nn,aa)}
    
    c(sum(table*data1a),
      sum(table*data1b))
    }
    write.csv(result[1,],paste0('stensola2mom_marg_from',nn,'to',n,'.csv'),row.names=F)
    write.csv(result[2,],paste0('stensola4mom_marg_from',nn,'to',n,'.csv'),row.names=F)
    

    data2a <- unlist(c(read.csv(paste0('riehledist2mom',nn,'.csv'),header=F)))
    data2b <- unlist(c(read.csv(paste0('riehledist4mom',nn,'.csv'),header=F)))
    n <- 159

    result <- foreach(a=0:n, .combine='cbind') %do% {
    table <- foreach(aa=0:nn, .combine='c', .export=c('condp','nn','n','a','data1a','data1b')) %dox% {condp(n,a,nn,aa)}
    
    c(sum(table*data2a),
      sum(table*data2b))
    }
    write.csv(result[1,],paste0('riehle2mom_marg_from',nn,'to',n,'.csv'),row.names=F)
    write.csv(result[2,],paste0('riehle4mom_marg_from',nn,'to',n,'.csv'),row.names=F)
}
    
if(cores>1){stopCluster(cl)}
