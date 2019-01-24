## Calculation of expected value, std, thetas, spread of conditional frequencies
## of each symptom given each SNP - extra data set

#### PREAMBLE
## libraries and colour-blind palette from http://www.sron.nl/~pault/
##runfunction <- function(aa=1e3L){
library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('png')
library('plot3D')
library('doParallel')
library('LaplacesDemon')
library('dplyr')
library('GA')
#
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mycolours <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mycolours)
barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
dev.off()
mmtoin <- 0.0393701
#### END PREAMBLE

mehomsolution <- function(popsize,popmean,popcouple){
    popsize2 <- (popsize^2-popsize)/2
    rpopmean <- popmean*popsize
    rpopcouple <- popcouple*popsize2

    zterm <- function(js){
        m <- 0:popsize
        log(sum(exp(lchoose(popsize,m) + js[1]*m + js[2]*(m^2-m)/2))) - sum(js*c(rpopmean,rpopcouple))
    }
    
    zterm <- function(js){
        m <- 0:popsize
        log(sum(choose(popsize,m)*exp( js[1]*m + js[2]*(m^2-m)/2))) - sum(js*c(rpopmean,rpopcouple))
    }
    
    result <- nlm(f=zterm,p=c(0,0),iterlim=1e6)

    result$estimate/c(popsize,popsize2)


    # this seems to lead to solutions closer to the maximum
    mehomsolution2 <- function(popsize,popmean,popcouple){
    popsize2 <- (popsize^2-popsize)/2
    rpopmean <- popmean*popsize
    rpopcouple <- popcouple*popsize2

    zterm <- function(js,cadd){
        m <- 0:popsize
    ##     log(sum(exp(lchoose(popsize,m)# -popsize*log(2)
    ##                 + js[1]*m/popsize + js[2]*(m^2-m)/2/popsize2))) - sum(js*c(popmean,popcouple))
    ## }
            log(sum(exp(lchoose(popsize,m) #
                        + js[1]*m/popsize + js[2]*(m^2-m)/2/popsize2))) - sum(js*c(popmean,popcouple))
        }

        zterm <- function(js){
        m <- 0:popsize
    ##     log(sum(exp(lchoose(popsize,m)# -popsize*log(2)
    ##                 + js[1]*m/popsize + js[2]*(m^2-m)/2/popsize2))) - sum(js*c(popmean,popcouple))
    ## }
            log(sum(exp(lchoose(popsize,m) #
                        + js[1]*(m/popsize-popmean) + js[2]*((m^2-m)/2/popsize2-popcouple)))) # - sum(js*c(popmean,popcouple))
        }

        result2b <- optim(fn=zterm,par=c(0,0))

    result2 <- nlm(f=zterm,p=c(-popsize,popsize),iterlim=1e6)

    meprob <- function(js,popsize){
        m <- 0:popsize
        up <-  exp(lchoose(popsize,m) -popsize*log(2) + js[1]*m/popsize + js[2]*(m^2-m)/2/popsize2)
        up/sum(up)
    }

    mecheck <- function(js,popsize){
        m <- 0:popsize
        popsize2 <- (popsize^2-popsize)/2
        up <-  exp(lchoose(popsize,m) -popsize*log(2) + js[1]*m/popsize + js[2]*(m^2-m)/2/popsize2)
        prob <- up/sum(up)
        c(sum(m*prob)/popsize, sum(((m^2-m)/2)*prob)/popsize2)
    }
    
    result$estimate
