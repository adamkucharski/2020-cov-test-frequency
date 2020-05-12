# - - - - - - - - - - - - - - - - - 
# Analysis functions for test sensitivity vs frequency
# Adam Kucharski (2020)
# https://github.com/adamkucharski/2020-cov-test-frequency
# - - - - - - - - - - - - - - - - - 

library(dplyr)
library(tibble)
library(readr)
library(magrittr)

col_def <-   list(col1="dark grey",col2=rgb(0.9,0.7,0),col3=rgb(0,0,0.8),col4=rgb(0.1,0.4,0.1),col5=rgb(1,0.4,1),col6=rgb(0.2,0,0.8),rgb(0,0.7,0.7))
col_def_F <- list(col1="grey",col2=rgb(0.9,0.7,0,0.5),col3=rgb(0,0,0.8,0.5),col5=rgb(0.1,0.4,0.1,0.5),col6=rgb(1,0.4,1,0.5),col7=rgb(0.2,0,0.8,0.5),rgb(0,0.7,0.7,0.2))


# Define functions ----------------------------------------------

incubation_p <- 5

scale_p <- 2
shift_p <- 0
#period_infectious <- function(x){dgamma(x+shift_p,shape=mean_p/scale_p,scale=scale_p)}

#xx1 <- seq(0,14,0.1)
#plot(xx1,sapply(xx1,period_infectious))

pos_test <- function(x){
  if(x<4){y=min(1,0.25*x)}
  if(x>=4 & x<10){y=1.0}
  if(x>=10){y=max(1-0.2*(x-10),0)}
  y
}


# Calculate probability of detection ----------------------------------------------

test_prop <- function(sensX,freqX){
  
  store_data <- NULL
  
  # Iterate over period
  for(ii in 1:freqX){
    
    testing_seq <- seq(ii,20,freqX)
    detect_prob_each <- sensX*sapply(testing_seq,pos_test)
    detect_prob_total <- 1-prod(1-detect_prob_each)
    
    pre_period <- testing_seq[testing_seq<incubation_p]
    if(length(pre_period)==0){
      detect_presymp_each <- 0
    }else{
      detect_presymp_each <- sensX*sapply(testing_seq[testing_seq<5],pos_test)
    }

    detect_presymp_total <- 1-prod(1-detect_presymp_each)
    
    store_data <- rbind(store_data,c(detect_presymp_total,detect_prob_total)/freqX)
  }
  
  colSums(store_data)

}


# Plot outputs ----------------------------------------------

xx <- seq(0,15,1)

par(mfrow=c(3,1),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=0)

# Baseline function
plot(xx,sapply(xx,pos_test),type="l",col="black",lwd=2,ylab="probability of testing positive",xlab="days post infection")
lines(c(5,5),c(0,10),lty=2)

title(LETTERS[1],adj=0)

range_sensitivity <- seq(0,1,0.01)

# Pre-symptomatic
plot(0,0,type="l",col="red",xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i",ylab="probability of detecting before symptoms",xlab="per-test sensitivity")
grid(ny = NULL, nx=NA, col = "lightgray")
day_list <- c(1,2,4,7,14)

for(ii in 1:length(day_list)){
  
  kk <- day_list[ii]
  yy_1 <- sapply(range_sensitivity,function(x){test_prop(x,kk)})
  
  lines(range_sensitivity,yy_1[1,],col=col_def[[ii]])
  text(labels=paste0("every ",kk," days"),x=0.6,y=0.9*yy_1[1,70],col=col_def[[ii]],pch=0.8,adj=0)

}

title(LETTERS[2],adj=0)

# Any time
plot(0,0,type="l",col="red",xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i",ylab="probability of detecting infection",xlab="per-test sensitivity")
grid(ny = NULL, nx=NA, col = "lightgray")

for(ii in 1:length(day_list)){
  
  kk <- day_list[ii]
  yy_1 <- sapply(range_sensitivity,function(x){test_prop(x,kk)})
  
  lines(range_sensitivity,yy_1[2,],col=col_def[[ii]])
  text(labels=paste0("every ",kk," days"),x=0.2,y=0.9*yy_1[2,30],col=col_def[[ii]],pch=0.8,adj=0)
  
}

title(LETTERS[3],adj=0)

dev.copy(png,paste0("outputs/test_freq.png"),units="cm",width=8,height=20,res=150)
dev.off()





