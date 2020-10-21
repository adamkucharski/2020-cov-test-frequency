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

# Define functions ----------------------------------------------

# Mean incubation period
incubation_p <- 5

# Probability of testing positive
pos_test <- function(x){
  if(x<4){y=min(0.8,0.25*x)}
  if(x>=4 & x<8){y=0.8}
  if(x>=8){y=max(0.8-0.2*(x-8),0)}
  y
}

# Incubation period  distribution
inc_cumulative <- function(x){pgamma(x,shape = incubation_p, scale = 1) }
inc_probability <- function(x){dgamma(x,shape = incubation_p, scale = 1) }


# Calculate probability of detection ----------------------------------------------

detection_prob <- function(freqX, # frequency of testing
                           detect_within = 7, # period aiming to detect within
                           delay_to_result = 2, # delay to result
                           prop_symp = 0.6 # proportional symptomatic (set to 0 for asymptomatic testing only)
                           ){
  
  
  #DEBUG:  detect_within = 7; delay_to_result = 2; freqX = 7; prop_symp = 0.6
  
  # Define effective period aiming for detection (accounting for delay to result)
  e_detect_within <- detect_within - delay_to_result
  
  # Define possible days at which could be detected
  seq_detect <- (1:e_detect_within)
  
  # - - -
  # Calculate probabilty detect symptomatic individuals - iterate over each possible day of onset
  store_data_symp <- NULL
  
  for(ii in 1:e_detect_within){
    
    # Full equation for onset on given day:
    # (1-[1-P(onset on day X)*P(test positive | onset)]*[1-P(detected by screening by day X)] )
  
    # P(onset by day X)*P(test positive | onset)
    prob1 <- inc_probability(ii)*pos_test(ii)
    
    # P(detected by screening by day X)
    prob2 <- 1-prod(1-sapply(1:ii,pos_test)/freqX)
    
    # Collect together:
    detect_if_symp <- (1-(1-prob1)*(1-prob2))
      
    # Store values
    store_data_symp <- c(store_data_symp,detect_if_symp)
  
  }
  
  detect_symp_total <- 1-prod(1-store_data_symp) # Estimate overall probability detect
  
  
  # - - -
  # Calculate probabilty detect asymptomatic individuals

  # P(detected by screening by day X)
  detect_apresymp_total <- 1-prod(1-sapply(1:e_detect_within,pos_test)/freqX)

  # - - -
  # Collate symptomatic and asymptomatic:
  probability_detect <- prop_symp*detect_symp_total + (1-prop_symp)*detect_apresymp_total
  
  probability_detect
  
}




# Calculate probability of detection ----------------------------------------------

test_prop <- function(sensX,freqX){
  
  store_data <- NULL
  
  # Iterate over period
  for(ii in 1:freqX){
    
    # Define testing intervals
    testing_seq <- seq(ii,20,freqX)
    
    detect_prob_each <- sensX*sapply(testing_seq,pos_test) # Estimate probability detect at each point
    detect_prob_total <- 1-prod(1-detect_prob_each) # Estimate overall probability detect
    
    pre_period <- testing_seq[testing_seq<incubation_p] # Pick tests during pre-symptomatic period
    if(length(pre_period)==0){
      detect_presymp_each <- 0 # If no tests during this period
    }else{
      detect_presymp_each <- sensX*sapply(testing_seq[testing_seq<5],pos_test) # If tests during period
    }

    detect_presymp_total <- 1-prod(1-detect_presymp_each) # Estimate overall probability detect during pre-symptomatic
    
    # Normalise by probability interval starts at given point (i.e. frequency) and store data
    store_data <- rbind(store_data,c(detect_presymp_total,detect_prob_total)/freqX) 
  }
  
  colSums(store_data)

}


# Plot outputs ----------------------------------------------


par(mfrow=c(3,1),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=0)
xx <- seq(0,15,1)

# Baseline function
plot(xx,sapply(xx,pos_test),type="l",col="black",lwd=2,ylab="probability of testing positive",xlab="days post infection")
lines(c(5,5),c(0,10),lty=2)
title(LETTERS[1],adj=0)

# Detect during pre-symptomatic period
plot(0,0,type="l",col="red",xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i",ylab="probability of detecting before symptoms",xlab="per-test sensitivity")
grid(ny = NULL, nx=NA, col = "lightgray")
day_list <- c(1,2,4,7,14)

range_sensitivity <- seq(0,1,0.01) # Pick range to iterate sensitivity over

for(ii in 1:length(day_list)){
  
  kk <- day_list[ii]
  yy_1 <- sapply(range_sensitivity,function(x){test_prop(x,kk)})
  
  lines(range_sensitivity,yy_1[1,],col=col_def[[ii]])
  text(labels=paste0("every ",kk," days"),x=0.6,y=0.9*yy_1[1,70],col=col_def[[ii]],pch=0.8,adj=0)

}

title(LETTERS[2],adj=0)

# Detect at any point in time
plot(0,0,type="l",col="red",xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i",ylab="probability of detecting infection",xlab="per-test sensitivity")
grid(ny = NULL, nx=NA, col = "lightgray")

for(ii in 1:length(day_list)){
  
  kk <- day_list[ii]
  yy_1 <- sapply(range_sensitivity,function(x){test_prop(x,kk)})
  
  lines(range_sensitivity,yy_1[2,],col=col_def[[ii]])
  text(labels=paste0("every ",kk," days"),x=0.2,y=0.9*yy_1[2,30],col=col_def[[ii]],pch=0.8,adj=0)
  
}

title(LETTERS[3],adj=0)

# Output figure
dev.copy(png,paste0("test_freq.png"),units="cm",width=8,height=20,res=150)
dev.off()





