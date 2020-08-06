# 2020-cov-test-frequency

Code for analysis of screening tests.

### Quick start guide for model

First, set local path in R to GitHub directory, e.g.:
`
setwd("~/Documents/GitHub/2020-cov-test-frequency/")
`
Main model run script is in `R/testing_model.r`.

### Schematic of testing and key delays

<img src="https://user-images.githubusercontent.com/8329046/89546705-cc726180-d7fc-11ea-8e60-6c50be29dad5.png" width="600px" />

*Timeline of transmission. Three key delays drive outbreak dynamics: time from onset of infectiousness to onset of symptoms (delay A); onset of symptoms to isolation of the index case (delay B); and time from onset of infectiousness in a contact to quarantine of that contact (delay C). Isolation of symptomatic cases & contact tracing aims to reduce delays B & C, while screening tests can reduce delay A & B.*