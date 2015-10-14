## Install the package and any dependencies from github
install.packages('devtools')
library(devtools)
install_github('ss3sim/ss3sim', dependencies = TRUE)
install_github('r4ss/r4ss')

## Make sure you can load all of these, if not the install them manually
## with install.packages()
library(ss3sim)
library(r4ss)
library(doParallel)
library(foreach)
library(ggplot2)

## Set your working directory to be in the same folder as this script
getwd()
setwd(dir=  )

## Run a simple model to see if the package is installed and working
## correctly. Don't worry about what this means, we'll cover that in the
## workshop.
run_ss3sim(iterations=1, scenarios='D0-F1-cod', case_folder='cases',
           case_files=list(F='F', D=c('index','lcomp','agecomp')),
           om_dir='model/om/', em_dir='model/em/')

## Checks to see if it worked
if(!file.exists('D0-F1-cod/1/om/Report.sso')){
    stop("OM failed to run, contact Cole or Kelli")
} else if(!file.exists('D0-F1-cod/1/em/Report.sso')){
    stop("EM failed to run, contact Cole or Kelli")
} else {
    message("OM and EM ran successfully: package installed and working!")
    unlink('D0-F1-cod', TRUE)
}
