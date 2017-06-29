## Check that a new enough R version is being used. R>= 3.2
if(version$major != "3" | as.numeric(version$minor) < 3.3)
  stop('R version needs to be > 3.3.3')

## Install the package and any dependencies from github
install.packages('devtools')
library(devtools)
install_github('r4ss/r4ss', dependencies = TRUE)
install_github('ss3sim/ss3sim', dependencies = TRUE)

## Make sure you can load all of these, if not the install them manually
## with install.packages()
library(r4ss)
library(ss3sim)
library(doParallel)
library(foreach)
library(ggplot2)

## Set your working directory to be in the same folder as this script
getwd()
setwd(dir= )

## Run a simple model to see if the package is installed and working
## correctly. Don't worry about what this means, we'll cover that in the
## workshop.
run_ss3sim(iterations=1, scenarios='D0-F1-cod', case_folder='cases',
           case_files=list(F='F', D=c('index','lcomp','agecomp')),
           om_dir='model/om/', em_dir='model/em/')

## Checks to see if it worked
if(!file.exists('D0-F1-cod/1/om/Report.sso')){
    stop("OM failed to run, contact Kelli")
} else if(!file.exists('D0-F1-cod/1/em/Report.sso')){
    stop("EM failed to run, contact Kelli")
} else {
    message("OM and EM ran successfully: package installed and working!")
    unlink('D0-F1-cod', recursive = TRUE)
}
