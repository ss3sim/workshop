## Install the package and any dependencies from github
if(!require(devtools)) install.packages('devtools')
devtools::install_github('ss3sim/ss3sim', dependencies = TRUE)
devtools::isntall_github('r4ss/r4ss')

## Make sure you can load all of these, if not the install them manually
## with install.packages()
library(ss3sim)
library(r4ss)
library(doParallel)
library(foreach)
library(ggplot2)

case_files <- list(F='F', D=c('index','lcomp','agecomp'))
## Path to the folders containing the models. <Stop and look at these>
om <- 'model/om/'; em <- 'model/em/'
run_ss3sim(iterations=1,                # vector of iterations to run
           scenarios='D0-F1-cod',       # vector of scenarios
           case_folder='cases',         # folder containing case files
           case_files=case_files,       # linking case files to cases
           om_dir=om, em_dir=em)        # where to find models to use

if(!file.exists('D0-F1-cod/1/om/Report.sso')){
    stop("OM failed to run, contact Cole or Kelli")
} else if(!file.exists('D0-F1-cod/1/em/Report.sso')){
    stop("EM failed to run, contact Cole or Kelli")
} else {
    message("OM and EM ran successfully: package installed and working!")
    unlink('D0-F1-cod', TRUE)
}
