### This file demonstrates the capabilities and workflow for the ss3sim
### stock assessment simulation framework.
### Contact the authors for any questions or
### issues. For more on the package see:
###
### Anderson, S.C., C.C. Monnahan, K.F. Johnson, K. Ono, and
### J.L. Valero. 2014. ss3sim: An R Package for Fisheries Stock Assessment
### Simulation with Stock Synthesis. Plos One 9:e92725.
### ------------------------------------------------------------

### ------------------------------------------------------------
### [Step 0]
### Startup the working environment

## Set your working directory to be in the same folder as this script
getwd()
setwd(dir= )

## Load the necessary libraries
library(r4ss)
library(ss3sim)
library(foreach)
library(doParallel, quietly = TRUE, verbose = FALSE)
library(ggplot2)

## Setup parallel cores and iterations to run
Sys.getenv("NUMBER_OF_PROCESSORS")
cores <- 2                              # cores for parallel
registerDoParallel(cores)

## Some checks before proceeding
<<<<<<< HEAD
packageVersion("ss3sim")                # should be 0.9.5.9000 (GitHub) or 0.9.6 (CRAN)
if(as.numeric(substr(packageVersion("r4ss"), 1, 4)) < 1.27) stop("r4ss should",
  " be version 1.27.0 or greater")
Sys.info()[5]                           # machine needs to be 64bit
## Check that a new enough R version is being used. R>= 3.3
if(version$major != "3" | as.numeric(version$minor) < 3)
  stop('R version needs to be > 3.3')

### End of [Step 0]
### ------------------------------------------------------------

### ------------------------------------------------------------
### [Step 1]
### Run a simple example of a scenario to demonstrate input syntax,
### output structure, and results.
### This is our 'cod' model with data-rich sampling history and
### a two-way trip fishing history.
### A full simulation for research purposes would have multiple scenarios.

## Tell ss3sim the specifications for the
## * operating model (OM),
## * sampling,
## * estimation model (EM).
## Case files (.txt files stored in a directory, here labeled 'cases')
## tell ss3sim what to do.
## Case files are linked with internal ss3sim functions.
## F == effort trajectory
## D == sampling protocol of index of abundance (index),
## length-composition data (lcomp), and age-composition data (agecomp).
## <Look at input folder structure.>
run_ss3sim(
  iterations=1,                # vector of iterations to run
  scenarios='D0-F1-cod',       # vector of scenarios
  case_folder='cases',         # folder containing case files
  case_files=list(             # linking case files to cases
    F='F',
    D=c('index','lcomp','agecomp')),
  om_dir='model/om/',          # location of the OM ss3 files
  em_dir='model/em/')          # location of the EM ss3 files
## Should see "Completed iterations: 1 for scenarios: D0-F1-cod" on console.

## <Look at the folder structure that was created.>
## We can also use r4ss to ## make standard plots.
## Plot the OM
r4ss::SS_plots(
  replist=r4ss::SS_output(
    'D0-F1-cod/1/om',          # Specify the OM folder
    forecast=FALSE, covar=FALSE,
    NoCompOK=TRUE, ncols=250, printstats=FALSE, verbose=FALSE),
  png=TRUE, uncertainty=FALSE, html=FALSE, verbose=FALSE)
## Plot the EM
r4ss::SS_plots(
  replist=r4ss::SS_output(
    'D0-F1-cod/1/em',          # Specify the EM folder
    forecast=FALSE, covar=FALSE,
    NoCompOK=TRUE, ncols=250, printstats=FALSE, verbose=FALSE),
  png=TRUE, uncertainty=FALSE, html=FALSE, verbose=FALSE)
dev.off()
## <Navigate to "D0-F1-cod/1/om/plots/catch9 harvest rate">
## <Navigate to "D0-F1-cod/1/em/plots/catch9 harvest rate">

## Now rerun with more replicates, in parallel (within scenario) this
## time. We recommend using serial execution during development since error
## messages in parallel can be more difficult to interpret.
## Specify arguments that will be recycled.
case_files <- list(F='F', D=c('index','lcomp','agecomp'))
om <- 'model/om/'; em <- 'model/em/'

run_ss3sim(
  iterations=2:4,              # Run iteration 2, 3, and 4 in parallel
  scenarios='D0-F1-cod',       # Same scenario as last time
  case_folder='cases',
  case_files=case_files,
  om_dir=om,
  em_dir=em,
  parallel=TRUE,               # Tells ss3sim to run things in parallel
  parallel_iterations=TRUE)    # Specifies iterations or scenarios
## Should see "Running iterations in parallel ..."

## Read in results of the runs.
## Use 'get_results_all', which uses r4ss to read in Report.sso files and
## write two .csv files for each scenario (stored in the scenario folders)
## and two .csv files combining the results from each scenario.
## Here, the files will be identical because we only have one scenario.
get_results_all(user='D0-F1-cod')
## <look at csv files>
results_sc_1 <- read.csv('ss3sim_scalar.csv') # Scalar quantities, e.g., M
results_ts_1 <- read.csv('ss3sim_ts.csv')     # Time-series quantities, e.g., F

## If ss3sim is not working on your computer or you want to look at the
## results without running the models, open the stored csv files in
## the 'results' folder.
## results_sc_1 <- read.csv(file.path('results', 'results_sc_1.csv'))
## results_ts_1 <- read.csv(file.path('results', 'results_ts_1.csv'))

## Add columns of relative error
results_sc_1 <- calculate_re(results_sc_1, add=TRUE)
results_ts_1 <- calculate_re(results_ts_1, add=TRUE)
## <see the additional columns>
head(results_ts_1)

## ss3sim has a set of functions for plotting scalar and time-series output:
## e.g., "plot_ts_lines", which will plot lines for each iteration of any
## column you wish in "results_ts_1" over time.
args(plot_ts_lines)
?plot_ts_lines
plot_ts_lines(results_ts_1, y='SpawnBio_om', vert='D')
plot_ts_lines(results_ts_1, y='SpawnBio_re', vert='D', re=TRUE)
## We can also look at scalar values, which will plot points for each iteration
plot_scalar_points(results_sc_1, x='D', y='depletion_re', re=TRUE)
plot_scalar_boxplot(results_sc_1, x='D', y='SSB_MSY_re', re=TRUE)

## If you had multiple scenarios, you could look at the results by scenario
plot_scalar_boxplot(results_sc_1, x='D', y='SSB_MSY_re', re=TRUE,
  horiz = "scenario")          # Boxplot for each scenario stacked vertically
plot_scalar_boxplot(results_sc_1, x='D', y='SSB_MSY_re', re=TRUE,
  horiz = "species")          # life histories stacked vertically

## Print a list of scalar values, for which relative error was calculated.
names(results_sc_1)[grep('_re', x=names(results_sc_1))]
### End of [Step 1]
### ------------------------------------------------------------

### ------------------------------------------------------------
### [Step 2].
## Create a new data case.
## Create three new text files manually,
## run the new scenario, and
## compare performance with the previous scenario from Step 1, "D0-F1-cod".
##
## The data case is referenced with letter D. Associated with this case are
## three files:
case_files$D
## Each case file contains arguments that are passed to sampling
## functions within ss3sim:
## an index of abundance (or CPUE),
## length compositions (lcomp), and
## age compositions (agecomp).

## ss3sim reads these text files in and parses the text to R.
## All of the arguments supplied in the text files are stored in a list
## within R and ss3sim passes the arguments to functions during the simulation
## at the appropriate time (e.g., F arguments before the OM runs).
## Internally, the parsing takes place with functions like the following:
caseargs <- get_caseargs(folder='cases', scenario='D0-F1-cod', case_files=case_files)
caseargs$F
## Plot the fishing trajectory passed to the OM.
plot(fvals~years, data=caseargs$F,
  type='l', ylab=expression(italic('F')), las = 1)
## Here are the agecomp case arguments:
(ageargs <- caseargs$agecomp)
## The arguments supplied within the text files specific to the D case
## lead to sampling of data generated by the OM.
## ss3sim passes these arguments after the OM is completed, but before the
## EM is ran. There are three case files associated with the arguments
## specific to sampling (agecomp0-cod.txt, lcomp0-cod.txt, and index0-cod.txt).
## For a given scenario, all D-case text files need to have the same number,
## whereas other cases contain a single text file.
##
## For example,
## (1) ss3sim extracts the expected data from the OM
## (2) ss3sim uses r4ss to read this data (e.g., 'data.ss_new') into memory
## (3) ss3sim samples from the data leading to data with observation error
## and writes it to the disk
ss3sim:::extract_expected_data('om_data.ss_new', data_out='expected.dat')
dat <- r4ss::SS_readdat(file='expected.dat', verbose=FALSE)
sample_agecomp(dat_list=dat, outfile='sampled.dat',
               fleets=ageargs$fleets,
               Nsamp=ageargs$Nsamp,
               years=ageargs$years,
               cpar=ageargs$cpar)
## <Look at differences between in text editor>.
r4ss::SS_writedat(dat, outfile = "expected_formatted.dat", verbose = FALSE)

## =====  EXERCISE: CREATE A NEW CASE 'D1' THAT HAS LESS DATA =====.

## [1] save a copy of agecomp0-cod.txt as agecomp1-cod.txt. Then reduce
## sample sizes from 125 for fleet 1 (the fishery) to 25 and from 500 to
## 100 for fleet 2 (the survey).
## [2] Repeat with the length comps.
## [3] Lastly, we need index1-cod too. Create this and and double the value
## for sds_obs to 0.4 (sd in log space of the index).

## Now make sure these can be parsed correctly:
caseargs <- get_caseargs(folder='cases', scenario='D1-F1-cod', case_files=case_files)
caseargs$agecomp
caseargs$lcomp
caseargs$index

## Now run the new scenario, "D1-F1-cod", where
## D0 was replaced with D1, so ss3sim uses the files you created.
run_ss3sim(
  iterations=1:cores,          # run iterations in parallel
  scenarios='D1-F1-cod',       # run the new scenario
  case_folder='cases', case_files=case_files, om_dir=om, em_dir=em,
  parallel=TRUE, parallel_iterations=TRUE)
## Compare D0 and D1 by reading in the ss3 files using
## ss3sim::get_results_all, which will create additional csv files for
## the new scenario and combine them with the old csv file for the original
## scenario.
get_results_all(user=c('D0-F1-cod', 'D1-F1-cod'))
results_sc_2 <- calculate_re(read.csv('ss3sim_scalar.csv'), add=TRUE)
results_ts_2 <- calculate_re(read.csv('ss3sim_ts.csv'), add=TRUE)
## In case the simulation did not work for you uncomment the following
## two lines to read in the stored result files.
## results_sc_2 <- read.csv('results/results_sc_2.csv')
## results_ts_2 <- read.csv('results/results_ts_2.csv')
plot_scalar_boxplot(results_sc_2, x='D', y='depletion_re', re=TRUE)
plot_scalar_boxplot(results_sc_2, x='D', y='SSB_MSY_re', re=TRUE)
plot_scalar_boxplot(results_sc_2, x='D', y='SizeSel_1P_1_Fishery_re', re=TRUE)
plot_ts_lines(results_ts_2, y='SpawnBio_re', vert='D', re=TRUE)
## Note: Each replicate number always has the same recdevs so the scenarios
## are more comparable:
plot_ts_lines(subset(results_ts_2, replicate==1), y='SpawnBio_om', vert='D')
### End of [Step 2]
### ------------------------------------------------------------

### ------------------------------------------------------------
### [Step 3]. In our last example we will explore the interaction of process
### error and data weighting.

## In addition to the sampling functions, there are a variety of other
## functions that manipulate the models during runtime. For instance
## 'change_e' is used to change starting parameters and phases for
## parameters in the EM.
?change_e
## We will use it to turn on and off estimation of M and steepness. Growth
## parameters will also be turned off (and fixed at their true values).
E.cases <- c(0,1,2)
E.df <- data.frame(
    E=paste0('E', E.cases),
    estimated=factor(c("M & h fixed", "h estimated", "M estimated"),
    levels=c("M & h fixed", "h estimated", "M estimated")))
E.df
## <look at E0, E1 and E2>

## 'change_tv' adds time-varying (tv) changes to OM parameters using
## environmental deviations. This function is very useful for adding
## process error in the OM.
?change_tv
## We will use it to add a trend to a selex parameter
selex.scalar.vec <- c(0,20)
S.cases <- seq_along(selex.scalar.vec)
S.df <- data.frame(selex=factor(c("No Process Error", "Process Error")),
                   S=paste0('S', S.cases))
S.df
## <look at S0 and S1>

## The sampling functions have an optional argument for the effective
## sample size (ESS). Here we specify it incorrectly to explore different
## data weightings on the age and length compositions.
ESS.scalar.vec <- sort(unique(c(1,exp(seq(log(.1),log(10), len=10)))))
D.cases <- seq_along(ESS.scalar.vec)
D.df <- data.frame(ess.ratio=ESS.scalar.vec, D=paste0('D', D.cases))
D.df
## <look at lcomp11-cod.txt>

## The case files are in the 'extra_cases' folder and ready to go. Notice
## that we've added new letters 'S' and 'E' to our case files. The package will now
## look for files named "S" to match up to our S cases (S0 and S1).
case_folder <- 'cases/extra_cases'
case_files <- list(F="F", D=c("index","lcomp","agecomp"), S='S', E='E')
scenarios <- expand_scenarios(cases=list(D=D.cases, F=1, S=S.cases, E=E.cases),
                     species='cod')
scenarios
## run_ss3sim(iterations=1:50, scenarios=scenarios, parallel=TRUE,
##            parallel_iterations=TRUE, case_folder=case_folder, om_dir=om,
##            em_dir=em, case_files=case_files)

## This simulation takes too long to run ,but here are the results,
## processed a bit
results_sc_3 <- readRDS('results/results_sc_3.Rdata')
results_ts_3 <- readRDS('results/results_ts_3.Rdata')

## Look at some performance measures
myylim <- ylim(-1,1)
g <- ggplot(results_sc_3, aes(ess.ratio, SSB_MSY_re, group=replicate))+
    geom_line(alpha=.5) + facet_grid(estimated~selex.scalar) +
    geom_vline(xintercept=1, col='blue') +
    geom_hline(yintercept=0, col='red') + myylim+ theme_bw()
g
ggsave('plots/SSB_MSY_re.png', g, width=9, height=6)
g <- ggplot(results_sc_3, aes(ess.ratio, depletion_re, group=replicate))+
    geom_line(alpha=.5) + facet_grid(estimated~selex.scalar) +
    geom_vline(xintercept=1, col='blue') +
    geom_hline(yintercept=0, col='red') + myylim+ theme_bw()
g
ggsave('plots/depletion_re.png', g, width=9, height=6)
### [End of step 3]
### ------------------------------------------------------------

### End of File
