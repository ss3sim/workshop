### This file demonstrates the capabilities and workflow for the ss3sim
### stock assessment simulation framework. There is an accompanying set of
### files necessary to run. Contact the authors for any questions or
### issues. For more on the package see:
###
### Anderson, S.C., C.C. Monnahan, K.F. Johnson, K. Ono, and
### J.L. Valero. 2014. ss3sim: An R Package for Fisheries Stock Assessment
### Simulation with Stock Synthesis. Plos One 9:e92725.
###
### File developed for the CAPAM workshop on data weighting, likelihoods,
### and process error. Cole Monnahan | monnahc@uw.edu.
### ------------------------------------------------------------

### ------------------------------------------------------------
### [Step 0] Startup the working environment

## Set your working directory to be in the same folder as this script
getwd()
setwd(dir= )

## Load the neccessary libraries
library(ss3sim)
library(r4ss)
library(doParallel)
library(foreach)
library(ggplot2)

## Setup parallel cores and iterations to run
Sys.getenv("NUMBER_OF_PROCESSORS")
cores <- 2                              # cores for parallel
registerDoParallel(cores)

## Some checks before proceeding
packageVersion("ss3sim")                # should be 0.8.9.9000
packageVersion("r4ss")                  # should be 1.23.5
Sys.info()[5]                           # machine needs to be 64bit
## Check that a new enough R version is being used. R>= 3.2
if(version$major != "3" | as.numeric(version$minor) < 2)
  stop('R version needs to be > 3.2')

### End of [Step 0]
### ------------------------------------------------------------

### ------------------------------------------------------------
### [Step 1] Run a very simple example to demonstrate syntax and output
### structure and results. This is our 'cod' model and a data rich case.
##
## Tell package how to link case files with internal functions (more on
## this later). F refers to the effort trajectory and D refers to the data
## used.
case_files <- list(F='F', D=c('index','lcomp','agecomp'))
## Path to the folders containing the models. <Stop and look at these>
om <- 'model/om/'; em <- 'model/em/'
run_ss3sim(iterations=1,                # vector of iterations to run
           scenarios='D0-F1-cod',       # vector of scenarios
           case_folder='cases',         # folder containing case files
           case_files=case_files,       # linking case files to cases
           om_dir=om, em_dir=em)        # where to find models to use

## <Look at the folder structure that was created.>
## We can also use r4ss to ## make standard plots.
out <- SS_output('D0-F1-cod/1/em', forecast=FALSE, covar=FALSE,
                 NoCompOK=TRUE, ncols=250, printstats=FALSE, verbose=FALSE)
SS_plots(replist=out, png=TRUE, uncertainty=FALSE, html=FALSE, verbose=FALSE)

## Now rerun with more replicates, in parallel (within scenario) this
## time. We recommend using serial execution during development since error
## messages in parallel can be more difficult to interpret.
run_ss3sim(iterations=2:(cores+1), scenarios='D0-F1-cod', case_folder='cases',
           case_files=case_files, om_dir=om, em_dir=em,
           parallel=TRUE, parallel_iterations=TRUE)

## Read in results of the runs. This function uses r4ss to read in and then
## write a csv for each scenario.
get_results_all(user='D0-F1-cod') ## <look at csv files>
results_sc_1 <- read.csv('ss3sim_scalar.csv')
results_ts_1 <- read.csv('ss3sim_ts.csv')
## Write and read in results (in case of failure or too slow)
## write.csv(results_sc_1, file='results/results_sc_1.csv')
## write.csv(results_ts_1, file='results/results_ts_1.csv')
## results_sc_1 <- read.csv('results/results_sc_1.csv')
## results_ts_1 <- read.csv('results/results_ts_1.csv')

## Add columns of relative error
results_sc_1 <- calculate_re(results_sc_1, add=TRUE)
results_ts_1 <- calculate_re(results_ts_1, add=TRUE)

## ss3sim has a set of functions for plotting scalar and ts output:
?plot_ts_lines
plot_ts_lines(results_ts_1, y='SpawnBio_om', vert='D')
plot_ts_lines(results_ts_1, y='SpawnBio_re', vert='D', re=TRUE)
## We can also look at scalar values instead of time series (year on x-axis)
plot_scalar_points(results_sc_1, x='D', y='depletion_re', re=TRUE)
plot_scalar_boxplot(results_sc_1, x='D', y='SSB_MSY_re', re=TRUE)
names(results_sc_1)[grep('_re', x=names(results_sc_1))]
### End of [Step 1]
### ------------------------------------------------------------

### ------------------------------------------------------------
### [Step 2]. In our second example we will add a new data case manually,
### run it, and compare peformance with the previous example.
##
## The data case is referenced with letter D. Associated with this case are
## three files:
case_files$D
## Each case file contains three arguments that are passed to sampling
## functions: an index of abundance (or CPUE), length comps (lcomp) or age
## comps (agecomp)

## ss3sim reads these files in and parses the text to R and then stores
## these in a list to be passed to the sampling function. We can see this
## parsing with this function
caseargs <- get_caseargs(folder='cases', scenario='D0-F1-cod', case_files=case_files)
caseargs$F
plot(fvals~years, data=caseargs$F, type='l', ylab='F')
## Here are the agecomp case arguments:
(ageargs <- caseargs$agecomp)
## Let's look at what happens with these arguments. First read in a
## 'data.ss_new' file output from an arbitrary run of the OM with this model.
ss3sim:::extract_expected_data('om_data.ss_new', data_out='expected.dat')
dat <- r4ss::SS_readdat(file='expected.dat', verbose=FALSE)
sample_agecomp(dat_list=dat, outfile='sampled.dat',
               fleets=ageargs$fleets,
               Nsamp=ageargs$Nsamp,
               years=ageargs$years,
               cpar=ageargs$cpar)
## <Look at differences in text editor>.
## These arguments are passed to this function **during** the simulation at
## the appropriate step (in this case after the OM runs but before the
## EM). In general this is how ss3sim works. We currently have three case
## files associated with 'D0'.

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

## Now we run the new scenario, replacing D0 with D1 so it looks for our
## new files and uses those arguments in the sampling functions.
run_ss3sim(iterations=1:Nsim, scenarios='D1-F1-cod', case_folder='cases',
           case_files=case_files, om_dir=om, em_dir=em,
           parallel=TRUE, parallel_iterations=TRUE)
## We want to compare D0 and D1, so read them both in to combine
## together. The old files for the D0 case still exist so it skips those
## but adds them to the final output files.
get_results_all(user=c('D0-F1-cod', 'D1-F1-cod'))
results_sc_2 <- calculate_re(read.csv('ss3sim_scalar.csv'), add=TRUE)
results_ts_2 <- calculate_re(read.csv('ss3sim_ts.csv'), add=TRUE)
## write.csv(results_sc_2, file='results/results_sc_2.csv')
## write.csv(results_ts_2, file='results/results_ts_2.csv')
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

## 'change_tv' addes time-varying (tv) changes to OM parameters using
## environmental deviations. This function is very useful for adding
## process error in the OM.
?change_tv
## We will use it to add a trend to a selex parameter
selex.scalar.vec <- c(0,20)
S.cases <- seq_along(selex.scalar.vec)
S.df <- data.frame(selex.scalar=factor(paste0('selex.scalar=',selex.scalar.vec),
                   levels=paste0('selex.scalar=',selex.scalar.vec)),
                   S=paste0('S', S.cases))
S.df <- data.frame(selex=factor(c("No Process Error", "Process Error")),
                   S=paste0('S', S.cases))
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
## that we've added a new letter S to our case files. The package will now
## look for files named "S" to match up to our S cases (S0 and S1).
case_folder <- 'cases/extra_cases'
case_files <- list(F="F", D=c("index","lcomp","agecomp"), S='S', E='E')
scenarios <- expand_scenarios(cases=list(D=D.cases, F=1, S=S.cases, E=E.cases),
                     species='cod')
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
ggsave('plots/SSB_MSY_re.png', g, width=9, height=6)
g <- ggplot(results_sc_3, aes(ess.ratio, depletion_re, group=replicate))+
    geom_line(alpha=.5) + facet_grid(estimated~selex.scalar) +
    geom_vline(xintercept=1, col='blue') +
    geom_hline(yintercept=0, col='red') + myylim+ theme_bw()
ggsave('plots/depletion_re.png', g, width=9, height=6)
### [End of step 3]
### ------------------------------------------------------------

### ------------------------------------------------------------
### [Some extra code for running step 3 above, DO NOT RUN]
###
## ## Read in results and process for plotting
## ## Look at trend in process error in selex
## S1.devs <- as.numeric(unlist(get_caseargs(case_folder, 'D1-E0-F1-S1-cod', case_files=case_files)$tv_params))
## S2.devs <- as.numeric(unlist(get_caseargs(case_folder, 'D1-E0-F1-S2-cod', case_files=case_files)$tv_params))
## devs.df <- data.frame(cbind(years=1:100, S1=S1.devs, S2=S2.devs))
## devs.df.long <- reshape2::melt(devs.df, 'years', variable.name='S')
## devs.df.long <- merge(devs.df.long, S.df, by='S')
## g <- ggplot(devs.df.long, aes(years, value+50.8, group=selex.scalar, color=selex.scalar))+
##     geom_line() + ylab('Sel_P1') + theme_bw()
## ggsave('plots/selex_patterns.png', g, width=6, height=3)
## ##
## ## Look at differences in OM biomass trajectories for an arbitrary iteration
## yy1 <- ddply(subset(results_ts_3, E=='E0' & D=='D1'), .(year, replicate), mutate,
##              SSB_re=(SpawnBio_om-SpawnBio_om[which(S=='S1')])/
##                  SpawnBio_om[which(S=='S1')])
## ggplot(subset(yy1,  replicate==1), aes(year, SpawnBio_om, group=selex.scalar, color=selex.scalar))+
##     geom_line() + ylim(0,5e9) + theme_bw()
## ## read in and process the results
## get_results_all(user=scenarios, parallel=TRUE)
## results_ts_3 <- read.csv("ss3sim_ts.csv"); results_ts_3 <- merge(results_ts_3, D.df, by="D")
## results_ts_3 <- merge(results_ts_3, S.df, by="S"); results_ts_3 <- merge(results_ts_3, E.df, by="E")
## ## create meta dataframe to merge into results
## results_sc_3 <- calculate_re(read.csv("ss3sim_scalar.csv"))
## results_sc_3 <- merge(results_sc_3, results_ts_3[results_ts_3$year==100, c('ID', 'year', 'SpawnBio_om',
##                    'SpawnBio_em')], by='ID')
## results_sc_3 <- within(results_sc_3, {
##              depletion_om = SpawnBio_om/SSB_Unfished_om
##              depletion_em = SpawnBio_em/SSB_Unfished_em
##          depletion_re=(depletion_om-depletion_em)/depletion_om})
## results_sc_3 <- merge(results_sc_3, D.df, by="D"); results_sc_3 <- merge(results_sc_3, S.df, by="S")
## results_sc_3 <- merge(results_sc_3, E.df, by="E")
## write.csv(results_sc_3, file='results/results_sc_3.csv')
## write.csv(results_ts_3, file='results/results_ts_3.csv')
## saveRDS(results_sc_3, file='results/results_sc_3.RData')
## saveRDS(results_ts_3, file='results/results_ts_3.RData')
## temp <- results_sc_3
## temp$depletion_re[results_sc_3$ess.ratio !=1] <- NA
## levels(temp$selex.scalar) <- c("No Process Error", "Process Error")
## g <- plot_scalar_points(temp, x='ess.ratio', y='depletion_re', rel=TRUE, vert='selex.scalar', horiz='estimated', print=FALSE)
## g <- g + geom_vline(xintercept=1, col='blue', alpha=.5, lwd=.5) + theme_bw()
## ggsave('plots/depletion_re_empty.png', g, width=6, height=4)
## rm(temp)
### ------------------------------------------------------------

### End of File
