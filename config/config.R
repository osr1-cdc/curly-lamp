# Update the following each run ----------------------------------------
svy.type <- "svyNEW" #variable specifying which survey design to use
ci.type <- "KG" #variable specifying which CI estimator to use for weighted proportion esimtates

#variable for whether or not to include the state tagged sequences
# Update the following each run ----------------------------------------
time_end <- as.Date("2021-10-23") #set end date for national and regional survey estimates  
state_time_end=c(as.Date("2021-09-18"),as.Date("2021-09-25"), as.Date("2021-10-02"), as.Date("2021-10-09"), as.Date("2021-10-16")) # set end dates for state-level estimates
data_date <- Sys.Date()  # set date for data creation, set to current date to allow more portability


#variable for whether or not to include the state tagged sequences
state_source <- "state_tag_included" #argument indicating whether to include state tagged data
tag <- paste0("_",state_source,"_Run", opts$run_number) #a tag for the filename to indicate which run from Lab TF request the results are for


#arguments to indicate whether sublineages should be aggregated to parent lineage
P.1_agg=TRUE
B.1.351_agg=TRUE
AY_agg=TRUE 
Q.1_3_agg=TRUE
B.1.621_agg=TRUE
B429_7_agg=TRUE

# List of variants to track (not just VOC or VOI):
# VOCs
# The current branch has the following custom defined lineages:
# - AY.4.2+ - AY.4.2 with Y145H and A222V
# - AY.35+ - AY.35 with E484Q
#
# All other lineages (including AY.4.2 and AY.35) are from default pangolin calls.
if(length(grep("Run1",tag))>0){
  
  voc=c("AY.1",
        "AY.2",
        "B.1.1.7",#*
        "B.1.617.2"#*
  )
  
} else if(length(grep("Run2",tag))>0) {
  voc=c("AY.1",
        "AY.2",
        "AY.3",
        "AY.3.1",
        "AY.10",
        "AY.14",
        "AY.16",
        "AY.20",
        "AY.24",
        "AY.25",
        "AY.26",
        "AY.35",
        "AY.39",
        "AY.4",
        "B.1.1.7",#*
        "B.1.617.2"#* (non-enumerated AYs aggregated)
  )
  
  #State-level run
} else if(length(grep("Run3",tag))>0) {
  voc=c("B.1.1.7",# with Q.1 to 8* 
        "B.1.351", #and B.1.351.* 
        "P.1", #and P.* 
        "B.1.617.2", #and AY.3-AY.25* 
        "AY.1", 
        "AY.2", 
        "B.1.427",#/B.1.429* 
        "B.1.525", 
        "B.1.526", 
        "B.1.617.1", 
        "B.1.617.3", 
        "P.2", 
        "B.1.621"# and B.1.621.1* 
  ) 
}


#Argument determining whether figures should be output as jpgs
fig_gen_run = TRUE

`%notin%` <- Negate(`%in%`)

# Some parameters defining what is modeled and displayed ---------------
n_top = 10 # Top by variant share that must be included in output
n_recent_weeks = 7 # Window for estimates (focus to top variants in this model)
model_weeks = 20 # Lookback for modeling (past 16 weeks of collection dates)
share_cutoff = 0.01 # Criterion for inclusion in model (i.e to be included in model weighted share must be at least 0.01 in the n_recent_weeks)
week0day1 = get0("week0day1", ifnotfound=as.Date("2020-01-05"))
current_week = as.numeric(as.Date(data_date) - week0day1)%/%7 #as.numeric(Sys.Date() - week0day1)%/%7


display_option = c("top7", "voc")[1] #define the display option for plotting
display_lookback = 8 # Lookback for display
mean_generation_time = 6/7 # weeks; CDC proposed modeling scenarios 2021-03-19