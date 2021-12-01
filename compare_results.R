################################################################################
# This script is intended to compare the .csv results output from multiple runs
# of "weekly_variant_report_nowcast.R"
# The files that it compares includes:
# - "state_weighted_roll4wk_KGCI_svyNEW_YYYY-MM-DD_state_tag_included_Run3.csv"
# - "updated_nowcast_fortnightly_YYYY-MM-DD_state_tag_included_Run1.csv"
# - "updated_nowcast_fortnightly_YYYY-MM-DD_state_tag_included_Run2.csv"
# - "updated_nowcast_weekly_YYYY-MM-DD_state_tag_included_Run1.csv"
# - "updated_nowcast_weekly_YYYY-MM-DD_state_tag_included_Run2.csv"
# - "variant_share_weekly_weighted_KGCI_svyNEW_YYYY-MM-DD_state_tag_included_Run1.csv"
# - "variant_share_weekly_weighted_KGCI_svyNEW_YYYY-MM-DD_state_tag_included_Run2.csv"
# - "variant_share_weighted_KGCI_svyNEW_YYYY-MM-DD_state_tag_included_Run1.csv"
# - "variant_share_weighted_KGCI_svyNEW_YYYY-MM-DD_state_tag_included_Run2.csv"
# - "wow_growth_variant_shareYYYY-MM-DD_state_tag_included_Run2.csv"

# currently script only compares the estimated "Share" (i.e. not CI's)
################################################################################

# Put results of each run into seperate folders & specify them here ------------
# setwd(dir = '/scicomp/home-pure/rsv4/labTF/proportion_model/git')
folder1 <- 'production_run/results'
folder2 <- 'results'
# if 1 run has incomplete results, put those in "folder1"

# Make the comparison ----------------------------------------------------------
library(tidyverse)

# vector of files in folder 1
lf1 = file.path(folder1,
                list.files(path = file.path(getwd(),
                                            folder1),
                           pattern = '*.csv'))

# vector of files in folder 2
lf2_all = file.path(folder2,
                list.files(path = file.path(getwd(),
                                            folder2),
                           pattern = '*.csv'))

# subset the file names from folder 2 to match the files in folder 1
lf2 = vector()
for(ff in seq(lf1)){
  # get 1 of the filenames from folder 1
  file = basename(lf1[ff])

  # replace the date with the grep pattern that matches a date
  pattern = sub(pattern = '[:0-9:]{4}-[:0-9:]{2}-[:0-9:]{2}',
                replacement = '[:0-9:]{4}-[:0-9:]{2}-[:0-9:]{2}',
                x = file)

  # get the matching file from the 2nd folder (file name can only differ by date)
  lf2[ff] <- lf2_all[grep(pattern = pattern, x = sapply(X = lf2_all, FUN = basename))]
}


# empty dataframe to hold results
diffs <- data.frame(
  file1 = NA,
  file2 = NA,
  same  = NA,
  comp_column = NA,
  max_diff = NA
)[0,]

# number of decimals to round before comparing values in files that do not match exactly
n_round = 3

# loop through files
for(ii in seq(lf1)){

  # read in the files
  df1 = read.csv(file = lf1[ii])
  df2 = read.csv(file = lf2[ii])
  names(df1) <- tolower(names(df1))
  names(df2) <- tolower(names(df2))

  # add filenames to the results
  diffs[ii,'file1'] <- lf1[ii]
  diffs[ii,'file2'] <- lf2[ii]

  # check if files are exactly equal
  diffs[ii,'same'] <- isTRUE(all.equal(df1, df2))

  # sort them the same way just in case the variants were listed in a different order
  columns_to_sort <- names(df1)[ (names(df1) %in% c('state',
                                                     'roll_fourweek_ending',
                                                     'variant',
                                                     'usa_or_hhsregion',
                                                     'fortnight_ending',
                                                     'week_end'))]
  df1 <- df1 %>% arrange(!!! rlang::syms(columns_to_sort))
  df2 <- df2 %>% arrange(!!! rlang::syms(columns_to_sort))

  # if files are not exactly equal, check if they're approximately equal
  if(!diffs[ii,'same']){

    # print out that they are not exactly equivalent
    print(paste0(basename(lf1[ii]), ' files are not exactly equal.'))

    # make sure the files contain the same columns
    if( !all(names(df1) == names(df2)) ){
      errorCondition(message = paste('Column names do not match!\n',
                                     lf1[ii],
                                     '\n',
                                     lf2[ii]))
    }
    # make sure the files contain the same variants
    # dplyr::setequal(df1$variant_1, df2$variant)
    if( !all.equal(target = df1$variant, current = df2$variant) ){
      errorCondition(message = paste('Variant names do not match!\n',
                                     lf1[ii],
                                     '\n',
                                     lf2[ii]))
    }

    # change the column names from one of the files so that names are unique
    df1_1 = df1
    names(df1_1) = paste0(names(df1), '_1')

    # Check if the files are approximately equal after rounding
    # (for files that have "Share" columns)
    if( grepl(pattern = '^updated_nowcast|^variant_share|^state_weighted_roll',
             x = basename(lf2[ii]) ) ){


      # calculate the difference in values between the two files
      diffs[ii,'max_diff'] <- cbind(df1_1, df2) %>%
        mutate(
          share_diff    = round(share_1 - share, n_round),
          share_lo_diff = round(share_lo_1 - share_lo, n_round),
          share_hi_diff = round(share_hi_1 - share_hi, n_round),
        ) %>%
        pull(share_diff) %>%
        max(abs(.), na.rm = T)

      diffs[ii, 'comp_column'] = "Share"

      print(paste0('maximum difference in predicted share = ', diffs[ii,'max_diff']))
    }

    # Check if the files are approximately equal after rounding
    # (for files that have "share" & "growth_rate" columns)
    if( grepl(pattern = '^wow_growth_variant',
              x = basename(lf2[ii]) ) ){

      # calculate the difference in values between the two files
      df_comp <- df1_1 %>% # dplyr::arrange(., variant_1) %>%
        cbind(.,
              df2) %>%  # %>% dplyr::arrange(., variant)) %>%
        mutate(
          share_diff = round(variant_share_1 - variant_share, n_round),
          growth_rate_diff = round(growth_rate_1 - growth_rate, n_round)
        )

      diffs[ii,'max_diff'] <- df_comp %>%
        select(share_diff, growth_rate_diff) %>%
        max(abs(.), na.rm = T)

      diffs[ii, 'comp_column'] = ifelse(
        test = diffs[ii,'max_diff'] %in% df_comp$share_diff,
        yes = 'Share',
        no = 'Growth Rate'
      )

      print(paste0('maximum difference in predicted share or growth rate = ', diffs[ii,'max_diff']))
    }
  }
}

# view similarities/differences
diffs
#