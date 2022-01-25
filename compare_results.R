library(tidyverse)

folder1 <- 'C:/ctemp/omicron_updates/results'
folder2 <- 'C:/ctemp/omicron_special_run_2021-12-23/results_2022-01-13_reducedvoc'


# compare variant share estimates
res1 <- read.csv(file = paste0(folder1, '/variant_share_weekly_weighted_KGCI_svyNEW_2022-01-13_state_tag_included_Run2_reduced_vocs.csv'))
res2 <- read.csv(file = paste0(folder2, '/variant_share_weekly_weighted_KGCI_svyNEW_2022-01-13_state_tag_included_Run2.csv'))

all.equal(res1, res2)

dat <- rbind(res1 %>% mutate(group = 'updated'),
             res2 %>% mutate(group = 'original')) %>% 
  mutate(region = USA_or_HHSRegion,
         date   = as.Date(WEEK_END))

ggplot(dat %>% filter(region == 'USA'), 
       aes(x = date, y = Share)) + 
  theme_bw() + 
  geom_ribbon(aes(ymin = Share_lo, 
                  ymax = Share_hi, 
                  fill = group),
              alpha = 0.2) + 
  geom_line(aes(color = group)) + 
  facet_grid(rows = vars(Variant))

diff = merge(res1[,c("USA_or_HHSRegion", 
                     "WEEK_END", 
                     "Variant", 
                     "Share")],
             res2[,c("USA_or_HHSRegion", 
                     "WEEK_END", 
                     "Variant", 
                     "Share")],
             all = TRUE, 
             by = c("USA_or_HHSRegion", 
                    "WEEK_END", 
                    "Variant")) %>% 
  mutate(share_diff = Share.x - Share.y) %>% 
  arrange(share_diff)

range(diff$share_diff)




# compare Nowcast estimates
nc1 <- read.csv(file = paste0(folder1, '/updated_nowcast_weekly_2022-01-13_state_tag_included_Run2_reduced_vocs_daily.csv'))
nc2 <- read.csv(file = paste0(folder2, '/updated_nowcast_weekly_2022-01-13_state_tag_included_Run2_daily.csv'))


all.equal(nc1, nc2)

dat <- rbind(nc1 %>% 
               mutate(group = 'updated',
                      week_number = NA) %>% 
               select(c(names(nc2), 'group')),
             nc2 %>% 
               mutate(group = 'original')) %>% 
  mutate(region = USA_or_HHSRegion,
         date = as.Date(date))

(ggplot(dat %>% filter(region == 'USA'), 
       aes(x = date, y = Share)) + 
  theme_bw() + 
  geom_ribbon(aes(ymin = Share_lo, 
                  ymax = Share_hi, 
                  fill = group),
              alpha = 0.2) + 
  geom_line(aes(color = group)) + 
  facet_grid(rows = vars(Variant))) %>% 
  plotly::ggplotly(p = ., 
                   # tooltip = 'text', 
                   width = 1600, 
                   height = 900) %>% 
  plotly::layout(legend = list( y =0.5))

diff = merge(nc1[,c("USA_or_HHSRegion", 
                     "date", 
                     "Variant", 
                     "Share")] %>% 
               mutate(date = as.Date(date)),
             nc2[,c("USA_or_HHSRegion", 
                     "date", 
                     "Variant", 
                     "Share")],
             all = FALSE, 
             by = c("USA_or_HHSRegion", 
                    "date", 
                    "Variant")) %>% 
  mutate(share_diff = Share.x - Share.y) %>% 
  arrange(share_diff)

range(diff$share_diff)

head(diff, n = 10)
tail(diff, n = 10)


sc1 <- readRDS(file = 'results_update1/src.moddat_2022-01-13_state_tag_included_Run2_reduced_vocs.RDS')
sc2 <- readRDS(file = 'results_original/src.moddat_2022-01-13_state_tag_included_Run2reduced_voc.RDS')

all.equal(sc1$SIMPLE_ADJ_WT, 
          sc2$SIMPLE_ADJ_WT)

all.equal(sc1$wts_trimmed, 
          sc2$wts_trimmed)

all.equal(sc1$weights, sc2$myweights)









# IS the difference just caused by weight trimming??
# if so, then the estimates should be the same if we run the original and updated code analyses 

nc1 <- read.csv(file = 'results_update_untrimmed2/updated_nowcast_weekly_2022-01-13_state_tag_included_Run2_reduced_vocs_daily.csv')
nc2 <- read.csv(file = 'results_original_untrimmed/updated_nowcast_weekly_2022-01-13_state_tag_included_Run2_daily.csv')


all.equal(nc1, nc2)

dat <- rbind(nc1 %>% 
               mutate(group = 'updated',
                      week_number = NA) %>% 
               select(c(names(nc2), 'group')),
             nc2 %>% 
               mutate(group = 'original')) %>% 
  mutate(region = USA_or_HHSRegion,
         date = as.Date(date))

(ggplot(dat %>% filter(region == 'USA'), 
        aes(x = date, y = Share)) + 
    theme_bw() + 
    geom_ribbon(aes(ymin = Share_lo, 
                    ymax = Share_hi, 
                    fill = group),
                alpha = 0.2) + 
    geom_line(aes(color = group)) + 
    facet_grid(rows = vars(Variant))) %>% 
  plotly::ggplotly(p = ., 
                   # tooltip = 'text', 
                   width = 1600, 
                   height = 2700) %>% 
  plotly::layout(legend = list( y =0.5))

diff = merge(nc1[,c("USA_or_HHSRegion", 
                    "date", 
                    "Variant", 
                    "Share")] %>% 
               mutate(date = as.Date(date)),
             nc2[,c("USA_or_HHSRegion", 
                    "date", 
                    "Variant", 
                    "Share")],
             all = FALSE, 
             by = c("USA_or_HHSRegion", 
                    "date", 
                    "Variant")) %>% 
  mutate(share_diff = Share.x - Share.y) %>% 
  arrange(share_diff)

range(diff$share_diff)

head(diff, n = 10)
tail(diff, n = 10)
