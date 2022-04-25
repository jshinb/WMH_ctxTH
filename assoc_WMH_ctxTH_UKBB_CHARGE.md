Creating forest and correlation plots for the association profiles of
associations between cortical thickness and WMH in UKBB and CHARGE
cohorts
================

\#1. Load packages

``` r
rm(list=ls())
library(tidyverse)
library(data.table)
library(patchwork)
```

\#2. Read association results

``` r
# read UKBB-results
ukbb = fread("../results/assoc_res_sex_specific_and_combined_2022-04-21.tsv")
ukbb_sex_spcific = subset(ukbb,sex!="sex-combined")
ukbb_sex_combined1 <- ukbb_sex_combined2 <- subset(ukbb,sex=="sex-combined")
ukbb_sex_combined1$sex = "F"
ukbb_sex_combined2$sex = "M"
ukbb2 = rbind(ukbb_sex_spcific,ukbb_sex_combined1,ukbb_sex_combined2)
head(ukbb2)
```

    ##    age_cat    Estimate         SE      roi sex  P       U95CI      L95CI
    ## 1: [45,53] -0.06081402 0.02490723 bankssts   F NA -0.01199585 -0.1096322
    ## 2: (53,56] -0.08820384 0.02746039 bankssts   F NA -0.03438147 -0.1420262
    ## 3: (56,59] -0.09116029 0.02570746 bankssts   F NA -0.04077367 -0.1415469
    ## 4: (59,62] -0.10806890 0.02381419 bankssts   F NA -0.06139309 -0.1547447
    ## 5: (62,64] -0.11182877 0.02933549 bankssts   F NA -0.05433120 -0.1693263
    ## 6: (64,67] -0.07732085 0.02351502 bankssts   F NA -0.03123142 -0.1234103

``` r
rm(list=ls(pattern="ukbb_sex_"))

# read CHARGE-results
charge = fread("../results/CHARGE_profile.csv")
charge[['age_cat']] = "CHARGE"
charge[['sex']] = 'F'
charge2 = charge
charge2$sex = 'M'
charge = rbind(charge,charge2)
charge = charge %>% 
  dplyr::rename(Estimate=beta_CHARGE,
                SE=se_CHARGE,
                P=p_CHARGE,
                roi=vars) %>%
  mutate(U95CI=Estimate+1.96*SE,L95CI=Estimate-1.96*SE,roi=str_remove_all(roi," "))

# merge two sets of results
ukbb_charge = rbind(ukbb2,subset(charge,select=names(ukbb2)))

rm(ukbb,ukbb2,charge,charge2)
```

\#3 Forest plot

\##3.1 Plot age_cat-by-WMH interaction estimates in all age bins and
UKBB_All and CHARGE
![](assoc_WMH_ctxTH_UKBB_CHARGE_files/figure-gfm/forestplot-1.png)<!-- -->

    ## quartz_off_screen 
    ##                 2

\##3.2 Plot age_cat-by-WMH interaction estimates in youngest and oldest
age bins, UKBB_All, and CHARGE
![](assoc_WMH_ctxTH_UKBB_CHARGE_files/figure-gfm/forestplot_old_young-1.png)<!-- -->

    ## quartz_off_screen 
    ##                 2

\##3.3 Plot age_cat-by-WMH interaction estimates in youngest and oldest
age bins, UKBB_All, and CHARGE (z-scored)

\#4. Correlation plot

    ## corrplot 0.92 loaded

![](assoc_WMH_ctxTH_UKBB_CHARGE_files/figure-gfm/corrplot-1.png)<!-- -->
\#5. Genotype

\#6. For VH:

``` r
to_write=subset(ukbb_charge,(age_cat=="UKBB_all" & sex=="F"), select=c(roi,Estimate))
write_tsv(to_write,"~/OneDrive - SickKids/4752955/profile_WMH_ctxTH_UKBB_all_assoc_2022-04-19.tsv")
```
