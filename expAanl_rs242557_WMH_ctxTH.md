Exploratory analysis of WMH, cortical thickness and rs242557
================

# 1. Load packages

``` r
rm(list=ls());gc()
```

    ##          used (Mb) gc trigger (Mb) limit (Mb) max used (Mb)
    ## Ncells 467069 25.0    1000788 53.5         NA   667388 35.7
    ## Vcells 870772  6.7    8388608 64.0      16384  1824350 14.0

``` r
library(tidyverse)
library(data.table)
library(stringr)
library(patchwork)
library(psych)
source("../scripts/add_rs242577.R")
source("../scripts/getFit_interaction_sex_combined.R")
source("../scripts/getFit_interaction_geno.R")
source("../scripts/getFit_interaction.R")
```

# 2. Read data files

``` r
# read and merge data files
d1 = fread('../data/hardcalled_genotypes_for_MAPT_SNP_ctxTH_WMH_for_tomas_2022-04.19.tsv') %>% dplyr::rename(eid=IID)
d2 = fread('../data/anal_dat_ctxTH_WMH_age_2022-04-25.txt')
d = merge(d1,d2)#n=31,082 (participants with genotype information)
rm(d1,d2); gc()
```

    ##           used (Mb) gc trigger  (Mb) limit (Mb) max used (Mb)
    ## Ncells 1179193 63.0    1844201  98.5         NA  1844201 98.5
    ## Vcells 6194281 47.3   15976050 121.9      16384 11660663 89.0

# 3. Wrangle data

    ##        eid age1 sex assessment_centre   WMH     ICV BrainVolume   logWMH
    ## 1: 1000147   65   F           Cheadle  2029 1447730     1173770 3.307496
    ## 2: 1000353   79   F         Newcastle 15610 1542390     1142510 4.193431
    ## 3: 1000430   61   F           Cheadle   952 1482940     1206850 2.979093
    ## 4: 1000508   51   F         Newcastle   832 1475470     1098740 2.920645
    ## 5: 1000636   51   F           Cheadle  2860 1516400     1128270 3.456518
    ## 6: 1000655   70   M           Cheadle 11853 1622660     1220270 4.073865
    ##          rntWMH      age.c     age.c2     age_cat age_cat2
    ## 1: -0.351031608  0.1629417 0.02655001   (63,66.6]  (64,67]
    ## 2:  1.581078404  2.0338525 4.13655601   (77.4,81]  (73,81]
    ## 3: -1.187517682 -0.3716042 0.13808967   (59.4,63]  (59,62]
    ## 4: -1.339049184 -1.7079690 2.91715814 (48.6,52.2]  [45,53]
    ## 5:  0.006774287 -1.7079690 2.91715814 (48.6,52.2]  [45,53]
    ## 6:  1.319634491  0.8311242 0.69076737 (66.6,70.2]  (69,71]

    ## [1] 34

    ## bankssts
    ## caudalanteriorcingulate
    ## caudalmiddlefrontal
    ## cuneus
    ## entorhinal
    ## fusiform
    ## inferiorparietal
    ## inferiortemporal
    ## isthmuscingulate
    ## lateraloccipital
    ## lateralorbitofrontal
    ## lingual
    ## medialorbitofrontal
    ## middletemporal
    ## parahippocampal
    ## paracentral
    ## parsopercularis
    ## parsorbitalis
    ## parstriangularis
    ## pericalcarine
    ## postcentral
    ## posteriorcingulate
    ## precentral
    ## precuneus
    ## rostralanteriorcingulate
    ## rostralmiddlefrontal
    ## superiorfrontal
    ## superiorparietal
    ## superiortemporal
    ## supramarginal
    ## frontalpole
    ## temporalpole
    ## transversetemporal
    ## insula

# 4. Associtions without the SNP

\#5. Explore: cortical thickness vs. genotypes, adjusting for assessment
centre, age, age2 and sex

-   The SNP appears not associated with insular thickness - no main
    effect -\> *no interaction*.

``` r
d_anal = merge(subset(d,select=c(eid,rs242557_chr17_G_A)),
               covdat)
i=34
roii = roi.34[i];print(roii)
```

    ## [1] "insula"

``` r
d_roii = merge(d_anal,subset(avg.th,select=c('eid',roii)))
names(d_roii)[names(d_roii)==roii] <- 'y'
d_roii = d_roii %>% mutate(y = ifelse(abs(scale(y)[,1]) > 4, NA, y))

# combined:
lm0 = lm(y~(assessment_centre)*sex, data=d_roii,
         na.action = na.exclude)
d_roii[['y.resid']] = resid(lm0)+coef(lm0)[1]
d_roii %>% ggplot(aes(x=rs242557_chr17_G_A,y=y.resid,fill=sex)) + 
  geom_boxplot()+
  ggtitle(roii) + 
  ylab('adjusted cortical thickness')
```

    ## Warning: Removed 12 rows containing non-finite values (stat_boxplot).

![](expAanl_rs242557_WMH_ctxTH_files/figure-gfm/explore-1.png)<!-- -->

``` r
lm.geno = lm(y.resid~rs242557_chr17_G_A,data=d_roii)
summary(lm.geno)
```

    ## 
    ## Call:
    ## lm(formula = y.resid ~ rs242557_chr17_G_A, data = d_roii)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.55956 -0.09253  0.00162  0.09410  0.53962 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error  t value Pr(>|t|)    
    ## (Intercept)          2.942567   0.002053 1433.135   <2e-16 ***
    ## rs242557_chr17_G_AAG 0.001813   0.002356    0.770    0.441    
    ## rs242557_chr17_G_AGG 0.001539   0.002426    0.634    0.526    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1399 on 31067 degrees of freedom
    ##   (12 observations deleted due to missingness)
    ## Multiple R-squared:  1.939e-05,  Adjusted R-squared:  -4.499e-05 
    ## F-statistic: 0.3011 on 2 and 31067 DF,  p-value: 0.74

``` r
# Female
d_roii.F = subset(d_roii,sex=="F")
lm0 = lm(y~(age.c + age.c2 + assessment_centre), 
         data=d_roii.F, na.action = na.exclude)
d_roii.F[['y.resid']] = resid(lm0)+coef(lm0)[1]
lm.geno.F = lm(y.resid~rs242557_chr17_G_A,data=d_roii.F)
sm.F = summary(lm.geno.F)

# Male
d_roii.M = subset(d_roii,sex=="M")
lm0 = lm(y~(age.c + age.c2 + assessment_centre), 
         data=d_roii.M, na.action = na.exclude)
d_roii.M[['y.resid']] = resid(lm0)+coef(lm0)[1]
lm.geno.M = lm(y.resid~rs242557_chr17_G_A,data=d_roii.M)
sm.M = summary(lm.geno.M)

d_roii.FM = rbind(d_roii.F,d_roii.M)

# boxplot
box.A = d_roii.FM %>% ggplot(aes(x=age_cat2,y=y.resid,
                                fill=rs242557_chr17_G_A)) + 
  geom_boxplot(alpha=0.5)+
  ggtitle(roii) + 
  ylab('adjusted cortical thickness')
box.A
```

    ## Warning: Removed 12 rows containing non-finite values (stat_boxplot).

![](expAanl_rs242557_WMH_ctxTH_files/figure-gfm/explore-2.png)<!-- -->

``` r
#  
abin1='[45,53]'
abin2='(73,81]'
lm.geno.F = lm(y.resid~rs242557_chr17_G_A,data=d_roii.F, subset=age_cat2==abin1)
lm.geno.M = lm(y.resid~rs242557_chr17_G_A,data=d_roii.M, subset=age_cat2==abin1)
lm.geno2.F = lm(y.resid~rs242557_chr17_G_A,data=d_roii.F, subset=age_cat2==abin2)
lm.geno2.M = lm(y.resid~rs242557_chr17_G_A,data=d_roii.M, subset=age_cat2==abin2)
```

\#6. Brain volume: age vs. brain-volume in age-bins

``` r
names(covdat)
```

    ##  [1] "eid"               "age1"              "sex"              
    ##  [4] "assessment_centre" "WMH"               "ICV"              
    ##  [7] "BrainVolume"       "logWMH"            "rntWMH"           
    ## [10] "age.c"             "age.c2"            "age_cat"          
    ## [13] "age_cat2"

``` r
covdat = covdat %>% 
  mutate(ICV_cm3=ICV/10^3,BrainVolume_cm3=BrainVolume/10^3) %>% 
  mutate(ICV.no_outliers=ifelse(abs(scale(ICV_cm3)[,1])>4,NA,ICV_cm3),
         BrainVolume.no_outliers=ifelse(abs(scale(BrainVolume_cm3)[,1])>4,NA,BrainVolume_cm3))

ICV.boxplot = covdat %>% ggplot(aes(y=ICV.no_outliers,color=sex)) +
  geom_boxplot() + theme(legend.position = 'none')
BV.boxplot = covdat %>% ggplot(aes(y=BrainVolume.no_outliers,color=sex)) + 
  geom_boxplot()
print(ICV.boxplot+BV.boxplot)
```

    ## Warning: Removed 254 rows containing non-finite values (stat_boxplot).

    ## Warning: Removed 248 rows containing non-finite values (stat_boxplot).

![](expAanl_rs242557_WMH_ctxTH_files/figure-gfm/brain_volume-1.png)<!-- -->

``` r
cat.levels = c("[45,53]", "(53,56]", "(56,59]", "(59,62]", "(62,64]", "(64,67]", "(67,69]", 
               "(69,71]", "(71,73]", "(73,81]")
covdat = covdat %>% mutate(age_cat2 = factor(age_cat2,levels=cat.levels))
L = length(cat.levels)
COLS <- colorRampPalette(c("black","orange"))(L)
names(COLS) <- cat.levels

getPlot = function(roi,covdat,avg.th){
  yname=roi
  mridat.roi = merge(covdat,subset(avg.th,select=c('eid',yname)))
  names(mridat.roi)[names(mridat.roi)==yname] <- 'y'
  p = mridat.roi %>% ggplot(aes(x=BrainVolume.no_outliers,y=y,color=age_cat2))+
  scale_color_manual(values=COLS) + 
  facet_grid(cols=vars(sex)) +
  # geom_point(alpha=0.5) + 
  ylab(yname)+
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=4),alpha=0.5)
  p
}
p.insula = getPlot('insula',covdat,avg.th) + theme(legend.position = 'n')
p.superiorparietal = getPlot('superiorparietal',covdat,avg.th)
p_all = p.insula + p.superiorparietal
p_all +  plot_layout(guides = 'collect')
```

    ## Warning: Removed 248 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 248 rows containing non-finite values (stat_smooth).

![](expAanl_rs242557_WMH_ctxTH_files/figure-gfm/brain_volume-2.png)<!-- -->
