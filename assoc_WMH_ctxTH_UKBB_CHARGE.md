Associations between cortical thickness and WMH in UKBB and CHARGE
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
ukbb = fread("../results/assoc_res_sex_specific_and_combined_2022-04-20.tsv")
ukbb_sex_spcific = subset(ukbb,sex!="sex-combined")
ukbb_sex_combined1 = subset(ukbb,sex=="sex-combined")
ukbb_sex_combined2 = subset(ukbb,sex=="sex-combined")
ukbb_sex_combined1$sex = "F"
ukbb_sex_combined2$sex = "M"
ukbb2 = rbind(ukbb_sex_spcific,ukbb_sex_combined1,ukbb_sex_combined2)
head(ukbb2)
```

    ##    age_cat    Estimate         SE      roi sex  P
    ## 1: [45,53] -0.06278199 0.02477923 bankssts   F NA
    ## 2: (53,56] -0.08781686 0.02743803 bankssts   F NA
    ## 3: (56,59] -0.09095723 0.02568741 bankssts   F NA
    ## 4: (59,62] -0.10890190 0.02377321 bankssts   F NA
    ## 5: (62,64] -0.11214055 0.02933650 bankssts   F NA
    ## 6: (64,67] -0.07770050 0.02350095 bankssts   F NA

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

``` r
cat.levels = c("[45,53]", "(53,56]", "(56,59]", "(59,62]", "(62,64]", "(64,67]", "(67,69]", 
               "(69,71]", "(71,73]", "(73,81]","UKBB_all",'CHARGE')
ukbb_charge = ukbb_charge %>% mutate(age_cat = factor(age_cat,levels=cat.levels))
L = length(cat.levels)
COLS <- c(colorRampPalette(c("black","orange"))(L-2),
          "darkred","dodgerblue4")
names(COLS) <- cat.levels

roi.34.ordered = c("superiorparietal", "postcentral", 
"precuneus", "cuneus", "caudalmiddlefrontal", "pericalcarine", 
"paracentral", "precentral", "inferiorparietal", "transversetemporal", 
"lateraloccipital", "superiorfrontal", "lingual", "bankssts", 
"parstriangularis", "temporalpole", "entorhinal", "middletemporal", 
"posteriorcingulate", "isthmuscingulate", "supramarginal", "frontalpole", 
"rostralmiddlefrontal", "parsopercularis", "superiortemporal", 
"parahippocampal", "inferiortemporal", "parsorbitalis", "fusiform", 
"caudalanteriorcingulate", "lateralorbitofrontal", "medialorbitofrontal", 
"rostralanteriorcingulate", "insula")
ukbb_charge = ukbb_charge %>% mutate(roi = factor(roi,levels=roi.34.ordered))

age_specific_coef = ukbb_charge
age_specific_coef = age_specific_coef %>% arrange(roi)
```

\##3.1 Plot age_cat-by-WMH interaction estimates in all age bins and
UKBB_All and CHARGE

``` r
pF = subset(age_specific_coef,sex=="F") %>%
  ggplot(aes(x=roi,y=Estimate,color=age_cat,group=age_cat)) +
  geom_point(alpha=0.5) + 
  geom_path(alpha=0.5) + 
  scale_color_manual(values=COLS) + 
  geom_hline(yintercept = 0)+
  xlab(NULL) + 
  #ylab("Estimate (adjusted for imaging centre)")+ 
  theme(legend.position = 'none')+
  # theme(axis.text.x = element_text(angle = 90,hjust =1,vjust=0.5)) + #when not coord_flip() 
  # coord_cartesian(ylim=(range(age_specific_coef$Estimate))) + #when not coord_flip()
  coord_flip() +
  ylim((range(age_specific_coef$Estimate)))+#only when coord_flip
  ggtitle("UKBB Female + CHARGE")
# male
pM = subset(age_specific_coef,sex=="M") %>%
  ggplot(aes(x=roi,y=Estimate,color=age_cat,group=age_cat)) +
  geom_point(alpha=0.5) + 
  geom_path(alpha=0.5) + 
  scale_color_manual(values=COLS) + 
  geom_hline(yintercept = 0)+
  # theme(axis.text.x = element_text(angle = 90,hjust =1,vjust=0.5)) + #when not coord_flip()
  # coord_cartesian(ylim=(range(age_specific_coef$Estimate))) + #when not coord_flip()
  # ylab(NULL) + 
  xlab(NULL) + 
  coord_flip() +
  ylim((range(age_specific_coef$Estimate)))+
  ggtitle("UKBB Male + CHARGE")

print(pF+pM)
```

![](assoc_WMH_ctxTH_UKBB_CHARGE_files/figure-gfm/forestplot-1.png)<!-- -->

``` r
ggsave("../results/age_specific_beta_F_M_wi_CHARGE_UKBB.png",
       width=12,height=6.5,units="in")

pdf("../results/age_specific_beta_F_M_wi_CHARGE_UKBB.pdf",
       width=12,height=6.5)
print(pF+pM)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

\##3.2 Plot age_cat-by-WMH interaction estimates in youngest and oldest
age bins, UKBB_All, and CHARGE

``` r
selelct_cat.levels = c("[45,53]", "(73,81]", "UKBB_all", "CHARGE")
pF = subset(age_specific_coef,sex=="F" & age_cat %in% selelct_cat.levels) %>%
  ggplot(aes(x=roi,y=Estimate,color=age_cat,group=age_cat)) +
  geom_point(alpha=0.5) + 
  geom_path(alpha=0.5) + 
  scale_color_manual(values=COLS) + 
  geom_hline(yintercept = 0)+
  xlab(NULL) + 
  #ylab("Estimate (adjusted for imaging centre)")+ 
  theme(legend.position = 'none')+
  # theme(axis.text.x = element_text(angle = 90,hjust =1,vjust=0.5)) + #when not coord_flip() 
  # coord_cartesian(ylim=(range(age_specific_coef$Estimate))) + #when not coord_flip()
  coord_flip() +
  ylim((range(age_specific_coef$Estimate)))+#only when coord_flip
  ggtitle("UKBB Female + CHARGE")
# male
pM = subset(age_specific_coef,sex=="M" & age_cat %in% selelct_cat.levels) %>%
  ggplot(aes(x=roi,y=Estimate,color=age_cat,group=age_cat)) +
  geom_point(alpha=0.5) + 
  geom_path(alpha=0.5) + 
  scale_color_manual(values=COLS) + 
  geom_hline(yintercept = 0)+
  # theme(axis.text.x = element_text(angle = 90,hjust =1,vjust=0.5)) + #when not coord_flip()
  # coord_cartesian(ylim=(range(age_specific_coef$Estimate))) + #when not coord_flip()
  # ylab(NULL) + 
  xlab(NULL) + 
  coord_flip() +
  ylim((range(age_specific_coef$Estimate)))+
  ggtitle("UKBB Male + CHARGE")

print(pF+pM)
```

![](assoc_WMH_ctxTH_UKBB_CHARGE_files/figure-gfm/forestplot_old_young-1.png)<!-- -->

``` r
ggsave("../results/age_specific_beta_F_M_wi_old_young_UKBB_CHARGE.png",
       width=12,height=6.5,units="in")
pdf("../results/age_specific_beta_F_M_wi_old_young_UKBB_CHARGE.pdf",
       width=12,height=6.5)
print(pF+pM)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

\##3.3 Plot age_cat-by-WMH interaction estimates in youngest and oldest
age bins, UKBB_All, and CHARGE (z-scored)

``` r
selelct_cat.levels = c("[45,53]", "(73,81]", "UKBB_all", "CHARGE")
dF= subset(age_specific_coef,sex=="F" & age_cat %in% selelct_cat.levels) 
dF$zEstimate = dF %>% 
  group_by(age_cat) %>%
  .$Estimate %>% scale()
pF = dF %>%
  ggplot(aes(x=roi,y=scale(Estimate),color=age_cat,group=age_cat)) +
  geom_point(alpha=0.5) + 
  geom_path(alpha=0.5) + 
  scale_color_manual(values=COLS) + 
  geom_hline(yintercept = 0)+
  xlab(NULL) + 
  #ylab("Estimate (adjusted for imaging centre)")+ 
  theme(legend.position = 'none')+
  # theme(axis.text.x = element_text(angle = 90,hjust =1,vjust=0.5)) + #when not coord_flip() 
  # coord_cartesian(ylim=(range(age_specific_coef$Estimate))) + #when not coord_flip()
  coord_flip() +
  #ylim((range(age_specific_coef$Estimate)))+#only when coord_flip
  ggtitle("UKBB Female + CHARGE")
# male
dM = subset(age_specific_coef,sex=="M" & age_cat %in% selelct_cat.levels)
dM$zEstimate = dM %>% 
  group_by(age_cat) %>%
  .$Estimate %>% scale()
pM = dM %>%
  ggplot(aes(x=roi,y=scale(Estimate),color=age_cat,group=age_cat)) +
  geom_point(alpha=0.5) + 
  geom_path(alpha=0.5) + 
  scale_color_manual(values=COLS) + 
  geom_hline(yintercept = 0)+
  # theme(axis.text.x = element_text(angle = 90,hjust =1,vjust=0.5)) + #when not coord_flip()
  # coord_cartesian(ylim=(range(age_specific_coef$Estimate))) + #when not coord_flip()
  # ylab(NULL) + 
  xlab(NULL) + 
  coord_flip() +
  #ylim((range(age_specific_coef$Estimate)))+
  ggtitle("UKBB Male + CHARGE")

print(pF+pM)
```

![](assoc_WMH_ctxTH_UKBB_CHARGE_files/figure-gfm/forestplot_old_young_zscored-1.png)<!-- -->

``` r
ggsave("../results/age_specific_beta_F_M_wi_old_young_UKBB_CHARGE_zscored.png",
       width=12,height=6.5,units="in")
pdf("../results/age_specific_beta_F_M_wi_old_young_UKBB_CHARGE_zscored.pdf",width=12,height=6.5)
print(pF+pM)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

\#4. Correlation plot

``` r
assoc_table_wide <- dcast(age_specific_coef, roi+sex ~ age_cat, 
                          value.var=c("Estimate","SE"))
assoc_table_wide2 <- dcast(age_specific_coef, roi ~ age_cat+sex, 
                          value.var=c("Estimate","SE"))

head(assoc_table_wide2,2)
```

    ##                 roi Estimate_[45,53]_F Estimate_[45,53]_M Estimate_(53,56]_F
    ## 1: superiorparietal        0.004033731       -0.013424000       -0.008112359
    ## 2:      postcentral       -0.084189814        0.003385468       -0.047403302
    ##    Estimate_(53,56]_M Estimate_(56,59]_F Estimate_(56,59]_M Estimate_(59,62]_F
    ## 1:        -0.01039685        -0.00421088         0.01273981         0.02162196
    ## 2:        -0.06148335        -0.04250030        -0.03694747        -0.04521201
    ##    Estimate_(59,62]_M Estimate_(62,64]_F Estimate_(62,64]_M Estimate_(64,67]_F
    ## 1:        0.062718487         0.11771436       0.0515872080         0.03114953
    ## 2:        0.003903627         0.05013833      -0.0005515082        -0.03834677
    ##    Estimate_(64,67]_M Estimate_(67,69]_F Estimate_(67,69]_M Estimate_(69,71]_F
    ## 1:         0.08982861         0.13869514        0.106645794         0.13796226
    ## 2:         0.01694470         0.03677877       -0.001512708         0.03752227
    ##    Estimate_(69,71]_M Estimate_(71,73]_F Estimate_(71,73]_M Estimate_(73,81]_F
    ## 1:         0.11097451         0.18674006         0.07498810         0.16185668
    ## 2:         0.03049615         0.06807955        -0.03196127         0.07996383
    ##    Estimate_(73,81]_M Estimate_UKBB_all_F Estimate_UKBB_all_M Estimate_CHARGE_F
    ## 1:         0.10668900         0.049498874         0.049498874       0.003342855
    ## 2:        -0.02108222        -0.004343942        -0.004343942      -0.018857657
    ##    Estimate_CHARGE_M SE_[45,53]_F SE_[45,53]_M SE_(53,56]_F SE_(53,56]_M
    ## 1:       0.003342855   0.02470729   0.02866763   0.02735837   0.03205866
    ## 2:      -0.018857657   0.02504831   0.02752369   0.02773598   0.03077940
    ##    SE_(56,59]_F SE_(56,59]_M SE_(59,62]_F SE_(59,62]_M SE_(62,64]_F
    ## 1:   0.02561283   0.03040230   0.02370419   0.02779088   0.02925132
    ## 2:   0.02596635   0.02918913   0.02403136   0.02668193   0.02965506
    ##    SE_(62,64]_M SE_(64,67]_F SE_(64,67]_M SE_(67,69]_F SE_(67,69]_M
    ## 1:   0.03151446   0.02343271   0.02453604   0.02822495   0.02919296
    ## 2:   0.03025692   0.02375614   0.02355697   0.02861452   0.02802805
    ##    SE_(69,71]_F SE_(69,71]_M SE_(71,73]_F SE_(71,73]_M SE_(73,81]_F
    ## 1:   0.02979344   0.02762985   0.03374499   0.03201495   0.03225255
    ## 2:   0.03020466   0.02652731   0.03421075   0.03073743   0.03269771
    ##    SE_(73,81]_M SE_UKBB_all_F SE_UKBB_all_M SE_CHARGE_F SE_CHARGE_M
    ## 1:   0.02729426   0.005400042   0.005400042 0.005831988 0.005831988
    ## 2:   0.02620512   0.005323613   0.005323613 0.007822213 0.007822213

``` r
cols_select = c(names(assoc_table_wide2)[str_detect(names(assoc_table_wide2),"Estimate_")& str_detect(names(assoc_table_wide2),"_F")],
                names(assoc_table_wide2)[str_detect(names(assoc_table_wide2),"Estimate_")& str_detect(names(assoc_table_wide2),"_M")])

library(corrplot)
```

    ## corrplot 0.92 loaded

``` r
corm = cor(subset(assoc_table_wide2,select=cols_select),method = 's')
rownames(corm) <- colnames(corm) <- str_remove(cols_select,"Estimate_")
rownames(corm) <- colnames(corm) <- str_split(rownames(corm),"_",simplify = T)[,1]

#ggcorrplot: https://github.com/kassambara/ggcorrplot--------------------------
library(ggcorrplot)
dim_corm = nrow(corm)
pF=ggcorrplot(corm[1:(dim_corm/2),(dim_corm/2):1],
              outline.col = "white",lab = TRUE,lab_size=2.5)+
  ggtitle("Female + whole UKBB + CHARGE") +
  theme(text = element_text(size = 10),
        legend.position = 'none',
        axis.text.y = element_text(size=8),
        # axis.text.x = element_blank()
        # axis.text.y = element_blank(),
        axis.text.x = element_text(size=8)
  )
pM=ggcorrplot(corm[(dim_corm/2+1):dim_corm,dim_corm:(dim_corm/2+1)],
              outline.col = "white",lab = TRUE,lab_size=2.5)+
  ggtitle("Male + whole UKBB + CHARGE")+
  theme(text = element_text(size = 11),
        legend.position = 'none',
        # axis.text.y = element_text(size=8),
        # axis.text.x = element_blank()
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=8)
  )
pFpM=ggcorrplot(corm[(dim_corm/2+1):dim_corm,(dim_corm/2):1],
                outline.col = "white",
                lab = TRUE,lab_size=2.5)+
  ggtitle("Female (rows) vs. Male (columns)") +
  theme(text = element_text(size = 11),
             legend.title=element_text(size=8), 
             legend.text=element_text(size=8),
        # axis.text.y = element_text(size=8),
        # axis.text.x = element_blank()
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=8)
        )
pCorr = pF+pM+pFpM
print(pCorr)
```

![](assoc_WMH_ctxTH_UKBB_CHARGE_files/figure-gfm/corrplot-1.png)<!-- -->

``` r
ggsave('../results/corr_plot_spearman_age_specific_slopes_wi_CHARGE_UKBB_all.png',
       width=12,height=5,units='in')
```

\#5. Genotype

\#6. For VH:

``` r
to_write=subset(ukbb_charge,(age_cat=="UKBB_all" & sex=="F"), select=c(roi,Estimate))
write_tsv(to_write,"~/OneDrive - SickKids/4752955/profile_WMH_ctxTH_UKBB_all_assoc_2022-04-19.tsv")
```
