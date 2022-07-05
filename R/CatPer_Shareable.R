#############################################################
# Categorical Perception.
#############################################################
library(MASS)
library(glmmTMB)
library(lme4)
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(eyetrackingR)
library(PupillometryR)
library(ggplot2); theme_set(theme_classic(base_size = 26))
############################################
#BabyCAT

# opt1 75G56Sim 75BG56Sim (background greenish, within is green, between blue)
# opt2 25B56Sim 25BG56Sim (background blueish, within is blue, between green)
# opt3 2G68 2B68 (far - within green between blue)


babs <- read_csv('data/BabyCat.csv', na = c("", "NA"))
#setup
babs$Condition[babs$`VIS WI STM` == '75G56Sim'] <- 'Green'
babs$Condition[babs$`VIS WI STM` == '25B56Sim'] <- 'Blue'
babs$Condition[is.na(babs$`VIS WI STM`)] <- 'Far'
#targets
babs$Within[babs$`VIS WI ISFXT` == 1 | babs$`VIS FWI ISFXT` == 1] <- 1
babs$Within[is.na(babs$Within)] <- 0
babs$Between[babs$`VIS BT ISFXT` == 1 | babs$`VIS FBT ISFXT` == 1] <- 1
babs$Between[is.na(babs$Between)] <- 0
#Age
babs$Age <- substr(babs$ID, 1, 2)
#out
babs$Other[babs$Within == 0 & babs$Between == 0] <- 1
babs$Other[is.na(babs$Other)] <- 0
#rename
names(babs)[3] <- 'Trial'

babs$Left[babs$`VIS WI STM POS` == 2|babs$`VIS FWI STM POS` == 2] <-  'Within'
babs$Left[babs$`VIS BT STM POS` == 2|babs$`VIS FBT STM POS` == 2] <-  'Between'

babs <- babs %>% 
  mutate(Between_Side = ifelse(Left == 'Between', 'Left', 'Right'))

baby <- babs %>% 
  select('ID', 'Age', 'Trial', 'Timestamp', 'Condition', 'Within', 'Between', 'Between_Side', 'Other')

#number
length(unique(baby$ID))
#N = 63
length(unique(subset(baby, Age ==19)$ID))
#N(19) = 30
#####################################################################
# Add CDI
#####################################################################
CDI <- read.csv('data/Participants_Shareable.csv')
#setup
CDI[CDI==""] <- 0 #coercion is not pretty but works here
CDI[is.na(CDI)] <- 0
CDI[CDI=="U/S"] <- 2
CDI[CDI=="U"] <- 1
CDI[CDI=="No"] <- 0


cdi <- subset(CDI, CDI$BaBYCat == '1')
names(cdi)[1] <- 'ID'
names(cdi)[7] <- 'Age'
#make a knowledge for EITHER colour word
cdi$Knowledge[cdi$Blue == 1 | cdi$Green == 1] <- 1
cdi$Knowledge[is.na(cdi$Knowledge)] <- 0
#select
babscdi <- cdi %>% 
  select('ID', 'Gender', 'AgeTest', 'Age', 'Knowledge')

#reformat both
babscdi$ID <- as.character(babscdi$ID)
baby$ID <- as.character(baby$ID)
babscdi$Age <- as.character(babscdi$Age)
#merge
babydata <- left_join(baby, babscdi, by = c('ID', 'Age'))

bbb <- babydata %>% 
  group_by(ID) %>% 
  summarize(Knowledge = first(Knowledge))
###########################################################
# ETR
###########################################################
baby_track <- make_eyetrackingr_data(babydata,
                                     participant_column = 'ID',
                                     trackloss_column = 'Other',
                                     time_column = 'Timestamp',
                                     trial_column = 'Trial',
                                     aoi_columns = c('Within','Between'),
                                     treat_non_aoi_looks_as_missing = T)

#trackloss
trackloss <- trackloss_analysis(data = baby_track)
trackloss_subjects <- unique(trackloss[,c('ID','TracklossForParticipant')])

#remove files with trackloss
clean_baby <- clean_by_trackloss(data = baby_track,
                                 trial_prop_thresh = 0.9) #was 6

#(removed 170 trials, removed 9 participants)
trackloss_clean <- trackloss_analysis(data = clean_baby)
trackloss_clean_subjects <- unique(trackloss_clean[, c('ID','TracklossForParticipant')])
#find out trackloss
trackloss_mean<-mean(1 - trackloss_clean_subjects$TracklossForParticipant)
trackloss_sd<-sd(1- trackloss_clean_subjects$TracklossForParticipant)
#See trials contributed by each participant
interim_summary <- describe_data(baby_track, 'Within', 'ID')
final_summary <- describe_data(clean_baby, 'Within', 'ID')
mean_num_trials <- mean(final_summary$NumTrials)
sd_sum_trials <- sd(final_summary$NumTrials)
#count
length(unique(final_summary$ID))
#63
length(unique(subset(clean_baby, Age ==19)$ID))
#N(19) = 30

#Ages
mean(subset(babscdi, Age == 19)$AgeTest) #19.29
sd(subset(babscdi, Age == 19)$AgeTest) #0.47
mean(subset(babscdi, Age == 12)$AgeTest) #12.23
sd(subset(babscdi, Age == 12)$AgeTest) #0.39

###################################################################
# Analysis
###################################################################
#including between side here
clean_baby_analyze <- subset_by_window(clean_baby,
                                       window_start_time = 400,
                                       window_end_time = 3000,
                                       remove = T,
                                       rezero = F)

babycat <- make_time_sequence_data(clean_baby_analyze,
                                   time_bin_size = 100,
                                   aois = 'Between',
                                   predictor_columns = c('Age', 'Condition', 
                                                         'Knowledge'),# 'Between_Side'
                                   summarize_by = 'ID')

plot(babycat, predictor_column = 'Condition') +
  facet_wrap(~Age)

#cut
babycat2 <- subset(babycat, SamplesTotal > 0)
#center
babycat2$Age_C <- ifelse(babycat2$Age == 19, 0.5, -0.5)
babycat2$Know <- ifelse(babycat2$Knowledge == 1, 'Known', 'Unknown')
babycat2$Know_C <- ifelse(babycat2$Know == 'Known', 0.5, -0.5)

#options(contrasts = c('contr.sum','contr.poly'))

model1.0 <- glmmTMB(cbind(SamplesInAOI, SamplesTotal - SamplesInAOI) ~
                    (ot1 + ot2 + ot3 + ot4) * Age_C * Condition +
                    (ot1 + ot2 + ot3 + ot4|ID:Condition), #maximalish random effect
                    contrasts = list(Condition = contr.sum),
                  family = 'binomial',
                  data = babycat2)

# model1.01 <- glmer(cbind(SamplesInAOI, SamplesTotal - SamplesInAOI) ~
#                       (ot1 + ot2 + ot3 + ot4 ) * Age_C * Condition +
#                       (ot1 + ot2 + ot3 + ot4 |ID:Condition), #maximalish random effect
#                     family = 'binomial',
#                    # hessianMethod = 'simple',
#                     data = babycat2)


summary(model1.0)
car::Anova(model1.0, type = 'III')
a <- DHARMa::simulateResiduals(model1.0, plot = T) #outlier test inflated positive
DHARMa::testOutliers(a, type = 'bootstrap') #NS should be bootstrap for binomial
DHARMa::testUniformity(a)
DHARMa::testDispersion(a)
DHARMa::plotResiduals(a, babycat2$Condition)
DHARMa::plotResiduals(a, quantreg = T)
drop1(model1.0, test = 'Chisq')


babycat2$model1.0 <- predict(model1.0, type = 'response')

babycat2 %>% 
  ggplot(aes(x = Time, y = Prop,
             colour = Condition, shape = Condition, linetype = Condition, fill = Condition)) +
  stat_summary(geom = 'ribbon', fun.data = 'mean_se', alpha = 0.3, colour = NA) +
  stat_summary(aes(y = model1.0), geom = 'line', fun = 'mean', size = 1) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = 'red') +
  labs(x = 'Time in trial (ms)',
       y = 'Prop looks to between colour') +
  scale_color_manual(values = c('blue', 'forestgreen', 'green')) +
  scale_fill_manual(values = c('blue', 'forestgreen', 'green')) +
  facet_wrap(~Age)

#knowledge (19mo only)
# babycat3 <- subset(babycat2, Age == '19')
bbct <- babycat2 %>% 
  group_by(ID) %>% 
  summarise(Knowledge = first(Knowledge)) %>% 
  ungroup()

babycat3 <- babycat2 %>% 
  filter(ID != 120323)

model2.0 <- glmmTMB(cbind(SamplesInAOI, SamplesTotal - SamplesInAOI) ~
                    (ot1 + ot2 + ot3 + ot4) * Know_C * Condition +
                    (ot1 + ot2 + ot3 + ot4 |ID:Condition),
                    contrasts = list(Condition = contr.sum),
                  family = 'binomial',
                  data = babycat3)

summary(model2.0)
car::Anova(model2.0, type = 'III')



babycat3$model2 <- predict(model2.0, type = 'response')

babycat3 %>% 
  ggplot(aes(x = Time, y = Prop,
             colour = Condition, shape = Condition, linetype = Condition, fill = Condition)) +
  stat_summary(geom = 'ribbon', fun.data = 'mean_se', alpha = 0.3, colour = NA) +
  stat_summary(aes(y = model2), geom = 'line', fun = 'mean', size = 1) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = 'red') +
  labs(x = 'Time in trial (ms)',
       y = 'Prop looks to between colour') +
  scale_color_manual(values = c('blue', 'forestgreen', 'green')) +
  scale_fill_manual(values = c('blue', 'forestgreen', 'green')) +
  facet_wrap(~Know)

babycat3 %>% 
  ggplot(aes(x = Time, y = Prop,
             colour = Know, fill = Know)) +
  stat_summary(geom = 'ribbon', fun.data = 'mean_se', alpha = 0.3, colour = NA) +
  stat_summary(aes(y = model2), geom = 'line', fun = 'mean', size = 1) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = 'red') +
  labs(x = 'Time in trial (ms)',
       y = 'Prop looks to between colour') 


babycat3 %>% 
  filter(Know == 'Unknown') %>% 
  ggplot(aes(x = Time, y = Prop,
             colour = Condition, shape = Condition, linetype = Condition, fill = Condition)) +
  stat_summary(geom = 'ribbon', fun.data = 'mean_se', alpha = 0.3, colour = NA) +
  # stat_summary(aes(y = model2), geom = 'line', fun = 'mean', size = 1) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = 'red') +
  labs(x = 'Time in trial (ms)',
       y = 'Prop looks to between colour') +
  scale_color_manual(values = c('blue', 'forestgreen', 'green')) +
  scale_fill_manual(values = c('blue', 'forestgreen', 'green')) +
  facet_wrap(~Age)


###########################################################------------
# CATPAT

kitten <- read_csv('data/CatPat_combined.csv',na = c("", "NA"))
info <- read_csv('data/Info_Final.csv')

names(kitten)[c(2:3)] <- c('Block', 'Trial')
#check
kitty <- inner_join(kitten, info, by = c('ID', 'Block', 'Trial'))

kitty$Condition[kitty$Pic1 == '75BG56' & kitty$Pic2 == '25B56'] <- 'Within_Blue'
kitty$Condition[kitty$Pic1 == '25B56' & kitty$Pic2 == '75BG56'] <- 'Within_Blue'
kitty$Condition[kitty$Pic1 == '25BG56' & kitty$Pic2 == '75G56'] <- 'Within_Green'
kitty$Condition[kitty$Pic1 == '75G56' & kitty$Pic2 == '25BG56'] <- 'Within_Green'
kitty$Condition[kitty$Pic1 == '25B56' & kitty$Pic2 == '25BG56'] <- 'Between'
kitty$Condition[kitty$Pic1 == '25BG56' & kitty$Pic2 == '25B56'] <- 'Between'

kitty$Type <- ifelse(kitty$Condition == 'Between', 'Between', 'Within')

kitty$Utrial <- kitty$Trial
kitty$Utrial[kitty$Block == 2 & kitty$Trial == 1] <- 7
kitty$Utrial[kitty$Block == 2 & kitty$Trial == 2] <- 8
kitty$Utrial[kitty$Block == 2 & kitty$Trial == 3] <- 9
kitty$Utrial[kitty$Block == 2 & kitty$Trial == 4] <- 10
kitty$Utrial[kitty$Block == 2 & kitty$Trial == 5] <- 11
kitty$Utrial[kitty$Block == 2 & kitty$Trial == 6] <- 12

kitty$Target <- ifelse(kitty$`VIS W2 ISFXT` == -1 | kitty$`VIS W2 ISFXT` == 1 |
                         kitty$`VIS W1 ISFXT` == -1 | kitty$`VIS W1 ISFXT` == 1, 1, 0)

# kitty$Target <- ifelse(kitty$`VIS W2 ISFXT` == -1 | kitty$`VIS W2 ISFXT` == 1 |
#                        kitty$`VIS W1 ISFXT` == -1 | kitty$`VIS W1 ISFXT` == 1 |
#                          kitty$`FXT outside AOI`, 1, 0)

kitty$Target[is.na(kitty$Target)] <- 0
kitty$Other <- ifelse(kitty$Target == 0, 1, 0)
kitty$Age <- substr(kitty$ID, 1, 2)
kitty$Age[kitty$Age == 16] <- 12
names(kitty)[21] <- 'Pupil'

# left right
kitty2 <- kitty %>%  
  group_by(ID, Block, Utrial) %>% 
  filter(!is.na(`VIS W2 ISFXT`) | !is.na(`VIS W1 ISFXT`)) %>% 
  summarise(mean = mean(`Gaze filtered x`)) %>% 
  ungroup() %>% 
  mutate(Side = ifelse(mean < 1000, 'Left', 'Right'))

kitty3 <- left_join(kitty, kitty2)

#make short
catty <- kitty3 %>% 
  select('ID', 'Utrial', 'Block', 'Age', 'Timestamp', 'Condition',
                              'Type', 'Side', 'Target', 'Other', 'Pupil') %>% 
  mutate(ID = as.character(ID))

catty <- catty %>% 
  filter(substr(ID, 1,2) != 16) #pilots
################################################################
#CDI
################################################################

CDI <- read.csv('data/Participants_Shareable.csv')
#setup
CDI[CDI==""] <- 0
CDI[CDI=="U/S"] <- 2
CDI[CDI=="U"] <- 1
CDI[CDI=="No"] <- 0

#subset
cdi <- subset(CDI, CDI$CatPat == '1')
names(cdi)[1] <- 'ID'
names(cdi)[7] <- 'Age'

#make a knowledge for EITHER colour word
cdi$Knowledge[cdi$Blue == 1 | cdi$Green == 1] <- 1
cdi$Knowledge[is.na(cdi$Knowledge)] <- 0

#select
babscdi <- select(cdi, one_of('ID', 'Gender', 'AgeTest', 'Age',
                              'Knowledge'))
#reformat both
babscdi$ID <- as.character(babscdi$ID)
catty$ID <- as.character(catty$ID)
babscdi$Age <- as.character(babscdi$Age)

#merge
catdata <- left_join(catty, babscdi, by = c('ID', 'Age'))  %>% 
  filter(substr(ID, 1,2) != 16)

cc <- catdata %>% 
  group_by(ID,
           Age) %>% 
  summarise(Knowledge = first(Knowledge),
            AgeTest = first(AgeTest)) %>% 
  group_by(Age) %>% 
  summarize(N = length(ID),
            MAge = mean(AgeTest, na.rm = T),
            SDAge = sd(AgeTest, na.rm = T))

a1 <- anti_join(bbb, cc)
a2 <- anti_join(cc, bbb)
######################################################
# ETR
catdata <- catdata %>% 
  mutate(Target = ifelse(is.na(Target), 0, Target))

cat_track <- make_eyetrackingr_data(catdata,
                                    participant_column = 'ID',
                                    trackloss_column = 'Other',
                                    time_column = 'Timestamp',
                                    trial_column = 'Utrial',
                                    aoi_columns = c('Target'),
                                    treat_non_aoi_looks_as_missing = F) # key

cat_track$Target[is.na(cat_track$Target)] <- FALSE #removes perfect prop bug

trackloss <- trackloss_analysis(data = cat_track)
trackloss_subjects <- unique(trackloss[,c('ID','TracklossForParticipant')])

#remove files with trackloss
clean_cat <- clean_by_trackloss(data = cat_track,
                                trial_prop_thresh = 0.9) #as above changed from .99
#(removed 106 trials)

trackloss_clean <- trackloss_analysis(data = clean_cat)
trackloss_clean_subjects <- unique(trackloss_clean[, c('ID','TracklossForParticipant')])
#find out trackloss

trackloss_mean<-mean(1 - trackloss_clean_subjects$TracklossForParticipant)
trackloss_sd<-sd(1- trackloss_clean_subjects$TracklossForParticipant)
#See trials contributed by each participant

final_summary <- describe_data(clean_cat, 'Target', 'ID')
mean_num_trials<-mean(final_summary$NumTrials)
sd_sum_trials<-sd(final_summary$NumTrials)
#count
length(unique(final_summary$ID))
#57
length(unique(subset(clean_cat, Age == 19)$ID))
#29

#Ages
mean(subset(babscdi, Age == 19)$AgeTest) #19.27
sd(subset(babscdi, Age == 19)$AgeTest) #0.47
mean(subset(babscdi, Age == 12)$AgeTest) #12.24
sd(subset(babscdi, Age == 12)$AgeTest) #0.39

####################################################################
# Prop
catall <- make_time_window_data(clean_cat,
                                aois = 'Target',
                                predictor_columns = c(
                                                      'Knowledge', #'Side',
                                                      'Type',
                                                      'Age'),
                                summarize_by = 'ID')

catall <- catall %>% 
  mutate(Know = ifelse(Knowledge == 1, 'Known', 'Unknown'))

#plot(catall, predictor_columns = 'Condition')
plot(catall, predictor_columns = 'Type') +
  facet_wrap(~Age)

#final plot
source('misc/geom_flat_violin.R')

y <- catall %>% 
  ggplot(aes(x = Type, y = Prop, colour = Type, fill = Type)) +
  geom_flat_violin(position = ggplot2::position_nudge(x = .2, y = 0), alpha = .5, colour = NA) +
  ggplot2::geom_boxplot(width = .2,  outlier.shape = NA, alpha = 0.2) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = .15), size = .8, alpha = 0.8) +
  facet_wrap(~Age) +
  theme(legend.position = 'none') +
  labs(y = 'Prop looking', x = 'Trial type')


catall_mod <- catall %>% 
  mutate(Type_C = ifelse(Type == 'Between', 0.5, -0.5),
         Age_C = ifelse(Age == 19, 0.5, -0.5),
         Know_C = ifelse(Knowledge == 1, 0.5, -0.5))

catmod <- lmer(Prop ~ Type_C * Age_C + (1|ID), data = catall_mod)
summary(catmod)
car::Anova(catmod, type = 'III')
drop1(catmod, ~., test = 'Chisq')

############################################################
# knowledge

cat19 <- catall %>% 
#  filter(Age == '19') %>% 
  filter(!is.na(Knowledge)) %>% 
  mutate(Type_C = ifelse(Type == 'Between', 0.5, -0.5),
         Age_C = ifelse(Age == 19, 0.5, -0.5),
         Know_C = ifelse(Know == 'Known', 0.5, -0.5))

length(unique(catall$ID))
length(unique(cat19$ID))

plot(cat19, predictor_columns = 'Know')

t.test(Prop ~ Know, data = cat19)

#rerun without missing
model2 <- lmerTest::lmer(Prop ~ Type_C * Age_C +
                           (1|ID),
                         data = cat19)

model3 <- lmerTest::lmer(Prop ~ Type_C * Age_C * Know_C +
                 (1|ID),
               data = cat19)
summary(model3)
car::Anova(model3, type = 'III')

anova(model2, model3)

############################################################ß
#Pupil
pcat <- make_pupillometryr_data(catdata,
                                subject = ID,
                                trial = Utrial,
                                time = Timestamp,
                                condition = Type,
                                other = Knowledge)

clean_pcat <- clean_missing_data(pcat,
                                 pupil = Pupil,
                                 trial_threshold = .9,
                                 subject_trial_threshold = 1)

filter_cat <- filter_data(clean_pcat,
                          pupil = Pupil,
                          filter = 'hanning')

downsample_cat <- downsample_time_data(filter_cat,
                                       pupil = Pupil,
                                       timebin_size = 50,
                                       option = 'mean')

trim_cat <- subset_data(downsample_cat,
                        start = 1000,
                        stop = 6000,
                        rezero = T,
                        remove = T)

plot(trim_cat, pupil = Pupil, group = 'condition')

base_cat <- baseline_data(trim_cat,
                          pupil = Pupil,
                          start = 0,
                          stop = 50)

aw <- plot(base_cat, pupil = Pupil, group = 'condition')

##########################################
#test
win_cat <- create_window_data(base_cat, pupil = Pupil)
win_cat2 <- win_cat %>% 
  mutate(Age = substr(ID, 1, 2),
         Age = ifelse(Age == '16', '12', Age))

plot(win_cat2, pupil = Pupil, windows = F, geom = 'raincloud') +
  facet_wrap(~Age) + 
  labs(x = 'Trial type') 

z <- win_cat2 %>% 
  ggplot(aes(x = Type, y = Pupil, colour = Type, fill = Type)) +
  geom_flat_violin(position = ggplot2::position_nudge(x = .2, y = 0), alpha = .5, colour = NA) +
  ggplot2::geom_boxplot(width = .2,  outlier.shape = NA, alpha = 0.2) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = .15), size = .8, alpha = 0.8) +
  facet_wrap(~Age) +
  theme(legend.position = 'none') +
  labs(y = 'Change in pupil size', x = 'Trial type')

library(gridExtra)
library(cowplot)
library(ggpubr)

source('misc/Multiplot.R')
multiplot(y,z)

yz <- arrangeGrob(y,z)
as_ggplot(yz) +
  cowplot::draw_plot_label(label = c('A', 'B'),
                           x = c(0,  0),
                           y = c(1, 0.5))

########################################
#get differences
diff_cat <- create_difference_data(base_cat,
                       pupil = Pupil)

ax <- plot(diff_cat, pupil = Pupil, geom = 'line')

fun_cat <- create_functional_data(diff_cat,
                                  pupil = Pupil,
                                  basis = 10,
                                  order = 4)

ay <- plot(fun_cat, pupil = Pupil, geom = 'line', colour = 'blue')

#test
cattest <- run_functional_t_test(fun_cat,
                                 pupil = Pupil)

az <- plot(cattest, colour = 'red', fill = 'orange', show_divergence = T)

multiplot(aw, ax, ay, az, cols = 2)

alll <- arrangeGrob(aw, ax, ay, az) 

as_ggplot(alll) +
  cowplot::draw_plot_label(label = c('A', 'B', 'C', 'D'),
                           x = c(0, 0.5, 0, 0.5),
                           y = c(1, 1, 0.5, 0.5))
###############################################
#GAM
otherdata <- catdata %>% 
  group_by(ID, Utrial) %>% 
  summarise(Age = first(Age)) %>% 
            #Side = first(Side)) %>% 
  ungroup()

base_cat2 <- base_cat %>% 
  left_join(otherdata)

bamdata <- base_cat2 %>% 
 group_by(ID, Age, Type, Knowledge, Timestamp) %>%
 summarise(Pupil = mean(Pupil, na.rm = T)) %>%
 ungroup() %>% #cut across trials
  mutate(Type = as_factor(Type),
         Age = as_factor(Age),
         # Side = as_factor(Side),
         ID = as_factor(ID),
         # Trial = as_factor(as.character(Utrial)),
         Knowledge = as_factor(Knowledge),
         Type_s = ifelse(Type == 'Between', 0.5, -0.5),
         Age_s = ifelse(Age == '19', 0.5, -0.5),
         # Side_s = ifelse(Side == 'Left', -0.5, 0.5),
         Know_s = ifelse(Knowledge == 1, 0.5, -0.5))# %>%
  # filter(!is.na(Side))

library(mgcv)
library(itsadug)

bamdata2 <- data.frame(bamdata)
bamdata2$Event <- interaction(bamdata2$ID,  drop = T)
bamdata2 <- start_event(bamdata2, column = 'Timestamp', event = 'Event')



bamdata3 <- bamdata2 %>% 
  filter(!is.na(Knowledge)) %>% 
  mutate(Know = ifelse(Knowledge == 1, 'Known', 'Unknown'))
bamdata3$KnowType <- interaction(bamdata3$Knowledge:bamdata3$Type)
bamdata3$IDType <- interaction(bamdata3$ID:bamdata3$Type)
bamdata3$KT2 <- as.factor(bamdata3$KnowType)

bamdata4 <- bamdata3 %>% 
  filter(!is.na(Pupil))

bamdata4 %>% 
  ggplot(aes(x = Pupil)) +
  geom_density()
#scaled t

model_base <- bam(Pupil ~ 1 +
                    s(Timestamp, IDType, bs = 'fs', m = 1),
                  data = bamdata4,
                  family = 'scat',
                  discrete = T)

summary(model_base)

r1 <- start_value_rho(model_base)
r1 

model_base2 <- bam(Pupil ~ 1 +
                    s(Timestamp, IDType, bs = 'fs', m = 1),
                  data = bamdata4,
                  AR.start = bamdata4$start.event, rho = r1, #update to start with AR1
                  family = 'scat',
                  discrete = T)

summary(model_base2)

model_b1 <- bam(Pupil ~ 
                  s(Timestamp, by = Know_s, k = 11) +
                  s(Timestamp, IDType, bs = 'fs', m = 1),
                data = bamdata4,
                AR.start = bamdata4$start.event, rho = r1, 
                family = 'scat',
                discrete = T)



model_b2 <- bam(Pupil ~ 
                  s(Timestamp, by = Type_s, k = 11) +
                  s(Timestamp, by = Know_s, k = 11) +
                  s(Timestamp, IDType, bs = 'fs', m = 1),
                data = bamdata4,
                AR.start = bamdata4$start.event, rho = r1, 
                family = 'scat',
                discrete = T)


# model_b3 <- bam(Pupil ~ 
#                   s(Timestamp, by = Type_s, k = 11) +
#                   s(Timestamp, by = Know_s, k = 11) +
#                   s(Timestamp, by = KnowType, bs = c('cr'), k = 11) +
#                   s(Timestamp, IDType, bs = 'fs', m = 1),
#                 data = bamdata4,
#                 AR.start = bamdata4$start.event, rho = r1, 
#                 family = 'scat',
#                 discrete = T)

# model_b4 <- bam(Pupil ~ KnowType +
#                   s(Timestamp, by = Type_s, k = 11) +
#                   s(Timestamp, by = Know_s, k = 11) +
#                   s(Timestamp, by = KnowType, bs = c('cr'), k = 11) +
#                   s(Timestamp, IDType, bs = 'fs', m = 1),
#                 data = bamdata4,
#                 AR.start = bamdata4$start.event, rho = r1, 
#                 family = 'scat',
#                 discrete = T)

model_b5 <- bam(Pupil ~ Know_s * Type_s +
                  s(Timestamp, by = Type_s, k = 11) +
                  s(Timestamp, by = Know_s, k = 11) +
                  s(Timestamp, IDType, bs = 'fs', m = 1),
                data = bamdata4,
                AR.start = bamdata4$start.event, rho = r1, #updated
                #na.action = na.exclude,
                family = 'scat',
                discrete = T)

model_b6 <- bam(Pupil ~ Know_s * Type_s +
                  s(Timestamp, by = Type_s, k = 11) +
                  s(Timestamp, by = Know_s, k = 11) +
                  s(Timestamp, by = KnowType, k = 11) + #seems not to add much
                  s(Timestamp, IDType, bs = 'fs', m = 1),
                data = bamdata4,
                AR.start = bamdata4$start.event, rho = r1, #updated
                #na.action = na.exclude,
                family = 'scat',
                discrete = T)

summary(model_b1)
summary(model_b2)
summary(model_b3)
summary(model_b4)
summary(model_b5)

AIC(model_base)
AIC(model_base2)
AIC(model_b1)
AIC(model_b2)
AIC(model_b3)

acf_resid(model_base)

compareML(model_base2, model_b1)
compareML(model_b1, model_b2)
compareML(model_b3, model_b2) #model b3 is a red herring as the parametric is needed, and b4 does not require an interaction smooth
compareML(model_b5, model_b2)
compareML(model_b6, model_b5)


qqnorm(resid(model_base))
qqnorm(resid(model_b5))
acf_resid(model_b5)

save(model_base, file = 'models/model_base')
save(model_base2, file = 'models/model_base2')
save(model_b1, file = 'models/model_b1')
save(model_b2, file = 'models/model_b2')
save(model_b3, file = 'models/model_b3')
save(model_b4, file = 'models/model_b4')
save(model_b5, file = 'models/model_b5')


  bamdata4$Predict = fitted(model_b5)

  bamdata4 %>% 
    filter(!is.na(Knowledge)) %>% 
    ggplot(aes(x = Timestamp, y = Pupil, colour = Type, fill = Type)) +
    stat_summary(geom = 'ribbon', fun.data = 'mean_se', alpha = 0.5, colour = NA) +
    stat_summary(aes(y = Predict), geom = 'line', fun = 'mean', size = 1) +
    facet_wrap( ~ Know)

