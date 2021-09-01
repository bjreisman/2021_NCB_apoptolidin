

library(tidyverse)
list.files()

mydata <- read_csv("nsgs-pk.csv")


mydata %>%
  separate(Sample, c("date", "mouse", "time", "replicate")) %>%
  mutate(time = as.numeric(time)) %>%
  group_by(time, mouse) %>%
  mutate(conc = mean(Concentration), 
         sd = sd(Concentration), 
         conc.max = conc+sd*1.96, 
         conc.min = conc-sd*1.96) %>%
  ggplot(aes(x=time, y = Concentration, color = mouse, group = mouse)) + 
  geom_point(position = position_dodge(width = 30), size = 0.7) + 
  theme_bw() +
  geom_hline(yintercept = 50, linetype = 2) + 
  labs(x = "Time (min)", 
       y = "Concentration (nM)", 
       title = "IP Ammocidin A Pharmacokinetics", 
       subtitle = "Serum conc. after 0.5mg/kg IP dose") + 
  theme(legend.position = "bottom", 
        text = element_text(size = 8))
  
ggsave('ammocidinIPPK.pdf', width = 3, height = 2.75, units = "in", dpi = 300)


dodge <- position_dodge(width = 20)
mydata %>%
  separate(Sample, c("date", "mouse", "time", "replicate")) %>%
  mutate(time = as.numeric(time)) %>%
  group_by(time, mouse) %>%
  mutate(conc = mean(Concentration), 
         sd = sd(Concentration), 
         conc.max = conc+sd*1.96, 
         conc.min = conc-sd*1.96) %>%
  ggplot(aes(x=time, y = conc, color = mouse, group = mouse)) + 
  geom_point(position = dodge) + 
  geom_errorbar(aes(ymin = conc.min, ymax = conc.max), position = dodge) + 
  theme_bw() +
  geom_hline(yintercept = 50, linetype = 2) + 
  labs(x = "Time (min)", 
       y = "Concentration (nM)", 
       title = "IP Ammocidin A Pharmacokinetics", 
       subtitle = "Serum conc. after 0.5mg/kg IP dose") + 
  theme(legend.position = "bottom")

?geom_errorbar
mydata.summarized <- mydata %>%
  separate(Sample, c("date", "mouse", "time", "replicate")) %>%
  mutate(time = as.numeric(time)) %>%
  group_by(time, mouse) %>%
  summarise(conc = mean(Concentration), 
         sd = sd(Concentration), 
         conc.max = conc+sd*1.96, 
         conc.min = conc-sd*1.96)# %>%
  filter(time!= 0)


install.packages("PK")
library(PK)
#PKNews()
#?PK::biexp
conc1 <-mydata.summarized[-c(1:3),'conc']
time1 <- mydata.summarized[-c(1:3),'time']
mydata.biexp <- PK::biexp(conc1$conc, time = time1$time)

plot(mydata.biexp)

mydata.biexp
?auc

mydata.summarized.0 <- mydata.summarized %>%
  filter(time == 0) %>%
  rename(conc.0 = conc) %>%
  ungroup() %>%
  select(mouse, conc.0)

mydata.summarized.0
mydata.i <- mydata.summarized %>%
  left_join(mydata.summarized.0) %>%
  mutate(conc = conc-conc.0) %>%
  filter(mouse == "A") %>%
  ungroup() %>%
  mutate(time = time/60)

mydata.i 
auc.i <- auc(conc = mydata.i$conc,
    time = mydata.i$time,
 #   group = mydata.summarized$mouse,
    design = "complete")

summary(auc.i)
plot(auc.i)
auc.i

auc(conc, time, group=NULL, method=c("t", "z", "boott"),
    alternative=c("two.sided", "less", "greater"),
    conf.level=0.95, strata=NULL, nsample=1000,
    design=c("ssd","batch","complete"), data)
