library(tidyverse)
list.files(pattern = "csv")

mydata <- read_csv("chimerism_day28_rep1.csv") %>%
  gather(compartment, chimerism, Blood, BM, SPL)

library(ggpubr)
mydata %>%
  ggplot(aes(x=Group, y = chimerism, group = Group)) + 
  geom_point() + 
  facet_grid(.~compartment)

mypal0 <- c("black","#bcbddc", "#756bb1")
mypal <- c("black","#edb5d3", "#d01d8a")

mydata
#?geom_errorbar
fig5b <- mydata %>%
  group_by(Group, compartment) %>%
  mutate(mean = mean(chimerism), 
         sem = sqrt(var(chimerism)/n()), 
         err.min = mean - sem, 
         err.max = mean + sem) %>%
  mutate(Group = as.factor(Group)) %>%
  mutate(Group = fct_relevel(Group, "Vehicle")) %>%
  ggplot(aes(y = chimerism, fill = Group, x = Group, color = Group)) + 
  geom_errorbar(aes(ymin = err.min, ymax = err.max, y = mean), width = 0.5, 
                color = "black") + 
  geom_errorbar(aes(ymin = mean, ymax = mean, y = mean), width = 1, 
                color = "black") + 
  geom_point(position = position_jitter(width = 0.2), 
             size  = 1) + 
  theme_classic()  + 
  geom_signif(map_signif_level = T,
              comparisons = list(c("Vehicle", "AmmoA_0.1"), 
                                 c("Vehicle", "AmmoA_0.03"), 
                                 c("AmmoA_0.1", "AmmoA_0.03")
                                 ), 
              y_position = c(80, 75, 70),
              textsize = 2, 
              test = "wilcox.test", 
              tip_length = 0) + 
  scale_color_manual(values = mypal) +   
  scale_fill_manual(values = mypal) +   
  facet_grid(.~compartment, switch  = "x") + 
  theme(legend.position = "bottom", axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        strip.placement = "outside", 
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(), 
        text = element_text(size = 7))


ggsave("rep1_d28.pdf", width = 4, height = 3, units = "in")


#fig5b <- 
  mydata %>%
  group_by(Group, compartment) %>%
  mutate(mean = mean(chimerism), 
         sem = sqrt(var(chimerism)/n()), 
         err.min = mean - sem, 
         err.max = mean + sem) %>%
  mutate(Group = as.factor(Group)) %>%
  mutate(Group = fct_relevel(Group, "Vehicle")) %>%
  ggplot(aes(y = chimerism, fill = Group, x = Group, color = Group)) + 
  geom_errorbar(aes(ymin = err.min, ymax = err.max, y = mean), width = 0.5, 
                color = "black") + 
  geom_errorbar(aes(ymin = mean, ymax = mean, y = mean), width = 1, 
                color = "black") + 
  geom_point(position = position_jitter(width = 0.2), 
             size  = 1) + 
  theme_classic()  + 
  stat_compare_means(
                comparisons = list(c("Vehicle", "AmmoA_0.1"),
                                   c("Vehicle", "AmmoA_0.03"),
                                   c("AmmoA_0.1", "AmmoA_0.03")), 
                hide.ns = TRUE, 
                method = "wilcox.test"
  ) + 
  # geom_signif(map_signif_level = T,
  #             comparisons = list(c("Vehicle", "AmmoA_0.1"), 
  #                                c("Vehicle", "AmmoA_0.03"), 
  #                                c("AmmoA_0.1", "AmmoA_0.03")
  #             ), 
  #             y_position = c(80, 75, 70),
  #             textsize = 2, 
  #             test = "wilcox.test", 
  #             tip_length = 0) + 
  scale_color_manual(values = mypal) +   
  scale_fill_manual(values = mypal) +   
  facet_grid(.~compartment, switch  = "x") + 
  theme(legend.position = "bottom", axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        strip.placement = "outside", 
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(), 
        text = element_text(size = 7))

library(rstatix)
mydata %>%
  group_by(Group, compartment) %>%
  mutate(mean = mean(chimerism), 
         sem = sqrt(var(chimerism)/n()), 
         err.min = mean - sem, 
         err.max = mean + sem) %>%
  mutate(Group = as.factor(Group)) %>%
  mutate(Group = fct_relevel(Group, "Vehicle")) %>%
  group_by(compartment) %>%
  wilcox_test(chimerism ~ Group,  comparisons = list(c("Vehicle", "AmmoA_0.1"),
                                                      c("Vehicle", "AmmoA_0.03"),
                                                      c("AmmoA_0.1", "AmmoA_0.03")), 
              p.adjust.method = "holm", 
              alternative = "two.sided") %>%
  write_csv("wilcox.csv")
