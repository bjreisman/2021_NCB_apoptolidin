list.files()
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

##------------------------------------------------------------------------------
rep(c("A", "B"), each = 2)
str(rep(LETTERS[1:8], each = 3))
str(rep(6:8, times = 8))
rep(c(100e-6 / c(10^c(0:6)), 0), times = 3)
lookuptable <- 
  tibble(row = rep(LETTERS[1:8], times = 12), 
         col = as.character(rep(1:12, each = 8))) %>%
  #mutate(col = as.numeric(col)) %>%
    mutate(conc = rep(c(1000 / 2^c(0:7)), times = 12), 
           condition = rep(c("Vehicle",
                             "Vehicle",
                             "ApoA - 1.6 µM",
                             "ApoA - 1.6 µM",
                             "ApoA - 0.8 µM",
                             "ApoA - 0.8 µM",
                             "AmmoA - 0.8 µM",
                             "AmmoA - 0.8 µM",
                             "AmmoA - 0.4 µM",
                             "AmmoA - 0.4 µM",
                             "Uncat",
                             "Uncat"), each = 8))

lookuptable
extinction.coef <- 6220
##------------------------------------------------------------------------------

parse_plate <- function(data.preproc.i){
  data.preproc.i %>%
    mutate(`Time(hh:mm:ss)` = as.numeric(.[1,"Time(hh:mm:ss)"]/60)) %>%
    rename("time" = `Time(hh:mm:ss)`) %>%
    mutate(row = LETTERS[1:8]) %>%
    gather(col,abs, -time, -timepoint, -row) %>%
    dplyr::select(-timepoint) #%>%
  #  mutate(col =str_pad(col, 2, "left", "0"))
  
}

##------------------------------------------------------------------------------
list.files()
data.raw <- read_tsv("20191019_exp7.txt", skip = 2,
                     col_types = cols("Time(hh:mm:ss)" = col_time()))

data.raw


data.preproc <- 
  data.raw[c(1:nrow(data.raw) %% 9) != 0,1:14] %>% 
  dplyr::select(-2) %>%
  mutate(timepoint = c(0:c(n()-1)) %/%8) %>%
  slice(1:(n()-10))
data.preproc

data.tidy <-
  data.preproc %>% 
  dplyr::group_by(timepoint) %>%
  do(parse_plate(.)) %>%
  ungroup() %>%
  dplyr::select(-timepoint) %>%
  mutate(abs = as.numeric(abs)) %>%
  dplyr::select(row, col, abs, time) %>%
  inner_join(lookuptable)

data.tidy
data.tidy %>%
  ggplot(aes(x=time, y = abs)) + 
  geom_point() + 
  facet_grid(row~col) + 
  scale_y_continuous(limits = c(-1.6,0.3))

data.norm <- data.tidy %>%
  left_join(data.tidy %>%
              filter(time == 0) %>%
              dplyr::select(-time) %>%
              rename(abs.0 = "abs")) %>%
  mutate(abs.norm = abs - abs.0)

data.norm
data.norm %>%
  #mutate(fraction = as.factor(fraction)) %>%
  ggplot(aes(x=time, y = abs.norm)) +
  geom_smooth() + 
  geom_point() +
  labs(y= "OD340", x = "[ADP] (uM)")  + 
  theme_bw() + 
  scale_y_continuous(limits = c(-1.6,0.2)) + 
  facet_grid(~condition) 

ggsave('adp_titration_abs_raw.png', height = 3, width = 6, units = "in", dpi = 300)

data.norm %>%
  group_by(condition, conc, row, col) %>%
  arrange(time) %>%
  mutate(velocity = (abs.norm - lag(abs.norm, 3))) %>%
#  filter(time > 100) %>%
 # filter(time < 300) %>%
  ggplot(aes(x=time, y = velocity)) +
  geom_smooth() +
  geom_point(shape = ".") +
  geom_vline(xintercept = 150) +
  geom_vline(xintercept = 250) + 
  #  geom_smooth() + 
  labs(y= "Velocity\n(OD340/min)", x = "time(s)")  + 
  theme_bw() + 
#   scale_y_continuous(limits= c(-0.01,0.3)) + 
  facet_grid(condition~conc) 

#ggsave('adp_titration_velocity.png', width = 8, height = 4, units = "in", dpi = 300)


volume <- 200e-6 # in L
protein.ugmL <- 1.5
protein.mg <- protein.ugmL*volume
protein.mg
path.length <- 0.5
extinction.coef

data.mm <- data.norm %>%
  #mutate(fraction = as.factor(fraction)) %>%
  group_by(conc, condition, row, col) %>%
  arrange(time) %>%
  mutate(d.od = (abs.norm - lag(abs.norm, 6)), 
         d.t = (time - lag(time, 6)), 
         velocity = d.od/d.t*60) %>% #od/min
  filter(time > 150, time < 250) %>%
  group_by(conc, condition, row, col) %>%
  summarise(V.od = mean(velocity)) %>%
  mutate(V = V.od/extinction.coef*volume/protein.mg*1e9/path.length/375000*60)

data.mm
data.0 <- data.mm %>%
  ungroup() %>%
  filter(condition == "Uncat") %>%
  dplyr::select(conc, V) %>%
  group_by(conc) %>%
  summarise(V.0 = mean(V))


data.0
data.mm.norm <- data.mm %>%
  left_join(data.0) %>%
  mutate(V = V - V.0) %>%
  mutate(V = V*-1)
# data.mm <- data.norm %>%
#   #mutate(fraction = as.factor(fraction)) %>%
#   group_by(conc, condition, row, col) %>%
#   arrange(time) %>%
#   mutate(d.od = (abs.norm - lag(abs.norm, 6)), 
#          d.t = (time - lag(time, 6)), 
#          velocity = d.od/d.t*60) %>% #od/min
#   filter(time > 60, time < 400) %>%
#   group_by(conc, condition, row, col) %>%
#   summarise(V.od = mean(velocity)) %>%
#   mutate(V = V.od/extinction.coef*1e-6/protein) %>% #od/E = Moles/min
#   left_join(data.0) %>%
#   mutate(V = V - V.0)
#   
# data.mm$d.t
# data.norm
# 
# data.mm
library(drc)


library(forcats)
model.drm.kd <- drm (V ~ conc, condition,
                  data = data.mm.norm %>%
                    filter(condition != "Uncat") %>%
                    ungroup() %>%
                    mutate(condition = as.factor(condition)) %>%
                    mutate(condition = fct_relevel(condition, "Vehicle")),
                  pmodels = data.frame(condition, 1),
                  fct = MM.2(names = c("Vmax", "Kd")))

plot(model.drm.kd)
model.drm.Vmax <- drm (V ~ conc, condition,
                     data = data.mm.norm %>%
                       filter(condition != "ApoA - 3 ÂµM") %>%
                       ungroup() %>%
                       mutate(condition = as.factor(condition)) %>%
                       mutate(condition = fct_relevel(condition, "Vehicle")),
                     
                     pmodels = data.frame(1, condition),
                     fct = MM.2(names = c("Vmax", "Kd")))
data.mm.norm
model.drm <- drm (V ~ conc, condition,
                       data = data.mm.norm %>%
                         filter(condition != "Uncat") %>%
                         ungroup() %>%
                         mutate(condition = as.factor(condition)) %>%
                         mutate(condition = fct_relevel(condition, "Vehicle")),
                       pmodels = data.frame(condition, condition),
                       fct = MM.2(names = c("Vmax", "Kd")))

model.drm.2 <- drm (V ~ conc, condition,
                  data = data.mm.norm %>%
                    filter(condition != "Uncat") %>%
                    ungroup() %>%
                    mutate(condition = as.factor(condition)) %>%
                    mutate(condition = fct_relevel(condition, "Vehicle")),
                  fct = MM.2(names = c("Vmax", "Kd")))

summary(model.drm.2)
modelFit(model.drm.2)
dev.off()

anova(model.drm, model.drm.kd) 

####----------------------------------------------------------------------------
dev.off()
pdf("MMplot_winsert.pdf", width = 6, height = 5)
#dev.off()
par(mar = c(5,5,4,4)) #bltr

uconds <- unique(model.drm$origData$condition)
par(fig = c(0,1,0,1))

ymax <-max(model.drm$data$V)*2

c('#bcbddc', '#756bb1', '#b2df8a', '#33a02c', '#000000')
#?plot.drc

insert.data.points <- as_tibble(model.drm$data[-3]) %>%
  mutate(conc.r = 1/conc, V.r = 1/V) %>%
  filter(conc != 0)

plot(model.drm,
     type = "all",
     log = "",
     pch = 1,
     xlab = "[ATP] (ÂµM)", 
     ylab =  expression("nMoles ATP sec"^-1),
     col = c('#bcbddc', '#756bb1', '#b2df8a', '#33a02c', '#000000'), 
     lty = 1, 
     ylim = c(0,ymax), 
     legend = F
     )


legend(x = 0, y = ymax,
      legend = levels(insert.data.points$orig.condition), 
      col = mypallete, pch = c(2, 1), lty = c(1,1), 
      bty = 'n')

mypallete <- c('#000000', '#bcbddc', '#756bb1', '#b2df8a', '#33a02c') 
names(mypallete) <-  levels(insert.data.points$orig.condition)

coefs <- coef(model.drm)
coefs
tibble(condition = substr(names(coefs[1:5]), 6, 100),
       Vmax = signif(coefs[1:5],3), 
       Kd = signif(coefs[6:10],3)) %>%
  mutate(Vmax_Kd = signif(Vmax/Kd,3)) %>%
  as.data.frame()
#km.y <- predict(model.drm, data.frame(conc = unname(coefs[c(3,4)]),
 #                             condition = c("ApoA - 1 ÂµM ", "Vehicle")))

library(tibble)
ki.data <- data.frame(Vmax = coefs[1:5],
                      Km = coefs[6:10]) %>%
  rownames_to_column('condition') %>%
  mutate(condition = substr(condition, 6, 100), 
         Vmax.r = 1/Vmax)

ki.data %>%
  mutate(Vmax.Km = Vmax/Km)
#ki.data[3,]
segments(ki.data[1,3], -20, y1 = ki.data[1,2]/2, lty = 2,
         col = "#000000")
segments(ki.data[2,3], -20, y1 = ki.data[2,2]/2, lty = 2,
         col = "#bcbddc")
segments(ki.data[3,3], -20, y1 = ki.data[3,2]/2, lty = 2,
         col = "#756bb1")
segments(ki.data[4,3], -20, y1 = ki.data[4,2]/2, lty = 2,
         col = "#b2df8a")
segments(ki.data[5,3], -20, y1 = ki.data[5,2]/2, lty = 2,
         col = "#33a02c")

segments(x0 = -20,
         y0 = ki.data[1,2]/2,
         x1 =ki.data[1,3],
         y1 = ki.data[1,2]/2, lty = 2,
         col = "#000000")

segments(x0 = -20,
         y0 = ki.data[2,2]/2,
         x1 =ki.data[2,3],
         y1 = ki.data[2,2]/2, lty = 2,
         col = "#bcbddc")
segments(x0 = -20,
         y0 = ki.data[3,2]/2,
         x1 =ki.data[3,3],
         y1 = ki.data[3,2]/2, lty = 2,
         col = "#756bb1")
segments(x0 = -20,
         y0 = ki.data[4,2]/2,
         x1 =ki.data[4,3],
         y1 = ki.data[4,2]/2, lty = 2,
         col = "#b2df8a")

segments(x0 = -20,
         y0 = ki.data[5,2]/2,
         x1 =ki.data[5,3],
         y1 = ki.data[5,2]/2, lty = 2,
         col = "#33a02c")


#dev.off()
#insert.data.points
#?par
#x1, x2, y1, y2
#dev.off()
par(fig = c(0.4,0.95, 0.35, 0.95), new = T)  
#insert.data.points$condition
#insert.data.points$condition
#insert.data.points$condition

insert.data.points$orig.condition
table(mypallete[insert.data.points$orig.condition])
insert.data.points$V.r
plot(V.r ~ conc.r, data = insert.data.points,
     col = mypallete[insert.data.points$orig.condition],
     #type = "p", 
#     col = c('#bcbddc', '#756bb1', '#b2df8a', '#33a02c', '#000000'), 
     ylab = expression("V"^-1), 
     xlab = expression("[ATP]"^-1), 
     pch = 1, 
     cex = 0.5,
     axes = F, 
     frame.plot = T, 
     mgp = c(0.3,1,1), 
     xlim = c(-0.03, 0.14), 
     ylim = c(-2, 600)
    )
#
#ki.data

#?plot
abline(v = 0, lwd = 2, col = "grey50")
abline(h = 0, lwd = 2, col = "grey50")
abline(a = 1/ki.data[1,2], b= ki.data[1,3]/ki.data[1,2], col = "#000000")
abline(a = 1/ki.data[2,2], b= ki.data[2,3]/ki.data[2,2], col = "#bcbddc")
abline(a = 1/ki.data[3,2], b= ki.data[3,3]/ki.data[3,2], col = "#756bb1")
abline(a = 1/ki.data[4,2], b= ki.data[4,3]/ki.data[4,2], col = "#b2df8a")
abline(a = 1/ki.data[5,2], b= ki.data[5,3]/ki.data[5,2], col = "#33a02c")

dev.off()

####----------------------------------------------------------------------------
png("MMplot_winsert_noncompetitive.png", width = 6, height = 5, units = "in", res = 300)
#dev.off()
par(mar = c(5,5,4,4)) #bltr

uconds <- unique(model.drm.kd$origData$condition)
par(fig = c(0,1,0,1))

ymax <-max(model.drm$data$V)*2
plot(model.drm.kd,
     log = "",
     xlab = "[ATP] (ÂµM)", 
     ylab =  expression("nMoles ATP sec"^-1),
     col = c('#bcbddc', '#756bb1', '#b2df8a', '#33a02c', '#000000'), 
     lty = 1, 
     ylim = c(0,ymax), 
     legend = F)


legend(x = 0, y = ymax,
       legend = levels(insert.data.points$condition), 
       col = mypallete, pch = c(2, 1), lty = c(1,1), 
       bty = 'n')

mypallete <- c('#000000', '#bcbddc', '#756bb1', '#b2df8a', '#33a02c') 
names(mypallete) <-  levels(insert.data.points$condition)

coefs <- coef(model.drm.kd)
#km.y <- predict(model.drm, data.frame(conc = unname(coefs[c(3,4)]),
#                             condition = c("ApoA - 1 ÂµM ", "Vehicle")))

library(tibble)
ki.data <- data.frame(Vmax = coefs[1:5],
                      Km = coefs[6:10]) %>%
  rownames_to_column('condition') %>%
  mutate(condition = substr(condition, 6, 100), 
         Vmax.r = 1/Vmax)
ki.data
segments(ki.data[1,3], -20, y1 = ki.data[1,2]/2, lty = 2,
         col = "#000000")
segments(ki.data[1,3], -20, y1 = ki.data[2,2]/2, lty = 2,
         col = "#bcbddc")
segments(ki.data[1,3], -20, y1 = ki.data[3,2]/2, lty = 2,
         col = "#756bb1")
segments(ki.data[1,3], -20, y1 = ki.data[4,2]/2, lty = 2,
         col = "#b2df8a")
segments(ki.data[1,3], -20, y1 = ki.data[5,2]/2, lty = 2,
         col = "#33a02c")

segments(x0 = -20,
         y0 = ki.data[1,2]/2,
         x1 =ki.data[1,3],
         y1 = ki.data[1,2]/2, lty = 2,
         col = "#000000")

segments(x0 = -20,
         y0 = ki.data[2,2]/2,
         x1 =ki.data[1,3],
         y1 = ki.data[2,2]/2, lty = 2,
         col = "#bcbddc")
segments(x0 = -20,
         y0 = ki.data[3,2]/2,
         x1 =ki.data[1,3],
         y1 = ki.data[3,2]/2, lty = 2,
         col = "#756bb1")
segments(x0 = -20,
         y0 = ki.data[4,2]/2,
         x1 =ki.data[1,3],
         y1 = ki.data[4,2]/2, lty = 2,
         col = "#b2df8a")

segments(x0 = -20,
         y0 = ki.data[5,2]/2,
         x1 =ki.data[1,3],
         y1 = ki.data[5,2]/2, lty = 2,
         col = "#33a02c")

insert.data.points <- as_tibble(model.drm$data[-3]) %>%
  mutate(conc.r = 1/conc, V.r = 1/V) %>%
  filter(conc != 0)

#dev.off()
#insert.data.points
#?par
#x1, x2, y1, y2
#dev.off()
par(fig = c(0.4,0.95, 0.35, 0.95), new = T)  
#insert.data.points$condition
#insert.data.points$condition
#insert.data.points$condition


plot(V.r ~ conc.r, data = insert.data.points,
     col = mypallete[insert.data.points$condition],
     #type = "p", 
     # col = c('#bcbddc', '#756bb1', '#b2df8a', '#33a02c', '#000000'), 
     ylab = expression("V"^-1), 
     xlab = expression("[ATP]"^-1), 
     pch = 16, 
     axes = F, 
     frame.plot = T, 
     mgp = c(0.3,1,1), 
     xlim = c(-0.03, 0.08), 
     ylim = c(-2, 15))
#
#ki.data

#?plot
abline(v = 0, lwd = 2, col = "grey50")
abline(h = 0, lwd = 2, col = "grey50")
#ki.data[2,2]
abline(a = 1/ki.data[1,2], b= ki.data[1,3]/ki.data[1,2], col = "#000000")
abline(a = 1/ki.data[2,2], b= ki.data[1,3]/ki.data[2,2], col = "#bcbddc")
abline(a = 1/ki.data[3,2], b= ki.data[1,3]/ki.data[3,2], col = "#756bb1")
abline(a = 1/ki.data[4,2], b= ki.data[1,3]/ki.data[4,2], col = "#b2df8a")
abline(a = 1/ki.data[5,2], b= ki.data[1,3]/ki.data[5,2], col = "#33a02c")

dev.off()

summary(model.drm)
summary(model.drm.kd)
summary(model.drm.Vmax)

modelFit(model.drm.kd)
modelFit(model.drm)
modelFit(model.drm.Vmax)

## Comparing the four-parameter log-logistic model 
##  to a one-way ANOVA model using an approximate F test
## in other words applying a lack-of-fit test
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
modelFit(ryegrass.m1)   

summary(model.drm)

?MM.2()
?plot.drc
#?MM.3
mml <- data.frame(conc = seq(0, max(data.mm$conc), length.out = 100))


mml$V <- predict(model.drm, newdata = mml)


?plot.drc

data.mm %>%
  ggplot(aes(y = V, x = conc, col = condition)) + 
  geom_point() + 
  theme_bw() + 
  #geom_line(data = mml, aes(y = V, x = conc, col = condition)) +
  #scale_x_log10() + 
  #geom_smooth(alpha = 0.2, se = F, formula = y~ x/(1 + x), method = "lm") + 
  scale_color_brewer(type = "qual", palette = 6)

 ggsave('adp_mm.png', width = 5, height = 3, units = "in", dpi = 300)

 
hist(rbinom(100, 10, 0.1))

ppois(0, 0.5)
(ppois(1, 0.2) - ppois(0, 0.2)) / (1 - ppois(0, 0.2))

(ppois(1, 0.2) - ppois(0, 0.2))
table(rpois(60, 1.875/5))
ppois(1, 1.875/5)
1.875/5

u.conds <- unique(data.mm$condition)
lb.mod.1 <- lm(V.r ~ conc.r,
             data = data.mm %>%
               filter(condition ==u.conds[1]) %>%
               filter(conc != 0) %>%
               mutate(V.r = 1/V, conc.r = 1/conc))
lb.mod.2 <- lm(V.r ~ conc.r,
               data = data.mm %>%
                 filter(condition ==u.conds[2]) %>%
                 filter(conc != 0) %>%
                 mutate(V.r = 1/V, conc.r = 1/conc))
lb.mod.3 <- lm(V.r ~ conc.r,
               data = data.mm %>%
                 filter(condition ==u.conds[3]) %>%
                 filter(conc != 0) %>%
                 mutate(V.r = 1/V, conc.r = 1/conc))



lb.x1 <- chemCal::inverse.predict(lb.mod.1, 0)
lb.x2 <- chemCal::inverse.predict(lb.mod.2, 0)
lb.x3 <- chemCal::inverse.predict(lb.mod.3, 0)

km <- 1/(-lb.x$Prediction)
Vmax <- 1/lb.mod$coefficients[1]

data.mm

data.mm %>%
  filter(conc > 0) %>%
  ggplot(aes(y = 1 / V, x = 1 / conc, col = condition)) +
  geom_point() +
  geom_abline(slope = lb.mod.1$coefficients[2],
              intercept = lb.mod.1$coefficients[1]) +
  geom_abline(slope = lb.mod.2$coefficients[2],
              intercept = lb.mod.2$coefficients[1]) +
  geom_abline(slope = lb.mod.3$coefficients[2],
              intercept = lb.mod.3$coefficients[1]) +
  geom_hline(yintercept = lb.mod$coefficients[1], linetype = 2) +
  geom_hline(yintercept = 0, size = 1)  + 
  geom_vline(xintercept = 0, size = 1)  + 
  geom_vline(xintercept = lb.x$Prediction, linetype = 2) + 
  theme_bw()

ggsave('adp_lb.png', width = 5, height = 5, units = "in", dpi = 300)
 
 
compound.i <-  unique(data.norm$compound)[[1]]

crc.data <-  data.norm %>%
   group_by(row, col, conc, compound) %>%
   arrange(time) %>%
   mutate(velocity = (abs.norm - lag(abs.norm, 5))) %>%
   filter(time > 300) %>%
   group_by(row, col, conc, compound) %>%
   summarise(Vmax.od = mean(velocity)) %>%
   mutate(Vmax = Vmax.od/extinction.coef*2*1e9*0.0002/c(50*0.2)*1000)
 
library(drc)
crc.data
?drm
crc.model <- drm(Vmax ~ conc, compound, data = crc.data, fct = LL.4())
png("yeastcrc.png", width = 6, height = 5, units = "in", res = 300)
plot(crc.model, broken = T,
     xlab = "[Compound]", ylab = "nMoles ATP / mg Protein / min", 
     main = "HK/G6PD Coupled ATP Synthase Assay\nin Isolated Yeast (USY006) Mitochondria", 
     col = TRUE, 
     legendPos = c(1e-8, 300))

dev.off()

ic50 <- drc::ED(crc.model, 50)
-
ic50[[1]]*1e9

data.norm %>%
  group_by(conc) %>%
  arrange(time) %>%
  mutate(velocity = (abs.norm - lag(abs.norm, 5))) %>%
  filter(time > 180) %>%
  group_by(conc, row) %>%
  summarise(Vmax.od = mean(velocity)) %>%
  mutate(Vmax = Vmax.od/extinction.coef*2*1e9*0.0002) %>%
  filter(row != "E") %>%
  filter(conc != 0) %>%
  ggplot(aes(y = 1/Vmax, x = 1/conc)) + 
  geom_point() + 
  theme_bw() + 
  geom_smooth(alpha = 0.2, se = F, method = "lm") + 
  scale_color_brewer(type = "qual", palette = 6) + 
  labs(title = "lineweaver-burk")

ggsave('adp_lineweaverburk.png', width = 5, height = 5, units = "in", dpi = 300)


data.norm %>%
  group_by(conc) %>%
  arrange(time) %>%
  mutate(velocity = (abs.norm - lag(abs.norm, 5))) %>%
  filter(time > 180) %>%
  group_by(conc, row) %>%
  summarise(Vmax.od = mean(velocity)) %>%
  mutate(Vmax = Vmax.od/extinction.coef*2*1e9*0.0002) %>%
  filter(row != "E") %>%
  filter(conc != 0) %>%
  ggplot(aes(y = conc/Vmax, x = conc)) + 
  geom_point() + 
  theme_bw() + 
  geom_smooth(alpha = 0.2, se = F, method = "lm") + 
  scale_color_brewer(type = "qual", palette = 6) + 
  labs(title = "HanesâWoolf plot")

ggsave('adp_hanes-woolf.png', width = 5, height = 5, units = "in", dpi = 300)

