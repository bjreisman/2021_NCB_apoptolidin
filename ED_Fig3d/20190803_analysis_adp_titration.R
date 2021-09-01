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
  tibble(row = rep(LETTERS[1:8], times = 3), 
         col = as.character(rep(7:9, each = 8))) %>%
    mutate(conc = rep(c(400 / 2^c(0:6),0), times = 3), 
           condition = rep(c("Vehicle", "ApoA - 3 µM", "ApoA - 1 µM "), each = 8))

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

data.raw <- read_tsv("20190809_apoA_adptitration.txt" , skip = 2,
                     col_types = cols("Time(hh:mm:ss)" = col_time()))

data.raw


data.preproc <- 
  data.raw[c(1:nrow(data.raw) %% 9) != 0,1:14] %>% 
  dplyr::select(-2) %>%
  mutate(timepoint = c(0:c(n()-1)) %/%8) %>%
  slice(1:(n()-2))
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
  scale_y_continuous(limits = c(0,2))

data.norm <- data.tidy %>%
  left_join(data.tidy %>%
              filter(time == 0) %>%
              dplyr::select(-time) %>%
              rename(abs.0 = "abs")) %>%
  mutate(abs.norm = abs - abs.0)

data.norm %>%
  #mutate(fraction = as.factor(fraction)) %>%
  ggplot(aes(x=time, y = abs.norm)) +
  geom_smooth() + 
  geom_point() +
  labs(y= "OD340", x = "[ADP] (uM)")  + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,0.75)) + 
  facet_grid(~condition) 

ggsave('adp_titration_abs_raw.png', height = 3, width = 6, units = "in", dpi = 300)

data.norm %>%
  group_by(condition, conc) %>%
  arrange(time) %>%
  mutate(velocity = (abs.norm - lag(abs.norm, 3))) %>%
  filter(time > 40) %>%
  ggplot(aes(x=time, y = velocity)) +
  geom_smooth() +
  geom_point(shape = ".") +
  #  geom_smooth() + 
  labs(y= "Velocity\n(OD340/min)", x = "time(s)")  + 
  theme_bw() + 
#   scale_y_continuous(limits= c(-0.01,0.3)) + 
  facet_grid(condition~conc) 

#ggsave('adp_titration_velocity.png', width = 8, height = 4, units = "in", dpi = 300)


volume <- 200e-6 # in L
protein.ugmL <- 50
protein.mg <- protein.ugmL*volume
protein <- protein.mg
path.length <- 0.5
extinction.coef

data.mm <- data.norm %>%
  #mutate(fraction = as.factor(fraction)) %>%
  group_by(conc, condition, row, col) %>%
  arrange(time) %>%
  mutate(d.od = (abs.norm - lag(abs.norm, 6)), 
         d.t = (time - lag(time, 6)), 
         velocity = d.od/d.t*60) %>% #od/min
  filter(time > 60, time < 400) %>%
  group_by(conc, condition, row, col) %>%
  summarise(V.od = mean(velocity)) %>%
  mutate(V = V.od/extinction.coef*volume/protein*1e9/path.length)

data.0 <- data.mm %>%
  ungroup() %>%
  filter(conc == 0) %>%
  dplyr::select(-conc, - row, -col, - V.od) %>%
  rename(V.0 = V)


data.mm.norm <- data.mm %>%
  left_join(data.0) %>%
  mutate(V = V - V.0)
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
                    filter(condition != "ApoA - 3 µM") %>%
                    ungroup() %>%
                    mutate(condition = as.factor(condition)) %>%
                    mutate(condition = fct_relevel(condition, "Vehicle")),
                  
                  pmodels = data.frame(condition, 1),
                  fct = MM.2(names = c("Vmax", "Kd")))

model.drm.Vmax <- drm (V ~ conc, condition,
                     data = data.mm.norm %>%
                       filter(condition != "ApoA - 3 µM") %>%
                       ungroup() %>%
                       mutate(condition = as.factor(condition)) %>%
                       mutate(condition = fct_relevel(condition, "Vehicle")),
                     
                     pmodels = data.frame(1, condition),
                     fct = MM.2(names = c("Vmax", "Kd")))

model.drm <- drm (V ~ conc, condition,
                       data = data.mm.norm %>%
                         filter(condition != "ApoA - 3 µM") %>%
                         ungroup() %>%
                         mutate(condition = as.factor(condition)) %>%
                         mutate(condition = fct_relevel(condition, "Vehicle")),
                       
                       pmodels = data.frame(condition, condition),
                       fct = MM.2(names = c("Vmax", "Kd")))
summary(model.drm)
dev.off()

anova(model.drm, model.drm.kd) 

####----------------------------------------------------------------------------
dev.off()
pdf("MMplot_winsert_uncompetitive.pdf", width = 6, height = 5)
par(mar = c(5,5,4,4)) #bltr

uconds <- unique(model.drm$origData$condition)
par(fig = c(0,1,0,1))

plot(model.drm,
     log = "",
     xlab = "[ADP] (�M)", 
     ylab =  expression("nMoles ATP min"^-1*"mg"^-1), 
     main = "Apoptolidin A is an uncompetitive inhibitor vs. ADP?", 
     col = c("red3", "black"),
     lty = c(1,1), 
     legend = F, 
     ylim = c(0,550))

legend(x = 0, y = 550, legend = c("Vehicle", "Apoptolidin A - 1µM"), 
      col = c("black", "red3"), pch = c(2, 1), lty = c(1,1), 
      bty = 'n')

coefs <- coef(model.drm)
#km.y <- predict(model.drm, data.frame(conc = unname(coefs[c(3,4)]),
 #                             condition = c("ApoA - 1 µM ", "Vehicle")))

ki.data <- data.frame(Vmax = coefs[1:2], 
           conc = c(0,1)) %>%
  mutate(Vmax.r = 1/Vmax) 

# ki.data %>%
#   ggplot(aes(x=conc, y = Vmax.r)) + 
#   geom_point() + 
#   geom_abline(slope = coefficients(ki.model)[2], 
#               intercept = coefficients(ki.model)[1]) + 
#   scale_x_continuous(limits =c(-2, 2)) + 
#   scale_y_continuous(limits = c(-0.001, 0.005))

ki.model <- lm(Vmax.r ~ conc + 1, data = ki.data)
coefficients(ki.model)[1]/coefficients(ki.model)[2]
segments(coefs[3], -20, y1 = coefs[1]/2, lty = 2,
         col = adjustcolor("black", 1))
segments(coefs[4], -20, y1 = coefs[2]/2, lty = 2,
         col = adjustcolor("red", 1))
segments(x0 = -20, y0 = coefs[1]/2, x1 =coefs[3], y1 = coefs[1]/2, lty = 2,
         col = adjustcolor("black", 1))
segments(x0 = -20, y0 = coefs[2]/2, x1 =coefs[4], y1 = coefs[2]/2, lty = 2,
         col = adjustcolor("red", 1))

insert.data.points <- as_tibble(model.drm$data[-3]) %>%
  mutate(conc.r = 1/conc, V.r = 1/V) %>%
  filter(conc != 0)


par(fig = c(0.4,0.95, 0.05, 0.65), new = T)  

plot(V.r ~ conc.r, data = insert.data.points, type = "p", 
     col = insert.data.points$orig.condition, 
     ylab = expression("V"^-1), 
     xlab = expression("[ADP]"^-1), 
     axes = F, 
     frame.plot = T, 
     mgp = c(0.3,1,1), 
     xlim = c(-0.03, 0.2), 
     ylim = c(-0.001, 0.02))
abline(v = 0, lwd = 2, col = "grey50")
abline(h = 0, lwd = 2, col = "grey50")

abline(a = 1/coefs[[1]], b= coefs[[3]]/coefs[[1]], col = "black")
abline(a = 1/coefs[[2]], b= coefs[[4]]/coefs[[2]], col = "red")


x <- rnorm(100) 
y <- rbinom(100, 1, 0.5)

dev.off()

####----------------------------------------------------------------------------
png("MMplot_winsert_noncompetitive.png", width = 6, height = 5, units = "in", res = 300)
par(mar = c(5,5,4,4)) #bltr

uconds <- unique(model.drm$origData$condition)
par(fig = c(0,1,0,1))

plot(model.drm.kd,
     log = "",
     xlab = "[ADP] (µM)", 
     ylab =  expression("nMoles ATP min"^-1*"mg"^-1), 
     main = "Apoptolidin A is an noncompetitive inhibitor vs. ADP?", 
     col = c("red3", "black"),
     lty = c(1,1), 
     legend = F, 
     ylim = c(0,550))

legend(x = 0, y = 550, legend = c("Vehicle", "Apoptolidin A - 1µM"), 
       col = c("black", "red3"), pch = c(2, 1), lty = c(1,1), 
       bty = 'n')

coefs <- coef(model.drm.kd)
km.y <- predict(model.drm, data.frame(conc = unname(coefs[c(3,3)]),
                                      condition = c("ApoA - 1 µM ", "Vehicle")))

segments(coefs[3], -20, y1 = km.y[1], lty = 2,
         col = adjustcolor("red3", 1))

segments(coefs[3], -20, y1 = km.y[2], lty = 2,
         col = adjustcolor("black", 1))
segments(x0 = -20, y0 = km.y[1], x1 =coefs[3], y1 = km.y[1], lty = 2,
         col = adjustcolor("red3", 1))
segments(x0 = -20, y0 = km.y[2], x1 =coefs[3], y1 = km.y[2], lty = 2,
         col = adjustcolor("black", 1))

insert.data.points <- as_tibble(model.drm.kd$data[-3]) %>%
  mutate(conc.r = 1/conc, V.r = 1/V) %>%
  filter(conc != 0)


par(fig = c(0.4,0.95, 0.05, 0.65), new = T)  

plot(V.r ~ conc.r, data = insert.data.points, type = "p", 
     col = insert.data.points$condition, 
     ylab = expression("V"^-1), 
     xlab = expression("[ADP]"^-1), 
     axes = F, 
     frame.plot = T, 
     mgp = c(0.3,1,1), 
     xlim = c(-0.03, 0.2), 
     ylim = c(-0.001, 0.02))
abline(v = 0, lwd = 2, col = "grey50")
abline(h = 0, lwd = 2, col = "grey50")

abline(a = 1/coefs[[1]], b= coefs[[3]]/coefs[[1]], col = "black")
abline(a = 1/coefs[[2]], b= coefs[[3]]/coefs[[2]], col = "red3")


x <- rnorm(100) 
y <- rbinom(100, 1, 0.5)

dev.off()

summary(model.drm)
summary(model.drm.kd)
summary(model.drm.Vmax)

modelFit(model.drm.kd)


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
  labs(title = "Hanes–Woolf plot")

ggsave('adp_hanes-woolf.png', width = 5, height = 5, units = "in", dpi = 300)

