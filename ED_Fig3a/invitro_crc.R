library(tidyverse)

mydata1 <-
  do.call("rbind", sapply(list.files(pattern = "*\\d.txt"), read.table,
                        col.names = c("conc", "rep_1", "rep_2", "rep_3"), 
                        simplify = FALSE)) %>%
  rownames_to_column("file") %>%
  as_tibble() %>%
  mutate(file = sapply(str_split(file, pattern = "\\."), `[[`, 1)) %>%
  separate(file, c("compound", "bio_rep")) %>%
  gather(tech_rep, activity, contains("rep_")) %>%
  separate(tech_rep, c('drop', 'tech_rep')) %>%
  dplyr::select(-drop) %>%
  filter(bio_rep == 1)

mydata2 <-
  do.call("rbind", sapply(list.files(pattern = "*\\d.txt"), read.table,
                          col.names = c("conc", "rep_1", "rep_2", "rep_3"), 
                          simplify = FALSE)) %>%
  rownames_to_column("file") %>%
  as_tibble() %>%
  mutate(file = sapply(str_split(file, pattern = "\\."), `[[`, 1)) %>%
  separate(file, c("compound", "bio_rep")) %>%
  gather(tech_rep, activity, contains("rep_")) %>%
  separate(tech_rep, c('drop', 'tech_rep')) %>%
  dplyr::select(-drop) %>%
  filter(bio_rep == 2) %>%
  mutate()


library(drc)

drc1 <- drm(activity ~ conc, compound,
            data = mydata1,
           fct = LL.4())

drc2 <- drm(activity ~ conc, compound,
            data = mydata2,
            fct = LL.4())

mypal <- c("Black", "DodgerBlue3", "#d95f02")
pdf(paste0("bio_rep1", ".pdf"), width = 3.5, height = 2.5, 
    pointsize = 8)  #size of text (10pt)

par(mfrow = c(1,1), mar = c(4.1, 2.6, 1.1, 1.1))

#?plot.drc
plot(drc1, col = mypal, 
     broken = T,
     bp = 0.01, 
     xlab = "Dose(nM)", 
     ylab = "ATP:ADP [PercevalHR]\n(Ex:488/Ex:405, Em:525)", 
     xtsty = "standard", 
     legendPos = c(20000, 0.9))

dev.off()
