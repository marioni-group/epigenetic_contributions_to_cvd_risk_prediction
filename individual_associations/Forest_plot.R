library(dplyr)
library(ggplot2)
setwd("/Volumes/marioni-lab/Ola/Shiny/CVD_EpiScores/data/")
path = '/Volumes/marioni-lab/Ola/Shiny/CVD_EpiScores/data/'
full_annots = readRDS(paste0(path,"EpiAnnots_prepped.RDS"))

f1 = 'HRs_assign_epi_wbc.csv'
f2 = 'HRs_assign_epi_cTnI_corrected_wbc.csv'

mod = read.csv(paste0(path, f1))
mod = left_join(mod, full_annots, by=c("protein" = "ID")) # 109
mod$Model = "ASSIGN + Protein EpiScore"
mod$Name = mod$Gene
mod = subset(mod, !duplicated(Name)) #103
mod = subset(mod, local > 0.05 & global > 0.05) #89
number_of_protein_episcores = nrow(mod)
threshold = 0.05/number_of_protein_episcores # 0.0005617978
bonf_corrected = subset(mod, p<threshold) #29
mod = subset(mod, p<0.05) # 50

troponin = read.csv(paste0(path, f2))
troponin = left_join(troponin, full_annots, by=c("protein" = "ID")) # 109
troponin$Model = "ASSIGN + cTnI + Protein EpiScore"
troponin$Name = troponin$Gene
troponin = subset(troponin, !duplicated(Name)) #103
troponin = subset(troponin, local > 0.05 & global > 0.05) #96
number_of_protein_episcores = nrow(troponin)
threshold = 0.05/number_of_protein_episcores # 0.0005208333
bonf_corrected = subset(troponin, p<threshold) #27
troponin = subset(troponin, troponin$protein %in% mod$protein) # 49

increase = subset(troponin, hr > 1) # 31
decrease = subset(troponin, hr < 1) # 18

x = rbind(mod, troponin)
episcore_names = x$Name
x = x %>% mutate(Name = reorder(Name, hr))
x$Outcome <- "Protein Episcores, p<0.05"
x = mutate(x, Name = reorder(Name, hr))
stacked = ggplot(x,aes(y=hr, x=Name, group=Model, colour=Model)) + 
  geom_point(size = 2, position = position_dodge(0.5))+
  geom_errorbar(aes(ymin = lci, ymax = uci),
                position = position_dodge(0.5), width = 0.1)+
  labs(x="", y="Hazard Ratio per SD [95% Confidence Interval]",
       col="Model covariates") +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 8, vjust = 0.5), 
        axis.text.y = element_text(size = 8), legend.position = "bottom",
        plot.title = element_text(size = 8))+ theme(legend.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#4E84C4", "#D73526")) +
  coord_flip()

pdf("/Users/shirin/Desktop/Forest.pdf", width=7, height=9)
stacked
dev.off()

