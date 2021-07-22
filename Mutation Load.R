# Preamble - load these packages#######
library(ggplot2)
library(dplyr)
library(lme4)
library(tidyr)
library(dplyr)
library(tidyverse)
library(glmmTMB)
library(bbmle)
library(gridExtra)
library(dabestr)
library(DHARMa)
library(BaSTA)
library(dplyr)
library(tidyverse)
library(doSNOW)
library(snowfall)
library(MuMIn)
library(effects)
library(survival)
library(survminer)
library(coxme)
library(emmeans)
library(ggpubr)
library(gt)
library(remotes)
library(devtools)
library(githubinstall)
library(lavaan)
library(PupillometryR)
library(ggbeeswarm)
library(gghalves)
library(car)
library(hablar)
library(sjPlot)
library(grid)
library(ggthemes)

theme_Publication <- function(base_size = 15, base_family = "Arial") {
  (theme_foundation(base_size = base_size, base_family = base_family)
  + theme(
      plot.title = element_text(
        face = "bold",
        size = rel(1.2), hjust = 0.5
      ),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.direction = "horizontal",
      # legend.key.size = unit(0.2, "cm"),
      legend.position = "bottom",
      # legend.justification = "bottom",
      # legend.box.just = "right",
      # legend.margin = margin(6, 6, 6, 6),
      legend.title = element_text(face = "italic"),
      # legend.spacing.x = unit(1.0, "cm"),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      legend.key.width = unit(2, "line"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    ))
}

# Data ########

#Load the data and create a new treatment variable
Data <- read_csv("~/OneDrive - University of East Anglia/Postdoc/MS/Mutation Load/Mutation_Load.csv") %>%
  convert(fct(Treatment.ID, Breeding, Temperature, Parent.ID, Egg.ID, ID, Hatch.Y.N, Sex.M.F, Parent.No)) %>%
  mutate(Treatment = paste(paste(Breeding, "-", Temperature))) %>%
  convert(fct(Treatment))

#check levels
Data$Temperature
Data$Breeding
Data$Treatment
Data <- within(Data, Temperature <- relevel(Temperature, ref = 2))
Data <- within(Data, Breeding <- relevel(Breeding, ref = 2))
Data <- within(Data, Treatment <- relevel(Treatment, ref = 3))

# Development Time #######
# Graph

#Use the dsitribueat function
D2<-Distribeaut(data = Data, xvar = "Breeding", yvar = "D.T.", fillvar = "Temperature") +
  xlab("Breeding Status") +
  ylab("Development Time (Days)")
  
#Save the plot
ggsave(file = "DT.png", plot = D2, width = 297, height = 250, units = "mm", dpi = 300)

# Model using glmmtmb
m1 <- glmmTMB(D.T. ~ Breeding + Temperature + Breeding:Temperature + (1 | Parent.ID), data = Data)
m2 <- glmmTMB(D.T. ~ Breeding + Temperature + (1 | Parent.ID), data = Data)

#check for overall effects
Anova(m1, type = "III")

#summarise
summary(mq)
summary(m2)

#create a summary table
tab_model(m1,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)
tab_model(m2,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#Emmeans at each level
em1<-emmeans(m1, specs = pairwise ~ Breeding*Temperature, adjust = "none",type = "response")
em2<-emmeans(m2, specs = pairwise ~ Breeding, adjust = "none",type = "response")


# Hatching Success #########

# ggplot Graph
g1 <- ggplot(Data, aes(x = Breeding, y = Hatch, colour = Temperature, fill = Temperature)) +
  geom_jitter(height = 0.02, alpha = 0.3) +
  #geom_smooth(aes(group = 1), data = Data %>% filter(Temperature == "Elevated"), method = "glm", method.args = list(family = "binomial")) +
  #geom_smooth(aes(group = 1), data = Data %>% filter(Temperature == "Control"), method = "glm", method.args = list(family = "binomial")) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.05) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 1.5) +
  #stat_summary(fun.data = "mean_cl_boot", geom = "line", size = 0.8) +
  xlab("Breeding Status") +
  ylab("Eclosion Success") +
  scale_x_discrete(limits = c("Inbred", "Outbred")) +
  theme_Publication() +
  coord_cartesian(y = c(0,1))

#view the ggplot
g1

#save he plot
ggsave(file = "Eclosion.png", plot = g1, width = 297, height = 250, units = "mm", dpi = 300)

# Model using glmmtmb
m1 <- glmmTMB(Hatch ~ Breeding  + Temperature + Breeding:Temperature + (1 | Parent.ID), data = Data, family = "binomial")
m2 <- glmmTMB(Hatch ~ Breeding + Temperature + (1 | Parent.ID), data = Data, family = "binomial")

#overall differences
Anova(m1, type = "III")

#summarise
summary(m1)
summary(m2)

#table of summary
tab_model(m1,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)
tab_model(m2,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#Emmeans at each level
em1<-emmeans(m1, specs = pairwise ~ Breeding*Temperature, adjust = "none",type = "response")
em2<-emmeans(m2, specs = pairwise ~ Breeding, adjust = "none",type = "response")

# Age-specific reproduction #####

#convert to longform from wideform and convert to factors
longrepo <- Data %>% 
  gather(Day, Offspring, 17:21) %>%
  convert(num(Day)) %>%
  mutate(Day2 = Day^2) %>%
  mutate(obs = seq_len(nrow(.))) %>%
  drop_na(Offspring) %>%
  convert(fct(obs))

#graph
g2 <- ggplot(longrepo, aes(x = Day, y = Offspring, colour = Temperature, linetype = Breeding, group = interaction(Temperature, Breeding))) +
  geom_jitter(alpha = 0.2) +
  # geom_boxplot() +
  # facet_grid(~Day, scales = "free", space = "free") +
  # theme(strip.text.y = element_text(angle = 0)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2) +
  stat_summary(fun.data = "mean_cl_boot", geom = "line", size = 1.1) +
  theme_Publication() +
  ylab("Offspring Count") +
  scale_linetype_manual(values = c("dotdash", "solid")) +
  guides(linetype=guide_legend(title="Breeding Status"))

#view graph
g2

#save graph
ggsave(file = "ASR.png", plot = g2, width = 297, height = 250, units = "mm", dpi = 300)

#create a list for comparison
ASRRepro <- list()

# Model checks for overdispersion
ASRRepro[[1]] <- glmmTMB(Offspring ~ Temperature * Breeding * Day + Temperature * Breeding * Day2 + (1 | Parent.ID / ID), data = longrepo, family = "poisson")
ASRRepro[[2]] <- glmmTMB(Offspring ~ Temperature * Breeding * Day + Temperature * Breeding * Day2 + (1|obs) + (1 | Parent.ID / ID), data = longrepo, family = "poisson")

#dharma simulation and zeroinflation
s1 <- simulateResiduals(ASRRepro[[1]])
plot(s1)
testZeroInflation(s1)

s1 <- simulateResiduals(ASRRepro[[2]])
plot(s1)
testZeroInflation(s1)

#create the two formulas (w/wo olre)
formula <- c(Offspring ~ Temperature * Breeding * Day + Temperature * Breeding * Day2 + (1 | Parent.ID / ID))
formula2 <- c(Offspring ~ Temperature * Breeding * Day + Temperature * Breeding * Day2 + (1|obs) + (1 | Parent.ID / ID))

#run with different zi components and dists
ASRRepro[[	1	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~0	)																																															
ASRRepro[[	2	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~1)																																																
ASRRepro[[	3	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~Breeding)																																																
ASRRepro[[	4	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~Temperature)																																																
ASRRepro[[	5	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~Day)																																																
ASRRepro[[	6	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day)																																													
ASRRepro[[	7	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Day2)																																													
ASRRepro[[	8	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Temperature)																																													
ASRRepro[[	9	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Temperature)																																													
ASRRepro[[	10	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2)																																											
ASRRepro[[	11	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature)																																											
ASRRepro[[	12	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Day2	+	Temperature)																																											
ASRRepro[[	13	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature)																																									
ASRRepro[[	14	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Breeding:Day)																																											
ASRRepro[[	15	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Temperature	+	Day	+	Temperature:Day)																																											
ASRRepro[[	16	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Breeding:Day)																																									
ASRRepro[[	17	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day)																																									
ASRRepro[[	18	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day2	+	Day	+	Breeding:Day2)																																									
ASRRepro[[	19	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day)																																							
ASRRepro[[	20	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2)																																							
ASRRepro[[	21	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Breeding:Day	+	Breeding:Day2)																																							
ASRRepro[[	22	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2)																																					
ASRRepro[[	23	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Temperature	+	Breeding:Temperature)																																											
ASRRepro[[	24	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Temperature)																																									
ASRRepro[[	25	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature)																																							
ASRRepro[[	26	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature)																																							
ASRRepro[[	27	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature)																																					
ASRRepro[[	28	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature)																																					
ASRRepro[[	29	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature)																																			
ASRRepro[[	30	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Day:Temperature)																																									
ASRRepro[[	31	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day:Temperature)																																									
ASRRepro[[	32	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day:Temperature)																																							
ASRRepro[[	33	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Day:Temperature)																																							
ASRRepro[[	34	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day:Temperature)																																					
ASRRepro[[	35	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day:Temperature)																																					
ASRRepro[[	36	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day:Temperature)																																			
ASRRepro[[	37	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Temperature	+	Day:Temperature)																																							
ASRRepro[[	38	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day:Temperature)																																					
ASRRepro[[	39	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature)																																					
ASRRepro[[	40	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature)																																			
ASRRepro[[	41	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature)																																			
ASRRepro[[	42	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature)																																	
ASRRepro[[	43	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day2:Temperature)																																									
ASRRepro[[	44	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day2:Temperature)																																							
ASRRepro[[	45	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day2:Temperature)																																					
ASRRepro[[	46	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day2:Temperature)																																					
ASRRepro[[	47	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day2:Temperature)																																			
ASRRepro[[	48	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day2:Temperature)																																					
ASRRepro[[	49	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	50	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	51	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	52	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day:Temperature	+	Day2:Temperature)																																							
ASRRepro[[	53	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day:Temperature	+	Day2:Temperature)																																					
ASRRepro[[	54	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	55	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	56	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	57	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	58	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	59	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	60	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																															
ASRRepro[[	61	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																																			
ASRRepro[[	62	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																																	
ASRRepro[[	63	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																															
ASRRepro[[	64	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature)																															
ASRRepro[[	65	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature)																													
ASRRepro[[	66	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																																	
ASRRepro[[	67	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																															
ASRRepro[[	68	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																															
ASRRepro[[	69	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																													
ASRRepro[[	70	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature	+	Breeding:Day2:Temperature)																											

ASRRepro[[	71	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~0	)																																															
ASRRepro[[	72	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~1)																																																
ASRRepro[[	73	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~Breeding)																																																
ASRRepro[[	74	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~Temperature)																																																
ASRRepro[[	75	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~Day)																																																
ASRRepro[[	76	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day)																																													
ASRRepro[[	77	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Day2)																																													
ASRRepro[[	78	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Temperature)																																													
ASRRepro[[	79	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Temperature)																																													
ASRRepro[[	80	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2)																																											
ASRRepro[[	81	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature)																																											
ASRRepro[[	82	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Day2	+	Temperature)																																											
ASRRepro[[	83	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature)																																									
ASRRepro[[	84	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Breeding:Day)																																											
ASRRepro[[	85	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Temperature	+	Day	+	Temperature:Day)																																											
ASRRepro[[	86	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Breeding:Day)																																									
ASRRepro[[	87	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day)																																									
ASRRepro[[	88	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day2	+	Day	+	Breeding:Day2)																																									
ASRRepro[[	89	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day)																																							
ASRRepro[[	90	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2)																																							
ASRRepro[[	91	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Breeding:Day	+	Breeding:Day2)																																							
ASRRepro[[	92	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2)																																					
ASRRepro[[	93	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Temperature	+	Breeding:Temperature)																																											
ASRRepro[[	94	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Temperature)																																									
ASRRepro[[	95	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature)																																							
ASRRepro[[	96	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature)																																							
ASRRepro[[	97	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature)																																					
ASRRepro[[	98	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature)																																					
ASRRepro[[	99	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature)																																			
ASRRepro[[	100	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Day:Temperature)																																									
ASRRepro[[	101	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day:Temperature)																																									
ASRRepro[[	102	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day:Temperature)																																							
ASRRepro[[	103	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Day:Temperature)																																							
ASRRepro[[	104	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day:Temperature)																																					
ASRRepro[[	105	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day:Temperature)																																					
ASRRepro[[	106	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day:Temperature)																																			
ASRRepro[[	107	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Temperature	+	Day:Temperature)																																							
ASRRepro[[	108	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day:Temperature)																																					
ASRRepro[[	109	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature)																																					
ASRRepro[[	110	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature)																																			
ASRRepro[[	111	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature)																																			
ASRRepro[[	112	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature)																																	
ASRRepro[[	113	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day2:Temperature)																																									
ASRRepro[[	114	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day2:Temperature)																																							
ASRRepro[[	115	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day2:Temperature)																																					
ASRRepro[[	116	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day2:Temperature)																																					
ASRRepro[[	117	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day2:Temperature)																																			
ASRRepro[[	118	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day2:Temperature)																																					
ASRRepro[[	119	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	120	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	121	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	122	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day:Temperature	+	Day2:Temperature)																																							
ASRRepro[[	123	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day:Temperature	+	Day2:Temperature)																																					
ASRRepro[[	124	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	125	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	126	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	127	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	128	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	129	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	130	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																															
ASRRepro[[	131	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																																			
ASRRepro[[	132	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																																	
ASRRepro[[	133	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																															
ASRRepro[[	134	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature)																															
ASRRepro[[	135	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature)																													
ASRRepro[[	136	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																																	
ASRRepro[[	137	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																															
ASRRepro[[	138	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																															
ASRRepro[[	139	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																													
ASRRepro[[	140	]]	<-	glmmTMB(as.formula(paste(formula2)),	data	=	longrepo,	family	=	poisson,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature	+	Breeding:Day2:Temperature)																											

ASRRepro[[	141	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~0	)																																															
ASRRepro[[	142	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~1)																																																
ASRRepro[[	143	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~Breeding)																																																
ASRRepro[[	144	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~Temperature)																																																
ASRRepro[[	145	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~Day)																																																
ASRRepro[[	146	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day)																																													
ASRRepro[[	147	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Day	+	Day2)																																													
ASRRepro[[	148	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Day	+	Temperature)																																													
ASRRepro[[	149	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Temperature)																																													
ASRRepro[[	150	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2)																																											
ASRRepro[[	151	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Temperature)																																											
ASRRepro[[	152	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Day	+	Day2	+	Temperature)																																											
ASRRepro[[	153	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature)																																									
ASRRepro[[	154	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Breeding:Day)																																											
ASRRepro[[	155	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Temperature	+	Day	+	Temperature:Day)																																											
ASRRepro[[	156	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Breeding:Day)																																									
ASRRepro[[	157	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day)																																									
ASRRepro[[	158	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day2	+	Day	+	Breeding:Day2)																																									
ASRRepro[[	159	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day)																																							
ASRRepro[[	160	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2)																																							
ASRRepro[[	161	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Breeding:Day	+	Breeding:Day2)																																							
ASRRepro[[	162	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2)																																					
ASRRepro[[	163	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Temperature	+	Breeding:Temperature)																																											
ASRRepro[[	164	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Temperature)																																									
ASRRepro[[	165	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature)																																							
ASRRepro[[	166	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature)																																							
ASRRepro[[	167	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature)																																					
ASRRepro[[	168	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature)																																					
ASRRepro[[	169	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature)																																			
ASRRepro[[	170	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Day:Temperature)																																									
ASRRepro[[	171	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day:Temperature)																																									
ASRRepro[[	172	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day:Temperature)																																							
ASRRepro[[	173	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Day:Temperature)																																							
ASRRepro[[	174	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day:Temperature)																																					
ASRRepro[[	175	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day:Temperature)																																					
ASRRepro[[	176	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day:Temperature)																																			
ASRRepro[[	177	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Temperature	+	Day:Temperature)																																							
ASRRepro[[	178	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day:Temperature)																																					
ASRRepro[[	179	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature)																																					
ASRRepro[[	180	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature)																																			
ASRRepro[[	181	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature)																																			
ASRRepro[[	182	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature)																																	
ASRRepro[[	183	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day2:Temperature)																																									
ASRRepro[[	184	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day2:Temperature)																																							
ASRRepro[[	185	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day2:Temperature)																																					
ASRRepro[[	186	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day2:Temperature)																																					
ASRRepro[[	187	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day2:Temperature)																																			
ASRRepro[[	188	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day2:Temperature)																																					
ASRRepro[[	189	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	190	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	191	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	192	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day:Temperature	+	Day2:Temperature)																																							
ASRRepro[[	193	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day:Temperature	+	Day2:Temperature)																																					
ASRRepro[[	194	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	195	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	196	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	197	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	198	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	199	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	200	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																															
ASRRepro[[	201	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																																			
ASRRepro[[	202	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																																	
ASRRepro[[	203	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																															
ASRRepro[[	204	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature)																															
ASRRepro[[	205	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature)																													
ASRRepro[[	206	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																																	
ASRRepro[[	207	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																															
ASRRepro[[	208	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																															
ASRRepro[[	209	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																													
ASRRepro[[	210	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	nbinom2,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature	+	Breeding:Day2:Temperature)																											

ASRRepro[[	211	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~0	)																																															
ASRRepro[[	212	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~1)																																																
ASRRepro[[	213	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~Breeding)																																																
ASRRepro[[	214	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~Temperature)																																																
ASRRepro[[	215	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~Day)																																																
ASRRepro[[	216	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day)																																													
ASRRepro[[	217	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Day	+	Day2)																																													
ASRRepro[[	218	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Day	+	Temperature)																																													
ASRRepro[[	219	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Temperature)																																													
ASRRepro[[	220	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2)																																											
ASRRepro[[	221	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Temperature)																																											
ASRRepro[[	222	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Day	+	Day2	+	Temperature)																																											
ASRRepro[[	223	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature)																																									
ASRRepro[[	224	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Breeding:Day)																																											
ASRRepro[[	225	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Temperature	+	Day	+	Temperature:Day)																																											
ASRRepro[[	226	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Breeding:Day)																																									
ASRRepro[[	227	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day)																																									
ASRRepro[[	228	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day2	+	Day	+	Breeding:Day2)																																									
ASRRepro[[	229	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day)																																							
ASRRepro[[	230	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2)																																							
ASRRepro[[	231	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Breeding:Day	+	Breeding:Day2)																																							
ASRRepro[[	232	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2)																																					
ASRRepro[[	233	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Temperature	+	Breeding:Temperature)																																											
ASRRepro[[	234	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Temperature)																																									
ASRRepro[[	235	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature)																																							
ASRRepro[[	236	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature)																																							
ASRRepro[[	237	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature)																																					
ASRRepro[[	238	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature)																																					
ASRRepro[[	239	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature)																																			
ASRRepro[[	240	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Day:Temperature)																																									
ASRRepro[[	241	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day:Temperature)																																									
ASRRepro[[	242	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day:Temperature)																																							
ASRRepro[[	243	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Day:Temperature)																																							
ASRRepro[[	244	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day:Temperature)																																					
ASRRepro[[	245	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day:Temperature)																																					
ASRRepro[[	246	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day:Temperature)																																			
ASRRepro[[	247	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Temperature	+	Day:Temperature)																																							
ASRRepro[[	248	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day:Temperature)																																					
ASRRepro[[	249	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature)																																					
ASRRepro[[	250	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature)																																			
ASRRepro[[	251	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature)																																			
ASRRepro[[	252	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature)																																	
ASRRepro[[	253	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day2:Temperature)																																									
ASRRepro[[	254	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day2:Temperature)																																							
ASRRepro[[	255	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day2:Temperature)																																					
ASRRepro[[	256	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day2:Temperature)																																					
ASRRepro[[	257	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day2:Temperature)																																			
ASRRepro[[	258	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day2:Temperature)																																					
ASRRepro[[	259	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	260	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	261	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	262	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Day	+	Day2	+	Temperature	+	Day:Temperature	+	Day2:Temperature)																																							
ASRRepro[[	263	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Day:Temperature	+	Day2:Temperature)																																					
ASRRepro[[	264	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	265	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	266	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	267	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																			
ASRRepro[[	268	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	269	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																																	
ASRRepro[[	270	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature)																															
ASRRepro[[	271	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																																			
ASRRepro[[	272	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																																	
ASRRepro[[	273	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Breeding:Day:Temperature)																															
ASRRepro[[	274	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature)																															
ASRRepro[[	275	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature)																													
ASRRepro[[	276	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																																	
ASRRepro[[	277	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																															
ASRRepro[[	278	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																															
ASRRepro[[	279	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day2:Temperature)																													
ASRRepro[[	280	]]	<-	glmmTMB(as.formula(paste(formula)),	data	=	longrepo,	family	=	genpois,	ziformula	=	~	Breeding	+	Day	+	Day2	+	Temperature	+	Breeding:Day	+	Breeding:Day2	+	Breeding:Temperature	+	Day:Temperature	+	Day2:Temperature	+	Breeding:Day:Temperature	+	Breeding:Day2:Temperature)																											


#check the best model (using dharmas dispersion and zeroninflation summaries)
modelchecker(check = ASRRepro)

#simulate the best model and check again (underdispersed)
s1b <- simulateResiduals(ASRRepro[[	237	]])
plot(s1b)
testZeroInflation(s1b)
testDispersion(s1b)

# summary table of the best model
tab_model(ASRRepro[[	237	]],transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)


#overall anova
Anova(ASRRepro[[	237	]], type = "III")


# Total Reproduction ######

# Remove unhatched indivudals and NA dta - create OLRE
DataNoHCen <- Data %>%
  filter(Hatch == "1") %>%
  filter(!is.na(Total.eclosed)) %>%
  mutate(obs = seq_len(nrow(.)))

#subset by temperature for graphs
DataElevated <- subset(DataNoHCen, Temperature == "Elevated")
DataControl <- subset(DataNoHCen, Temperature == "Control")


#dabest KRS
analysis1a <- DataElevated %>%
  drop_na(Total.eclosed) %>%
  dabest(Breeding, Total.eclosed,
    idx = list(c("Inbred", "Outbred"))
  ) %>%
  mean_diff()

analysis2a <- DataControl %>%
  drop_na(Total.eclosed) %>%
  dabest(Breeding, Total.eclosed,
    idx = list(c("Inbred", "Outbred"))
  ) %>%
  mean_diff()

#create titles
title1 <- textGrob("Elevated Temperature", gp = gpar(fontface = "italic", fontfamily = "Arial"))
title2 <- textGrob("Control Temperature", gp = gpar(fontface = "italic", fontfamily = "Arial"))

#arrange and plot
p1 <- grid.arrange(arrangeGrob(plot(analysis1a, theme = theme_Publication(), rawplot.ylim = c(0, 150), palette = "Dark2", effsize.ylabel = "", rawplot.ylabel = "Total Reproduction"), top = title1),
  arrangeGrob(plot(analysis2a, theme = theme_Publication(), rawplot.ylim = c(0, 150), rawplot.ylabel = ""), top = title2),
  ncol = 2
)

#save
ggsave(file = "LRS.png", plot = p1, width = 297, height = 250, units = "mm", dpi = 300)

# Create model list
LRSRepro <- list()

#check for overdispersion (w/wo olre)
LRSRepro[[1]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), data = DataNoHCen, family = "poisson")
LRSRepro[[2]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1|obs) + (1 | Parent.ID), data = DataNoHCen, family = "poisson")

#simulate residulas with dharma
summary(LRSRepro[[1]])
s1 <- simulateResiduals(LRSRepro[[1]])
plot(s1)
testZeroInflation(s1)

summary(LRSRepro[[2]])
s1 <- simulateResiduals(LRSRepro[[2]])
plot(s1)
testZeroInflation(s1)

#run models with different zero inflation and error
LRSRepro[[1]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~0, data = DataNoHCen, family = "poisson" )
LRSRepro[[2]]<- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~1, data = DataNoHCen, family = "poisson" )
LRSRepro[[3]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~Temperature, data = DataNoHCen, family = "poisson" )
LRSRepro[[4]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~Breeding, data = DataNoHCen, family = "poisson" )
LRSRepro[[5]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~ Temperature + Breeding, data = DataNoHCen, family = "poisson" )
LRSRepro[[6]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~ Temperature * Breeding, data = DataNoHCen, family = "poisson" )

LRSRepro[[19]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1|obs) +(1 | Parent.ID), zi = ~0, data = DataNoHCen, family = "poisson" )
LRSRepro[[20]]<- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1|obs) + (1 | Parent.ID), zi = ~1, data = DataNoHCen, family = "poisson" )
LRSRepro[[21]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1|obs) + (1 | Parent.ID), zi = ~Temperature, data = DataNoHCen, family = "poisson" )
LRSRepro[[22]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1|obs) + (1 | Parent.ID), zi = ~Breeding, data = DataNoHCen, family = "poisson" )
LRSRepro[[23]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1|obs) + (1 | Parent.ID), zi = ~ Temperature + Breeding, data = DataNoHCen, family = "poisson" )
LRSRepro[[24]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1|obs) + (1 | Parent.ID), zi = ~ Temperature * Breeding, data = DataNoHCen, family = "poisson" )

LRSRepro[[7]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~0, data = DataNoHCen, family = "nbinom2" )
LRSRepro[[8]]<- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~1, data = DataNoHCen, family = "nbinom2" )
LRSRepro[[9]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~Temperature, data = DataNoHCen, family = "nbinom2" )
LRSRepro[[10]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~Breeding, data = DataNoHCen, family = "nbinom2" )
LRSRepro[[11]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~ Temperature + Breeding, data = DataNoHCen, family = "nbinom2" )
LRSRepro[[12]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~ Temperature * Breeding, data = DataNoHCen, family = "nbinom2" )

LRSRepro[[13]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~0, data = DataNoHCen, family = "genpois" )
LRSRepro[[14]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~1, data = DataNoHCen, family = "genpois" )
LRSRepro[[15]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~Temperature, data = DataNoHCen, family = "genpois" )
LRSRepro[[16]]<- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~Breeding, data = DataNoHCen, family = "genpois" )
LRSRepro[[17]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~ Temperature + Breeding, data = DataNoHCen, family = "genpois" )
LRSRepro[[18]] <- glmmTMB(Total.eclosed ~ Temperature * Breeding + (1 | Parent.ID), zi = ~ Temperature * Breeding, data = DataNoHCen, family = "genpois" )

#check the best model (using dharmas dispersion and zeroninflation summaries)
modelchecker(check = LRSRepro)

# Anova summary
Anova(LRSRepro[[15]], type = "III")

#table of best model
tab_model(LRSRepro[[15]],transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

# check without interaction for overall effects
g1 <- glmmTMB(Total.eclosed ~ Temperature + Breeding + (1 | Parent.ID), zi = ~Temperature, data = DataNoHCen, family = "genpois" )

#again, table of best model
tab_model(g1,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#Emmeans of differences
em1<-emmeans(LRSRepro[[15]], specs = pairwise ~ Breeding*Temperature, adjust = "none",type = "response")
em2<-emmeans(g1, specs = pairwise ~ Temperature, adjust = "none",type = "response")


# Lambda #####


#Subset by removing total eclosion that are NA
LamDat <- subset(Data, Total.eclosed != "NA")

#remove row names for loop
row.names(LamDat) <- NULL


#fitness loop to calculate lambda
L <- matrix(nrow = nrow(LamDat), ncol = 7)
for (i in 1:nrow(LamDat)) {
  Les <- matrix(0, ncol = LamDat$D.T.[i] + 5, nrow = LamDat$D.T.[i] + 5)
  diag(Les[-1, ]) <- rep(1, LamDat$D.T.[i] + 4)
  Fert <- c(rep(0, LamDat$D.T.[i]), as.numeric(as.vector(LamDat[i, ][17:21])))
  Fert[is.na(Fert)] <- 0
  Les[1, ] <- c(Fert)
  class(Les) <- "leslie.matrix"
  Lambda <- popbio::eigen.analysis(Les)$lambda1
  eigen(Les)
  L[i, 1:7] <- c(paste0(LamDat$ID[i]), paste0(LamDat$Parent.ID[i]), paste0(LamDat$Breeding[i]), paste0(LamDat$Temperature[i]), Lambda, paste0(LamDat$D.T.[i]), paste0(LamDat$Treatment[i]))
}

#create dataframe and change names
L <- as.data.frame(L)
names(L) <- c("ID", "ParentID", "Breeding", "Temperature", "Lambda", "DT", "Treatment")

#convert
L %<>% convert(num(Lambda, DT),
               fct(ID, ParentID, Breeding, Temperature, Treatment))


#subset by temperature
LElevated <- subset(L, Temperature == "Elevated")
LControl <- subset(L, Temperature == "Control")

#dabestr 
analysis1aL <- LElevated %>%
  dabest(Breeding, Lambda, idx = list(c("Inbred", "Outbred"))) %>%
  mean_diff()

analysis2aL <- LControl %>%
  dabest(Breeding, Lambda, idx = list(c("Inbred", "Outbred"))) %>%
  mean_diff()

#title
title1 <- textGrob("Elevated Temperature", gp = gpar(fontface = "italic", fontfamily = "Arial"))
title2 <- textGrob("Control Temperature", gp = gpar(fontface = "italic", fontfamily = "Arial"))

#arrange
p2 <- grid.arrange(arrangeGrob(plot(analysis1aL, theme = theme_Publication(), rawplot.ylim = c(0, 1.5), palette = "Dark2", effsize.ylabel = "", rawplot.ylabel = "Lambda"), top = title1),
                   arrangeGrob(plot(analysis2aL, theme = theme_Publication(), rawplot.ylim = c(0, 1.5), rawplot.ylabel = ""), top = title2),
                   ncol = 2
)

#save
ggsave(file = "Lambda.png", plot = p2, width = 297, height = 250, units = "mm", dpi = 300)

#change level of model for convenience
L$Breeding
L <- within(L, Breeding <- relevel(Breeding, ref = 2))

#model with and without interaction
m1 <- glmmTMB(Lambda ~ Breeding + Temperature + Breeding:Temperature + (1 | ParentID), data = L, family = "gaussian")
m2 <- glmmTMB(Lambda ~ Breeding + Temperature + (1 | ParentID), data = L, family = "gaussian")

# anova of overall effects
Anova(m1, type = "III")

#summary table
tab_model(m1,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)
tab_model(m2,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#Emmeans at different levels
em1<-emmeans(m1, specs = pairwise ~ Breeding*Temperature, adjust = "none",type = "response")
em2<-emmeans(m2, specs = pairwise ~ Temperature, adjust = "none",type = "response")


# Lifespan ######

#remove indvidausl that havent hatched r are censored
DataNoH <- Data %>%
  filter(Hatch == "1") %>%
  filter(!is.na(Death))

# Cox Model 
F0 <- coxme(Surv(Death, Observed) ~ Breeding*Temperature + (1 | Parent.No), data = DataNoH)
summary(F0)

# look at overall effects
Anova(F0, type = "III")

# summary table
tab_model(F0,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#Emmeans of hazard ratio
em1<-emmeans(F0, specs = pairwise ~ Breeding*Temperature, adjust = "none")

# use treatment for graph
cox <- coxme(Surv(Death, Observed) ~ Treatment + (1 | Parent.No), data = DataNoH)

#look at summary
summary(cox)

#meforest plot of hazard ratios
s1 <- meforestplot(cox, "Outbred - Control") + 
  scale_colour_discrete(labels = c( "Outbred - Control", "Inbred - Control", "Inbred - Elevated", "Outbred - Elevated"))

#save
ggsave(file = "Cox.png", plot = s1, width = 297, height = 250, units = "mm", dpi = 300)

#distribution of lifespan
D2<-Distribeaut(data = DataNoH, xvar = "Breeding", yvar = "Death", fillvar = "Temperature") +
  xlab("Breeding Status") +
  ylab("Lifespan (Days)") +
  scale_x_discrete(limits = rev)

#save
ggsave(file = "Life.png", plot = D2, width = 297, height = 250, units = "mm", dpi = 300)

# Inbreeding depression calc (optional) ########

#create subsets for DT
InbredC <- Data %>%
  filter(Breeding == "Inbred" & Temperature == "Control") %>%
  drop_na(D.T.)

OutbredC <- Data %>%
  filter(Breeding == "Outbred" & Temperature == "Control") %>%
  drop_na(D.T.)

InbredE <- Data %>%
  filter(Breeding == "Inbred" & Temperature == "Elevated") %>%
  drop_na(D.T.)

OutbredE <- Data %>%
  filter(Breeding == "Outbred" & Temperature == "Elevated") %>%
  drop_na(D.T.)


#create inbreeding measures for DT
CDT<-(mean(InbredC$D.T.)- mean(OutbredC$D.T.))/mean(OutbredC$D.T.)

EDT<-(mean(InbredE$D.T.)- mean(OutbredE$D.T.))/mean(OutbredE$D.T.)

ECDT<-(mean(InbredE$D.T.)- mean(OutbredC$D.T.))/mean(OutbredC$D.T.)

#create subsets for Hatch
InbredC <- Data %>%
  filter(Breeding == "Inbred" & Temperature == "Control") %>%
  drop_na(Hatch)

OutbredC <- Data %>%
  filter(Breeding == "Outbred" & Temperature == "Control") %>%
  drop_na(Hatch)

InbredE <- Data %>%
  filter(Breeding == "Inbred" & Temperature == "Elevated") %>%
  drop_na(Hatch)

OutbredE <- Data %>%
  filter(Breeding == "Outbred" & Temperature == "Elevated") %>%
  drop_na(Hatch)


#create inbreeding dep for hatch
CHatch<-(mean(OutbredC$Hatch)- mean(InbredC$Hatch))/mean(OutbredC$Hatch)

EHatch<-(mean(OutbredE$Hatch)- mean(InbredE$Hatch))/mean(OutbredE$Hatch)

CEHatch<-(mean(OutbredC$Hatch) - mean(InbredE$Hatch))/mean(OutbredC$Hatch)


#create subsets for longevity
InbredC <- Data %>%
  filter(Breeding == "Inbred" & Temperature == "Control") %>%
  drop_na(Death)

OutbredC <- Data %>%
  filter(Breeding == "Outbred" & Temperature == "Control") %>%
  drop_na(Death)

InbredE <- Data %>%
  filter(Breeding == "Inbred" & Temperature == "Elevated") %>%
  drop_na(Death)

OutbredE <- Data %>%
  filter(Breeding == "Outbred" & Temperature == "Elevated") %>%
  drop_na(Death)


#calculate inbreeeding dep for longevity
CSurv<-(mean(OutbredC$Death)- mean(InbredC$Death))/mean(OutbredC$Death)

ESurv<-(mean(OutbredE$Death)- mean(InbredE$Death))/mean(OutbredE$Death)

ECSurv<-(mean(OutbredC$Death) - mean(InbredE$Death))/mean(OutbredC$Death)


#subset for lambda
InbredC <- L %>%
  filter(Breeding == "Inbred" & Temperature == "Control") %>%
  drop_na(Lambda)

OutbredC <- L %>%
  filter(Breeding == "Outbred" & Temperature == "Control") %>%
  drop_na(Lambda)

InbredE <- L %>%
  filter(Breeding == "Inbred" & Temperature == "Elevated") %>%
  drop_na(Lambda)

OutbredE <- L %>%
  filter(Breeding == "Outbred" & Temperature == "Elevated") %>%
  drop_na(Lambda)

#calculate inbreeding dep for lambda
CLam<-(mean(OutbredC$Lambda)- mean(InbredC$Lambda))/mean(OutbredC$Lambda)

ELam<-(mean(OutbredE$Lambda)- mean(InbredE$Lambda))/mean(OutbredE$Lambda)

ECLam<-(mean(OutbredC$Lambda) - mean(InbredE$Lambda))/mean(OutbredC$Lambda)

#subset for LRS
InbredC <- Data %>%
  filter(Breeding == "Inbred" & Temperature == "Control") %>%
  drop_na(Total.eclosed)

OutbredC <- Data %>%
  filter(Breeding == "Outbred" & Temperature == "Control") %>%
  drop_na(Total.eclosed)

InbredE <- Data %>%
  filter(Breeding == "Inbred" & Temperature == "Elevated") %>%
  drop_na(Total.eclosed)

OutbredE <- Data %>%
  filter(Breeding == "Outbred" & Temperature == "Elevated") %>%
  drop_na(Total.eclosed)

#calculate inbreeding dep for LRS
CTotal<- (mean(OutbredC$Total.eclosed)- mean(InbredC$Total.eclosed))/mean(OutbredC$Total.eclosed)

ETotal<-(mean(OutbredE$Total.eclosed)- mean(InbredE$Total.eclosed))/mean(OutbredE$Total.eclosed)

CETotal<-(mean(OutbredC$Total.eclosed) - mean(InbredE$Total.eclosed))/mean(OutbredC$Total.eclosed)

#create a nice table
IBred <- rbind(CDT, EDT, ECDT, CHatch, EHatch,CEHatch, CTotal, ETotal, CETotal, CLam, ELam, ECLam, CSurv, ESurv, ECSurv)
Trait <- rbind("Development Time", "Development Time", "Development Time", "Hatching Success", "Hatching Success", "Hatching Success", "LRS", "LRS", "LRS", "Lambda", "Lambda", "Lambda", "Survival", "Survival", "Survival")
Temp <- rbind("C-C", "E-E", "C-E", "C-C", "E-E", "C-E","C-C", "E-E", "C-E","C-C", "E-E", "C-E","C-C", "E-E", "C-E")

IBreed <- tibble(delta = IBred,Trait = Trait,Temperature = Temp) %>%
  mutate(delta = sprintf("%0.3f", delta))

#Inbreeding depresion tables
IBreed %>%
  gt(groupname_col = "Trait", rowname_col = "delta") %>%
  cols_label(Temperature = " ") %>%
  tab_header(
    title = "Inbreeding Depression (\u03b4)") %>%
  tab_options(row_group.background.color = "#FFEFDB") %>%
  tab_footnote(
    footnote = md("Outbred Control - Inbred Control"),
    locations = cells_body(
      columns = vars(Temperature),
      rows = Temperature == "C-C")) %>%
  tab_footnote(
    footnote = md("Outbred Elevated - Inbred Elevated"),
    locations = cells_body(
      columns = vars(Temperature),
      rows = Temperature == "E-E")) %>%
  tab_footnote(
    footnote = md("Outbred Control - Inbred Elevated"),
    locations = cells_body(
      columns = vars(Temperature),
      rows = Temperature == "C-E"))


