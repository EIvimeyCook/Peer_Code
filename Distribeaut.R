#this function will create a nice alternative to a boxplot, twovariants give different alternatives.

Distribeaut_Shape <- function(data, xvar, yvar, fillvar, shapevar){
  require(grid)
  require(ggthemes)
  require(ggplot2)
  require(PupillometryR)
  require(gghalves)
  require(ggrepel)
  require(rlang)
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
    data_mean <- data %>%
    group_by(!!as.symbol((xvar)), !!as.symbol((fillvar)),!!as.symbol((shapevar))) %>%
    summarise(meanyvar = mean(!!as.symbol((yvar))))

    dodge <- position_dodge(0.2)    
    
Dis<-ggplot(data, aes_string(x = xvar, y = yvar, colour = fillvar, fill = fillvar, shape = shapevar)) +
  geom_flat_violin(position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5)+
  geom_jitter(position = position_jitter(width = .1), size = 1, shape = 20, alpha = 0.05)+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position=dodge) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2, position=dodge) +
  theme_Publication() +
  geom_label_repel(data = data_mean,
                   mapping = aes(x = !!as.symbol((xvar)), y = meanyvar, label = sprintf("%0.2f", round(meanyvar,3))),
                   fontface = 'bold',
                   color = 'black',
                   max.iter = 3e2,
                   box.padding = 0.35,
                   point.padding = 0.5,
                   position = dodge,
                   segment.color = 'grey50',
                   segment.linetype = 6,
                   segment.curvature = -1e-20,
                   force = 2, show.legend = F)
  

return(Dis)
}


Distribeaut <- function(data, xvar, yvar, fillvar){
  require(grid)
  require(ggthemes)
  require(ggplot2)
  require(PupillometryR)
  require(gghalves)
  require(ggrepel)
  require(rlang)
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
  
  data_mean <- data %>%
    group_by(!!as.symbol((xvar)), !!as.symbol((fillvar))) %>%
    summarise(meanyvar = mean(!!as.symbol((yvar))))
  
  dodge <- position_dodge(0.2)    
  
  Dis<-ggplot(data, aes_string(x = xvar, y = yvar, colour = fillvar, fill = fillvar)) +
    geom_flat_violin(position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5)+
    geom_jitter(position = position_jitter(width = .1), size = 1, shape = 20, alpha = 0.05)+
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position=dodge) +
    stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2, position=dodge) +
    theme_Publication() +
    geom_hline(data = data_mean, mapping = aes(yintercept = meanyvar[1]), color = 'black', linetype = "dotted") +
    geom_hline(data = data_mean, mapping = aes(yintercept = meanyvar[1]), color = 'black', linetype = "dotted") +
    geom_hline(data = data_mean, mapping = aes(yintercept = meanyvar[1]), color = 'black', linetype = "dotted") +
    geom_label_repel(data = data_mean,
                     mapping = aes(x = !!as.symbol((xvar)), y = meanyvar, label = sprintf("%0.3f", round(meanyvar,3))),
                     fontface = 'bold',
                     color = 'black',
                     max.iter = 3e2,
                     box.padding = 0.35,
                     point.padding = 0.5,
                     position = dodge,
                     segment.color = 'grey50',
                     segment.linetype = 6,
                     segment.curvature = -1e-20,
                     force = 2, show.legend = F)
  
  
  return(Dis)
}


