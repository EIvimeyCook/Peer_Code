#create a forest plot from coxme values
meforestplot <-function(cox, reference){
  require(AICcmodavg)
  require(ggplot2)
  store <- matrix(nrow = length(cox$coefficients) + 1, ncol = 4)
  ref <- c(paste(reference), 0,NA, NA)
  store[1,] <- ref
  for (x in 1:length(cox$coefficients)){
    y = x+1
    mean<-cox$coefficients[x]
    CIL <- cox$coefficients[x] - (1.96*extractSE(cox)[x])
    CUL <- cox$coefficients[x] + (1.96*extractSE(cox)[x])
    store[y, 1:4] <- c(names(cox$coefficients[x]), mean, CIL, CUL)}
  colnames(store) <- c("Treatment", "mean", "CIL","CUL")
  store <- as.data.frame(store)
  store$mean <- as.numeric(as.character(store$mean))
  store$CIL <- as.numeric(as.character(store$CIL))
  store$CUL <- as.numeric(as.character(store$CUL))
  
  forest<-ggplot(store, aes(x = mean, xmax = CUL, xmin = CIL, y = Treatment, colour = Treatment)) +
    geom_point(size = 4) +
    geom_errorbarh(height = 0.3, size = 2) +
    theme_Publication() +
    xlab("Coefficient Estimate") +
    geom_vline(xintercept=0, linetype="dotted") +
    scale_y_discrete(limits = rev) +
    #scale_colour_discrete(labels = c( "Inbred-High", "Inbred-Low", "Outbred-High", paste(c.name))) +
    ylab(NULL) +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank())
    #theme(legend.position="none") 
    
  return(forest)
  
}


