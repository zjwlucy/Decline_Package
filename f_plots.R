## functions to 
## generate plots 

library(ggplot2)
#library(gridExtra)


f_plots <- function(dat){
        print("----------Plotting the full dataset----------")
        dat <- dat[order(dat$FID, dat$IID, dat$timefactor_spiro),]

   
        p1 <- ggplot(data = dat, aes(x = timefactor_spiro,  y = pre_fev1, col = factor(IID)) ) +
                     viridis::scale_color_viridis(discrete = TRUE,option = "viridis") + 
                     geom_line(alpha = 0.32) + theme_minimal() + 
                     theme(legend.position="none", axis.text = element_text(size=30), axis.title = element_text(size=35) ) + 
                     labs(title = "", col = "ID") + xlab("Time (Years)") + ylab("FEV1") + ylim(0, 6000)
        png('PLOTSP_CohortName_2023.png', width=1000, height=850)
            print(p1)
        dev.off()


      # random sample 100 individuals to plot
        set.seed(2023)
        r_id <- sample(unique(dat$IID), 100, replace = FALSE)
        r100 <- dat[which(dat$IID %in% r_id), ]
        
        p2 <- ggplot(data = r100, aes(x = timefactor_spiro,  y = pre_fev1, col = factor(IID)) ) +
                     viridis::scale_color_viridis(discrete = TRUE,option = "viridis") + 
                     geom_line(alpha = 0.32) + theme_minimal() + 
                     theme(legend.position="none", axis.text = element_text(size=30), axis.title = element_text(size=35) ) + 
                     labs(title = "", col = "ID") + xlab("Time (Years)") + ylab("FEV1") + ylim(0, 6000)
        
        png('PLOTSP_CohortName_sample100_2023.png', width=1000, height=850)
            print(p2)
        dev.off()


       }

