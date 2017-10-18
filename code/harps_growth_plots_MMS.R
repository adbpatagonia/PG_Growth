setwd('D:/Buren_files/DFO/jobs/RES 2016/ppt')
load('harps_growth_plots.Rdata')
library(cowplot)

febg <- subset(grmoddataw,month == 2 & decade == 1980 | decade == 1990  | decade == 2000)
febg90 <- subset(grmoddataw,month == 2 & decade == 1990)
febg00 <- subset(grmoddataw,month == 2 & decade == 2000)



vb <- ggplot(febg, aes(x = age, y = bw, colour = as.factor(decade)))
vb <- vb + geom_point(alpha = 0.5) 
vb <- vb + theme_set(theme_cowplot())
vb <- vb + labs(x = 'Age', y = 'Weight (kg)')
vb <- vb + scale_colour_discrete(name = "",
                               breaks = c("1980", "1990", "2000"),
                               labels = c("1980s", "1990s", "2000s"))
vb <- vb + theme(legend.justification = c(1,1), legend.position = c(.98,0.35))
vb <- vb + geom_line(data = MyData2.feb80, aes(x = age, y = mean), col = 'red', size = 1.5)
vb <- vb + geom_line(data = MyData2.feb90, aes(x = age, y = mean), col = 'darkgreen', size = 1.5)
vb <- vb + geom_line(data = MyData2.feb00, aes(x = age, y = mean), col = 'blue', size = 1.5)
vb <- vb + geom_ribbon(data = MyData2.feb00, aes(x = age, y = mean, ymin = lb,  ymax = ub), alpha = 0.2, fill = 'blue', color = NA)
vb <- vb + geom_ribbon(data = MyData2.feb80, aes(x = age, y = mean, ymin = lb,  ymax = ub), alpha = 0.2, fill = 'red', color = NA)
vb <- vb + geom_ribbon(data = MyData2.feb90, aes(x = age, y = mean, ymin = lb,  ymax = ub), alpha = 0.2, fill = 'darkgreen', color = NA)
vb


aw <- ggplot(data = as.data.frame(aw.feb80.mcmc), aes(x = V1))
aw <- aw + geom_density(alpha = 0.1, col = 'red', size = 1.5)
aw <- aw + geom_density(data = as.data.frame(aw.feb90.mcmc), aes(x = V1), alpha = 0.1, col = 'darkgreen', size = 1.5)
aw <- aw + geom_density(data = as.data.frame(aw.feb00.mcmc), aes(x = V1), alpha = 0.1, col = 'blue', size = 1.5)
aw <- aw + labs(x = 'Asymptotic weight (kg)', y = 'Density')
aw


gr <- ggplot(data = as.data.frame(gr.feb80.mcmc), aes(x = V1))
gr <- gr + geom_density(alpha = 0.1, col = 'red', size = 1.5)
gr <- gr + geom_density(data = as.data.frame(gr.feb90.mcmc), aes(x = V1), alpha = 0.1, col = 'darkgreen', size = 1.5)
gr <- gr + geom_density(data = as.data.frame(gr.feb00.mcmc), aes(x = V1), alpha = 0.1, col = 'blue', size = 1.5)
gr <- gr + labs(x = 'Growth rate', y = 'Density')
gr


bw <- ggplot(data = as.data.frame(bw.feb80.mcmc), aes(x = V1))
bw <- bw + geom_density(alpha = 0.1, col = 'red', size = 1.5)
bw <- bw + geom_density(data = as.data.frame(bw.feb90.mcmc), aes(x = V1), alpha = 0.1, col = 'darkgreen', size = 1.5)
bw <- bw + geom_density(data = as.data.frame(bw.feb00.mcmc), aes(x = V1), alpha = 0.1, col = 'blue', size = 1.5)
bw <- bw + labs(x = 'Weight at weaning (kg)', y = 'Density')
bw

setwd('D:/Buren_files/GitHub/PG_Growth/')

save_plot("output/toGarry/growth-feb-decades-trans.png", vb, base_aspect_ratio = 1.4,  base_width = 6, bg = "transparent") # make room for figure legend)
save_plot("output/toGarry/growth-feb-decades-white.png", vb, base_aspect_ratio = 1.4,  base_width = 6, bg = "white") # make room for figure legend)

save_plot("output/toGarry/growthrate-feb-decades-trans.png", gr, base_height = 6,  base_width = 6, bg = "transparent") # make room for figure legend)
save_plot("output/toGarry/growthrate-feb-decades-white.png", gr, base_height = 6,  base_width = 6, bg = "white") # make room for figure legend)

save_plot("output/toGarry/asympweight-feb-decades-trans.png", aw, base_aspect_ratio = 1,  base_width = 6, bg = "transparent") # make room for figure legend)
save_plot("output/toGarry/asympweight-feb-decades-white.png", aw, base_aspect_ratio = 1,  base_width = 6, bg = "white") # make room for figure legend)

save_plot("output/toGarry/weightwean-feb-decades-trans.png", bw, base_aspect_ratio = 1,  base_width = 6, bg = "transparent") # make room for figure legend)
save_plot("output/toGarry/weightwean-feb-decades-white.png", bw, base_aspect_ratio = 1,  base_width = 6, bg = "white") # make room for figure legend)
