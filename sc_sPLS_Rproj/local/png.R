# Add the plot ----
library(ggplot2)
library(data.table)


gg_box(met_all[1:3e5,], aesX="anno", aesY="rhat", x="", y="Methylation Rate", global=global)
ggsave('methylation-rate-global-1.png')
