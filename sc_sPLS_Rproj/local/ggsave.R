library(ggplot2)

# Print ----

Sys.sleep(3)

# Add the plot ----

gg_box(met_all, aesX="anno", aesY="rhat", x="", y="Methylation Rate", global=global)
ggsave('../img/ggsave-coordflip.png', dpi=600)