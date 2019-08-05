ui <- fluidPage(
  fluidRow(
    column(4,
           selectInput("x", "X-axis",
                       choices = colnames(sample_stats), selected = "PC1"),
           selectInput("y", "Y-axis",
                       choices = colnames(sample_stats), selected = "PC2"),
           selectInput("col", "Color",
                       choices = colnames(sample_stats), selected = "cpg_mean"),
           selectInput("size", "Size",
                       choices = colnames(sample_stats), selected = "rnaseq_libsize" ),
           textInput("caption", "Caption", ""),
           textInput("legend_col", "Color Legend", "Mean CpG Methylation (%)"),
           textInput("legend_size", "Size Legend", "scRNAseq Library Size"),
           textInput("xlab", "X Label", "Transcriptome PC 1"),
           textInput("ylab", "Y Label", "Transcriptome PC 2"),
           textInput("col_low", "Color Low:", "red"),
           textInput("col_high", "Color High", "green")
    ),
    column(8,
           plotOutput("ggplot", height=600)
    )
  )
)

server = function(input, output) {
  output$ggplot = renderPlot({
    ggplot(sample_stats, aes_string(x=input$x, y=input$y, col=input$col, size=input$size)) + 
      geom_point() + scale_color_gradient(low=input$col_low, high=input$col_high, 
                                          guide = guide_legend(title=input$legend_col)) +
      guides(size=guide_legend(title = input$legend_size)) + theme_bw() +
      labs(x=input$xlab, y=input$ylab, title=input$caption)  +
      geom_text(aes(label=rownames(sample_stats)), vjust=-1, col="grey50", hjust=1, size=3, show.legend = FALSE)
    
  })
}

shinyApp(ui, server)