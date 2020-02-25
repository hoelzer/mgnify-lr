process ideel {
    publishDir "${params.output}/${name}/ideel", mode: 'copy', pattern: "${name}_ideel.pdf"
    label 'ggplot2'

  input:
    tuple val(name), file(results) 

  output:
	  tuple val(name), file("${name}_ideel.pdf") 
      
  script:
    """
    #!/usr/bin/Rscript
    library(ggplot2)  
    data <- read.table("${results}", header=FALSE, sep='\\t')
    names(data) <- c('qlen', 'slen')
    pseudogenes <- sum(data\$qlen / data\$slen < 0.9)
    print(paste0('Encountered genes < 0.9 reference length: ', pseudogenes))
    theme_min = function (
      size=10, font=NA, face='plain',
      panelColor=backgroundColor, axisColor='#999999',
      gridColor=gridLinesColor, textColor='black')
    {
      theme_text = function(...)
      ggplot2::theme_text(family=font, face=face, colour=textColor, size=size, ...)
      opts(
          axis.text.x = theme_text(),
          axis.text.y = theme_text(),
          axis.ticks = theme_segment(colour=axisColor, size=0.25),
          panel.border = theme_rect(colour=backgroundColor),
          legend.background = theme_blank(),
          legend.key = theme_blank(),
          legend.key.size = unit(1.5, 'lines'),
          legend.text = theme_text(hjust=0),
          legend.title = theme_text(hjust=0),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = theme_rect(fill=NA, colour=NA),
          strip.text.x = theme_text(hjust=0),
          strip.text.y = theme_text(angle=-90),
          plot.title = theme_text(hjust=0),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'lines'))
      }
      # PLOT
        p <- ggplot(data, aes(x=qlen/slen)) +
          geom_histogram(fill='white', color='grey25', bins=20) +
          xlab('query length / hit length') +
          ylab('frequency') +
          scale_x_continuous(limits=c(0, 1.3)) +
          theme_minimal()
      ggsave("${name}_ideel.pdf", p, height=7, width=7, units='cm')
    """
}