# chord
# Written by Nanbo Sun, Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# Edited by Sid Chopra


make_chord <- function(data=NULL,
                       link_colorscale="bkr",
                       min_thre=NULL,
                       max_thre=NULL,
                       output_png=FALSE,
                       outname=NULL) {

  if (!require(circlize)){
    install.packages('circlize')
  }
  if (!require(igraph)){
    install.packages('igraph')
  }
  if (!require(ComplexHeatmap)){
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install("ComplexHeatmap")
  }
  if (!require(gridExtra)){
    install.packages('gridExtra')
  }


  library(circlize)
  library(igraph)
  library(ComplexHeatmap)
  library(gridExtra)

  # load data
  #data <- read.csv(mat18, header=FALSE)
  if (is.null(min_thre)) {min_thre <- min(data)}
  if (is.null(max_thre)) {max_thre <- max(data)}


  networks = c('1','2','3','4','5','6','7','8','9','10')
  colnames(data) <- as.numeric(networks)
  rownames(data) <- as.numeric(networks)
  d <- data.matrix(data)
  g <- graph.adjacency(d, mode="undirected",weighted=TRUE)
  df <- get.data.frame(g)

  ###
  # define parameters
  ###
  # outer track parameters

  outer_labels = c('Visual', 'SomMot','DorsAttn','VentAttn', 'Limbic', 'Frontoparietal', 'Default', 'Subcortex')
  #outer_color = c('#0C30FF', '#cd3e4e','#e69422','#dcf8a4','#c43afa','#00760e','#4682b4','#781286','#BF9959')
  outer_color = rep('#000000', 9)
  outer_textcol = '#7B7D7D'
  outer_textcol = '#000000'
  outer_fontsize = 1.5

  # inner track parameters
  inner_labels = c('','','','','','','','MTL','Thal','Stri')
  inner_color = c("#781286", "#4682b4", "#00760e", "#c43afa",
                  "#dcf8a4",  "#e69422", "#cd3e4e", "#5B7C96",
                  "#CADBC8", "#46706B")
  inner_textcol = rep('#000000', 10)
  inner_textcol[2] = '#FFFFFF'
  inner_textcol[6] = '#FFFFFF'
  inner_textcol[17] = '#FFFFFF'
  inner_fontsize = 1

  ###
  # output the figure as eps
  ###
  #setEPS()
  #postscript(paste(outname,'.eps',sep=''), width=10, height=10)

  # start circos plot
  lg = 5
  gap_vec = c(lg,lg,lg,lg,lg,lg,lg,lg,lg,lg)
  circos.clear()
  circos.par("start.degree"= 215,"gap.after"= gap_vec, unit.circle.segments= 10000, canvas.ylim=c(-1.2,1.2))

  f = factor(networks, levels = networks)
  circos.initialize(factors = f, xlim = c(0, 1))
  # 1st track
  circos.track(factors = f, ylim = c(0, 1), "track.height"= 0.02,
               track.margin = c(0.01, 0), bg.border='white', cell.padding=c(0,1,0,1))
  highlight.sector(c('1'), track.index = 1, col = outer_color[1], text = outer_labels[1],
                   text.col = outer_textcol, text.vjust = "7mm", font = 2, cex=outer_fontsize, niceFacing=TRUE,  facing = "bending.inside")
  highlight.sector(c('2'), track.index = 1, col = outer_color[2], text = outer_labels[2],
                   text.col = outer_textcol, text.vjust = "7mm", font = 2, cex=outer_fontsize, niceFacing=TRUE,facing = "bending.inside")
  highlight.sector(c('3'), track.index = 1, col = outer_color[3], text = outer_labels[3],
                   text.col = outer_textcol, text.vjust = "7mm", font = 2, cex=outer_fontsize, niceFacing=TRUE,facing = "bending.inside")
  highlight.sector(c('4'), track.index = 1, col = outer_color[4], text = outer_labels[4],
                   text.col = outer_textcol, text.vjust = "7mm", font = 2, cex=outer_fontsize, niceFacing=TRUE,facing = "bending.inside")
  highlight.sector(c('5'), track.index = 1, col = outer_color[5], text = outer_labels[5],
                   text.col = outer_textcol, text.vjust = "7mm", font = 2, cex=outer_fontsize, niceFacing=TRUE,facing = "bending.inside")
  highlight.sector(c('6'), track.index = 1, col = outer_color[6], text = outer_labels[6],
                   text.col = outer_textcol, text.vjust = "7mm", font = 2, cex=outer_fontsize, niceFacing=TRUE,facing = "bending.inside")
  highlight.sector(c('7'), track.index = 1, col = outer_color[7], text = outer_labels[7],
                   text.col = outer_textcol, text.vjust = "7mm", font = 2, cex=outer_fontsize, niceFacing=TRUE,facing = "bending.inside")
  highlight.sector(c('8', '9', '10'), track.index = 1, col = outer_color[8], text = outer_labels[8],
                   text.col = outer_textcol, text.vjust = "7mm", font = 2, cex=outer_fontsize, niceFacing=TRUE,facing = "bending.inside")

  # 2nd track
  circos.track(factors = f, ylim = c(0, 1), "track.height"= 0.1, bg.col = inner_color, bg.lwd = rep(2, 10))
  circos.trackText(f, rep(0.5, 10), rep(0.5, 10), inner_labels, track.index = 2, facing="bending",
                   niceFacing = TRUE, font=2, cex=inner_fontsize, col=inner_textcol)

  # links
  ind = !(df$weight < min_thre & df$weight > -min_thre)
  df = df[ind,]
  df$weight[df$weight > max_thre] = max_thre
  df$weight[df$weight < -max_thre] = -max_thre+1e-10
  resolution <- 80
  if (min_thre == 0){
    breaks = c(seq(-max_thre, max_thre, length.out = resolution))
  } else {
    breaks = c(seq(-max_thre,-min_thre,length.out = resolution/2),
               seq(min_thre, max_thre, length.out = resolution/2))
  }
  df$color_level <- cut(df$weight, breaks, labels=FALSE)
  # generate bwr color
  # bwr <- colorRampPalette(c("blue","white","red"))
  # bwr_colors = bwr(resolution)

  if (link_colorscale == 'bwr'){
    bwr <- colorRamp2(c(-1, 0, 1), c("blue","white","red"), space='RGB')
    link_colors = bwr(seq(-1, 1, length.out=resolution))
  }

  if  (link_colorscale == 'bkr'){
    link_colors = scan("~/Dropbox/Sid/python_files/PredictingCognition/scripts/visualisation/functions/bkr_colorscale.txt", what="character", sep="\n")
    bkr = function(vector){
      return(link_colors[vector])
    }
  }

  if (link_colorscale == 'whitered'){
    bwr <- colorRamp2(c(0, 1), c("white","red"), space='RGB')
    link_colors = bwr(seq(0, 1, length.out=resolution))
  }

  if(link_colorscale != 'bkr' & link_colorscale != 'bwr' &  link_colorscale != 'whitered'){
    hcl_scale = colorRamp2(c(0, 1), space='RGB',hcl_palette = link_colorscale)
    link_colors = hcl_scale(seq(0, 1, length.out=resolution))
  }

  len = max(table(c(df$from,df$to))) + 1.1
  link_breaks <- seq(0, 1, length.out=len)
  count = rep(1, 10)
  for (row in 1:nrow(df)) {
    from <- df[row, "from"]
    to <- df[row, "to"]
    color_level <- df[row, "color_level"]
    if (from == to) {
      from_num = as.numeric(from)
      circos.link(from, link_breaks[1:2], to, link_breaks[(len-1):len], col=link_colors[color_level])
    } else {
      from_num = as.numeric(from)
      to_num = as.numeric(to)
      count[from_num] <- count[from_num] + 1
      count[to_num] <- count[to_num] + 1
      #print(link_breaks[count[from_num]:(count[from_num]+1)])
      circos.link(from, link_breaks[count[from_num]:(count[from_num]+1)], to,
                  link_breaks[count[to_num]:(count[to_num]+1)], col=link_colors[color_level])
    }
  }
  grid.arrange(nullGrob(), newpage = FALSE)
  ## add color bar
  #if (link_colorscale == 'bwr'){
  #  color_bar = Legend(at=c(-1,-0.2,0.2,1), labels=c(-max_thre,-min_thre,min_thre,max_thre),
  #                     col_fun=bwr, grid_height = unit(4, "mm"), grid_width = unit(6, "mm"),
  #                     labels_gp = gpar(fontsize = 20), direction = "horizontal")
  #}
  #if (link_colorscale == 'bkr'){
  #  color_bar = Legend(at=c(1,30,51,80), labels=c(-max_thre,-min_thre,min_thre,max_thre),
  #                     col_fun=bkr, grid_height = unit(4, "mm"), grid_width = unit(6, "mm"),
  #                     labels_gp = gpar(fontsize = 20), direction = "horizontal")
  #}

  #grid.arrange(nullGrob(), bottom=color_bar, newpage = FALSE)

}
#dev.off()

###
# output the figure as png
###
#png(paste(outname,'.png',sep=''), width=10, height=10, units="in", res=300)
