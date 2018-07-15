library(tidyverse)
library(Seurat)


PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

DotPlot2 <- function(
  object,
  genes.plot,
  color_scaling = "zero-one", # or "mean-var"
  size_scaling = "area", # or "radius"
  cols.use = c("grey90", "red"),
  horizontal = FALSE,
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  axis.label.size = 15,
  group.by,
  legend.position = "none",
  do.return = FALSE
) {
  color_scaling_types <- c("mean-var", "zero-one")
  if (! color_scaling %in% color_scaling_types) {
      stop(paste("color scaling must be one of: ", paste(color_scaling_types, collapse = ", ")))
  }
  size_scaling_types <- c("area", "radius")
  if (! size_scaling %in% size_scaling_types) {
      stop(paste("size scaling must be one of: ", paste(size_scaling_types, collapse = ", ")))
  }

  if (! missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }

  if (horizontal) {
      genes.plot <- rev(genes.plot)
  }
  genes.plot <- genes.plot %>% unique()
  genes.plot.exist <- genes.plot[genes.plot %in% rownames(object@data)]
  genes.plot.nonexist <- genes.plot[!genes.plot %in% rownames(object@data)]
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot.exist))
  data.to.plot[sub("-", ".", genes.plot.nonexist)] <- 0.0
  data.to.plot <- data.to.plot[sub("-", ".", genes.plot)]

  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  data.to.plot %>% gather(
    key = genes.plot,
    value = expression,
    -c(cell, id)
  ) -> data.to.plot

  data.to.plot <- data.to.plot %>%
    group_by(id, genes.plot) %>%
    summarize(
      avg.exp.raw = mean(expm1(x = expression)),
      pct.exp = PercentAbove(x = expression, threshold = 0)
    )
  data.to.plot <- data.to.plot %>%
    ungroup() %>%
    group_by(genes.plot)

  if (color_scaling == "mean-var") {
      data.to.plot <- data.to.plot %>%
        mutate(avg.exp = scale(x = avg.exp.raw)) %>%
        mutate(avg.exp = MinMax(
          data = avg.exp,
          max = col.max,
          min = col.min
        ))
  } else if (color_scaling == "zero-one") {
      data.to.plot <- data.to.plot %>%
        mutate(avg.exp = avg.exp.raw / max(avg.exp.raw))
  }


  data.to.plot$genes.plot <- factor(
    x = data.to.plot$genes.plot,
    levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
  )
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA

  if (size_scaling == "area") {
      scale_func <- scale_size
  } else if (size_scaling == "radius"){
      scale_func <- scale_radius
  }

  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
      geom_point(mapping = aes(size = pct.exp, color = avg.exp)) +
      scale_func(limits = c(0.001,1), range = c(0, dot.scale))

  if (color_scaling == "mean-var") {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  } else if (color_scaling == "zero-one") {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2], limits = c(0,1))
  } else {
      stop("color scaling error.")
  }
  p <- p +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(face="bold", size=axis.label.size),
          axis.text.y = element_text(face="bold", size=axis.label.size),
          legend.position = legend.position)

  if (horizontal) {
    p <- p + scale_y_discrete(limits = rev(levels(data.to.plot$id)))
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  } else {
    p <- p + coord_flip()
  }

  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}
