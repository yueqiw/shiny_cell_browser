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
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by,
  plot.legend = FALSE,
  do.return = FALSE,
  x.lab.rot = FALSE
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
      avg.exp = mean(expm1(x = expression)),
      pct.exp = PercentAbove(x = expression, threshold = 0)
    )
  data.to.plot <- data.to.plot %>%
    ungroup() %>%
    group_by(genes.plot) %>%
    mutate(avg.exp.scale = scale(x = avg.exp)) %>%
    mutate(avg.exp.scale = MinMax(
      data = avg.exp.scale,
      max = col.max,
      min = col.min
    )) %>%
    mutate(avg.exp.scale.one = avg.exp / max(avg.exp))


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

  if (color_scaling == "mean-var") {
    p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
        geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
        scale_func(limits = c(0.001,1), range = c(0, dot.scale)) +
        scale_color_gradient(low = cols.use[1], high = cols.use[2])
  } else if (color_scaling == "zero-one") {
    p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
        geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale.one)) +
        scale_func(limits = c(0.001,1), range = c(0, dot.scale)) +
        scale_color_gradient(low = cols.use[1], high = cols.use[2], limits = c(0,1))
  } else {
      stop("color scaling error.")
  }
  p <- p +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())

  if (! plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}
