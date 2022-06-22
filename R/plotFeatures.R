.testplot <- function(object, x, set = NULL, x.lab = NULL, x.lim = NULL,
                      y.type = "count", y.lab = NULL,
                      group.by = NULL, group.pos = "dodge"){
    if(grepl("%", x)){
        split <- str_split(x, "%")[[1]]
        set <- split[1]
        x <- split[2]
    }
  dat <- features(object, set = set)
  if(is.null(group.by)){
      group.by = x
  }


  # base plot
  plot <- ggplot(dat, aes_string(x=x, fill = group.by)) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")

  # plot y-axis
  if(y.type == "percent"){
      plot <- plot +
          geom_bar(position = "fill") +
          scale_y_continuous(labels = scales::percent) +
          ylab("percent")
  } else{
      plot <- plot +
          geom_bar(position = group.pos)
  }

  # change axis labels
  if(!is.null(x.lab)){
      plot <- plot + xlab(x.lab)
  }
  if(!is.null(y.lab)){
      plot <- plot + ylab(y.lab)
  }

  # change x-axis limits
  if(!is.null(x.lim) & !class(dat[[x]]) %in% c("character","factor")){
      plot <- plot + xlim(x.lim)
  }

      plot
}
