.testplot <- function(object, x, y = "count",
                      x.set = NULL, y.set = NULL,
                      geom = "auto",
                      x.lab = NULL, y.lab = NULL,
                      x.lim = NULL, y.lim = NULL,
                      group.by = NULL, group.pos = "dodge",
                      min.exp = 5, y.trans = NULL, x.rot = NULL){

  # check if a set var is given for x
    if(grepl("\\$", x)){
        split <- str_split(x, "\\$")[[1]]
        x.set <- split[1]
        x <- split[2]
    }
    if(grepl("\\$", y)){
        split <- str_split(y, "\\$")[[1]]
        y.set <- split[1]
        y <- split[2]
    }

  # get data
  if(x== "sample"){
      if(is.null(y.set)){
          y.set = object@active.set
      }
    
    var <- paste0(y.set, "_id")

      # get expression
      counts <- object@sets[[y.set]]@data %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var) %>%
          tidyr::pivot_longer(-var, names_to = "sample", values_to = "exp")

      
      # get features
      dat <- features(object, set = y.set)
      
      # join expression to features and filter
      # for min expression
      dat <- counts %>% 
        dplyr::left_join(dat, by = var) %>%
        filter(exp >= min.exp)
  } else {
      dat <- features(object, set = x.set)
  }


  # if no grouping is given, set x as groupings
  if(is.null(group.by)){
      group.by = x
  }


  # base plot
  plot <- ggplot(dat, aes_string(x=x, fill = group.by, color = group.by)) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2")

  # plot y-axis
  ## get geometry
  if(geom == "auto"){
    if(y %in% colnames(dat)){
      if(class(dat[[y]]) %in% c("character","factor")) {
        geom.plot <- "geom_count(aes_string(y=y))"
      } else{
        geom.plot <- "geom_violin(aes_string(y=y))"
      }
    }
  } else {
    geom.plot <- sprintf("geom_%s(aes_string(y=y))", geom)
  }

  if(y == "percent"){
      plot <- plot +
          geom_bar(position = "fill") +
          scale_y_continuous(labels = scales::percent) +
          ylab("percent")
  } else if(y == "count"){
      plot <- plot +
          geom_bar(position = group.pos)
  } else if(y %in% colnames(dat)){
    plot <- plot +
      eval(parse(text = geom.plot))


    # # plot cont variables from dat
    # if(!class(dat[[y]]) %in% c("character","factor")){
    #
    #     # geom_bar(aes_string(y = y), stat = "summary", fun = "mean") +
    #     # stat_summary(aes_string(y=y),
    #     #              geom = "errorbar",
    #     #              width = 1,
    #     #              fun = mean,
    #     #              fun.min = function(x) mean(x) - (sd(x)/sqrt(length(x))),
    #     #              fun.max = function(x) mean(x) + (sd(x)/sqrt(length(x))))
    # }


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
  if(!is.null(y.trans)){
    plot <- plot + coord_trans(y = y.trans)
  }

      plot
}
