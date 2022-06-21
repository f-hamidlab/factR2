.testplot <- function(object, x, group.by = NULL, set = NULL){
  dat <- features(object, set = set)
  ggplot(dat, aes_string(x=x, fill = group.by)) + 
    geom_bar() +
    theme_minimal() + 
    scale_fill_brewer(palette = "Set2")
}