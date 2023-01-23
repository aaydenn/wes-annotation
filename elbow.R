# get elbow----
get_elbow_points_indices <- function(x, y, jump) {
  d1 <- diff(y) / diff(x) # first derivative
  d2 <- diff(d1) / diff(x[-1]) # second derivative
  indices <- which(abs(d2) > jump)
  return(indices)
}

get_elbow_points <- function(x, y, iter, jump) {
  ap <- approx(x,
               y,
               n = iter,
               yleft = min(y),
               yright = max(y))
  x_i <- ap$x
  y_i <- ap$y
  indices <- get_elbow_points_indices(x_i, y_i, jump)
  return(tibble(x = x_i[indices], y = y_i[indices]))
}
#----

# get density basad elbow----
get_delbow_points <- function(x, iter, probs = 0.05) {
  counts <- hist(x,breaks = iter, plot = FALSE)$counts
  breaks <- hist(x, breaks = iter, plot = FALSE)$breaks
  names(counts) <- breaks[-length(breaks)]
  peak_indx <- c(F,diff(sign(c(diff(counts))))==-2,F) %>% which()
  topcounts <- counts[peak_indx]
  return(quantile(as.numeric(names(topcounts)), probs = probs))
}
nrow(dt1)/10
#----

freq_table <- function(table) {
  tibble(x = round(sort(as.numeric(
    table$Lab_frequency
  )), 2),
  y = 1:nrow(table))
}

freq_list <- lapply(input_list, freq_table)

points_list <- lapply(freq_list, function(tbl) {get_elbow_points(tbl$x, tbl$y, 10000, 10000)} )

thresholds <- lapply(points_list, function(tbl) {min(tbl$x) + 1} )

p1 <- ggplot(data = freq_list$`MG22-5563-5564`, aes(x, y)) +

  geom_point() +
  geom_point(data = points_list$`MG22-5563-5564`, aes(x, y), colour = "red") +
  geom_text_repel(data = points_list$`MG22-5563-5564`, aes(x, y, label = x)) +
  labs(title = "MG22-5563-5564",
       x = "Frequency",
       y = "Number of variants")

p2 <- ggplot(data = freq_list$`MG22-5795`, aes(x, y)) +
  
  geom_point() +
  geom_point(data = points_list$`MG22-5795`, aes(x, y), colour = "red") +
  geom_text_repel(data = points_list$`MG22-5795`, aes(x, y, label = x)) +
  labs(title = "MG22-5795",
       x = "Frequency",
       y = "Number of variants")

p1 |> ggsave(filename = "MG22-5563-5564_elbow.png")
p2 |> ggsave(filename = "MG22-5795_elbow.png")


ggplot(freq_list$`MG22-5563-5564`, aes(x,1)) + geom_violin() + geom_vline(xintercept = 3.1)
ggplot(freq_list$`MG22-5795`, aes(x,1)) + geom_violin() + geom_vline(xintercept = 3.5)
