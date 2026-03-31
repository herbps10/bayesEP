assign_groups_to_sites <- function(C, K) {
  groups_per_site <- ceiling(C / K)
  lapply(1:K, function(w) {
    start <- (w - 1) * groups_per_site + 1
    end <- min(w * groups_per_site, C)
    if (start <= C) start:end else integer(0)
  })
}
