write_and_return_path <- function(x) {
  path <- paste0("output/", deparse(substitute(x)), ".csv")
  write_csv(x, path)
  return(path)
}

write_gt_and_return_path <- function(x) {
	path <- paste0("output/", deparse(substitute(x)), ".html")
	gtsave(x, path)
	return(path)
}

write_svg_and_return_path <- function(x) {
	path <- paste0("output/", deparse(substitute(x)), ".svg")
	ggsave(plot=x, filename=path, device="svg")
	return(path)
}

write_png_and_return_path <- function(x) {
	path <- paste0("output/", deparse(substitute(x)), ".png")
	ggsave(plot=x, filename=path, device="png")
	return(path)
}

write_ggsurvplot_svg_and_return_path <- function(x) {
	#see https://github.com/kassambara/survminer/issues/152
	path <- paste0("output/", deparse(substitute(x)), ".svg")
	ggsave(plot=survminer:::.build_ggsurvplot(x), filename=path, device="svg")
	return(path)
}

write_ggsurvplot_png_and_return_path <- function(x) {
	#see https://github.com/kassambara/survminer/issues/152
	path <- paste0("output/", deparse(substitute(x)), ".png")
	ggsave(plot=survminer:::.build_ggsurvplot(x), filename=path, device="png")
	return(path)
}

write_parquet_and_return_path <- function(x) {
	path <- paste0("output/", deparse(substitute(x)), ".csv")
	write_parquet(x, path)
	return(path)
}

write_lines_and_return_path <- function(x) {
	path <- paste0("output/", deparse(substitute(x)), ".txt")
	write_lines(x, path)
	return(path)
}

save_and_return_path <- function(x) {
	path <- paste0("output/", deparse(substitute(x)), ".RDS")
	saveRDS(x, file=path)
	return(path)
}
