#' Generate a named color palette for discrete classes
#'
#' @param classes Character string.
#' The discrete class names or factor levels to map to colors.
#' @param others_label Character, (default = "Others").
#'
#' @returns A named character vector where names correspond to
#' @export
#'
#' @examples NULL
#'
get_palette <- function(classes, others_label = "Others") {
  base_colors <- c(
    '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
    '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
    '#ccebc5','#ffed6f','#a6cee3','#b2df8a','#fb9a99',
    '#bdbdbd','#1f78b4','#33a02c','#e31a1c','#fdbf6f',
    '#ff7f00','#6a3d9a','#ffff99','#b15928'
  )

  classes <- as.character(classes)
  uniq_classes <- unique(stats::na.omit(classes))

  has_others <- others_label %in% uniq_classes
  uniq_main  <- setdiff(uniq_classes, others_label)

  color_map <- setNames(
    rep(base_colors, length.out = length(uniq_main)),
    uniq_main
  )

  if (has_others) {
    color_map <- c(color_map, setNames("#bdbdbd", others_label))
  }
  color_map
}


#' Custom discrete color scale for ggNetView
#'
#' @param classes Character string.
#' The discrete class names or factor levels to map to colors.
#' @param ... Additional arguments passed to `ggplot2::scale_color_manual()`.
#' @param others_label Character, (default = "Others").
#' @param na_value Color for missing values. Default `"#e0e0e0"`.
#' @param drop Logical, passed to `ggplot2::scale_color_manual()`.
#'
#' @returns A `ggplot2` scale object.
#' @export
#'
#' @examples NULL
#'
scale_color_ggnetview <- function(classes,
                                  ...,
                                  others_label = "Others",
                                  na_value = "#e0e0e0",
                                  drop = FALSE) {
  pal <- get_palette(classes, others_label = others_label)
  ggplot2::scale_color_manual(values = pal, na.value = na_value, drop = drop, ...)
}

#' Custom discrete fill scale for ggNetView
#'
#' @param classes Character string.
#' The discrete class names or factor levels to map to colors.
#' @param ... Additional arguments passed to `ggplot2::scale_color_manual()`.
#' @param others_label Character, (default = "Others").
#' @param na_value Color for missing values. Default `"#e0e0e0"`.
#' @param drop Logical, passed to `ggplot2::scale_color_manual()`.
#'
#' @returns A `ggplot2` scale object.
#' @export
#'
#' @examples NULL
scale_fill_ggnetview <- function(classes,
                                 ...,
                                 others_label = "Others",
                                 na_value = "#e0e0e0",
                                 drop = FALSE) {
  pal <- get_palette(classes, others_label = others_label)
  ggplot2::scale_fill_manual(values = pal, na.value = na_value, drop = drop, ...)
}
