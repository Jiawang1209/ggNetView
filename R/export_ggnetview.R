#' Save a ggNetView with sensible defaults
#'
#' @param p Plot to save, defaults to last plot displayed.
#' @param filename 	File name to create on disk.
#' @param height Plot size in units expressed by the units argument. If not supplied, uses the size of the current graphics device.
#' @param width Plot size in units expressed by the units argument. If not supplied, uses the size of the current graphics device.
#' @param limitsize When TRUE (the default), export_ggnetview() will not save images larger than 50x50 inches, to prevent the common error of specifying dimensions in pixels.
#' @param dpi Plot resolution. Also accepts a string input: "retina" (320), "print" (300), or "screen" (72). Only applies when converting pixel units, as is typical for raster output types.
#'
#' @returns a graph output
#' @export
#'
#' @examples NULL
export_ggnetview <- function(p, filename, height, width, limitsize = F, dpi = 300){
  ggplot2::ggsave(filename = filename,
                  plot = p,
                  height = height,
                  width = width,
                  limitsize = limitsize,
                  dpi = dpi
                  )
}
