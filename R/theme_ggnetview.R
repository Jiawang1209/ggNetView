#' ggNetView Theme
#'
#' @param base_size Base text size. Default 12.
#' @param base_family Base font family. Default NULL.
#' @param title_face Title font face (e.g., "bold"). Default "bold".
#' @param title_size Title font size. Default 18.
#' @param title_hjust Title horizontal justification (0–1). Default 0.5 (center).
#' @param subtitle_size Subtitle size. Default 12.
#' @param caption_size Caption size. Default 9.
#' @param background Plot background color; NULL for transparent. Default "white".
#' @param foreground Foreground color for strip background/border; NULL to skip.
#' @param border Logical; draw panel border with `foreground` color if TRUE. Default FALSE.
#' @param plot_margin Outer plot margin (ggplot2::margin). Default margin(10,10,10,10).
#' @param grid "none","x","y","both" to toggle major grid lines. Default "none".
#'
#' @returns A ggplot2 theme object.
#' @export
#'
#' @examples
#' library(ggplot2)
#' library(ggNetView)
#' ggplot(mtcars, aes(wt, mpg, color = hp)) +
#' geom_point(size = 2) +
#' ggtitle("Example of ggNetView Theme") +
#' theme_ggnetview()
#'
theme_ggnetview <- function(
    base_size    = 12,
    base_family  = NULL,
    title_face   = "bold",
    title_size   = 18,
    title_hjust  = 0.5,
    subtitle_size= 12,
    caption_size = 9,
    background   = "white",
    foreground   = NULL,
    border       = FALSE,
    plot_margin  = ggplot2::margin(10, 10, 10, 10),
    grid         = c("none","x","y","both")
){
  grid <- match.arg(grid)

  th <- ggraph::theme_graph(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(face = title_face, size = title_size, hjust = title_hjust),
      plot.subtitle= ggplot2::element_text(size = subtitle_size),
      plot.caption = ggplot2::element_text(size = caption_size),
      plot.margin  = plot_margin
    )

  # 背景：NULL 透明；否则上色
  if (is.null(background)) {
    th <- th + ggplot2::theme(plot.background = ggplot2::element_blank())
  }else{
    th <- th + ggplot2::theme(plot.background = ggplot2::element_rect(fill = background, colour = NA))
  }

  # 条带背景（分面）与边框
  if (!is.null(foreground)) {
    th <- th + ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = foreground, colour = foreground)
    )
  }
  if (isTRUE(border)) {
    th <- th + ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, colour = if (is.null(foreground)) "grey70" else foreground, linewidth = 0.5)
    )
  }

  # 网格线简易开关（theme_graph 默认很干净，这里按需补）
  th <- th + switch(
    grid,
    none = ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                          panel.grid.minor = ggplot2::element_blank()),
    x    = ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour = "grey85", linewidth = 0.25),
                          panel.grid.major.y = ggplot2::element_blank(),
                          panel.grid.minor   = ggplot2::element_blank()),
    y    = ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "grey85", linewidth = 0.25),
                          panel.grid.major.x = ggplot2::element_blank(),
                          panel.grid.minor   = ggplot2::element_blank()),
    both = ggplot2::theme(panel.grid.major   = ggplot2::element_line(colour = "grey85", linewidth = 0.25),
                          panel.grid.minor   = ggplot2::element_blank())
  )

  return(th)
}
