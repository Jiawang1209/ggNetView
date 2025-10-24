#' Translate TOM matrin to create graph object
#'
#' @param TOM matrix
#' TOM matrix from WGCNA result
#'
#' @param mat matrix
#' matrox to WGCNA analysis
#'
#' @returns A Data frame contain from, to and weight
#' @export
#'
#' @examples NULL
trans_TOM_in_WGCNA <- function(TOM, mat){
  TOM_mat <- as.matrix(TOM) %>%
    as.data.frame() %>%
    magrittr::set_colnames(colnames(mat)) %>%
    magrittr::set_rownames(colnames(mat)) %>%
    tibble::rownames_to_column(var = "from") %>%
    tidyr::pivot_longer(cols = -from, names_to = "to", values_to = "weight") %>%
    dplyr::filter(from != to) %>%
    dplyr::mutate(tmp = ifelse(from > to, str_c(from, to), str_c(to, from))) %>%
    dplyr::distinct(tmp, .keep_all = T) %>%
    dplyr::select(-tmp)

  return(TOM_mat)
}
