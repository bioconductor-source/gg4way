#' @rdname gg4wayExtractors
#' @name extractors
#' @title Helper Functions for \code{gg4way}
#' @param p1 The plot from \link{gg4way}
#' @description These helper functions provide data used in the plot:
#'  \tabular{ll}{
#'   \code{getCor} \tab Get the correlation of the logFC of all genes \cr
#'   \tab \cr
#'   \code{getShared} \tab Get only the shared genes that pass the thresholds \cr
#'   \tab \cr
#'   \code{getTotals} \tab Get the totals of overlap categories \cr
#'   }
#' @return Each function returns a different result:
#'  \tabular{ll}{
#'   \code{getCor} \tab A numeric \cr
#'   \tab \cr
#'   \code{getShared} \tab A \link[tibble:tbl_df-class]{tibble} \cr
#'   \tab \cr
#'   \code{getTotals} \tab A \link[janitor]{tabyl} \cr
#'   }
#' @examples
#' data("airwayFit")
#' p1 <- airwayFit |>
#'     gg4way(x = "N61311 vs N052611",
#'            y = "N061011 vs N052611")
#'
#' ## Correlation
#' getCor(p1)
#'
#' ## Shared
#' getShared(p1)
#'
#' ## Totals
#' getTotals(p1)
#'
NULL

#' @rdname gg4wayExtractors
#' @importFrom methods is
#' @export
#'
getCor <- function(p1){
    stopifnot(is(p1, "gg4way"))

    p1$data |>
        .testCor() |>
        round(digits = 2)
}

#' @rdname gg4wayExtractors
#' @importFrom methods is
#' @importFrom dplyr select
#' @export
#'
getShared <- function(p1){
    stopifnot(is(p1, "gg4way"))

    p1$data |>
        dplyr::filter(Significant == "Significant in Both")
}

#' @rdname gg4wayExtractors
#' @importFrom methods is
#' @importFrom dplyr select
#' @export
#'
getTotals <- function(p1){
    stopifnot(is(p1, "gg4way"))

    p1[["layers"]][[7]][["data"]] |>
        dplyr::select(countGroup, n, color, Significant)
}
