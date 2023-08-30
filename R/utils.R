#' @title Tidy DGE data for plotting
#' @description Tidy DGE data for \link[gg4way]{gg4way}
#' @param fit A \link[limma:MArrayLM-class]{MArrayLM} object from
#'  \link[limma]{eBayes}, a \link[DESeq2:DESeqDataSet-class]{DESeqDataSet} from
#'  \link[DESeq2]{DESeq}, or a list of
#'  \link[DESeq2:DESeqResults-class]{DESeqResults} from \link[DESeq2]{results}
#' @return A named list of \link[tibble:tbl_df-class]{tibbles}
#'  with at least the following columns
#'  from a \link[limma:MArrayLM-class]{MArrayLM} object:
#' \itemize{
#' \item \code{symbol} Gene symbol
#' \item \code{ID} Gene ID
#' \item \code{logFC} Log2FC
#' \item \code{adj.P.Val} Adjusted p-value
#' }
#' @importFrom methods is
#' @importFrom purrr set_names map
#' @importFrom limma topTable
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom stringr str_replace str_subset
#'  str_replace_all str_remove str_trim
#' @importFrom DESeq2 resultsNames
#' @importFrom magrittr %>%
#' @examples
#' data("airwayFit")
#' airwayFit |>
#'     tidyDGE()
#' @export
#'
tidyDGE <- function(fit) {
    ## magrittr pipe used due to need for unnamed placeholder
    if (methods::is(fit, "MArrayLM")) {
        res <- fit$contrasts %>%
            colnames() %>%
            purrr::set_names() %>%
            purrr::map(~ fit %>%
                         limma::topTable(sort.by = "P", n = Inf, coef = .x) %>%
                         tibble::rownames_to_column() %>%
                         tibble::as_tibble()) %>%
          purrr::set_names(names(.) %>% stringr::str_replace("[-]", "vs"))
    } else if (methods::is(fit, "DESeqDataSet")) {
        res <- fit %>%
            DESeq2::resultsNames() %>%
            stringr::str_subset("Intercept", negate = TRUE) %>%
            purrr::set_names() %>%
            purrr::map(~ fit %>%
                         DESeq2::results(name = .x)) %>%
            purrr::set_names(names(.) %>%
                               stringr::str_replace_all("_", " ") %>%
                               stringr::str_remove(stringr::word(., 1)) %>%
                               stringr::str_trim()) %>%
            purrr::map(~ .x %>%
                         as.data.frame() %>%
                         tibble::rownames_to_column("ID") %>%
                         tibble::as_tibble())
    } else if (all(unlist(lapply(fit, methods::is, "DESeqResults")))) {
        res <- fit %>%
          purrr::map(~ .x %>%
                       as.data.frame() %>%
                       tibble::rownames_to_column("ID") %>%
                       tibble::as_tibble())
    } else {
        stop("Format not suppourted")
    }

    return(res)
}

utils::globalVariables(c(
    ".", "contrast", "LogFC", "FDR", "FDRpass", "Direction",
    "countGroup", "Significant", "Direction1", "n",
    "color", "hjust", "vjust", "symbol", "atop",
    "N052611", "N061011", "N61311", "airway"))
