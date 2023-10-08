#' @title Create example data
#' @description Load and process example data for \link[gg4way]{tidyDGE}
#' @keywords internal
#' @return A \link[limma:MArrayLM-class]{MArrayLM}
#' @import airway
#' @import SummarizedExperiment
#' @import org.Hs.eg.db
#' @importFrom edgeR SE2DGEList filterByExpr calcNormFactors
#' @importFrom limma makeContrasts voomWithQualityWeights
#'  lmFit contrasts.fit eBayes
#' @importFrom stats model.matrix
#' @import utils
#' @importFrom AnnotationDbi mapIds
#'
gg4way_example_data <- function() {
    data_env <- new.env(parent = emptyenv())
    data("airway", package = "airway", envir = data_env)
    se <- data_env[["airway"]]

    geneSymbols <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                         keys = rownames(se),
                                         column = "SYMBOL",
                                         keytype = "ENSEMBL")

    SummarizedExperiment::rowData(se)$symbol <- geneSymbols

    SummarizedExperiment::rowData(se)$ID <- rownames(se)

    se <- se[!is.na(SummarizedExperiment::rowData(se)$symbol)]

    dge <- se |>
        edgeR::SE2DGEList()

    design <- model.matrix(~ 0 + cell + dex, data = dge$samples)
    colnames(design) <- gsub("cell", "", colnames(design))

    contr.matrix <- limma::makeContrasts(N61311 - N052611,
                                         N061011 - N052611,
                                         levels = c("N052611", "N061011",
                                                    "N080611", "N61311",
                                                    "dexuntrt"))

    keep <- edgeR::filterByExpr(dge, design)
    dge <- dge[keep, ]

    efit <- dge |>
        edgeR::calcNormFactors() |>
        limma::voom(design) |>
        limma::lmFit(design) |>
        limma::contrasts.fit(contrasts = contr.matrix) |>
        limma::eBayes()

    return(efit)
}

airwayFit <- gg4way_example_data()

usethis::use_data(airwayFit, overwrite = TRUE)
