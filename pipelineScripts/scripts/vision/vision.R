library(feather)
library(data.table)
library(VISION)

REPO_BASE <- system2("git", c("rev-parse",  "--show-toplevel"), stdout = TRUE)

exp_file <- snakemake@input[["exp"]]
out_file <- snakemake@output[["out"]]

dir.create(dirname(out_file), showWarnings = FALSE)

sigs <- c(
    paste0(REPO_BASE, "/Signatures/Yoseflab/signatures_NY_private.gmt"),
    paste0(REPO_BASE, "/Signatures/Yoseflab/netPath.gmt"),
    paste0(REPO_BASE, "/Signatures/MSIGDB/H_Hallmark.gmt"),
    paste0(REPO_BASE, "/Signatures/MSIGDB/c7_subsets/TCells.gmt"),
    paste0(REPO_BASE, "/Signatures/CellType/ciberSortSigs.gmt"),
    paste0(REPO_BASE, "/Signatures/CellCycle/whitfield_cell_cycle.gmt")
)

read_feather_matrix <- function(feather_file) {
    data <- read_feather(feather_file)
    data <- as.data.frame(data)
    rownames(data) <- data$index
    data$index <- NULL
    data <- data.matrix(data)
    return(data)
}

if (endsWith(exp_file, ".feather")){
    expression <- read_feather_matrix(exp_file)
} else {
    fr_input <- exp_file
    if (endsWith(fr_input, ".gz")){
        fr_input <- paste0("zcat ", fr_input)
    }
    expression <- fread(input = fr_input, sep = "\t",
                        header = TRUE, data.table = FALSE)

    rownames(expression) <- expression[, 1]
    expression[, 1] <- NULL
    expression <- as.matrix(expression)
}

if ("unnorm" %in% names(snakemake@input)){
    unnorm_file <- snakemake@input[["unnorm"]]
    if (endsWith(unnorm_file, ".feather")){
        unnorm <- read_feather_matrix(unnorm_file)
    } else {
        fr_input <- unnorm_file
        if (endsWith(fr_input, ".gz")){
            fr_input <- paste0("zcat ", fr_input)
        }
        unnorm <- fread(input = fr_input, sep = "\t",
                            header = TRUE, data.table = FALSE)

        rownames(unnorm) <- unnorm[, 1]
        unnorm[, 1] <- NULL
        unnorm <- as.matrix(unnorm)
    }
} else {
    unnorm <- NULL
}

pc <- data.frame(row.names = colnames(expression))

if ("meta" %in% names(snakemake@input)) {
    meta_file <- snakemake@input[["meta"]]
    meta <- read.delim(meta_file,
                       check.names = FALSE,
                       row.names = 1)

    pc <- merge(pc, meta, by = "row.names")
    row.names(pc) <- pc$Row.names
    pc$Row.names <- NULL
}

if ("qc" %in% names(snakemake@input)) {
    qc_file <- snakemake@input[["qc"]]
    qc <- read.delim(qc_file,
                    check.names = FALSE,
                    row.names = 1)

    pc <- merge(pc, qc, by = "row.names")
    row.names(pc) <- pc$Row.names
    pc$Row.names <- NULL
}

# Load in the clustering
if ("clusters" %in% names(snakemake@input)) {
    clusters_file <- snakemake@input[["clusters"]]
    clusters <- read.delim(clusters_file,
                    check.names = FALSE,
                    row.names = 1)

    pc <- merge(pc, clusters, by = "row.names")
    row.names(pc) <- pc$Row.names
    pc$Row.names <- NULL
}


# Load latent space coordinates
if ("latent" %in% names(snakemake@input)){
    latent_file <- snakemake@input[["latent"]]
    latentSpace <- read.delim(latent_file, check.names = FALSE, row.names = 1)
} else {
    latentSpace <- NULL
}

# Name the results
wd <- getwd()
name <- basename(dirname(wd))

projection_methods <- c("tSNE30")
if ("tsne" %in% names(snakemake@input)) {
    projection_methods <- character()
}

fp <- Vision(expression, signatures = sigs, meta = pc,
            projection_genes = rownames(expression),
            name = name,
            unnormalizedData = unnorm,
            projection_methods = projection_methods,
            latentSpace = latentSpace, pool = FALSE)

# Load projection coordinates
if ("tsne" %in% names(snakemake@input)) {
    proj_file <- snakemake@input[["tsne"]]
    tsne <- read.delim(proj_file, check.names = FALSE, row.names = 1)
    fp <- addProjection(fp, "xTSNE", tsne)
}

options(mc.cores = 10)
fp <- analyze(fp)

saveRDS(fp, out_file)
