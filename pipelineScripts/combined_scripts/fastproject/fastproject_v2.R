library(feather)
library(data.table)
library(VISION)

REPO_BASE <- system2("git", c("rev-parse",  "--show-toplevel"), stdout = TRUE)

exp_file <- snakemake@input[["exp"]]
meta_file <- snakemake@input[["meta"]]
out_file <- snakemake@output[["out"]]

out_dir <- dirname(out_file)
dir.create(out_dir, showWarnings = FALSE)

sigs <- c(
    paste0(REPO_BASE, "/Signatures/Yoseflab/misc.gmt"),
    paste0(REPO_BASE, "/Signatures/Yoseflab/exhaustion.gmt"),
    paste0(REPO_BASE, "/Signatures/MSIGDB/H_Hallmark.gmt"),
    paste0(REPO_BASE, "/Signatures/MSIGDB/c7_subsets/TCells_noCD8.gmt"),
    paste0(REPO_BASE, "/Signatures/CellType/ciberSortSigs.gmt"),
    paste0(REPO_BASE, "/Signatures/CellCycle/whitfield_cell_cycle.gmt"),
    paste0(REPO_BASE, "/Signatures/IBD_Colitis/Venema_Voskuil_2018.gmt"),
    paste0(REPO_BASE, "/Signatures/Smillie2019/InflamedSignatures.gmt")
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

pc <- read.delim(meta_file,
                check.names = FALSE,
                row.names = 1)

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

# Load any file that starts with 'meta_' as more metadata
for (input_name in names(snakemake@input)) {
    if (startsWith(input_name, "meta_")) {
        meta_file <- snakemake@input[[input_name]]
        newMeta <- read.delim(meta_file, check.names = FALSE, row.names = 1)
        pc <- merge(pc, newMeta, by = "row.names")
        row.names(pc) <- pc$Row.names
        pc$Row.names <- NULL
    }
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
name <- basename(wd)

projection_methods <- c("tSNE30", "UMAP")
if ("tsne" %in% names(snakemake@input) || "umap" %in% names(snakemake@input) ||
    any(startsWith(names(snakemake@input), "proj_"))) {
    projection_methods <- character()
}

fp <- Vision(expression, signatures = sigs, meta = pc,
            projection_genes = rownames(expression),
            name = name,
            unnormalizedData = unnorm,
            projection_methods = projection_methods,
            latentSpace = latentSpace, pool = FALSE,
            #sig_norm_method = "znorm_rows_then_columns"
            sig_norm_method = "znorm_columns"
)

# Load projection coordinates
if ("tsne" %in% names(snakemake@input)) {
    proj_file <- snakemake@input[["tsne"]]
    tsne <- read.delim(proj_file, check.names = FALSE, row.names = 1)
    fp <- addProjection(fp, "TSNE", tsne)
}

if ("umap" %in% names(snakemake@input)) {
    proj_file <- snakemake@input[["umap"]]
    umap <- read.delim(proj_file, check.names = FALSE, row.names = 1)
    fp <- addProjection(fp, "UMAP", umap)
}

# Load any file that starts with 'proj_' as a projection
for (input_name in names(snakemake@input)) {
    if (startsWith(input_name, "proj_")) {
        name <- gsub("proj_", "", input_name)
        proj_file <- snakemake@input[[input_name]]
        coords <- read.delim(proj_file, check.names = FALSE, row.names = 1)
        fp <- addProjection(fp, name, coords)
    }
}

options(mc.cores = 10)
fp <- analyze(fp)

saveRDS(fp, out_file)
