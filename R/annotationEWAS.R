#' @title Download and cache chip annotation tables
#' @description Download CpG annotation tables for supported Illumina chip types and
#' store them in a local cache directory. The downloaded files are saved as `.rds`
#' and reused in later analyses.
#' @usage downloadAnnotEWAS(chipType = c("EPICV2", "EPICV1", "450K", "27K", "MSA"),
#' cache_dir = NULL, force = FALSE, base_url = getOption("easyEWAS.annotation_base_url",
#' "https://github.com/ytwangZero/easyEWAS_materials/raw/main/annotation"), quiet = FALSE)
#' @param chipType One or more chip types. Supported values: `"EPICV2"`, `"EPICV1"`,
#' `"450K"`, `"27K"`, and `"MSA"`.
#' @param cache_dir Directory to store downloaded annotation files. If `NULL`, the
#' default cache path `tools::R_user_dir("easyEWAS", "cache")/annotation` is used.
#' @param force Logical. If `TRUE`, re-download and overwrite existing cached files.
#' @param base_url Base URL where annotation `.rds` files are hosted.
#' @param quiet Logical. Passed to `utils::download.file()`.
#' @return A named character vector of local file paths (invisibly).
#' @export
downloadAnnotEWAS <- function(
    chipType = c("EPICV2", "EPICV1", "450K", "27K", "MSA"),
    cache_dir = NULL,
    force = FALSE,
    base_url = getOption(
      "easyEWAS.annotation_base_url",
      "https://github.com/ytwangZero/easyEWAS_materials/raw/main/annotation"
    ),
    quiet = FALSE
) {
  valid_chip <- c("EPICV2", "EPICV1", "450K", "27K", "MSA")
  chipType <- unique(chipType)
  invalid_chip <- setdiff(chipType, valid_chip)
  if (length(invalid_chip) > 0) {
    stop("Unsupported chipType: ", paste(invalid_chip, collapse = ", "),
         ". Supported values are: ", paste(valid_chip, collapse = ", "), ".")
  }

  cache_dir <- .easyEWAS_annotation_cache_dir(cache_dir)
  manifest <- .easyEWAS_annotation_manifest(chipType = chipType, base_url = base_url)
  out_paths <- character(nrow(manifest))
  names(out_paths) <- manifest$chipType

  for (i in seq_len(nrow(manifest))) {
    chip <- manifest$chipType[i]
    url <- manifest$url[i]
    dest <- file.path(cache_dir, manifest$filename[i])

    if (file.exists(dest) && !isTRUE(force)) {
      message("Annotation already cached for ", chip, ": ", dest)
      out_paths[i] <- dest
      next
    }

    tmp <- tempfile(fileext = ".rds")
    tryCatch(
      utils::download.file(url = url, destfile = tmp, mode = "wb", quiet = quiet),
      error = function(e) {
        unlink(tmp)
        stop(
          "Failed to download annotation for ", chip, " from: ", url, "\n",
          "You can set a custom URL with `base_url=` or `options(easyEWAS.annotation_base_url = ...)`.\n",
          "Original error: ", e$message
        )
      }
    )

    annotation <- tryCatch(
      readRDS(tmp),
      error = function(e) {
        unlink(tmp)
        stop("Downloaded file for ", chip, " is not a valid .rds file. Error: ", e$message)
      }
    )

    annotation <- .easyEWAS_normalize_annotation(annotation, chipType = chip)
    saveRDS(annotation, file = dest, version = 2)
    unlink(tmp)
    out_paths[i] <- dest
    message("Annotation downloaded for ", chip, ": ", dest)
  }

  invisible(out_paths)
}

.easyEWAS_annotation_manifest <- function(chipType = NULL, base_url = getOption(
  "easyEWAS.annotation_base_url",
  "https://github.com/ytwangZero/easyEWAS_materials/raw/main/annotation"
)) {
  base_url <- sub("/+$", "", base_url)
  manifest <- data.frame(
    chipType = c("EPICV2", "EPICV1", "450K", "27K", "MSA"),
    filename = c("annotationV2.rds", "annotationV1.rds", "annotation450K.rds",
                 "annotation27K.rds", "annotationMSA.rds"),
    stringsAsFactors = FALSE
  )
  manifest$url <- paste0(base_url, "/", manifest$filename)
  if (is.null(chipType)) {
    return(manifest)
  }
  manifest[manifest$chipType %in% chipType, , drop = FALSE]
}

.easyEWAS_annotation_cache_dir <- function(cache_dir = NULL) {
  if (is.null(cache_dir) || !nzchar(cache_dir)) {
    cache_dir <- file.path(tools::R_user_dir("easyEWAS", which = "cache"), "annotation")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir
}

.easyEWAS_normalize_annotation <- function(annotation, chipType) {
  annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
  required_cols <- if (identical(chipType, "27K")) {
    c("probe", "chr", "pos", "gene")
  } else {
    c("probe", "chr", "pos", "relation_to_island", "gene", "location")
  }

  if (ncol(annotation) == length(required_cols)) {
    colnames(annotation) <- required_cols
  }

  alias <- list(
    probe = c("probe", "Probe", "Name", "IlmnID", "ID"),
    chr = c("chr", "CHR", "Chromosome", "chromosome"),
    pos = c("pos", "MAPINFO", "Position", "start"),
    relation_to_island = c("relation_to_island", "Relation_to_Island", "CpG_Island"),
    gene = c("gene", "UCSC_RefGene_Name", "Symbol", "Gene"),
    location = c("location", "UCSC_RefGene_Group", "Gene_Group", "Location")
  )

  for (nm in names(alias)) {
    if (nm %in% colnames(annotation)) {
      next
    }
    hit <- intersect(alias[[nm]], colnames(annotation))
    if (length(hit) > 0) {
      colnames(annotation)[match(hit[1], colnames(annotation))] <- nm
    }
  }

  missing_cols <- setdiff(required_cols, colnames(annotation))
  if (length(missing_cols) > 0) {
    stop(
      "Annotation data for chipType '", chipType, "' is missing required columns: ",
      paste(missing_cols, collapse = ", "), "."
    )
  }

  annotation[, required_cols, drop = FALSE]
}

.easyEWAS_load_annotation <- function(
    chipType,
    cache_dir = NULL,
    auto_download = FALSE,
    base_url = getOption(
      "easyEWAS.annotation_base_url",
      "https://github.com/ytwangZero/easyEWAS_materials/raw/main/annotation"
    )
) {
  manifest <- .easyEWAS_annotation_manifest(chipType = chipType, base_url = base_url)
  if (nrow(manifest) != 1) {
    stop("Unsupported chipType: ", chipType)
  }

  cache_dir <- .easyEWAS_annotation_cache_dir(cache_dir)
  cache_file <- file.path(cache_dir, manifest$filename)
  if (!file.exists(cache_file) && isTRUE(auto_download)) {
    downloadAnnotEWAS(
      chipType = chipType,
      cache_dir = cache_dir,
      force = FALSE,
      base_url = base_url,
      quiet = TRUE
    )
  }

  if (!file.exists(cache_file)) {
    return(NULL)
  }

  annotation <- tryCatch(
    readRDS(cache_file),
    error = function(e) {
      stop("Failed to read cached annotation file: ", cache_file, ". Error: ", e$message)
    }
  )
  .easyEWAS_normalize_annotation(annotation, chipType = chipType)
}
