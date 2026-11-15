# =============================================================================
# config.R - the only place in this repository where paths are defined
#
# Run every script from this directory (the one containing config.R), e.g.
#     cd code
#     Rscript 05_paf_attribution.R
# The first line of each script loads this file, so paths do not have to be
# edited inside the analysis code.
#
# Data are not distributed with this repository (see README.md, section
# "Data").  Point GBDSTROKE_DATA_ROOT at a folder that mirrors the layout
# described below, or edit the default here.
# =============================================================================

# Directory holding this config file, i.e. the repository root.
CODE_DIR <- normalizePath(".", mustWork = FALSE)

# Parent of the repository, used only to build a convenient default data path.
PROJECT_ROOT <- normalizePath(file.path(CODE_DIR, ".."), mustWork = FALSE)

# -----------------------------------------------------------------------------
# Data root
#   Expected layout (mirrors the authors' working directory):
#     <DATA_ROOT>/rawdat/   GBD extracts, CHARLS data, intermediate results
#     <DATA_ROOT>/GWAS/     GWAS summary statistics used by the MR scripts
#     <DATA_ROOT>/B5_outputs/  outputs of 50_b5_integration.R
#   Change it with the environment variable GBDSTROKE_DATA_ROOT, or edit the
#   default below.
# -----------------------------------------------------------------------------
DATA_ROOT <- Sys.getenv("GBDSTROKE_DATA_ROOT",
                        unset = file.path(PROJECT_ROOT, "data"))
RAWDAT    <- file.path(DATA_ROOT, "rawdat")
GWAS_DIR  <- file.path(DATA_ROOT, "GWAS")

# OpenGWAS needs a personal token.  Set it in your environment before running
# the MR scripts (never commit a token):
#     Sys.setenv(OPENGWAS_JWT = "<your token>")   # or export OPENGWAS_JWT=...
if (!nzchar(Sys.getenv("OPENGWAS_JWT")))
  message("Note: OPENGWAS_JWT is not set; the MR scripts will not be able to ",
          "query OpenGWAS until a token is provided.")

if (!dir.exists(DATA_ROOT))
  message("Data root not found: ", DATA_ROOT,
          "\nSet GBDSTROKE_DATA_ROOT to the folder that contains rawdat/ and ",
          "GWAS/ (see README.md).")
