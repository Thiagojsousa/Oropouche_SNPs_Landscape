# SNP CATALOGUE + CODON EFFECT (Aligned FASTA → SNP table + Syn/Nonsyn/Stop) ----
#
# OBJECTIVE
#   Given a multiple sequence alignment (MSA) FASTA containing one reference
#   sequence and multiple samples (e.g., Seq_01, Seq_02, ...), this script:
#     1) catalogues SNPs for each sample relative to the reference (1-based
#        coordinates on the gap-free reference; reference gaps are ignored);
#     2) annotates coding impact assuming the entire segment is a single CDS
#        (one gene per segment: L, M, S), classifying the effect as:
#        synonymous / nonsynonymous / stop_gained / stop_lost / unknown;
#     3) exports one table per segment (L, M, S) and a merged (LMS) file.
#
# EXPECTED INPUTS
#   - Aligned (DNA) FASTA files located in PARAM$workdir:
#       L_aln.fasta
#       M_aln.fasta
#       S_aln_N.fasta
#     Each file must contain:
#       - 1 reference sequence (header must match PARAM$segments$*$ref_id exactly)
#       - ≥1 sample sequence(s) (free headers, e.g., Seq, Seq_02, ...)
#     Note: the alignment must have all sequences with the same length
#     (valid MSA).
#
# OUTPUTS
#   - One CSV per segment (PARAM$segments$*$out_csv) with columns:
#       CHROM, POS, REF, ALT, REGION, STRAND, CDS_START, CDS_END,
#       CDS_POS, CODON_POS, CODON_REF, CODON_ALT, AA_POS, AA_REF, AA_ALT,
#       EFFECT, SAMPLE
#
#   - A consolidated CSV (PARAM$merged_out_csv) concatenating L+M+S.
#
# KEY DEFINITIONS
#   - CHROM: segment label (L/M/S).
#   - POS: 1-based position on the gap-free reference (alignment column is
#          converted to reference coordinates).
#   - REF/ALT: reference and sample bases observed at POS
#              (A/C/G/T only, by default).
#   - CDS_START/CDS_END: CDS coordinates on the reference (default: 1..len).
#   - CDS_POS: SNP position within the CDS (1-based, in the translation frame).
#   - CODON_POS: SNP position within the codon (1, 2, 3).
#   - CODON_REF/CODON_ALT: codon triplets before/after the SNP (DNA).
#   - AA_POS: affected amino-acid position (1-based).
#   - AA_REF/AA_ALT: amino acids (one-letter code) before/after the SNP.
#   - EFFECT:
#       synonymous     → AA_REF == AA_ALT
#       nonsynonymous  → AA_REF != AA_ALT (and not STOP)
#       stop_gained    → AA_ALT == "*"
#       stop_lost      → AA_REF == "*" and AA_ALT != "*"
#       unknown        → ambiguous cases (N, gaps, outside CDS, etc.)
#
# FILTERS / LIMITATIONS (default)
#   - Only considers SNPs where REF and ALT are in {A, C, G, T}.
#   - Ignores indels/gaps ("-") if PARAM$ignore_indels = TRUE.
#   - Does not interpret ambiguous IUPAC codes (N, R, Y, etc.) as variants.
#   - Assumes each segment is a single gene (one CDS). If UTRs are present, set
#     PARAM$cds_mode = "manual" and adjust cds_start/cds_end per segment.
#
# 1) PARAMETERS / OPTIONS (EDIT ONLY HERE) ----

PARAM <- list(
  
  ## 1.1) Working directory (where FASTAs are located and where outputs will be written) ----
  workdir = "/home/user/folder_path",
  output_dir = "/home/user/folder_path/output_files",
  
  ## 1.2) Segments: aligned FASTA, reference header, strand, and output ----
  segments = list(
    L = list(
      CHROM     = "L",
      fasta_aln = "L_aln.fasta",
      ref_id    = "PP154172.1",
      STRAND    = "+",
      out_csv   = "L_snps_catalogue.csv"
    ),
    M = list(
      CHROM     = "M",
      fasta_aln = "M_aln.fasta",
      ref_id    = "PP154171.1",
      STRAND    = "+",
      out_csv   = "M_snps_catalogue.csv"
    ),
    S = list(
      CHROM     = "S",
      fasta_aln = "S_aln.fasta",
      ref_id    = "PP154170.1",
      STRAND    = "+",
      out_csv   = "S_snps_catalogue.csv"
    )
  ),
  
  ## 1.3) CDS definition (each segment treated as a single gene): by default, CDS = full segment ----
  cds_mode = "full_segment",   # "full_segment" or "manual"
  cds_manual = list(
    L = list(cds_start = 1L, cds_end = NA_integer_),  
    M = list(cds_start = 1L, cds_end = NA_integer_),
    S = list(cds_start = 1L, cds_end = NA_integer_)
  ),
  
  ## 1.4) Variant filters ----
  valid_bases   = c("A","C","G","T"),
  ignore_indels = TRUE,   # TRUE: ignore '-' in the alignment (SNPs only)
  
  # 1.5) Effect labels (standardized output naming)
  effect_labels = list(
    synonymous     = "synonymous",
    nonsynonymous  = "nonsynonymous",
    stop_gained    = "stop_gained",
    stop_lost      = "stop_lost",
    unknown        = "unknown",
    noncoding      = "noncoding"
  ),
  
 ## 1.6) Merged output (L+M+S) ----
  merged_out_csv = "LMS_snps_catalogue.csv",
  
  ## 1.7) Verbose ----
  verbose = TRUE
)

# 2) PACKAGES AND CODON DICTIONARY ----

if (!requireNamespace("seqinr", quietly = TRUE)) install.packages("seqinr")
if (!requireNamespace("readr",  quietly = TRUE)) install.packages("readr")
if (!requireNamespace("dplyr",  quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")

library(seqinr)
library(readr)
library(dplyr)
library(tibble)

setwd(PARAM$workdir)
if (!dir.exists(PARAM$output_dir)) dir.create(PARAM$output_dir, recursive = TRUE)

# Build the segment_jobs list from PARAM (do not edit below)
segment_jobs <- lapply(PARAM$segments, function(s) {
  list(
    CHROM     = s$CHROM,
    fasta_aln = s$fasta_aln,
    ref_id    = s$ref_id,
    out_csv   = s$out_csv,
    STRAND    = s$STRAND
  )
})

# Shortcuts for options (do not edit below)
valid_bases   <- PARAM$valid_bases
ignore_indels <- PARAM$ignore_indels
verbose       <- PARAM$verbose


## 2.1) INTERNAL CODON DICTIONARY (DNA codon -> 1-letter amino acid) ----

CODON_AA <- c(
  # Phenylalanine / Leucine / Serine / Tyrosine / STOP / Cysteine / Tryptophan
  "TTT"="F","TTC"="F","TTA"="L","TTG"="L",
  "TCT"="S","TCC"="S","TCA"="S","TCG"="S",
  "TAT"="Y","TAC"="Y",
  "TAA"="*","TAG"="*","TGA"="*",
  "TGT"="C","TGC"="C","TGG"="W",
  
  # Leucine / Proline / Histidine / Glutamine / Arginine
  "CTT"="L","CTC"="L","CTA"="L","CTG"="L",
  "CCT"="P","CCC"="P","CCA"="P","CCG"="P",
  "CAT"="H","CAC"="H",
  "CAA"="Q","CAG"="Q",
  "CGT"="R","CGC"="R","CGA"="R","CGG"="R",
  
  # Isoleucine / Methionine / Threonine / Asparagine / Lysine / Serine / Arginine
  "ATT"="I","ATC"="I","ATA"="I",
  "ATG"="M",
  "ACT"="T","ACC"="T","ACA"="T","ACG"="T",
  "AAT"="N","AAC"="N",
  "AAA"="K","AAG"="K",
  "AGT"="S","AGC"="S",
  "AGA"="R","AGG"="R",
  
  # Valine / Alanine / Aspartic acid / Glutamic acid / Glycine
  "GTT"="V","GTC"="V","GTA"="V","GTG"="V",
  "GCT"="A","GCC"="A","GCA"="A","GCG"="A",
  "GAT"="D","GAC"="D",
  "GAA"="E","GAG"="E",
  "GGT"="G","GGC"="G","GGA"="G","GGG"="G"
)

## 3) HELPER FUNCTIONS ----

# 3.1) Reverse-complement base
comp_base <- function(x) {
  dplyr::recode(x,
                "A"="T","T"="A","C"="G","G"="C",
                .default = NA_character_
  )
}

## 3.2) Read aligned FASTA (returns a named character vector: name -> string) ----
read_aln_fasta_as_strings <- function(fasta_file) {
  fa <- seqinr::read.fasta(fasta_file, as.string = FALSE, seqtype = "DNA", forceDNAtolower = FALSE)
  seq_names <- names(fa)
  seq_str <- vapply(fa, function(x) paste0(toupper(x), collapse = ""), character(1))
  names(seq_str) <- seq_names
  seq_str
}

## 3.3) Alignment-column -> reference POS map (1-based, gap-free reference coordinates) ----
build_ref_pos_map <- function(ref_aln_string) {
  ref_chars <- strsplit(ref_aln_string, "")[[1]]
  pos_map <- rep(NA_integer_, length(ref_chars))
  pos <- 0L
  for (i in seq_along(ref_chars)) {
    if (ref_chars[i] != "-") {
      pos <- pos + 1L
      pos_map[i] <- pos
    }
  }
  list(ref_chars = ref_chars, pos_map = pos_map)
}

## 3.4) Codon dictionary (internal only) ----
# Uses the standard genetic code (DNA; STOP="*"). To modify the code table,
# edit CODON_AA directly in section 2.1.
CODON_AA[c("TAA","TAG","TGA")] <- "*"
if (any(is.na(CODON_AA[c("TAA","TAG","TGA")]))) {
  stop("STOP codons are not defined in CODON_AA.")
}

## 3.5) Annotate codon effect for a SNP within a single CDS ----
annotate_codon_effect <- function(POS, ALT, cds_start, cds_end, strand,
                                  ref_ungapped_chars, CODON_AA) {
  
  # CDS position in translation direction
  if (strand == "+") {
    cds_pos <- POS - cds_start + 1L
  } else {
    cds_pos <- cds_end - POS + 1L
  }
  
  codon_pos <- ((cds_pos - 1L) %% 3L) + 1L
  codon_start_cdspos <- cds_pos - (codon_pos - 1L)
  
  if (strand == "+") {
    g1 <- cds_start + (codon_start_cdspos - 1L)
    g2 <- g1 + 1L
    g3 <- g1 + 2L
    if (g1 < cds_start || g3 > cds_end) {
      return(list(CDS_POS=cds_pos, CODON_POS=codon_pos,
                  CODON_REF=NA, CODON_ALT=NA,
                  AA_POS=ceiling(cds_pos/3), AA_REF=NA, AA_ALT=NA,
                  EFFECT="unknown"))
    }
    trip <- c(ref_ungapped_chars[g1], ref_ungapped_chars[g2], ref_ungapped_chars[g3])
    codon_ref <- paste0(trip, collapse = "")
    
    trip_alt <- trip
    trip_alt[codon_pos] <- ALT
    codon_alt <- paste0(trip_alt, collapse = "")
    
  } else { # strand == "-"
    g1 <- cds_end - (codon_start_cdspos - 1L)
    g2 <- cds_end - (codon_start_cdspos)
    g3 <- cds_end - (codon_start_cdspos + 1L)
    
    if (g3 < cds_start || g1 > cds_end) {
      return(list(CDS_POS=cds_pos, CODON_POS=codon_pos,
                  CODON_REF=NA, CODON_ALT=NA,
                  AA_POS=ceiling(cds_pos/3), AA_REF=NA, AA_ALT=NA,
                  EFFECT="unknown"))
    }
    
    trip_gen <- c(ref_ungapped_chars[g1], ref_ungapped_chars[g2], ref_ungapped_chars[g3])
    trip_cds <- comp_base(trip_gen)
    codon_ref <- paste0(trip_cds, collapse = "")
    
    trip_alt <- trip_cds
    trip_alt[codon_pos] <- comp_base(ALT)  # genomic ALT -> CDS base
    codon_alt <- paste0(trip_alt, collapse = "")
  }
  
  codon_ref <- toupper(gsub("U","T",codon_ref))
  codon_alt <- toupper(gsub("U","T",codon_alt))
  
  aa_ref <- unname(CODON_AA[codon_ref])
  aa_alt <- unname(CODON_AA[codon_alt])
  
  aa_pos <- ceiling(cds_pos / 3)
  
  if (!is.na(aa_ref) && !is.na(aa_alt)) {
    if (aa_ref == aa_alt) {
      effect <- PARAM$effect_labels$synonymous
    } else if (aa_alt == "*") {
      effect <- PARAM$effect_labels$stop_gained
    } else if (aa_ref == "*" && aa_alt != "*") {
      effect <- PARAM$effect_labels$stop_lost
    } else {
      effect <- PARAM$effect_labels$nonsynonymous
    }
  } else {
    effect <- PARAM$effect_labels$unknown
  }
  
  list(
    CDS_POS   = cds_pos,
    CODON_POS = codon_pos,
    CODON_REF = codon_ref,
    CODON_ALT = codon_alt,
    AA_POS    = aa_pos,
    AA_REF    = aa_ref,
    AA_ALT    = aa_alt,
    EFFECT    = effect
  )
}

## 4) RUN PER SEGMENT (L/M/S) ----

all_results <- list()

for (job in segment_jobs) {
  
  chrom     <- job$CHROM
  fasta_aln <- job$fasta_aln
  ref_id    <- job$ref_id
  out_csv   <- job$out_csv
  strand    <- job$STRAND
  
  if (verbose) message("\n=== Segment ", chrom, " | FASTA: ", fasta_aln, " | Ref: ", ref_id, " ===")
  
  seqs_aln <- read_aln_fasta_as_strings(fasta_aln)
  
  if (!(ref_id %in% names(seqs_aln))) {
    stop("ref_id '", ref_id, "' was not found in ", fasta_aln, ". Available headers:\n",
         paste(names(seqs_aln), collapse = "\n"))
  }
  
  # Check MSA validity
  aln_lens <- nchar(seqs_aln)
  if (length(unique(aln_lens)) != 1) stop("Invalid MSA in ", fasta_aln, ": sequences have different lengths.")
  
  ref_aln <- seqs_aln[[ref_id]]
  map <- build_ref_pos_map(ref_aln)
  
  # Gap-free reference defines the coordinate system
  ref_ungapped <- gsub("-", "", ref_aln)
  ref_ungapped_chars <- strsplit(ref_ungapped, "")[[1]]
  ref_len <- length(ref_ungapped_chars)
  
  # Single CDS: full segment by default
  if (PARAM$cds_mode == "full_segment") {
    cds_start <- 1L
    cds_end   <- ref_len
  } else if (PARAM$cds_mode == "manual") {
    cfg <- PARAM$cds_manual[[chrom]]
    cds_start <- as.integer(cfg$cds_start)
    cds_end   <- as.integer(cfg$cds_end)
    if (is.na(cds_end)) cds_end <- ref_len
  } else {
    stop("Invalid cds_mode: ", PARAM$cds_mode)
  }
  
  # Samples = all sequences except the reference
  sample_ids <- setdiff(names(seqs_aln), ref_id)
  
  extract_snps_for_sample <- function(sample_id) {
    sample_aln <- seqs_aln[[sample_id]]
    
    ref_chars <- map$ref_chars
    sam_chars <- strsplit(sample_aln, "")[[1]]
    pos_map   <- map$pos_map
    
    keep <- !is.na(pos_map) &
      ref_chars %in% valid_bases &
      sam_chars %in% valid_bases &
      (sam_chars != ref_chars)
    
    if (ignore_indels) keep <- keep & (ref_chars != "-") & (sam_chars != "-")
    
    tibble(
      CHROM  = chrom,
      POS    = pos_map[keep],
      REF    = ref_chars[keep],
      ALT    = sam_chars[keep],
      SAMPLE = sample_id
    ) %>% arrange(POS)
  }
  
  snps_all <- bind_rows(lapply(sample_ids, extract_snps_for_sample))
  
  if (nrow(snps_all) == 0) {
    if (verbose) message("[", chrom, "] No SNPs found. Writing empty file: ", out_csv)
    out_csv_path <- file.path(PARAM$output_dir, out_csv)
    readr::write_csv(snps_all, out_csv_path)
    next
  }
  
  annotated <- snps_all %>%
    rowwise() %>%
    do({
      x <- .
      ann <- annotate_codon_effect(
        POS = x$POS,
        ALT = x$ALT,
        cds_start = cds_start,
        cds_end = cds_end,
        strand = strand,
        ref_ungapped_chars = ref_ungapped_chars,
        CODON_AA = CODON_AA
      )
      
      tibble(
        CHROM     = x$CHROM,
        POS       = x$POS,
        REF       = x$REF,
        ALT       = x$ALT,
        REGION    = "CDS",
        STRAND    = strand,
        CDS_START = cds_start,
        CDS_END   = cds_end,
        CDS_POS   = ann$CDS_POS,
        CODON_POS = ann$CODON_POS,
        CODON_REF = ann$CODON_REF,
        CODON_ALT = ann$CODON_ALT,
        AA_POS    = ann$AA_POS,
        AA_REF    = ann$AA_REF,
        AA_ALT    = ann$AA_ALT,
        EFFECT    = ann$EFFECT,
        SAMPLE    = x$SAMPLE
      )
    }) %>%
    ungroup() %>%
    arrange(SAMPLE, POS)
  
  all_results[[chrom]] <- annotated
  out_csv_path <- file.path(PARAM$output_dir, out_csv)
  readr::write_csv(annotated, out_csv_path)
  
  if (verbose) message("[", chrom, "] OK: ", out_csv_path)
  
  if (verbose) {
    message("[", chrom, "] OK: ", out_csv)
    summary_tab <- annotated %>%
      dplyr::group_by(SAMPLE, EFFECT) %>%
      dplyr::summarise(N = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(SAMPLE, dplyr::desc(N))
    
    print(summary_tab)
  }
}

## 5) MERGE L + M + S ----

merged_LMS <- dplyr::bind_rows(all_results) %>%
  dplyr::mutate(SEGMENT = CHROM) %>%
  dplyr::arrange(SAMPLE, SEGMENT, POS)

merged_out_path <- file.path(PARAM$output_dir, PARAM$merged_out_csv)
readr::write_csv(merged_LMS, merged_out_path)

summary_out_csv  <- sub("\\.csv$", "_summary.csv", PARAM$merged_out_csv)
summary_out_path <- file.path(PARAM$output_dir, summary_out_csv)
readr::write_csv(summary_LMS, summary_out_path)

if (verbose) {
  message("Merged catalogue: ", merged_out_path)
  message("Summary table:   ", summary_out_path)
}
