#!/usr/bin/env Rscript
# =====================================================================
# Pipeline para visualização de SNPs de Oropouche (LMS) (Parte 1 + Parte 2)
# Parte 1) Catálogo de SNPs + efeito em códons (sinônimo/não-sinônimo/STOP)

#     OBJETIVO: Dado um alinhamento de múltiplas sequências (MSA) FASTA contendo uma sequência de referência
# e múltiplas amostras (por exemplo, Seq_01, Seq_02, ...).
#     1) cataloga SNPs para cada amostra em relação à referência (coordenadas baseadas em 1
#        na referência sem lacunas; as lacunas de referência são ignoradas);
#     2) anota o impacto em região codificante assumindo que todo o segmento é um único CDS
#        (um gene por segmento: L, M, S), classificando o efeito como:
#        sinônimo / não sinônimo / stop_gained / stop_lost / desconhecido;
#     3) exporta uma tabela por segmento (L, M, S) and a merged (LMS) file.
#
# Part 2) Plot "SNPs landscape" (L/M/S) a partir do LMS_snps_catalogue.csv
#
# Uso (RStudio):
#   1) Abra este arquivo e ajuste apenas o BLOCO 1 (parâmetros).
#   2) Execute por blocos ou rode tudo de uma vez.
# =====================================================================

##### 1) PARÂMETROS (EDITE APENAS AQUI) ################################

##### 0) FLAGS (controle de execução) ###################################
# - RUN_PART1: gera/atualiza LMS_snps_catalogue.csv (Parte 1)
# - RUN_PART2: gera a figura a partir do LMS_snps_catalogue.csv (Parte 2)
# - FORCE_PART1: se TRUE, recalcula a Parte 1 mesmo se o output já existir

RUN_PART1  <- TRUE
RUN_PART2  <- TRUE
FORCE_PART1 <- FALSE

## 1.1) Parte 1 (catálogo de SNPs + efeitos em códons)
PARAM <- list(
  
  # Diretório de entrada contendo FASTAs alinhados (MSA) onde o script vai procurar:
  #   L_aln.fasta, M_aln.fasta, S_aln.fasta
  input_dir  = "/home/Users/workdir/examples/input",
  output_dir = "/home/Users/workdir/outputs",
  
  # Segmentos (FASTA alinhado + header exato da referência no FASTA)
  segments = list(
    L = list(CHROM="L", fasta_aln="L_aln.fasta", ref_id="PP154172.1", STRAND="+", out_csv="L_snps_catalogue.csv"),
    M = list(CHROM="M", fasta_aln="M_aln.fasta", ref_id="PP154171.1", STRAND="+", out_csv="M_snps_catalogue.csv"),
    S = list(CHROM="S", fasta_aln="S_aln.fasta", ref_id="PP154170.1", STRAND="+", out_csv="S_snps_catalogue.csv")
  ),
  
  # CDS: por padrão, trata o segmento inteiro como um único CDS (um gene por segmento)
  cds_mode   = "full_segment",   # "full_segment" ou "manual"
  cds_manual = list(
    L = list(cds_start = 1L, cds_end = NA_integer_),
    M = list(cds_start = 1L, cds_end = NA_integer_),
    S = list(cds_start = 1L, cds_end = NA_integer_)
  ),
  
  # filtros de variantes
  valid_bases   = c("A","C","G","T"),
  ignore_indels = TRUE,
  
  # rótulos de efeito
  effect_labels = list(
    synonymous     = "synonymous",
    nonsynonymous  = "nonsynonymous",
    stop_gained    = "stop_gained",
    stop_lost      = "stop_lost",
    unknown        = "unknown",
    noncoding      = "noncoding"
  ),
  
  # saída juntas L+M+S
  merged_out_csv = "LMS_snps_catalogue.csv",
  
  verbose = TRUE
)

## 1.2) Parte 2 (Plot)
PARAM_PLOT <- list(
  # IMPORTANTE: por padrão, a Parte 2 lê o arquivo produzido na Parte 1 dentro de output_dir.
  workdir          = PARAM$output_dir,
  datafr_file  = PARAM$merged_out_csv,   # "LMS_snps_catalogue.csv"
  out_dir_fig      = "figures",
  ncol_facets      = 2,
  width_in         = 10,
  height_in        = 4,
  dpi_png          = 600,
  
  # Prefixos esperados nos headers para ordenação (ex.: Node_01, Node_02, Seq_01 ...)
  header_prefixes  = c("Node_", "Seq_")
)

##### 2) PACOTES ########################################################

pkgs <- c("seqinr", "readr", "dplyr", "tibble", "ggplot2", "stringr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(seqinr)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)

if (!dir.exists(PARAM$output_dir)) dir.create(PARAM$output_dir, recursive = TRUE)

##### 3) PARTE 1 — CATÁLOGO DE SNPs + EFEITOS EM CÓDONS #####################

# (auto-skip) Se o arquivo final existir e FORCE_PART1=FALSE, a Parte 1 não é recalculada.
merged_out_path <- file.path(PARAM$output_dir, PARAM$merged_out_csv)

if (!RUN_PART1) {
  message("RUN_PART1=FALSE: pulando Parte 1 (catálogo de SNPs).")
} else if (file.exists(merged_out_path) && !FORCE_PART1) {
  message("Parte 1: encontrado '", merged_out_path, "'. Pulando recálculo (FORCE_PART1=FALSE).")
} else {
  
  setwd(PARAM$input_dir)
  
  # 3.1) Tabela de códons (DNA codon -> AA 1-letra; STOP="*")
  CODON_AA <- c(
    "TTT"="F","TTC"="F","TTA"="L","TTG"="L",
    "TCT"="S","TCC"="S","TCA"="S","TCG"="S",
    "TAT"="Y","TAC"="Y",
    "TAA"="*","TAG"="*","TGA"="*",
    "TGT"="C","TGC"="C","TGG"="W",
    
    "CTT"="L","CTC"="L","CTA"="L","CTG"="L",
    "CCT"="P","CCC"="P","CCA"="P","CCG"="P",
    "CAT"="H","CAC"="H",
    "CAA"="Q","CAG"="Q",
    "CGT"="R","CGC"="R","CGA"="R","CGG"="R",
    
    "ATT"="I","ATC"="I","ATA"="I",
    "ATG"="M",
    "ACT"="T","ACC"="T","ACA"="T","ACG"="T",
    "AAT"="N","AAC"="N",
    "AAA"="K","AAG"="K",
    "AGT"="S","AGC"="S",
    "AGA"="R","AGG"="R",
    
    "GTT"="V","GTC"="V","GTA"="V","GTG"="V",
    "GCT"="A","GCC"="A","GCA"="A","GCG"="A",
    "GAT"="D","GAC"="D",
    "GAA"="E","GAG"="E",
    "GGT"="G","GGC"="G","GGA"="G","GGG"="G"
  )
  
  # 3.2) Funções auxiliares
  comp_base <- function(x) {
    dplyr::recode(x, "A"="T","T"="A","C"="G","G"="C", .default = NA_character_)
  }
  
  read_aln_fasta_as_strings <- function(fasta_file) {
    fa <- seqinr::read.fasta(fasta_file, as.string = FALSE, seqtype = "DNA", forceDNAtolower = FALSE)
    seq_names <- names(fa)
    seq_str <- vapply(fa, function(x) paste0(toupper(x), collapse = ""), character(1))
    names(seq_str) <- seq_names
    seq_str
  }
  
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
  
  annotate_codon_effect <- function(POS, ALT, cds_start, cds_end, strand,
                                    ref_ungapped_chars, CODON_AA) {
    
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
                    EFFECT=PARAM$effect_labels$unknown))
      }
      trip <- c(ref_ungapped_chars[g1], ref_ungapped_chars[g2], ref_ungapped_chars[g3])
      codon_ref <- paste0(trip, collapse = "")
      
      trip_alt <- trip
      trip_alt[codon_pos] <- ALT
      codon_alt <- paste0(trip_alt, collapse = "")
      
    } else { # "-"
      g1 <- cds_end - (codon_start_cdspos - 1L)
      g2 <- cds_end - (codon_start_cdspos)
      g3 <- cds_end - (codon_start_cdspos + 1L)
      
      if (g3 < cds_start || g1 > cds_end) {
        return(list(CDS_POS=cds_pos, CODON_POS=codon_pos,
                    CODON_REF=NA, CODON_ALT=NA,
                    AA_POS=ceiling(cds_pos/3), AA_REF=NA, AA_ALT=NA,
                    EFFECT=PARAM$effect_labels$unknown))
      }
      
      trip_gen <- c(ref_ungapped_chars[g1], ref_ungapped_chars[g2], ref_ungapped_chars[g3])
      trip_cds <- comp_base(trip_gen)
      codon_ref <- paste0(trip_cds, collapse = "")
      
      trip_alt <- trip_cds
      trip_alt[codon_pos] <- comp_base(ALT)
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
  
  # 3.3) Montar segment_jobs
  segment_jobs <- lapply(PARAM$segments, function(s) {
    list(
      CHROM     = s$CHROM,
      fasta_aln = s$fasta_aln,
      ref_id    = s$ref_id,
      out_csv   = s$out_csv,
      STRAND    = s$STRAND
    )
  })
  
  valid_bases   <- PARAM$valid_bases
  ignore_indels <- PARAM$ignore_indels
  verbose       <- PARAM$verbose
  
  all_results <- list()
  
  for (job in segment_jobs) {
    
    chrom     <- job$CHROM
    fasta_aln <- job$fasta_aln
    ref_id    <- job$ref_id
    out_csv   <- job$out_csv
    strand    <- job$STRAND
    
    if (verbose) message("\n=== Segmento ", chrom, " | FASTA: ", fasta_aln, " | Ref: ", ref_id, " ===")
    
    fasta_path <- file.path(PARAM$input_dir, fasta_aln)
    seqs_aln <- read_aln_fasta_as_strings(fasta_path)
    
    if (!(ref_id %in% names(seqs_aln))) {
      stop("ref_id '", ref_id, "' was not found in ", fasta_path, ". Headers disponíveis:\n",
           paste(names(seqs_aln), collapse = "\n"))
    }
    
    aln_lens <- nchar(seqs_aln)
    if (length(unique(aln_lens)) != 1) stop("MSA inválido em ", fasta_path, ": sequences have different lengths.")
    
    ref_aln <- seqs_aln[[ref_id]]
    map <- build_ref_pos_map(ref_aln)
    
    ref_ungapped <- gsub("-", "", ref_aln)
    ref_ungapped_chars <- strsplit(ref_ungapped, "")[[1]]
    ref_len <- length(ref_ungapped_chars)
    
    if (PARAM$cds_mode == "full_segment") {
      cds_start <- 1L
      cds_end   <- ref_len
    } else if (PARAM$cds_mode == "manual") {
      cfg <- PARAM$cds_manual[[chrom]]
      cds_start <- as.integer(cfg$cds_start)
      cds_end   <- as.integer(cfg$cds_end)
      if (is.na(cds_end)) cds_end <- ref_len
    } else {
      stop("cds_mode inválido: ", PARAM$cds_mode)
    }
    
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
    
    out_csv_path <- file.path(PARAM$output_dir, out_csv)
    
    if (nrow(snps_all) == 0) {
      if (verbose) message("[", chrom, "] Nenhum SNP encontrado. Escrevendo arquivo vazio: ", out_csv_path)
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
    readr::write_csv(annotated, out_csv_path)
    
    if (verbose) message("[", chrom, "] OK: ", out_csv_path)
    
    if (verbose) {
      summary_tab <- annotated %>%
        dplyr::group_by(SAMPLE, EFFECT) %>%
        dplyr::summarise(N = dplyr::n(), .groups = "drop") %>%
        dplyr::arrange(SAMPLE, dplyr::desc(N))
      print(summary_tab)
    }
  }
  
  # 3.4) Merge L+M+S
  merged_LMS <- dplyr::bind_rows(all_results) %>%
    dplyr::mutate(SEGMENT = CHROM) %>%
    dplyr::arrange(SAMPLE, SEGMENT, POS)
  
  merged_out_path <- file.path(PARAM$output_dir, PARAM$merged_out_csv)
  readr::write_csv(merged_LMS, merged_out_path)
  
  summary_LMS <- merged_LMS %>%
    dplyr::group_by(SAMPLE, SEGMENT, EFFECT) %>%
    dplyr::summarise(N = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(SAMPLE, SEGMENT, dplyr::desc(N))
  
  summary_out_csv  <- sub("\\.csv$", "_summary.csv", PARAM$merged_out_csv)
  summary_out_path <- file.path(PARAM$output_dir, summary_out_csv)
  readr::write_csv(summary_LMS, summary_out_path)
  
  if (verbose) {
    message("Catálogo mesclado: ", merged_out_path)
    message("Tabela de resumo:   ", summary_out_path)
  }
  
}  # fim Parte 1

##### 4) PARTE 2 — PLOT DO SNPs LANDSCAPE ##################################

if (!RUN_PART2) {
  message("RUN_PART2=FALSE: pulando Parte 2 (plot).")
} else {
  input_path <- file.path(PARAM_PLOT$workdir, PARAM_PLOT$datafr_file)
  if (!file.exists(input_path)) {
    stop("Parte 2: arquivo de entrada não encontrado: ", input_path,
         "\nDica: rode a Parte 1 (RUN_PART1=TRUE) ou ajuste PARAM_PLOT$workdir/datafr_file.")
  }
  
  
  setwd(PARAM_PLOT$workdir)
  if (!dir.exists(PARAM_PLOT$out_dir_fig)) dir.create(PARAM_PLOT$out_dir_fig, recursive = TRUE)
  
  dados_raw <- readr::read_csv(PARAM_PLOT$datafr_file, show_col_types = FALSE)
  
  ugene_colors <- c(
    "Basic/Positively (RHK)"  = "#0000FF",
    "Acidic/Negatively (DE)"  = "#FF0000",
    "Polar/Uncharged (STNQ)"  = "#00FF00",
    "Hydrophobic (AVILM)"     = "#FFFF00",
    "Aromatic (FYW)"          = "#FFA500",
    "Small (CGP)"             = "#808080",
    "Stop"                    = "#000000"
  )
  
  ugene_aminoacid_class <- c(
    K="Basic/Positively (RHK)",
    R="Basic/Positively (RHK)",
    H="Basic/Positively (RHK)",
    D="Acidic/Negatively (DE)",
    E="Acidic/Negatively (DE)",
    F="Aromatic (FYW)",
    W="Aromatic (FYW)",
    Y="Aromatic (FYW)",
    A="Hydrophobic (AVILM)",
    V="Hydrophobic (AVILM)",
    I="Hydrophobic (AVILM)",
    L="Hydrophobic (AVILM)",
    M="Hydrophobic (AVILM)",
    C="Small (CGP)",
    G="Small (CGP)",
    P="Small (CGP)",
    S="Polar/Uncharged (STNQ)",
    T="Polar/Uncharged (STNQ)",
    N="Polar/Uncharged (STNQ)",
    Q="Polar/Uncharged (STNQ)",
    `*`="Stop"
  )
  
  dados <- dados_raw %>%
    mutate(
      Segment        = factor(if ("SEGMENT" %in% names(.)) SEGMENT else CHROM, levels = c("L","M","S")),
      Header         = SAMPLE,
      Amino_position = suppressWarnings(as.numeric(AA_POS)),
      Amino_change   = paste0(AA_REF, ">", AA_ALT),
      Base_position  = POS,
      Base_change    = paste0(REF, ">", ALT),
      Codon_change   = paste0(CODON_REF, ">", CODON_ALT),
      
      Mutation_type = case_when(
        EFFECT == "synonymous" ~ "Syn",
        EFFECT %in% c("nonsynonymous", "stop_gained", "stop_lost") ~ "Nonsyn",
        TRUE ~ "Unk"
      ),
      
      aa_ref = AA_REF,
      aa_alt = AA_ALT,
      AA_class_alt = unname(ugene_aminoacid_class[aa_alt])
    ) %>%
    filter(!is.na(Amino_position))
  
  order_headers_by_prefix <- function(headers, prefixes) {
    prefixes <- prefixes[!is.na(prefixes) & nzchar(prefixes)]
    if (length(prefixes) == 0) return(headers)
    
    prefix_regex <- paste0("^(", paste0(stringr::str_replace_all(prefixes,
                                                                 "([\\^\\$\\.|\\(\\)\\[\\]\\*\\+\\?\\\\])", "\\\\\\1"),
                                        collapse = "|"), ")")
    full_regex   <- paste0(prefix_regex, "\\d+$")
    
    if (!all(grepl(full_regex, headers))) return(headers)
    
    nums <- as.integer(sub(prefix_regex, "", headers))
    headers[order(nums)]
  }
  
  headers <- unique(dados$Header)
  headers_sorted <- order_headers_by_prefix(headers, PARAM_PLOT$header_prefixes)
  dados$Header <- factor(dados$Header, levels = headers_sorted)
  
  dados_plot <- dados %>%
    group_by(Header, Segment, Amino_position, Mutation_type, AA_class_alt) %>%
    summarise(
      Amino_label = paste(unique(gsub(">", ":", Amino_change)), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(
      seg_id   = as.numeric(Segment),
      y_bottom = seg_id - 0.25,
      y_top    = seg_id + ifelse(Mutation_type == "Syn", 0.15, 0.30),
      y_point  = y_top
    ) %>%
    arrange(Header, Segment, Amino_position) %>%
    group_by(Header, Segment) %>%
    mutate(
      label_jitter_step = row_number() %% 2,
      y_label = y_point + ifelse(label_jitter_step == 0, 0.15, 0.28)
    ) %>%
    ungroup()
  
  p <- ggplot(dados_plot, aes(x = Amino_position, y = y_point)) +
    
    geom_segment(aes(xend = Amino_position, y = y_bottom, yend = y_top),
                 color = "grey80", linewidth = 0.4) +
    
    geom_point(
      data = subset(dados_plot, Mutation_type != "Nonsyn"),
      aes(color = Mutation_type),
      size = 2
    ) +
    
    geom_point(
      data = subset(dados_plot, Mutation_type == "Nonsyn"),
      aes(color = Mutation_type, fill = AA_class_alt),
      shape = 21, size = 3.2, stroke = 1
    ) +
    
    geom_text(
      data = subset(dados_plot, Mutation_type == "Nonsyn"),
      aes(y = y_label, label = Amino_label),
      size = 2, vjust = 0, fontface = "plain"
    ) +
    
    scale_y_continuous(
      breaks = sort(unique(dados_plot$seg_id)),
      labels = levels(dados_plot$Segment),
      name   = "Segments"
    ) +
    
    scale_x_continuous(
      name   = "Amino acid position",
      breaks = seq(0, max(dados_plot$Amino_position, na.rm = TRUE), by = 500)
    ) +
    
    scale_color_manual(
      name   = "SNPs landscape of L, M and S segments",
      breaks = c("Syn", "Nonsyn", "Unk"),
      values = c(Syn = "#1f78b4", Nonsyn = "#e31a1c", Unk = "grey50"),
      labels = c(
        Syn    = "Synonymous",
        Nonsyn = "Non-synonymous (incl. stop)",
        Unk    = "Unknown/ambiguous"
      )
    ) +
    
    scale_fill_manual(
      name     = NULL,
      breaks   = names(ugene_colors),
      values   = ugene_colors,
      na.value = "white"
    ) +
    
    guides(
      color = guide_legend(order = 1, title.position = "top"),
      fill  = guide_legend(order = 2, title = NULL)
    ) +
    
    facet_wrap(~ Header, ncol = PARAM_PLOT$ncol_facets) +
    
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.background   = element_rect(fill = "grey90"),
      strip.text         = element_text(face = "bold"),
      
      legend.position    = "top",
      legend.box         = "vertical",
      legend.box.spacing = grid::unit(0, "lines"),
      legend.spacing.y   = grid::unit(0, "lines"),
      legend.margin      = margin(0, 0, 0, 0),
      legend.box.margin  = margin(t = -2, r = 0, b = 0, l = 0)
    )
  
  print(p)
  
  ggsave(file.path(PARAM_PLOT$out_dir_fig, "mutational_landscape_L_M_S.svg"),
         plot = p, width = PARAM_PLOT$width_in, height = PARAM_PLOT$height_in, units = "in")
  ggsave(file.path(PARAM_PLOT$out_dir_fig, "mutational_landscape_L_M_S.pdf"),
         plot = p, width = PARAM_PLOT$width_in, height = PARAM_PLOT$height_in, units = "in")
  ggsave(file.path(PARAM_PLOT$out_dir_fig, "mutational_landscape_L_M_S.png"),
         plot = p, width = PARAM_PLOT$width_in, height = PARAM_PLOT$height_in, units = "in", dpi = PARAM_PLOT$dpi_png)
  
  message("Figura final salva em: ",
          file.path(PARAM_PLOT$out_dir_fig, "mutational_landscape_L_M_S.png"))
  
}  # fim Parte 2

