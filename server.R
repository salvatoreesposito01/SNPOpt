library(shiny)
library(GA)
library(data.table)
library(ggplot2)

### ---------- Utils: detection & readers ----------

detect_file_type <- function(path) {
  first <- tryCatch(readLines(path, n = 50, warn = FALSE), error = function(e) character())
  first <- first[nzchar(first)]
  if (any(startsWith(first, "##fileformat=VCF") | startsWith(first, "#CHROM\t"))) return("vcf")
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("vcf","gz")) return("vcf")
  "hapmap"
}

read_hapmap <- function(path) {
  dt <- tryCatch(
    suppressWarnings(fread(path, header = TRUE, sep = "\t", data.table = FALSE,
                           na.strings = c("NA",".","","N"))),
    error = function(e) NULL
  )
  if (is.null(dt)) {
    dt <- tryCatch(
      suppressWarnings(fread(path, header = TRUE, sep = " ", data.table = FALSE,
                             na.strings = c("NA",".","","N"))),
      error = function(e) NULL
    )
  }
  if (is.null(dt) || ncol(dt) < 12) return(NULL)
  list(
    meta = dt[, 1:11, drop = FALSE],
    geno = as.data.frame(dt[, 12:ncol(dt), drop = FALSE], stringsAsFactors = FALSE)
  )
}

iupac_hetero <- function(a, b) {
  if (a == b) return(a)
  key <- paste(sort(c(a,b)), collapse = "")
  switch(key,
         "AC"="M", "AG"="R", "AT"="W",
         "CG"="S", "CT"="Y", "GT"="K",
         NA_character_)
}

read_vcf_minimal <- function(path, max_records = Inf) {
  con <- if (grepl("\\.gz$", tolower(path))) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)

  header <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "##")) next
    if (startsWith(line, "#CHROM")) { header <- line; break }
  }
  if (!length(header)) return(NULL)

  hdr <- strsplit(sub("^#", "", header), "\t", fixed = TRUE)[[1]]
  if (length(hdr) < 10) return(NULL)
  sample_names <- hdr[10:length(hdr)]

  chrom <- character(); pos <- integer(); refv <- character(); altv <- character()
  rows_chr <- list(); n_kept <- 0L

  while (length(line <- readLines(con, n = 1, warn = FALSE))) {
    if (n_kept >= max_records) break
    if (!nzchar(line) || startsWith(line, "#")) next
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(f) < 10) next
    CHROM <- f[1]; POS <- as.integer(f[2]); REF <- f[4]; ALT <- f[5]
    if (nchar(REF) != 1) next
    if (grepl(",", ALT, fixed = TRUE)) next
    if (nchar(ALT) != 1) next

    fmt <- f[9]; fmt_fields <- strsplit(fmt, ":", fixed = TRUE)[[1]]
    gt_idx <- match("GT", fmt_fields); if (is.na(gt_idx)) next
    alleles <- c(REF, ALT)
    samp_gt <- f[10:length(f)]
    row_vals <- character(length(sample_names))

    for (i in seq_along(samp_gt)) {
      tok <- strsplit(samp_gt[i], ":", fixed = TRUE)[[1]][gt_idx]
      if (!length(tok) || tok %in% c(".", "./.", ".|.")) { row_vals[i] <- NA_character_; next }
      tok <- gsub("\\|", "/", tok)
      sp <- strsplit(tok, "/", fixed = TRUE)[[1]]
      if (length(sp) != 2 || any(sp %in% c(".", ""))) { row_vals[i] <- NA_character_; next }
      a1 <- suppressWarnings(as.integer(sp[1])); a2 <- suppressWarnings(as.integer(sp[2]))
      if (any(is.na(c(a1,a2))) || a1 > 1 || a2 > 1 || a1 < 0 || a2 < 0) { row_vals[i] <- NA_character_; next }
      b1 <- alleles[a1 + 1]; b2 <- alleles[a2 + 1]
      if (b1 %in% c("A","C","G","T") && b2 %in% c("A","C","G","T")) {
        row_vals[i] <- if (b1 == b2) b1 else iupac_hetero(b1, b2)
      } else row_vals[i] <- NA_character_
    }

    rows_chr[[length(rows_chr) + 1L]] <- row_vals
    chrom <- c(chrom, CHROM); pos <- c(pos, POS); refv <- c(refv, REF); altv <- c(altv, ALT)
    n_kept <- n_kept + 1L
  }

  if (!length(rows_chr)) return(NULL)
  m <- do.call(rbind, rows_chr)
  colnames(m) <- sample_names
  rownames(m) <- paste0(chrom, ":", pos)

  meta <- data.frame(chrom = chrom, pos = pos, REF = refv, ALT = altv, stringsAsFactors = FALSE)
  list(meta = meta, geno = as.data.frame(m, stringsAsFactors = FALSE))
}

encode_numeric <- function(m_char) {
  map <- list("A"=0, "C"=0, "G"=2, "T"=2,
              "H"=1, "R"=1, "Y"=1, "S"=1, "W"=1, "K"=1, "M"=1)
  Mnum <- matrix(NA_real_, nrow = nrow(m_char), ncol = ncol(m_char),
                 dimnames = dimnames(m_char))
  for (sym in names(map)) Mnum[m_char == sym] <- map[[sym]]
  Mnum
}

unique_profile_count <- function(subM) {
  if (is.null(dim(subM)) || any(dim(subM) == 0)) return(0L)
  prof <- apply(subM, 2, function(col) paste0(col, collapse = "|"))
  length(unique(prof))
}

unique_profile_fraction <- function(subM) {
  uniq <- unique_profile_count(subM)
  if (is.null(dim(subM)) || any(dim(subM) == 0)) return(0)
  uniq / ncol(subM)
}

simpson_div <- function(subM) {
  if (is.null(dim(subM)) || any(dim(subM) == 0)) return(0)
  prof <- apply(subM, 2, function(col) paste0(col, collapse = "|"))
  N <- length(prof); if (N <= 1) return(0)
  freq <- as.integer(table(prof))
  1 - (sum(freq * (freq - 1)) / (N * (N - 1)))
}

dup_groups_list <- function(subM) {
  if (is.null(dim(subM)) || any(dim(subM) == 0)) return(list())
  prof <- apply(subM, 2, function(col) paste0(col, collapse = "|"))
  tbl <- split(colnames(subM), prof)
  Filter(function(v) length(v) >= 2, tbl)
}

fitness_multi <- function(m_char, Mnum, idx,
                          w_unique = 0.6, w_simpson = 0.25,
                          w_missing = 0.1, w_ld = 0.05) {
  if (length(idx) == 0) return(0)
  sub_chr <- m_char[idx, , drop = FALSE]
  sub_num <- Mnum[idx, , drop = FALSE]
  unique_comp  <- unique_profile_fraction(sub_chr)
  simpson_comp <- simpson_div(sub_chr)
  if (is.null(dim(sub_num))) sub_num <- matrix(sub_num, nrow = 1)
  missing_comp <- mean(colMeans(is.na(sub_num)))
  ld_comp <- 0
  if (nrow(sub_num) >= 2) {
    suppressWarnings({
      C <- cor(t(sub_num), use = "pairwise.complete.obs")
    })
    if (!all(is.na(C))) {
      offdiag <- C[upper.tri(C)]
      ld_comp <- mean(abs(offdiag), na.rm = TRUE)
      if (is.nan(ld_comp)) ld_comp <- 0
    }
  }
  score <- (w_unique*unique_comp) + (w_simpson*simpson_comp) -
           (w_missing*missing_comp) - (w_ld*ld_comp)
  if (!is.finite(score)) score <- 0
  max(min(score, 1), -1)
}

repair_k <- function(chrom, k) {
  chrom <- as.integer(round(chrom))
  ones <- which(chrom == 1L)
  if (length(ones) > k) {
    idx_off <- sample(ones, length(ones) - k)
    chrom[idx_off] <- 0L
  } else if (length(ones) < k) {
    zeros <- which(chrom == 0L)
    if (length(zeros) >= (k - length(ones))) {
      idx_on <- sample(zeros, k - length(ones))
      chrom[idx_on] <- 1L
    }
  }
  chrom
}

# -------- Discrimination curve (greedy) --------
discrimination_curve <- function(m_char, sel_idx) {
  if (length(sel_idx) == 0) return(list(x = numeric(0), y = numeric(0), order = integer(0)))
  remaining <- sel_idx
  used <- integer(0)
  y <- numeric(0)
  for (step in seq_along(sel_idx)) {
    best_gain <- -Inf
    best_snp <- NA_integer_
    base_count <- if (length(used) == 0) 0 else unique_profile_count(m_char[used, , drop = FALSE])
    for (cand in remaining) {
      test_idx <- c(used, cand)
      cnt <- unique_profile_count(m_char[test_idx, , drop = FALSE])
      gain <- cnt - base_count
      if (gain > best_gain) { best_gain <- gain; best_snp <- cand }
    }
    used <- c(used, best_snp)
    remaining <- setdiff(remaining, best_snp)
    y <- c(y, unique_profile_count(m_char[used, , drop = FALSE]))
  }
  list(x = seq_along(y), y = y, order = used)
}

read_fasta_simple <- function(path) {
  lines <- readLines(path)
  hdr_idx <- grep("^>", lines)
  hdr_idx <- c(hdr_idx, length(lines) + 1)
  seqs <- character(length(hdr_idx) - 1)
  names(seqs) <- character(length(hdr_idx) - 1)
  for (i in seq_along(seqs)) {
    nm <- sub("^>", "", lines[hdr_idx[i]])
    nm <- strsplit(nm, "\\s+")[[1]][1]
    names(seqs)[i] <- nm
    start <- hdr_idx[i] + 1
    end <- hdr_idx[i+1] - 1
    seqs[i] <- if (start <= end) paste0(lines[start:end], collapse = "") else ""
  }
  seqs
}

extract_flanks <- function(chr, pos, genome, flankL = 200) {
  if (!(chr %in% names(genome))) return(list(left=NA, ref=NA, right=NA))
  s <- genome[[chr]]
  n <- nchar(s); p <- as.integer(pos)
  if (is.na(p) || p < 1 || p > n) return(list(left=NA, ref=NA, right=NA))
  left_start  <- max(1, p - flankL); left_end <- p - 1
  right_start <- p + 1; right_end <- min(n, p + flankL)
  left_seq  <- if (left_end >= left_start)  substr(s, left_start, left_end) else ""
  ref_base  <- substr(s, p, p)
  right_seq <- if (right_end >= right_start) substr(s, right_start, right_end) else ""
  list(left = left_seq, ref = ref_base, right = right_seq)
}

make_primer3_block <- function(id, flank) {
  template <- paste0(flank$left, flank$ref, flank$right)
  paste0("SEQUENCE_ID=", id, "\n",
         "SEQUENCE_TEMPLATE=", template, "\n",
         "=\n")
}

### --------------- Server ---------------

server <- function(input, output, session) {

  r <- reactiveValues(
    meta = NULL,
    m_char = NULL,
    Mnum = NULL,
    best_idx = NULL,
    best_score = NULL,
    dup_groups = NULL,
    curve_xy = NULL,        # <- punti della curva
    curve_order = NULL,     # <- ordine greedy degli SNP
    flanks_fa = NULL,
    primer3_txt = NULL
  )

  # -------- Load HapMap/VCF --------
  observeEvent(input$hapmap, {
    req(input$hapmap$datapath)
    withProgress(message = "Loading dataset", value = 0, {
      incProgress(0.1, detail = "Detecting format")
      path <- input$hapmap$datapath
      ftype <- detect_file_type(path)

      hm <- NULL
      if (identical(ftype, "vcf")) {
        incProgress(0.2, detail = "Parsing VCF (biallelic SNPs)")
        hm <- read_vcf_minimal(path)
      } else {
        incProgress(0.2, detail = "Parsing HapMap")
        hm <- read_hapmap(path)
      }
      validate(need(!is.null(hm), "Failed to read file. Supported: HapMap (>=12 cols) or VCF (biallelic SNPs with GT)."))

      incProgress(0.2, detail = "Building genotype matrix")
      m <- as.matrix(hm$geno)
      if (is.null(dim(m))) {
        m <- matrix(m, nrow = length(m), ncol = 1)
        colnames(m) <- colnames(hm$geno)
      }
      validate(need(ncol(m) >= 1, "No sample columns found in the dataset."))
      validate(need(nrow(m) >= 1, "No SNP rows found in the dataset."))

      if (!is.null(hm$meta$pos) && !is.null(hm$meta$chrom)) {
        rn <- paste0(hm$meta$chrom, ":", hm$meta$pos)
        if (length(rn) == nrow(m)) rownames(m) <- rn
      }

      incProgress(0.2, detail = "Filtering by heterozygosity")
      het_mat <- m %in% c("H","R","Y","S","W","K","M")
      if (is.null(dim(het_mat))) het_mat <- matrix(het_mat, nrow = nrow(m), ncol = ncol(m))
      het_prop <- rowMeans(het_mat, na.rm = TRUE)
      keep <- het_prop <= input$het_thresh

      if (!any(keep)) {
        showNotification("All SNPs were filtered out by heterozygosity threshold. Increase the threshold.", type = "warning", duration = 6)
        return(NULL)
      }

      m <- m[keep, , drop = FALSE]
      meta <- hm$meta[keep, , drop = FALSE]

      incProgress(0.2, detail = "Encoding numeric")
      r$meta <- meta
      r$m_char <- m
      r$Mnum <- encode_numeric(m)
      r$best_idx <- NULL
      r$best_score <- NULL
      r$dup_groups <- NULL
      r$curve_xy <- NULL
      r$curve_order <- NULL
      incProgress(0.1, detail = "Ready")
    })
  })

  # -------- Discovery (GA) --------
  observeEvent(input$run, {
    req(r$m_char)
    withProgress(message = "Running discovery", value = 0, {
      if (input$useSeed) set.seed(input$seedValue)

      P <- nrow(r$m_char); k <- input$k
      validate(need(P >= k, sprintf("Not enough SNPs after filtering: %d available, k = %d", P, k)))

      incProgress(0.1, detail = "Preparing fitness")
      fitness_wrapper <- function(chrom) {
        chrom <- repair_k(chrom, k)
        idx <- which(as.integer(round(chrom)) == 1L)
        fitness_multi(r$m_char, r$Mnum, idx,
                      w_unique = 0.6, w_simpson = 0.25,
                      w_missing = 0.1, w_ld = 0.05)
      }

      incProgress(0.2, detail = "Evolving population")
      ga_res <- ga(type = "binary",
                   fitness = fitness_wrapper,
                   nBits   = P,
                   popSize = input$popSize,
                   maxiter = input$maxiter,
                   pmutation = 0.2,
                   pcrossover = 0.8,
                   elitism = 2,
                   run = max(30, round(input$maxiter * 0.2)),
                   monitor = FALSE)

      incProgress(0.2, detail = "Extracting best solution")
      sol <- ga_res@solution[1, ]
      sol <- repair_k(sol, k)
      idx <- which(as.integer(round(sol)) == 1L)
      r$best_idx <- idx
      r$best_score <- fitness_wrapper(sol)

      incProgress(0.2, detail = "Duplicates & curve")
      subM <- r$m_char[idx, , drop = FALSE]
      r$dup_groups <- dup_groups_list(subM)

      # Build discrimination curve (greedy)
      dc <- discrimination_curve(r$m_char, idx)
      r$curve_xy <- data.frame(SNPs = dc$x, Unique = dc$y)
      r$curve_order <- dc$order

      incProgress(0.3, detail = "Done")
    })
  })

  # -------- KPIs --------
  output$discoveryStats <- renderUI({
    req(r$m_char)
    if (is.null(r$best_idx)) {
      return(tags$p(class="muted", "No discovery run yet. Upload data and click 'Run discovery'."))
    }
    subM <- r$m_char[r$best_idx, , drop = FALSE]
    tot_snps <- nrow(r$m_char)
    n_sel <- length(r$best_idx)
    n_samples <- ncol(r$m_char)
    uniq_frac <- sprintf("%.3f", unique_profile_fraction(subM))
    simp <- sprintf("%.3f", simpson_div(subM))
    fitness <- sprintf("%.3f", r$best_score)
    n_dup_groups <- length(r$dup_groups)
    tagList(
      span(class = "kpi", strong("SNPs after filter: "), tot_snps),
      span(class = "kpi", strong("Selected SNPs: "), n_sel),
      span(class = "kpi", strong("Samples: "), n_samples),
      br(),
      span(class = "kpi", strong("Unique profile fraction: "), uniq_frac),
      span(class = "kpi", strong("Simpson diversity: "), simp),
      span(class = "kpi", strong("Overall fitness: "), fitness),
      br(),
      span(class = "kpi", strong("Duplicate groups: "), n_dup_groups)
    )
  })

  # -------- Discrimination curve plot --------
  output$discCurve <- renderPlot({
    req(r$curve_xy, r$m_char)
    df <- r$curve_xy
    total_samples <- ncol(r$m_char)
    ggplot(df, aes(x = SNPs, y = Unique)) +
      geom_line() +
      geom_point(size = 2) +
      scale_y_continuous(limits = c(0, total_samples)) +
      labs(x = "Number of SNPs", y = "Unique genotypes",
           title = "SNP Discrimination Curve",
           subtitle = paste0("Unique genotypes: ", tail(df$Unique, 1), " of ", total_samples)) +
      theme_minimal(base_size = 13)
  })

  # -------- Table + Downloads --------
  output$selectedSNPTable <- renderTable({
    req(r$meta, r$best_idx)
    sel <- r$meta[r$best_idx, , drop = FALSE]
    cols <- intersect(c("chrom","pos","rs#","alleles","REF","ALT"), colnames(sel))
    if (length(cols) == 0) sel else sel[, cols, drop = FALSE]
  }, rownames = TRUE)

  output$dlSelectedMeta <- downloadHandler(
    filename = function() sprintf("selected_snps_meta_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      req(r$meta, r$best_idx)
      sel <- r$meta[r$best_idx, , drop = FALSE]
      write.table(sel, file, sep = "\t", quote = FALSE, col.names = NA)
    }
  )

  output$dlSelectedMatrix <- downloadHandler(
    filename = function() sprintf("selected_genotype_matrix_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      req(r$m_char, r$best_idx)
      subM <- r$m_char[r$best_idx, , drop = FALSE]
      write.table(subM, file, sep = "\t", quote = FALSE, col.names = NA)
    }
  )

  output$dlDiscoverySummary <- downloadHandler(
    filename = function() sprintf("discovery_summary_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      req(r$m_char, r$best_idx, r$best_score)
      subM <- r$m_char[r$best_idx, , drop = FALSE]
      lines <- c(
        sprintf("SNPs after filter: %d", nrow(r$m_char)),
        sprintf("Selected SNPs: %d", length(r$best_idx)),
        sprintf("Samples: %d", ncol(r$m_char)),
        sprintf("Unique profile fraction: %.4f", unique_profile_fraction(subM)),
        sprintf("Simpson diversity: %.4f", simpson_div(subM)),
        sprintf("Overall fitness: %.4f", r$best_score),
        sprintf("Duplicate groups: %d", length(r$dup_groups)),
        "",
        "Selected SNP IDs (greedy order for curve):",
        paste(rownames(r$m_char)[r$curve_order], collapse = ", ")
      )
      writeLines(lines, file)
    }
  )

  output$dlPlot <- downloadHandler(
    filename = function() sprintf("SNP_discrimination_curve_%s.png", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      req(r$curve_xy, r$m_char)
      df <- r$curve_xy
      total_samples <- ncol(r$m_char)
      p <- ggplot(df, aes(x = SNPs, y = Unique)) +
        geom_line() + geom_point(size = 2) +
        scale_y_continuous(limits = c(0, total_samples)) +
        labs(x = "Number of SNPs", y = "Unique genotypes",
             title = "SNP Discrimination Curve",
             subtitle = paste0("Unique genotypes: ", tail(df$Unique, 1), " of ", total_samples)) +
        theme_minimal(base_size = 13)
      ggsave(file, p, width = 7, height = 5, dpi = 150)
    }
  )

  output$dlDupSamples <- downloadHandler(
    filename = function() sprintf("duplicated_samples_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      if (is.null(r$dup_groups) || !length(r$dup_groups)) {
        writeLines("No duplicate profiles remain.", file); return()
      }
      lines <- unlist(lapply(seq_along(r$dup_groups), function(i) {
        paste0("[Group ", i, "] ", paste(r$dup_groups[[i]], collapse = ", "))
      }))
      writeLines(lines, file)
    }
  )

  output$dupGroupsTxt <- renderText({
    if (is.null(r$dup_groups) || !length(r$dup_groups)) return("No duplicate profiles remain.")
    paste(vapply(seq_along(r$dup_groups),
                 function(i) paste0("[Group ", i, "] ", paste(r$dup_groups[[i]], collapse = ", ")),
                 FUN.VALUE = ""), collapse = "\n")
  })

  # -------- Flanking/Primer3 --------
  observeEvent(input$makeFlanks, {
    req(r$best_idx, input$ref_fasta)
    withProgress(message = "Generating flanking sequences", value = 0, {
      incProgress(0.2, detail = "Reading FASTA")
      genome <- read_fasta_simple(input$ref_fasta$datapath)

      sel_meta <- r$meta[r$best_idx, , drop = FALSE]
      has_chr <- "chrom" %in% colnames(sel_meta)
      has_pos <- "pos" %in% colnames(sel_meta)

      ids <- chr <- character(length(r$best_idx))
      pos <- integer(length(r$best_idx))

      incProgress(0.2, detail = "Collecting coordinates")
      for (i in seq_along(r$best_idx)) {
        ridx <- r$best_idx[i]
        if (has_chr && has_pos) {
          chr[i] <- as.character(sel_meta$chrom[i])
          pos[i] <- as.integer(sel_meta$pos[i])
          ids[i] <- if ("rs#" %in% colnames(sel_meta)) as.character(sel_meta$`rs#`[i]) else paste0(chr[i], ":", pos[i])
        } else {
          nm <- rownames(r$m_char)[ridx]
          ids[i] <- nm
          sp <- strsplit(nm, ":", fixed = TRUE)[[1]]
          chr[i] <- sp[1]
          pos[i] <- as.integer(sp[2])
        }
      }

      incProgress(0.3, detail = "Extracting flanks")
      flankL <- input$flankL
      fa_lines <- c(); p3_blocks <- c()
      for (i in seq_along(ids)) {
        flk <- extract_flanks(chr[i], pos[i], genome, flankL = flankL)
        if (is.na(flk$ref)) next
        fa_lines <- c(fa_lines, paste0(">", ids[i]), paste0(flk$left, flk$ref, flk$right))
        if (isTRUE(input$primer3_export)) p3_blocks <- c(p3_blocks, make_primer3_block(ids[i], flk))
      }

      incProgress(0.3, detail = "Finalizing")
      r$flanks_fa <- paste0(fa_lines, collapse = "\n")
      r$primer3_txt <- if (isTRUE(input$primer3_export)) paste0(p3_blocks, collapse = "") else NULL
    })
  })

  output$dlFlanksFA <- downloadHandler(
    filename = function() sprintf("flanking_%s.fa", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      validate(need(!is.null(r$flanks_fa), "Generate flanking sequences first."))
      writeLines(r$flanks_fa, file)
    }
  )

  output$dlPrimer3 <- downloadHandler(
    filename = function() sprintf("primer3_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      validate(need(!is.null(r$primer3_txt), "Primer3 export not generated."))
      writeLines(r$primer3_txt, file)
    }
  )

  # -------- Evaluation (unchanged logic) --------
  parse_snp_list <- function(text) {
    if (is.null(text) || !nzchar(text)) return(NULL)
    lines <- strsplit(text, "\n")[[1]]
    lines <- trimws(lines); lines <- lines[nzchar(lines)]
    chr <- c(); pos <- c()
    for (ln in lines) {
      if (grepl(":", ln)) {
        sp <- strsplit(ln, ":", fixed = TRUE)[[1]]
        chr <- c(chr, sp[1]); pos <- c(pos, as.integer(sp[2]))
      } else {
        sp <- strsplit(ln, "\\s+")[[1]]
        chr <- c(chr, sp[1]); pos <- c(pos, as.integer(sp[2]))
      }
    }
    data.frame(CHR = chr, POS = pos, stringsAsFactors = FALSE)
  }

  eval_data <- reactive({
    if (!is.null(input$eval_dataset)) {
      path <- input$eval_dataset$datapath
      ftype <- detect_file_type(path)
      hm <- if (identical(ftype, "vcf")) read_vcf_minimal(path) else read_hapmap(path)
      validate(need(!is.null(hm), "Failed to read evaluation dataset (HapMap/VCF)."))
      geno <- as.matrix(hm$geno); if (is.null(dim(geno))) geno <- matrix(geno, nrow = length(geno), ncol = 1)
      list(meta = hm$meta, geno = geno)
    } else if (!is.null(r$m_char)) {
      list(meta = r$meta, geno = r$m_char)
    } else NULL
  })

  observeEvent(input$evalBtn, {
    req(eval_data())
    showNotification("Evaluation started", type = "message", duration = 2)
  })

  output$evalPanel <- renderUI({
    ed <- eval_data()
    if (is.null(ed)) return(tagList(p("No dataset available yet. Upload one or run Discovery first.")))

    df_list <- NULL
    if (!is.null(input$snp_list_file)) {
      df_list <- tryCatch(
        suppressWarnings(fread(input$snp_list_file$datapath, header = TRUE, data.table = FALSE)),
        error = function(e) NULL
      )
      if (!is.null(df_list) && all(c("CHR","POS") %in% colnames(df_list))) {
        df_list <- df_list[, c("CHR","POS")]
      } else df_list <- NULL
    }
    if (is.null(df_list) && nzchar(input$snp_list_text)) df_list <- parse_snp_list(input$snp_list_text)

    if (is.null(df_list) || nrow(df_list) == 0) {
      return(tagList(h3("Evaluation"), p("Provide a SNP list (CHR, POS) to evaluate.")))
    }

    meta <- ed$meta; geno <- ed$geno
    has_chr <- "chrom" %in% colnames(meta)
    has_pos <- "pos" %in% colnames(meta)
    if (!(has_chr && has_pos)) {
      return(tagList(h3("Evaluation"),
                     p("CHR/POS not available in metadata. Ensure your dataset includes 'chrom' and 'pos' columns.")))
    }

    withProgress(message = "Subsetting matrix", value = 0, {
      incProgress(0.4)
      key_meta <- paste0(meta$chrom, ":", meta$pos)
      key_req  <- paste0(df_list$CHR, ":", df_list$POS)
      idx <- match(key_req, key_meta)
      found <- !is.na(idx)
      incProgress(0.4)
      subM <- geno[idx[found], , drop = FALSE]
      rn <- paste0(df_list$CHR[found], ":", df_list$POS[found])
      rownames(subM) <- rn
      incProgress(0.2)

      tagList(
        h3("Evaluation"),
        p(strong("Requested SNPs:"), nrow(df_list),
          " — ", strong("Found in dataset:"), sum(found)),
        if (any(!found)) {
          tagList(p("Not found:"), verbatimTextOutput("notFoundList"))
        },
        tags$hr(),
        h4("Subset matrix (first 50 rows)"),
        tableOutput("evalSubsetTable"),
        tags$hr(),
        downloadButton("dlEvalMatrix", "Download subset matrix")
      )
    })
  })

  output$evalSubsetTable <- renderTable({
    ed <- eval_data(); req(ed)
    meta <- ed$meta; geno <- ed$geno
    has_chr <- "chrom" %in% colnames(meta); has_pos <- "pos" %in% colnames(meta)
    if (!(has_chr && has_pos)) return(NULL)

    df_list <- NULL
    if (!is.null(input$snp_list_file)) {
      df_list <- tryCatch(
        suppressWarnings(fread(input$snp_list_file$datapath, header = TRUE, data.table = FALSE)),
        error = function(e) NULL
      )
      if (!is.null(df_list) && all(c("CHR","POS") %in% colnames(df_list))) df_list <- df_list[, c("CHR","POS")] else df_list <- NULL
    }
    if (is.null(df_list) && nzchar(input$snp_list_text)) {
      df_list <- parse_snp_list(input$snp_list_text)
    }
    if (is.null(df_list)) return(NULL)

    key_meta <- paste0(meta$chrom, ":", meta$pos)
    key_req  <- paste0(df_list$CHR, ":", df_list$POS)
    idx <- match(key_req, key_meta)
    found <- !is.na(idx)
    subM <- geno[idx[found], , drop = FALSE]
    rn <- paste0(df_list$CHR[found], ":", df_list$POS[found])
    rownames(subM) <- rn

    head(as.data.frame(subM), 50)
  }, rownames = TRUE)

  output$notFoundList <- renderText({
    ed <- eval_data(); req(ed)
    meta <- ed$meta
    has_chr <- "chrom" %in% colnames(meta); has_pos <- "pos" %in% colnames(meta)
    if (!(has_chr && has_pos)) return("")
    df_list <- NULL
    if (!is.null(input$snp_list_file)) {
      df_list <- tryCatch(
        suppressWarnings(fread(input$snp_list_file$datapath, header = TRUE, data.table = FALSE)),
        error = function(e) NULL
      )
      if (!is.null(df_list) && all(c("CHR","POS") %in% colnames(df_list))) df_list <- df_list[, c("CHR","POS")] else df_list <- NULL
    }
    if (is.null(df_list) && nzchar(input$snp_list_text)) df_list <- parse_snp_list(input$snp_list_text)
    if (is.null(df_list)) return("")
    key_meta <- paste0(meta$chrom, ":", meta$pos)
    key_req  <- paste0(df_list$CHR, ":", df_list$POS)
    paste(key_req[is.na(match(key_req, key_meta))], collapse = "\n")
  })

  output$dlEvalMatrix <- downloadHandler(
    filename = function() sprintf("subset_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      ed <- eval_data(); req(ed)
      meta <- ed$meta; geno <- ed$geno
      has_chr <- "chrom" %in% colnames(meta); has_pos <- "pos" %in% colnames(meta)
      validate(need(has_chr && has_pos, "CHR/POS not available in metadata."))

      df_list <- NULL
      if (!is.null(input$snp_list_file)) {
        df_list <- tryCatch(
          suppressWarnings(fread(input$snp_list_file$datapath, header = TRUE, data.table = FALSE)),
          error = function(e) NULL
        )
        if (!is.null(df_list) && all(c("CHR","POS") %in% colnames(df_list))) df_list <- df_list[, c("CHR","POS")] else df_list <- NULL
      }
      if (is.null(df_list) && nzchar(input$snp_list_text)) df_list <- parse_snp_list(input$snp_list_text)
      validate(need(!is.null(df_list) && nrow(df_list) > 0, "Provide a SNP list first."))

      key_meta <- paste0(meta$chrom, ":", meta$pos)
      key_req  <- paste0(df_list$CHR, ":", df_list$POS)
      idx <- match(key_req, key_meta)
      found <- !is.na(idx)
      subM <- geno[idx[found], , drop = FALSE]
      rn <- paste0(df_list$CHR[found], ":", df_list$POS[found])
      rownames(subM) <- rn

      write.table(subM, file, sep = "\t", quote = FALSE, col.names = NA)
    }
  )
}
