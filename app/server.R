library(shiny)
options(shiny.maxRequestSize = 1024 * 1024^2) # 1 GB
library(GA)
library(data.table)
library(ggplot2)

### ---------- Utils: detection & readers ----------

detect_file_type <- function(path) {
  first <- tryCatch({
    con <- if (grepl("\\.gz$", tolower(path))) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con), add = TRUE)
    readLines(con, n = 50, warn = FALSE)
  }, error = function(e) character())
  first <- first[nzchar(first)]
  if (any(startsWith(first, "##fileformat=VCF") | startsWith(first, "#CHROM\t"))) return("vcf")
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("vcf")) return("vcf")
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

score_subset <- function(m_char, Mnum, idx, sample_idx = NULL,
                         w_unique = 0.6, w_simpson = 0.25,
                         w_missing = 0.1, w_ld = 0.05) {
  if (is.null(sample_idx)) {
    return(fitness_multi(m_char, Mnum, idx, w_unique, w_simpson, w_missing, w_ld))
  }
  fitness_multi(m_char[, sample_idx, drop = FALSE],
                Mnum[, sample_idx, drop = FALSE],
                idx, w_unique, w_simpson, w_missing, w_ld)
}

greedy_select_snps <- function(m_char, Mnum, candidates, k,
                               base_idx = integer(0), sample_idx = NULL) {
  if (k <= 0 || length(candidates) == 0) return(sort(unique(base_idx)))
  selected <- sort(unique(base_idx))
  remaining <- setdiff(candidates, selected)
  target_k <- min(k, length(remaining))
  if (target_k <= 0) return(selected)

  for (step in seq_len(target_k)) {
    max_unique <- if (is.null(sample_idx)) ncol(m_char) else length(sample_idx)
    uniq_vec <- numeric(length(remaining))
    for (j in seq_along(remaining)) {
      cand <- remaining[j]
      idx <- c(selected, cand)
      chr_m <- if (is.null(sample_idx)) m_char[idx, , drop = FALSE] else m_char[idx, sample_idx, drop = FALSE]
      uniq_vec[j] <- unique_profile_count(chr_m)
    }
    best_uniq <- max(uniq_vec, na.rm = TRUE)
    tie_idx <- which(uniq_vec == best_uniq)
    tie_cands <- remaining[tie_idx]
    if (length(tie_cands) == 1) {
      best_cand <- tie_cands[1]
      selected <- c(selected, best_cand)
      remaining <- setdiff(remaining, best_cand)
      next
    }

    best <- list(idx = NA_integer_, uniq = best_uniq, simp = -Inf, miss = Inf, ld = Inf)
    for (cand in tie_cands) {
      idx <- c(selected, cand)
      chr_m <- if (is.null(sample_idx)) m_char[idx, , drop = FALSE] else m_char[idx, sample_idx, drop = FALSE]
      num_m <- if (is.null(sample_idx)) Mnum[idx, , drop = FALSE] else Mnum[idx, sample_idx, drop = FALSE]
      simp <- simpson_div(chr_m)
      miss <- mean(colMeans(is.na(num_m)))
      ld <- 0
      # LD is expensive; only evaluate it when all other criteria tie.
      if (nrow(num_m) >= 2 && simp == best$simp && miss == best$miss) {
        suppressWarnings({
          C <- cor(t(num_m), use = "pairwise.complete.obs")
        })
        if (!all(is.na(C))) {
          off <- C[upper.tri(C)]
          ld <- mean(abs(off), na.rm = TRUE)
          if (is.nan(ld)) ld <- 0
        }
      }
      better <- (simp > best$simp) ||
        (simp == best$simp && miss < best$miss) ||
        (simp == best$simp && miss == best$miss && ld < best$ld)
      if (better) best <- list(idx = cand, uniq = best_uniq, simp = simp, miss = miss, ld = ld)
    }
    if (is.na(best$idx)) best$idx <- tie_cands[1]
    selected <- c(selected, best$idx)
    remaining <- setdiff(remaining, best$idx)
    if (best_uniq >= max_unique && length(selected) >= length(base_idx) + 1 && length(selected) >= target_k) break
  }

  sort(unique(selected))
}

local_swap_optimize <- function(m_char, Mnum, selected, candidate_pool,
                                sample_idx = NULL, max_pass = 1, max_swap_candidates = 120) {
  selected <- sort(unique(selected))
  if (length(selected) == 0) return(selected)
  if (is.null(sample_idx)) {
    sub_chr <- m_char[selected, , drop = FALSE]
    if (unique_profile_count(sub_chr) >= ncol(m_char)) return(selected)
  }
  pool <- sort(unique(candidate_pool))
  others <- setdiff(pool, selected)
  if (length(others) == 0) return(selected)
  if (length(others) > max_swap_candidates) others <- head(others, max_swap_candidates)

  cur_score <- score_subset(m_char, Mnum, selected, sample_idx = sample_idx)
  eps <- 1e-12
  for (pass in seq_len(max_pass)) {
    improved <- FALSE
    for (i in seq_along(selected)) {
      best_score <- cur_score
      best_cand <- NA_integer_
      for (cand in others) {
        trial <- selected
        trial[i] <- cand
        trial <- sort(unique(trial))
        if (length(trial) != length(selected)) next
        sc <- score_subset(m_char, Mnum, trial, sample_idx = sample_idx)
        if (sc > best_score + eps) {
          best_score <- sc
          best_cand <- cand
        }
      }
              if (!is.na(best_cand)) {
                  selected[i] <- best_cand
                  selected <- sort(selected)
        others <- setdiff(pool, selected)
        if (length(others) > max_swap_candidates) others <- head(others, max_swap_candidates)
        cur_score <- best_score
        improved <- TRUE
      }
    }
    if (!improved) break
  }
  selected
}

repair_k <- function(chrom, k) {
  chrom <- as.integer(round(chrom))
  # Keep repair reproducible per chromosome to avoid noisy GA fitness.
  seed <- 17
  for (i in seq_along(chrom)) {
    seed <- (seed * 1103515245 + chrom[i] + i * 12345) %% 2147483647
  }
  seed <- as.integer(max(1L, seed + as.integer(k)))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  set.seed(seed)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

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
  con <- if (grepl("\\.gz$", tolower(path))) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  hdr_idx <- grep("^>", lines)
  if (!length(hdr_idx)) return(character())
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

normalize_chr_id <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- sub("^chromosome", "", x)
  x <- sub("^chr", "", x)
  x <- sub("^ch", "", x)
  x <- gsub("[^a-z0-9]", "", x)
  x <- sub("^0+([0-9]+)$", "\\1", x)
  x
}

header_chr_id <- function(hdr_token) {
  tok <- strsplit(as.character(hdr_token), "[|:]", perl = TRUE)[[1]][1]
  normalize_chr_id(tok)
}

is_chr_match <- function(hdr_norm, target_norm) {
  chr_aliases <- function(x) {
    if (!nzchar(x)) return(character(0))
    out <- c(x)
    # Handle headers like sl40ch01 -> 1
    m <- sub(".*ch0*([0-9]+)$", "\\1", x, perl = TRUE)
    if (!identical(m, x) && nzchar(m)) out <- c(out, m)
    # Generic trailing numeric suffix alias (e.g. contig001 -> 1)
    n <- sub(".*?([0-9]+)$", "\\1", x, perl = TRUE)
    if (!identical(n, x) && nzchar(n)) out <- c(out, n)
    out <- unique(out)
    out
  }

  if (!nzchar(hdr_norm) || !nzchar(target_norm)) return(FALSE)
  h <- chr_aliases(hdr_norm)
  t <- chr_aliases(target_norm)
  if (length(intersect(h, t)) > 0) return(TRUE)
  # Avoid numeric prefix ambiguity (e.g. 1 vs 11): allow prefix fallback only on non-numeric aliases.
  h_txt <- h[!grepl("^[0-9]+$", h)]
  t_txt <- t[!grepl("^[0-9]+$", t)]
  if (!length(h_txt) || !length(t_txt)) return(FALSE)
  any(vapply(h_txt, function(hx) any(startsWith(t_txt, hx) | startsWith(hx, t_txt)), logical(1)))
}

read_fasta_header_ids <- function(path, max_n = 2000) {
  con <- if (grepl("\\.gz$", tolower(path))) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  ids <- character(0)
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) break
    if (!startsWith(line, ">")) next
    nm <- sub("^>", "", line)
    nm <- strsplit(nm, "\\s+")[[1]][1]
    ids <- c(ids, nm)
    if (length(ids) >= max_n) break
  }
  ids
}

resolve_chr_header <- function(query_chr, fasta_headers) {
  qn <- normalize_chr_id(query_chr)
  hits <- fasta_headers[vapply(fasta_headers, function(h) {
    is_chr_match(header_chr_id(h), qn)
  }, logical(1))]
  hits <- unique(hits)
  if (length(hits) == 1) return(hits[1])
  NA_character_
}

ensure_fai_index <- function(path) {
  sam <- Sys.which("samtools")
  if (!nzchar(sam)) return(FALSE)
  fai <- paste0(path, ".fai")
  if (file.exists(fai)) return(TRUE)
  status <- suppressWarnings(system2(sam, c("faidx", path), stdout = TRUE, stderr = TRUE))
  file.exists(fai)
}

read_fai_index <- function(path) {
  fai <- paste0(path, ".fai")
  if (!file.exists(fai)) return(NULL)
  dt <- tryCatch(
    read.table(fai, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""),
    error = function(e) NULL
  )
  if (is.null(dt) || ncol(dt) < 2) return(NULL)
  data.frame(name = dt[[1]], len = as.integer(dt[[2]]), stringsAsFactors = FALSE)
}

extract_flanks_faidx <- function(path, chr_header, pos, chr_len, flankL = 200) {
  p <- as.integer(pos); n <- as.integer(chr_len)
  if (is.na(p) || is.na(n) || p < 1 || p > n) return(list(left=NA, ref=NA, right=NA))
  start <- max(1, p - flankL)
  end <- min(n, p + flankL)
  region <- paste0(chr_header, ":", start, "-", end)
  sam <- Sys.which("samtools")
  if (!nzchar(sam)) return(list(left=NA, ref=NA, right=NA))

  out <- suppressWarnings(system2(sam, c("faidx", path, region), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) return(list(left=NA, ref=NA, right=NA))
  seq_lines <- out[!startsWith(out, ">")]
  if (!length(seq_lines)) return(list(left=NA, ref=NA, right=NA))
  seq <- toupper(paste0(seq_lines, collapse = ""))
  ref_i <- p - start + 1L
  if (ref_i < 1 || ref_i > nchar(seq)) return(list(left=NA, ref=NA, right=NA))
  list(
    left = if (ref_i > 1) substr(seq, 1, ref_i - 1) else "",
    ref = substr(seq, ref_i, ref_i),
    right = if (ref_i < nchar(seq)) substr(seq, ref_i + 1, nchar(seq)) else ""
  )
}

is_chr_in_headers <- function(query_chr, headers) {
  qn <- normalize_chr_id(query_chr)
  any(vapply(headers, function(h) is_chr_match(header_chr_id(h), qn), logical(1)))
}

read_fasta_targets <- function(path, targets = NULL) {
  con <- if (grepl("\\.gz$", tolower(path))) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)

  keep_all <- is.null(targets) || length(targets) == 0
  target_norm <- if (keep_all) character(0) else unique(normalize_chr_id(targets))

  seq_buf <- list()
  cur_name <- NULL
  cur_keep <- FALSE

  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) break
    if (!nzchar(line)) next

    if (startsWith(line, ">")) {
      nm <- sub("^>", "", line)
      nm <- strsplit(nm, "\\s+")[[1]][1]
      nm_norm <- header_chr_id(nm)
      cur_name <- nm
      cur_keep <- keep_all || any(vapply(target_norm, function(tn) is_chr_match(nm_norm, tn), logical(1)))
      if (cur_keep && is.null(seq_buf[[cur_name]])) seq_buf[[cur_name]] <- character(0)
      next
    }

    if (cur_keep && !is.null(cur_name)) {
      seq_buf[[cur_name]] <- c(seq_buf[[cur_name]], toupper(trimws(line)))
    }
  }

  if (!length(seq_buf)) return(character(0))
  out <- vapply(seq_buf, paste0, collapse = "", FUN.VALUE = character(1), USE.NAMES = TRUE)
  out
}

get_chr_sequence <- function(genome, chr) {
  if (!length(genome)) return(NULL)
  chr <- as.character(chr)
  if (chr %in% names(genome)) return(genome[[chr]])
  chr_norm <- normalize_chr_id(chr)
  nm_norm <- vapply(names(genome), header_chr_id, FUN.VALUE = character(1))
  hits <- which(vapply(nm_norm, function(nm) is_chr_match(nm, chr_norm), logical(1)))
  if (length(hits) == 1) return(genome[[hits]])
  NULL
}

extract_flanks <- function(chr, pos, genome, flankL = 200) {
  s <- get_chr_sequence(genome, chr)
  if (is.null(s)) return(list(left=NA, ref=NA, right=NA))
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
    method_used = NULL,
    dup_groups = NULL,
    dup_groups_before_additional = NULL,
    additional_added_idx = integer(0),
    curve_xy = NULL,        # <- punti della curva
    curve_order = NULL,     # <- ordine greedy degli SNP
    fasta_cache_path = NULL,
    fasta_cache = NULL,
    flanks_fa = NULL,
    primer3_txt = NULL,
    flank_skipped = NULL
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
      r$method_used <- NULL
      r$dup_groups <- NULL
      r$dup_groups_before_additional <- NULL
      r$additional_added_idx <- integer(0)
      r$curve_xy <- NULL
      r$curve_order <- NULL
      r$fasta_cache_path <- NULL
      r$fasta_cache <- NULL
      r$flanks_fa <- NULL
      r$primer3_txt <- NULL
      r$flank_skipped <- NULL
      incProgress(0.1, detail = "Ready")
    })
  })

  # -------- Discovery (GA) --------
  observeEvent(input$run, {
    req(r$m_char)
    withProgress(message = "Running discovery", value = 0, {
      method <- if (is.null(input$search_method)) "Hybrid (greedy + local search)" else input$search_method
      r$method_used <- method
      if (isTRUE(input$useSeed)) set.seed(input$seedValue)

      P <- nrow(r$m_char); k <- input$k
      validate(need(P >= k, sprintf("Not enough SNPs after filtering: %d available, k = %d", P, k)))
      candidates_all <- seq_len(P)
      final_idx <- integer(0)
      r$additional_added_idx <- integer(0)

      if (identical(method, "GA")) {
        incProgress(0.1, detail = "Preparing GA fitness")
        fitness_wrapper <- function(chrom) {
          chrom <- repair_k(chrom, k)
          idx <- which(as.integer(round(chrom)) == 1L)
          fitness_multi(r$m_char, r$Mnum, idx,
                        w_unique = 0.6, w_simpson = 0.25,
                        w_missing = 0.1, w_ld = 0.05)
        }
        incProgress(0.25, detail = "Evolving GA population")
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

        sol <- ga_res@solution[1, ]
        sol <- repair_k(sol, k)
        final_idx <- which(as.integer(round(sol)) == 1L)
      } else {
        incProgress(0.25, detail = "Greedy selection")
        final_idx <- greedy_select_snps(r$m_char, r$Mnum, candidates_all, k)
      }

      incProgress(0.2, detail = "Checking duplicates")
      subM <- r$m_char[final_idx, , drop = FALSE]
      dup_groups <- dup_groups_list(subM)
      if (!identical(method, "GA") && length(dup_groups) > 0) {
        incProgress(0.15, detail = "Local refinement")
        final_idx <- local_swap_optimize(r$m_char, r$Mnum, final_idx, candidates_all,
                                         sample_idx = NULL, max_pass = 1, max_swap_candidates = 120)
        subM <- r$m_char[final_idx, , drop = FALSE]
        dup_groups <- dup_groups_list(subM)
      }
      r$dup_groups_before_additional <- dup_groups

      if (isTRUE(input$doAdditional) && length(dup_groups) > 0) {
        candidates <- setdiff(seq_len(P), final_idx)
        if (length(candidates) > 0) {
          k2 <- min(as.integer(input$n_additional), length(candidates))
          if (k2 > 0) {
            dup_samples <- unique(unlist(dup_groups, use.names = FALSE))
            incProgress(0.2, detail = "Second run on duplicate samples")
            if (identical(method, "GA")) {
              fitness_wrapper2 <- function(chrom) {
                chrom <- repair_k(chrom, k2)
                add_local <- which(as.integer(round(chrom)) == 1L)
                add_idx <- candidates[add_local]
                idx2 <- c(final_idx, add_idx)
                fitness_multi(r$m_char[, dup_samples, drop = FALSE],
                              r$Mnum[, dup_samples, drop = FALSE],
                              idx2,
                              w_unique = 0.6, w_simpson = 0.25,
                              w_missing = 0.1, w_ld = 0.05)
              }
              ga_res2 <- ga(type = "binary",
                            fitness = fitness_wrapper2,
                            nBits = length(candidates),
                            popSize = max(30, round(input$popSize * 0.7)),
                            maxiter = max(40, round(input$maxiter * 0.5)),
                            pmutation = 0.2,
                            pcrossover = 0.8,
                            elitism = 2,
                            run = max(20, round(input$maxiter * 0.1)),
                            monitor = FALSE)
              sol2 <- ga_res2@solution[1, ]
              sol2 <- repair_k(sol2, k2)
              add_local <- which(as.integer(round(sol2)) == 1L)
              added_idx <- if (length(add_local) > 0) candidates[add_local] else integer(0)
              if (length(added_idx) > 0) {
                final_idx <- sort(unique(c(final_idx, added_idx)))
                r$additional_added_idx <- sort(unique(added_idx))
              }
            } else {
              expanded <- greedy_select_snps(r$m_char, r$Mnum, candidates, k2,
                                             base_idx = final_idx, sample_idx = dup_samples)
              added_idx <- setdiff(expanded, final_idx)
              if (length(added_idx) > 0) {
                final_idx <- sort(unique(expanded))
                r$additional_added_idx <- sort(unique(added_idx))
              }
            }
          }
        }
      }

      final_idx <- sort(unique(final_idx))
      r$best_idx <- final_idx
      r$best_score <- fitness_multi(r$m_char, r$Mnum, final_idx,
                                    w_unique = 0.6, w_simpson = 0.25,
                                    w_missing = 0.1, w_ld = 0.05)

      incProgress(0.2, detail = "Final duplicates & curve")
      subM <- r$m_char[final_idx, , drop = FALSE]
      r$dup_groups <- dup_groups_list(subM)

      # Build discrimination curve (greedy)
      dc <- discrimination_curve(r$m_char, final_idx)
      r$curve_xy <- data.frame(SNPs = dc$x, Unique = dc$y)
      r$curve_order <- dc$order

      incProgress(0.1, detail = "Done")
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
    method_used <- if (is.null(r$method_used)) "NA" else r$method_used
    uniq_frac <- sprintf("%.3f", unique_profile_fraction(subM))
    simp <- sprintf("%.3f", simpson_div(subM))
    fitness <- sprintf("%.3f", r$best_score)
    n_dup_groups <- length(r$dup_groups)
    n_dup_before <- if (is.null(r$dup_groups_before_additional)) NA_integer_ else length(r$dup_groups_before_additional)
    n_added <- length(r$additional_added_idx)
    tagList(
      span(class = "kpi", strong("SNPs after filter: "), tot_snps),
      span(class = "kpi", strong("Selected SNPs: "), n_sel),
      span(class = "kpi", strong("Samples: "), n_samples),
      span(class = "kpi", strong("Method: "), method_used),
      br(),
      span(class = "kpi", strong("Unique profile fraction: "), uniq_frac),
      span(class = "kpi", strong("Simpson diversity: "), simp),
      span(class = "kpi", strong("Overall fitness: "), fitness),
      br(),
      span(class = "kpi", strong("Duplicate groups: "), n_dup_groups),
      if (!is.na(n_dup_before)) span(class = "kpi", strong("Duplicate groups before additional run: "), n_dup_before),
      if (n_added > 0) span(class = "kpi", strong("Additional SNPs added: "), n_added)
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
      n_dup_before <- if (is.null(r$dup_groups_before_additional)) NA_integer_ else length(r$dup_groups_before_additional)
      n_added <- length(r$additional_added_idx)
      add_ids <- if (n_added > 0) paste(rownames(r$m_char)[r$additional_added_idx], collapse = ", ") else "None"
      lines <- c(
        sprintf("Method: %s", if (is.null(r$method_used)) "NA" else r$method_used),
        sprintf("SNPs after filter: %d", nrow(r$m_char)),
        sprintf("Selected SNPs: %d", length(r$best_idx)),
        sprintf("Samples: %d", ncol(r$m_char)),
        sprintf("Unique profile fraction: %.4f", unique_profile_fraction(subM)),
        sprintf("Simpson diversity: %.4f", simpson_div(subM)),
        sprintf("Overall fitness: %.4f", r$best_score),
        sprintf("Additional run enabled: %s", if (isTRUE(input$doAdditional)) "yes" else "no"),
        if (!is.na(n_dup_before)) sprintf("Duplicate groups before additional run: %d", n_dup_before) else "Duplicate groups before additional run: NA",
        sprintf("Duplicate groups: %d", length(r$dup_groups)),
        sprintf("Additional SNPs added: %d", n_added),
        sprintf("Additional SNP IDs: %s", add_ids),
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
      r$flank_skipped <- data.frame(
        snp_id = character(0),
        chr_requested = character(0),
        pos = integer(0),
        reason = character(0),
        used_faidx = logical(0),
        stringsAsFactors = FALSE
      )
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

      req_chr <- unique(chr[nzchar(chr)])
      req_chr <- req_chr[!is.na(req_chr)]

      fasta_path <- input$ref_fasta$datapath
      can_faidx <- nzchar(Sys.which("samtools")) && !grepl("\\.gz$", tolower(fasta_path))
      used_faidx <- FALSE

      if (can_faidx) {
        incProgress(0.2, detail = "Indexing FASTA (samtools faidx)")
        if (ensure_fai_index(fasta_path)) {
          fai <- read_fai_index(fasta_path)
          if (!is.null(fai) && nrow(fai) > 0) {
            header_map <- setNames(character(length(req_chr)), req_chr)
            unresolved <- character(0)
            for (q in req_chr) {
              h <- resolve_chr_header(q, fai$name)
              if (is.na(h)) {
                unresolved <- c(unresolved, q)
              } else {
                header_map[q] <- h
              }
            }
            if (!length(unresolved)) {
              used_faidx <- TRUE
              incProgress(0.3, detail = "Extracting flanks with indexed FASTA")
              len_map <- setNames(fai$len, fai$name)
              flankL <- input$flankL
              fa_lines <- c(); p3_blocks <- c()
              skipped <- r$flank_skipped
              for (i in seq_along(ids)) {
                h <- header_map[chr[i]]
                p_i <- suppressWarnings(as.integer(pos[i]))
                n_i <- suppressWarnings(as.integer(len_map[[h]]))
                reason <- NA_character_
                if (is.na(p_i)) reason <- "pos_missing"
                if (is.na(reason) && (is.na(n_i) || n_i < 1L)) reason <- "chr_len_missing"
                if (is.na(reason) && (p_i < 1L || p_i > n_i)) reason <- "pos_out_of_range"
                if (!is.na(reason)) {
                  skipped <- rbind(skipped, data.frame(
                    snp_id = ids[i],
                    chr_requested = chr[i],
                    pos = if (is.na(p_i)) NA_integer_ else p_i,
                    reason = reason,
                    used_faidx = TRUE,
                    stringsAsFactors = FALSE
                  ))
                  next
                }
                flk <- extract_flanks_faidx(fasta_path, h, p_i, n_i, flankL = flankL)
                if (is.na(flk$ref)) {
                  skipped <- rbind(skipped, data.frame(
                    snp_id = ids[i],
                    chr_requested = chr[i],
                    pos = p_i,
                    reason = "faidx_extraction_failed",
                    used_faidx = TRUE,
                    stringsAsFactors = FALSE
                  ))
                  next
                }
                fa_lines <- c(fa_lines, paste0(">", ids[i]), paste0(flk$left, flk$ref, flk$right))
                if (isTRUE(input$primer3_export)) p3_blocks <- c(p3_blocks, make_primer3_block(ids[i], flk))
              }
              r$flank_skipped <- skipped
              incProgress(0.3, detail = "Finalizing")
              r$flanks_fa <- paste0(fa_lines, collapse = "\n")
              r$primer3_txt <- if (isTRUE(input$primer3_export)) paste0(p3_blocks, collapse = "") else NULL
            } else {
              fasta_ids <- head(fai$name, 10)
              msg <- paste0(
                "FASTA index ready, but CHR unresolved: ",
                paste(head(unique(unresolved), 10), collapse = ", "),
                " | FASTA headers (first): ",
                paste(fasta_ids, collapse = ", ")
              )
              showNotification(msg, type = "warning", duration = 12)
            }
          }
        }
      }

      if (!used_faidx) {
        incProgress(0.2, detail = "Loading required FASTA chromosomes")
        if (is.null(r$fasta_cache_path) || !identical(r$fasta_cache_path, fasta_path)) {
          r$fasta_cache_path <- fasta_path
          r$fasta_cache <- character(0)
        }
        cached <- if (is.null(r$fasta_cache)) character(0) else r$fasta_cache
        cached_hdrs <- names(cached)
        missing_chr <- req_chr[!vapply(req_chr, function(q) is_chr_in_headers(q, cached_hdrs), logical(1))]
        if (length(missing_chr) > 0) {
          loaded <- read_fasta_targets(fasta_path, missing_chr)
          if (!length(loaded)) {
            fasta_ids <- read_fasta_header_ids(fasta_path, max_n = 20)
            req_preview <- paste(head(unique(req_chr), 10), collapse = ", ")
            fa_preview <- if (length(fasta_ids)) paste(head(fasta_ids, 10), collapse = ", ") else "(none)"
            msg <- paste0(
              "No requested chromosomes found. Requested: ", req_preview,
              " | FASTA headers (first): ", fa_preview
            )
            showNotification(msg, type = "error", duration = 12)
            return(NULL)
          }
          r$fasta_cache <- c(cached, loaded)
        }
        genome <- r$fasta_cache

        incProgress(0.3, detail = "Extracting flanks")
        flankL <- input$flankL
        fa_lines <- c(); p3_blocks <- c()
        skipped <- r$flank_skipped
        for (i in seq_along(ids)) {
          p_i <- suppressWarnings(as.integer(pos[i]))
          s_i <- get_chr_sequence(genome, chr[i])
          reason <- NA_character_
          if (is.null(s_i)) reason <- "chr_not_found"
          if (is.na(reason) && is.na(p_i)) reason <- "pos_missing"
          if (is.na(reason) && !is.null(s_i) && (p_i < 1L || p_i > nchar(s_i))) reason <- "pos_out_of_range"
          if (!is.na(reason)) {
            skipped <- rbind(skipped, data.frame(
              snp_id = ids[i],
              chr_requested = chr[i],
              pos = if (is.na(p_i)) NA_integer_ else p_i,
              reason = reason,
              used_faidx = FALSE,
              stringsAsFactors = FALSE
            ))
            next
          }
          flk <- extract_flanks(chr[i], pos[i], genome, flankL = flankL)
          if (is.na(flk$ref)) {
            skipped <- rbind(skipped, data.frame(
              snp_id = ids[i],
              chr_requested = chr[i],
              pos = p_i,
              reason = "extraction_failed",
              used_faidx = FALSE,
              stringsAsFactors = FALSE
            ))
            next
          }
          fa_lines <- c(fa_lines, paste0(">", ids[i]), paste0(flk$left, flk$ref, flk$right))
          if (isTRUE(input$primer3_export)) p3_blocks <- c(p3_blocks, make_primer3_block(ids[i], flk))
        }
        r$flank_skipped <- skipped
        incProgress(0.3, detail = "Finalizing")
        r$flanks_fa <- paste0(fa_lines, collapse = "\n")
        r$primer3_txt <- if (isTRUE(input$primer3_export)) paste0(p3_blocks, collapse = "") else NULL
      }

      if (!is.null(r$flank_skipped) && nrow(r$flank_skipped) > 0) {
        prev <- paste(head(r$flank_skipped$snp_id, 5), collapse = ", ")
        msg <- sprintf(
          "Generated flanks with %d skipped SNP(s): %s%s. Download skipped report for details.",
          nrow(r$flank_skipped),
          prev,
          if (nrow(r$flank_skipped) > 5) " ..." else ""
        )
        showNotification(msg, type = "warning", duration = 10)
      }
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

  output$dlFlankSkipped <- downloadHandler(
    filename = function() sprintf("flanking_skipped_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      validate(need(!is.null(r$flank_skipped), "No flanking run data found."))
      validate(need(nrow(r$flank_skipped) > 0, "No SNPs were skipped in the latest flanking run."))
      write.table(r$flank_skipped, file, sep = "\t", quote = FALSE, row.names = FALSE)
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
