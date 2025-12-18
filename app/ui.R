library(shiny)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .well { background:#f8f9fa; border-radius:12px; }
      .btn-primary { font-weight:600; }
      .code { font-family: monospace; }
      .muted { color:#6c757d; }
      .kpi { font-size: 16px; margin-right: 18px; display:inline-block; }
    "))
  ),
  titlePanel("SNPoptimizer — Minimal SNP set for genotype discrimination"),

  sidebarLayout(
    sidebarPanel(
      selectInput(
        "mode", "Section",
        choices = c("Discovery", "Evaluation (given SNPs)"),
        selected = "Discovery"
      ),
      tags$hr(),

      # ============ DISCOVERY ============
      conditionalPanel(
        condition = "input.mode == 'Discovery'",
        h4("Input & Parameters"),
        fileInput(
          "hapmap",
          "Upload HapMap or VCF",
          accept = c(".txt", ".tsv", ".hmp", ".vcf", ".vcf.gz"),
          buttonLabel = "Browse..."
        ),
        tags$small("Supported: HapMap (genotypes from column 12) and VCF (biallelic SNPs; GT required)."),
        br(), br(),

        sliderInput(
          "het_thresh",
          label = "Max heterozygosity per SNP (proportion)",
          min = 0, max = 0.50, value = 0.05, step = 0.01
        ),
        tags$small("Example: 0.05 = 5%"),
        tags$hr(),

        numericInput("k", "Number of SNPs to select (k)", value = 5, min = 1, step = 1),

        tags$hr(),
        h5("Genetic Algorithm (main run)"),
        numericInput("popSize", "Population size", value = 60, min = 10, step = 5),
        numericInput("maxiter", "Max generations (iterations)", value = 300, min = 10, step = 10),

        tags$hr(),
        h5("Reproducibility"),
        checkboxInput("useSeed", "Use fixed random seed", value = FALSE),
        conditionalPanel(
          condition = "input.useSeed == true",
          numericInput("seedValue", "Random seed value", value = 123, min = 1, step = 1)
        ),

        tags$hr(),
        h5("Resolve remaining duplicates (optional)"),
        checkboxInput("doAdditional", "Run a second discovery on duplicates", value = FALSE),
        conditionalPanel(
          condition = "input.doAdditional == true",
          numericInput("n_additional", "Additional SNPs to search", value = 3, min = 1, step = 1)
        ),

        tags$hr(),
        actionButton("run", "Run discovery", class = "btn btn-primary btn-lg", width = "100%"),

        tags$hr(),
        h5("Post-discovery: flanking & primers"),
        fileInput("ref_fasta", "Reference genome (FASTA or .fa.gz)", accept = c(".fa", ".fasta", ".fa.gz")),
        numericInput("flankL", "Flanking length (bp, each side)", value = 200, min = 20, step = 10),
        checkboxInput("primer3_export", "Prepare Primer3 input", value = TRUE),
        actionButton("makeFlanks", "Generate flanking sequences", class = "btn btn-primary", width = "100%"),
        br(), br(),
        downloadButton("dlFlanksFA", "Download flanking FASTA"),
        downloadButton("dlPrimer3",  "Download Primer3 file")
      ),

      # ============ EVALUATION ============
      conditionalPanel(
        condition = "input.mode == 'Evaluation (given SNPs)'",
        h4("Evaluate SNP list (CHR, POS)"),
        helpText("Upload a dataset (HapMap/VCF) for evaluation, or leave empty to reuse the dataset loaded in 'Discovery'."),
        fileInput(
          "eval_dataset",
          "Dataset for evaluation (HapMap/VCF)",
          accept = c(".txt", ".tsv", ".hmp", ".vcf", ".vcf.gz")
        ),
        tags$hr(),
        helpText("Upload or paste a list of SNPs present in the dataset. Accepted formats:"),
        tags$ul(
          tags$li("TSV/TXT/CSV file with columns: CHR and POS"),
          tags$li("Or lines like: 1 12345  or  chr1:12345")
        ),
        fileInput("snp_list_file", "Upload SNP list", accept = c(".txt", ".tsv", ".csv")),
        textAreaInput("snp_list_text", "…or paste here", rows = 6,
                      placeholder = "1 12345\n1:67890\nchr2 34567"),
        actionButton("evalBtn", "Evaluate SNP list", class = "btn btn-primary", width = "100%"),
        tags$hr(),
        downloadButton("dlEvalMatrix", "Download subset matrix")
      )
    ),

    mainPanel(
      # ------------- DISCOVERY PANEL -------------
      conditionalPanel(
        condition = "input.mode == 'Discovery'",
        h3("Discovery results"),

        uiOutput("discoveryStats"),
        tags$hr(),

        h4("SNP Discrimination Curve"),
        plotOutput("discCurve", height = "420px"),
        tags$small(class = "muted",
                   "Y: unique genotypes (profiles) — X: number of SNPs (greedy order)"),
        tags$hr(),

        h4("Selected SNP table"),
        div(
          style="margin-bottom:10px",
          downloadButton("dlSelectedMeta", "Download selected SNP metadata"),
          downloadButton("dlSelectedMatrix", "Download selected genotype matrix"),
          downloadButton("dlDiscoverySummary", "Download discovery summary"),
          downloadButton("dlPlot", "Download plot"),
          downloadButton("dlDupSamples", "Download duplicated samples")
        ),
        tableOutput("selectedSNPTable"),
        tags$hr(),

        h5("Remaining duplicate groups"),
        verbatimTextOutput("dupGroupsTxt")
      ),

      # ------------- EVALUATION PANEL -------------
      conditionalPanel(
        condition = "input.mode == 'Evaluation (given SNPs)'",
        uiOutput("evalPanel")
      )
    )
  )
)
