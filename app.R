### Transcription Unit Annotation
### Rui Shao
### Adapt annotation scripts from Michael Lidschreiber, TT-seq 2016 Science paper

### part1 binning genome and count bam file, generate tu Granges
### part2 merge tu upon indicated features and annotate non-coding tu genomic locations
### part3 normalise TU expression, and cutoff low expressed ncRNAs

if (!require("pacman")) install.packages("pacman")
pacman::p_load(shiny, shinyjs, shinydashboard, bsplus, shinythemes, shinyFiles,
               dplyr, S4Vectors,
               Rsamtools, rtracklayer,
               IRanges, GenomicRanges, GenomeInfoDb, GenomicAlignments,
               STAN,
               foreach, doParallel,
               ggplot2, cowplot)

#setting----------------------------------------------------------------------

options(shiny.maxRequestSize=30*1024^15)
options("scipen"=999,"digits"=4)
# options(repos = BiocInstaller::biocinstallRepos())
options(warn=-1)

#ui----------------------------------------------------------------------

# Define UI for application that draws a histogram
ui <- dashboardPage(skin = "black",
                    dashboardHeader(title = "TU filter"),
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Input", tabName = "filesupload", icon = icon("file")),
                        menuItem("Setting", tabName = "setting", icon = icon("cog")),
                        menuItem("GenoSTAN", tabName = "STAN", icon = icon("align-justify")),
                        menuItem("Output", tabName = "stats", icon = icon("th")),
                        menuItem("Help", tabName = "help", icon = icon("question-circle"))
                      )
                    ),
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "filesupload",
                                fluidRow(
                                  box(title = 'Submit',background = "teal", solidHeader = TRUE,width = 4,
                                      uiOutput('buttonsUI4'),
                                      uiOutput('buttonsUI2'),
                                      collapsible = TRUE, collapsed = FALSE),
                                  
                                  box(title='Reference file',status = "info",solidHeader = TRUE,width = 4,
                                      shinyFilesButton("gene_ref", "Ensembl or GENCODE" ,multiple = F,
                                                       title = "Please select a reference file:",
                                                       buttonType = "default",class=NULL,tags$span(style="color:red", "*")),
                                      verbatimTextOutput('refInput'),br(),actionButton('load','Load'),
                                      collapsible = TRUE, collapsed = FALSE
                                  ),
                                  
                                  box(title = "Input bam files",status = "info",solidHeader = TRUE,width = 4,
                                      
                                      shinyFilesButton("bamfiles1", "Sample 1, bam file(s)", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default",class = NULL),
                                      verbatimTextOutput('filepaths1'),br(),
                                      
                                      shinyFilesButton("bamfiles2", "Sample 2, bam file(s)", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('filepaths2'),br(),
                                      
                                      shinyFilesButton("bamfiles3", "Sample 3, bam file(s)", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('filepaths3'),br(),
                                      
                                      shinyFilesButton("bamfiles4", "Sample 4, bam file(s)", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('filepaths4'),br(),
                                      
                                      shinyFilesButton("bamfiles5", "Sample 5, bam file(s)", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('filepaths5'),br(),
                                      
                                      shinyFilesButton("bamfiles6", "Sample 6, bam file(s)", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default",class = NULL),
                                      verbatimTextOutput('filepaths6'),br(),
                                      
                                      shinyFilesButton("bamfiles7", "Sample 7, bam file(s)", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('filepaths7'),br(),
                                      
                                      shinyFilesButton("bamfiles8", "Sample 8, bam file(s)", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('filepaths8'),br(),
                                      
                                      shinyFilesButton("bamfiles9", "Sample 9, bam file(s)", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('filepaths9'),br(),
                                      
                                      shinyFilesButton("bamfiles10", "Sample 10, bam file(s)" , multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('filepaths10'),br(),
                                      collapsible = TRUE, collapsed = FALSE
                                  ),
                                  box(title = 'Location preview', plotOutput("counts_plot"), collapsible = TRUE,
                                      width = 8, collapsed = FALSE,
                                      background = "teal", solidHeader = TRUE #status = "warning",solidHeader = TRUE
                                  )
                                )
                        ),
                        tabItem(tabName = "setting",
                                fluidRow(
                                  box(title='Option2: Annotate Transcribed Regions',status = "info",solidHeader = TRUE,
                                      fileInput("tx_file", "Upload a Transcript File",
                                                accept = c("text",".bed",".gtf"),placeholder = 'BED/GTF'
                                      )%>%
                                        shinyInput_label_embed(
                                          icon("question-circle") %>%
                                            bs_embed_tooltip(title = "e.g.Peak_calling_from_Homer")
                                        ),collapsible = T, collapsed = FALSE),
                                  
                                  box(title='Zooms',status = "warning",solidHeader = TRUE,
                                      
                                      sliderInput( "L_Promoter",
                                                   label = "Length of Promoter:",
                                                   min = 1,
                                                   max = 3000,
                                                   value = 1000),
                                      sliderInput("L_TSS",
                                                  label = "Length of TSS:",
                                                  min = 1,
                                                  max = 3000,
                                                  value = 1000),
                                      sliderInput( "L_TTS",
                                                   label = "Length of TTS:",
                                                   min = 1,
                                                   max = 3000,
                                                   value = 1000),
                                      numericInput(inputId = "binning.size", label = "Bin size", value = 200,
                                                   min = 50, max = 2000,step = 100,width = "60%")
                                  ),
                                  
                                  box(title = 'Annotation tuning', status = "warning",solidHeader = TRUE,
                                      
                                      fileInput("unmappable", HTML(paste('Join TU gaps' )),
                                                multiple = FALSE,
                                                accept = c("text/csv",".bed",".gtf"), placeholder = 'BED'
                                      ) %>%
                                        shinyInput_label_embed(
                                          icon("question-circle") %>%
                                            bs_embed_tooltip(title = "Unmappable_regions")
                                        ),
                                      
                                      checkboxInput("precBound", "Precise Boundary", value = FALSE, width = NULL)%>%
                                        shinyInput_label_embed(
                                          icon("question-circle") %>%
                                            bs_embed_tooltip(title = "Refine non-coding TU boundary by coverage
                                                             abrupt changes. Slow warning.")
                                        ),
                                      
                                      checkboxInput("jaccardCutoff", "Filter Common TUs", value = FALSE, width = NULL) %>%
                                        shinyInput_label_embed(
                                          icon("question-circle") %>%
                                            bs_embed_tooltip(title = "Cutoff expression below the threshold which 
                                            keeps the max jaccard similarity among all samples when lower expressed
                                                             TUs are removed.")
                                        ),
                                      
                                      collapsible = TRUE, collapsed = FALSE
                                  ),
                                  
                                  box(title='TU and Gene Cutoff', status = "warning",solidHeader = TRUE,
                                      checkboxGroupInput('merging_TU', label = 'Merge TU from genes ',
                                                         choices = c('Minimum TU overlap'='TU',
                                                                     'Minimum gene overlap'='Gene')),
                                      uiOutput(outputId = "overlap_cutoff_"),
                                      uiOutput(outputId = "jaccard_cutoff_"),
                                      useShinyjs(),
                                      fileInput(inputId = "mask_list", label = "Upload a mask list",
                                                placeholder = 'gene_id to remove'
                                      ) %>%
                                        shinyInput_label_embed(
                                          icon("question-circle") %>%
                                            bs_embed_tooltip(title = "Extremely long genes")
                                        ),collapsible = TRUE, collapsed = FALSE
                                  )
                                )),
                        tabItem(tabName = "STAN",
                                fluidRow(
                                  box(title = "ChIP-seq bam files",status = "info",solidHeader = TRUE, width = 6,
                                      
                                      shinyFilesButton("C_bamfiles1", "Total bam files 1", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default",class = NULL),
                                      verbatimTextOutput('C.filepaths1'),br(),
                                      
                                      shinyFilesButton("C_bamfiles2", "Total bam files 2", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('C.filepaths2'),br(),
                                      
                                      shinyFilesButton("C_bamfiles3", "Total bam files 3", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('C.filepaths3'),br(),
                                      
                                      shinyFilesButton("C_bamfiles4", "Total bam files 4", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('C.filepaths4'),br(),
                                      
                                      shinyFilesButton("C_bamfiles5", "Total bam files 5", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('C.filepaths5'),br(),
                                      
                                      shinyFilesButton("C_bamfiles6", "Total bam files 6", multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('C.filepaths6'),br(),
                                      
                                      shinyFilesButton("C_bamfiles7", "Total bam files 7" ,multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('C.filepaths7'),br(),
                                      
                                      shinyFilesButton("C_bamfiles8", "Total bam files 8" ,multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('C.filepaths8'),br(),
                                      
                                      shinyFilesButton("C_bamfiles9", "Total bam files 9" ,multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('C.filepaths9'),br(),
                                      
                                      shinyFilesButton("C_bamfiles10", "Total bam files 10" ,multiple = TRUE,
                                                       title = "Please select at least one bam file:",
                                                       buttonType = "default", class = NULL),
                                      verbatimTextOutput('C.filepaths10'),br(),
                                      collapsible = TRUE, collapsed = FALSE
                                  )
                                )
                        ),
                        tabItem(tabName = "stats",
                                fluidRow(
                                  box(title = 'Plot',uiOutput('buttonsUI3'), status = "primary",solidHeader = TRUE,
                                      collapsible = TRUE, collapsed = FALSE),
                                  
                                  box(title='Download', uiOutput('downloadUI'),status = "primary",solidHeader = TRUE,
                                      collapsible = TRUE, collapsed = FALSE),
                                  
                                  box(title = "Location table", tableOutput("table_output"),status = "primary",
                                      solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE)
                                )
                        ),
                        tabItem(tabName = "help",
                                fluidRow(
                                  box(title = "File requirments",status = "success",
                                      p('This app organizes the annotatation steps of TT-seq paper, combines TU 
                                        (transcription unit) calling and ncRNA location assignment for auto-annotation and parallelization.'),
                                      a(href="http://science.sciencemag.org/content/352/6290/1225.long",
                                        h5("TTseq paper")),
                                      p('Input bam files will be subject to binning genome, calling the active state by GenoSTAN HMM, and
                                        merging intervals upon the provided references and features.'),
                                      p('Step 1: bin genome, count reads, call TU ranges'),
                                      p('Step 2: merge scatter TUs upon indicated features, and assign non-coding TU locations'),
                                      p('Step 3: normalise TU expression if multiple samples are inputs, and cutoff lowly expressed ncRNAs'),
                                      br(),
                                      
                                      img(src = 'flow.png',  width = 500),br(),
                                      a(href="www.bioconductor.org/packages/release/bioc/html/STAN.html",
                                        h5("Bioconductor: STAN")),
                                      br(),
                                      
                                      p("Requirements:"),
                                      p('Bam files: paired-end and strand specific. First / second read in pair will be auto-corrected
                                        to match the forward strand.'),
                                      
                                      br(),
                                      p('Transcripts file: a list of BED / GTF intervals to be annotated, e.g. peaks from HOMER'),
                                      
                                      br(),
                                      p("Gene reference: each gene reference has the 'gene' and 'exon' features, and 'protein_coding'
                                        gene_type, e.g. GENCODE and ensembl. "),
                                      br(),
                                      img(src = 'example.png',  width = 700),br(),
                                      width = 600,
                                      collapsible = TRUE
                                      
                                  ),
                                  box(title = "Annotation methods",status = "success",#solidHeader = TRUE,
                                      p('After calling / uploading TUs, the scattered interval will be merged upon
                                        the indicated type of references:'),br(),
                                      strong('TU assembling method: '),
                                      br(),
                                      img(src = 'steps.png',  width = 700),br(),
                                      strong('Name the overlapped features.'),
                                      p("TUs overlapping the chosen features, 'Skip' merging, or 'Join TU' (keeps
                                        contineous interval), 'Join exons' (splits at gene boundary, up / down 
                                        stream parts will be named separately) "), 
                                      br(),
                                      
                                      p("(Note: redundant reference genes will interfere the reults, e.g. assigning
                                        TU to the name of embeded pseudogene.)"),
                                      br(),
                                      
                                      p('Non-coding RNA locations:'),
                                      img(src='abbre.png',  width = 700),br(),
                                      width = 600,
                                      collapsible = TRUE
                                      
                                  ),
                                  box(title = "Parameter setting", status = "success",#solidHeader = TRUE,
                                      p('Lowly expressed genes usually generate scattered TUs, parameters can help to merge 
                                        the intervals of high certainty.'),
                                      br(),
                                      
                                      strong('Joining methods: '),
                                      img(src='selection.png',  width = 700),br(),
                                      p('Transcript isoforms will not be accounted.'),
                                      br(),
                                      strong('Join TU gaps: '),
                                      p('Scattered intervals gapped by unmappable regions can be joined if a
                                        genome mappability file is supplied.'),
                                      br(),
                                      
                                      strong('Minimum overlap: '),
                                      p('To remove the lowly expressed TUs on gene reference.'),
                                      br(),
                                      strong('Mask list: '),
                                      p('To remove extraordinarily wide reference genes that covers many other genes.'),
                                      br(),
                                      width = 700,
                                      collapsible = TRUE, collapsed = FALSE
                                      
                                  ),
                                  
                                  box(title = "Info",status = "success",#solidHeader = TRUE,
                                      p('Michael Lidschreiber wrote the original annotation scripts, Rui Shao
                                        adapted and compiled to the shiny app.'),
                                      br(),
                                      p('Contact rui.shao@scilifelab.se if there is any question.'),
                                      width = 700,
                                      collapsible = TRUE, collapsed = FALSE)
                                )
                        )
                      )
                    )
)
############################################
###
###  server
###
############################################
# Define server logic required to draw a histogram
server <- shinyServer(function(input, output,session) {
  # variables to control the sequence of processes
  controlVar <<- reactiveValues(bamReady = FALSE,
                                fileReady1 = FALSE,
                                fileReady2 = FALSE,
                                exprcutoff = FALSE,
                                positionReady = FALSE,
                                plotReady = FALSE,
                                tableReady = FALSE,
                                heatmapReady = FALSE,
                                heatmapChoice = FALSE,
                                submit_trigger = FALSE,
                                nc_reduced = FALSE,
                                nc_given = FALSE,
                                gene_given = FALSE,
                                gene_merged = FALSE)
  # prepare bam inputs
  volumes <- getVolumes()
  # labeled RNA
  shinyFileChoose(input, "bamfiles1", roots=volumes, filetypes=c('bam'), session = session)
  shinyFileChoose(input, "bamfiles2", roots=volumes, filetypes=c('bam'), session = session)
  shinyFileChoose(input, "bamfiles3", roots=volumes, filetypes=c('bam'), session = session)
  shinyFileChoose(input, "bamfiles4", roots=volumes, filetypes=c('bam'), session = session)
  shinyFileChoose(input, "bamfiles5", roots=volumes, filetypes=c('bam'), session = session)
  shinyFileChoose(input, "bamfiles6", roots=volumes, filetypes=c('bam'), session = session)
  shinyFileChoose(input, "bamfiles7", roots=volumes, filetypes=c('bam'), session = session)
  shinyFileChoose(input, "bamfiles8", roots=volumes, filetypes=c('bam'), session = session)
  shinyFileChoose(input, "bamfiles9", roots=volumes, filetypes=c('bam'), session = session)
  shinyFileChoose(input, "bamfiles10", roots=volumes, filetypes=c('bam'), session = session)
  
  # reference
  shinyFileChoose(input, "gene_ref", roots=volumes, filetypes=c('bed','gtf'), session = session)
  
  # labeled RNA input
  L.sample.list <<- list()
  observeEvent(input$bamfiles1, {
    L.sample.list$bamfiles1 <<- parseFilePaths(roots = volumes, input$bamfiles1)$datapath
    output$filepaths1 <- renderPrint({(L.sample.list$bamfiles1)})
  })
  observeEvent(input$bamfiles2, {
    L.sample.list$bamfiles2 <<- parseFilePaths(roots = volumes, input$bamfiles2)$datapath
    output$filepaths2 <<- renderPrint({L.sample.list$bamfiles2})
  })
  observeEvent(input$bamfiles3, {
    L.sample.list$bamfiles3 <<- parseFilePaths(roots = volumes, input$bamfiles3)$datapath
    output$filepaths3 <<- renderPrint({L.sample.list$bamfiles3})
  })
  observeEvent(input$bamfiles4, {
    L.sample.list$bamfiles4 <<- parseFilePaths(roots = volumes, input$bamfiles4)$datapath
    output$filepaths4 <<- renderPrint({L.sample.list$bamfiles4})
  })
  observeEvent(input$bamfiles5, {
    L.sample.list$bamfiles5 <<- parseFilePaths(roots = volumes, input$bamfiles5)$datapath
    output$filepaths5 <<- renderPrint({L.sample.list$bamfiles5})
  })
  L.sample.list<<-list()
  observeEvent(input$bamfiles6, {
    L.sample.list$bamfiles6 <<- parseFilePaths(roots = volumes, input$bamfiles6)$datapath
    output$filepaths6 <- renderPrint({(L.sample.list$bamfiles6)})
  })
  observeEvent(input$bamfiles7, {
    L.sample.list$bamfiles7 <<- parseFilePaths(roots = volumes, input$bamfiles7)$datapath
    output$filepaths7 <<- renderPrint({L.sample.list$bamfiles7})
  })
  observeEvent(input$bamfiles8, {
    L.sample.list$bamfiles8 <<- parseFilePaths(roots = volumes, input$bamfiles8)$datapath
    output$filepaths8 <<- renderPrint({L.sample.list$bamfiles8})
  })
  observeEvent(input$bamfiles9, {
    L.sample.list$bamfiles9 <<- parseFilePaths(roots = volumes, input$bamfiles9)$datapath
    output$filepaths9 <<- renderPrint({L.sample.list$bamfiles9})
  })
  observeEvent(input$bamfiles10, {
    L.sample.list$bamfiles10 <<- parseFilePaths(roots = volumes, input$bamfiles10)$datapath
    output$filepaths10 <<- renderPrint({L.sample.list$bamfiles10})
  })
  
  # load reference file 
  observeEvent(input$gene_ref, {
    controlVar$fileReady1 <- FALSE
    ref.path <<- parseFilePaths(roots = volumes, input$gene_ref)$datapath
    output$refInput <- renderPrint({((ref.path))})
  })
  
  observeEvent(input$tx_file, {
    controlVar$fileReady2 <- FALSE
    TU_input <<- import.input.ranges(input$tx_file, TRUE)
    controlVar$fileReady2 <- TRUE
  })
  
  # process reference
  observeEvent(input$load, {
    if(length(ref.path)>0)
    {
      showNotification("Loading reference ...", duration = 30)
      Gene_input<<-import.input.ranges(ref.path,F)
      
      if(!is.null(Gene_input$gene_biotype))
      {
        colnames(mcols(Gene_input))[ colnames(mcols(Gene_input)) == 'gene_biotype' ] <<- 'gene_type'
      }
      
      if(!is.null(Gene_input$source))
      {
        colnames(mcols(Gene_input))[ colnames(mcols(Gene_input)) == 'source'] <<- 'gene_source'
      }
      
      removeNotification("Loading reference ...")
      controlVar$fileReady1 <<- TRUE
    }
  })
  
  # pop up setting if box is checked
  output$overlap_cutoff_<-renderUI({
    #Merging settings
    overlap_types <- c('TU','Gene')
    ov_indx <- paste0("overlap_cutoff_", overlap_types)
    
    Cutoff_inputs <<- lapply(1:2, function(i){
      numericInput(inputId = ov_indx[i], label = paste(overlap_types[i],"overlap cutoff"), value = 0, 
                   min = 0, max = 1,step = 0.01,width = "40%")
    }) %>% tagList()
    
    shinyjs::hidden(Cutoff_inputs)
    
  })
  
  # pop up which selection
  observe({
    ovCutoff <<- input$overlap_cutoff_TU
    annoCutoff <<- input$overlap_cutoff_Gene
    overlap_types <- c('TU','Gene')
    div(
      for(i in overlap_types){
        if (i %in% input$merging_TU) {
          shinyjs::show(id = paste0("overlap_cutoff_", i))
        } else {
          shinyjs::hide(id = paste0("overlap_cutoff_", i))
        }
      }
    )
  })
  
  #UI buttons----------------------------------------------------------------------
  # show buttons only when file is uploaded
  output$buttonsUI2 <- renderUI({
    if (controlVar$fileReady1  & (controlVar$fileReady2 | !is.null(input$bamfiles1)))
      div(
        selectInput(inputId = "stan_mode",
                    multiple = TRUE,
                    HTML(paste("GenoSTAN model")),
                    choices = c("NegativeBinomial" = "NB", "PoissonLogNormal" = "Poilog",
                                "ZINegativeBinomial" = "ZINB","IndependentGaussian" = "Gauss"),
                    width = '90%', selected = 'Poilog'),
        selectInput(inputId = "mergeMethod",
                    multiple = FALSE,
                    HTML(paste("TU merge method",tags$span(style = "color:red", "*"))),
                    choices = c("Skip joining TU" = "skip","Join TUs" = "tu","Join exons" = "exon"),
                    width = '90%', selected = 'skip'),
        
        collapsible = TRUE, collapsed = FALSE,
        # numericInput(inputId = "cutoff_quantile", label = "Cutoff Expression Quantile", value = 0.1, min = 0, max = 5,step = 0.01,width = "80%"),
        actionButton('go','Run')
      )
  })
  
  output$buttonsUI4 <- renderUI({
    if (controlVar$fileReady1)
      div(
        selectInput(inputId = "merge_features",
                    multiple = TRUE, selected="protein_coding",
                    HTML(paste("Multi-exon transcript types",tags$span(style="color:red", "*"))),
                    unique(Gene_input$gene_type),width = '90%')
      )
  })
  
  #Event go----------------------------------------------------------------------
  ### Part1: binning and counting
  
  observeEvent(input$go, {
    #initialize controlVars
    
    #conditional controls
    if (!exists('all.binned.counts.list'))controlVar$bamReady <<- FALSE #won't count bam file twice
    controlVar$exprcutoff <<- FALSE
    
    if (!exists('previous_mode')) {
      controlVar$stanset <<- FALSE
    } else if (previous_mode!=input$stan_mode) {
      controlVar$stanset <<- FALSE
    } else {controlVar$stanset <<- TRUE}
    
    if(!exists('previous_cutfactor')){
      controlVar$expcutoff<<- FALSE
    } else if (previous_cutfactor!=c(minSizeFactor,minRPKM)) {
      controlVar$expcutoff <<- FALSE
    } else {controlVar$expcutoff<<-TRUE}
    
    if (!exists('anno.gr')) controlVar$TUgrReady <<- FALSE
    if (exists('anno.gr') & exists('TU_input')) rm(TU_input)
    if (!exists('merge_method')) {
      controlVar$gene_given <<- FALSE #gene_merged won't re-initialize if nothing changes
    } else if (mergeMethod!=input$mergeMethod) {
      controlVar$gene_given <<- FALSE}
    
    # renewable controls
    controlVar$nc_reduced <<- FALSE
    controlVar$nc_given <<- FALSE
    controlVar$location <<- FALSE
    controlVar$positionReady <<- FALSE
    controlVar$tableReady <<- FALSE
    
    # set params--------
    if (T) {
      binning <<- input$binning.size
      nChrs <<-length(seqlevels(Gene_input))
      ncores <<- ifelse(detectCores() > nChrs, 20, detectCores() - 1)
      registerDoParallel(cores = ncores)
      
      outCoverage <<- FALSE
      mergeFeatures <<- input$merge_features
      fit_method <<- input$stan_mode
      nStates <<- 2
      mergeMethod <<- input$mergeMethod
      if (is.null(input$merging_TU)) {
        ovCutoff <<-0
        annoCutoff <<-0
      }
      txInput <<- input$tx_file
      precBound <<- input$precBound
      
      TSS_len <<- input$L_TSS
      Promoter_len <<- input$L_Promoter
      TTS_len <<- input$L_TTS
      
      jaccardCutoff <<- input$jaccardCutoff
      
      if (is.null(input$mask_list)) {
        mask_list <<- NULL
      } else { mask_list <<- input$mask_list}
      
      if (is.null(input$unmappable)) {
        unmappable <<- NULL
      } else { unmappable <<- input$unmappable }
    }
    
    if (xor(controlVar$fileReady2, controlVar$TUgrReady)) {
      TU_input$id = 1:length(TU_input)
      TU_input$gene_id = NA
      TU_input$transcript_id = NA
      TU_input$exonOV = FALSE
      TU_input$type = "ncRNA"
      TU_input$gene_type = NA
      TU_input <<- TU_input
    }
    
    # read bam files (or load transcript intervals)
    if (!is.null(txInput)) {
      TU_input <- import.input.ranges(txInput,is.datapath = TRUE)
      chr.names <<- paste0('chr', c(1:100,'X','Y'))[paste0('chr', c(1:100,'X','Y')) %in% seqlevels(TU_input)]
      gene.gr <<- gene_prep(Gene_input, TU_input)
      TU.gr <<- TU_prep(Gene_input, gene.gr, TU_input)
    } else if (length(L.sample.list)>0) 
    {
      start_time = Sys.time()
      sample.names <<- lapply(L.sample.list, function(x){
        strsplit(tail(unlist(strsplit(x[1],'\\/')),1),'\\.bam')%>%unlist()%>%head(1)
      })
      showNotification("Step1 Binning genome.", duration = 2000)
      
      all.TU.gr.list = list()
      
      # index bam files
      bais = lapply(L.sample.list, function(x) sapply(paste0(x,'.bai'), function(y) 
        if(!file.exists(y)) Rsamtools::indexBam(x)))
      
      all.readPaired <<- invisible(lapply(L.sample.list, function(x) sapply(x, testPairedEndBam)) )
      
      all.strandFlipped <<- setNames(rep(FALSE, length(L.sample.list)), names(L.sample.list)) # initiate
      
      all.binned.counts.list <<- bam_bin_list(L.sample.list, FUN=bamBinCount, all.readPaired)
      
      ### Part2: calling TU -----------------------------------------------------------------------------
      showNotification("Step2 Calling TU.", duration = 2000)
      for (bams in names(L.sample.list)) {
        strandFlipped <- FALSE
        TU.gr <- bin_transcribed_TU(binned.RNA.counts.list = all.binned.counts.list[[bams]],
                                    Gene_input = Gene_input,
                                    strandFlipped = FALSE,
                                    bam.files = L.sample.list[[bams]])
        # strandFlipped from function "TU_prep"; if any experiment prepared first-in-pair read on reverse strand, then "strandFlipped" is TRUE
        if (strandFlipped)
        {
          all.binned.counts.list[[bams]] <- lapply(all.binned.counts.list[[bams]],
                                                   function(x) x[, rev(seq_len(ncol(x)))])
          all.strandFlipped[bams] <- TRUE
        }
        all.TU.gr.list <- c(all.TU.gr.list, list(TU.gr))
      }
      
      # removeNotification("Step2 Calling TU ...")
      showNotification(paste("TU calling cost",
                             round(difftime(Sys.time(), start_time, units='mins'),
                                   digits = 2),"min"),
                       duration = 2000)
      
      names(all.TU.gr.list) <- names(L.sample.list)
      all.TU.gr.list <<- all.TU.gr.list
      gc()
      
      previous_mode <<- input$stan_mode
      controlVar$bamReady <<- TRUE
    } else {
      stop('Error: please input bam file or transcript. Stop.\n')
    }
    
    ### part3: transript unit expression cutoff by exon counts -----------------------------------------
    
    if (controlVar$bamReady) {
      showNotification("Step3 Expression cutoff ...", duration = 2000)
      
      all.TU.expr.list <<- cutoff_expression(L.sample.list, all.TU.gr.list, all.binned.counts.list)
      
      removeNotification("Step3 Expression cutoff.")
    }
    
    controlVar$TUgrReady <<- TRUE
    showNotification("Done! [^_^]",duration = 2000)
    registerDoSEQ()
  })
  # The end of annotation
  
  # prepare output ----------------------------------------------------------------------
  observeEvent(input$go, {
    if(controlVar$TUgrReady){
      location.table <- list()
      if (!is.null(txInput)) {
        location.table <<- list(table(TU.gr$location))
      } else {
        for (i in seq_along(all.TU.expr.list)) {
          temp.table <- table(all.TU.expr.list[[i]]$location)
          location.table <<- c(location.table, 
                               list(temp.table[order(temp.table, decreasing = TRUE)]) )
        }
      }
      controlVar$tableReady <<- TRUE
    }
  })
  
  # table_output
  output$table_output <- renderTable({
    if(controlVar$tableReady){
      cbind("Location" = names(location.table[[1]]), 
            "Number" = unname(location.table[[1]]))
    }
  })
  
  output$counts_plot <- renderPlot({
    controlVar$plotReady <<- FALSE
    if(controlVar$tableReady){
      
      glist <- list()
      
      for (i in seq_along(location.table)) {
        new.table <- as.data.frame(location.table[[i]])
        colnames(new.table) <- c("Location", "Number")
        new.table$Location = factor(new.table$Location, 
                                    levels = new.table$Location[order(new.table$Number, decreasing = T)])
        
        g <- ggplot(new.table, aes(x = Location, y = Number, fill = Location)) +
          geom_bar(stat = 'identity', width = 0.8) +
          theme_minimal() +
          theme(axis.title = element_text(size = 14),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                panel.grid.major = element_blank(),
                plot.title = element_text(size = 14, face="bold", vjust = 0.5),
                axis.text = element_text(size = 12, angle = 90),
                axis.text.y = element_text(size = 10, angle = 0)) +
          labs(title = ifelse(exists("sample.names"), sample.names[[i]], "TU by location"),
               y = "Annotated TU number", x = " ")
        
        glist <- c(glist, list(g))
      }
      
      print(cowplot::plot_grid(plotlist = glist, ncol = 1, align = "v"))
      
      controlVar$plotReady <<- TRUE
    }
  })
  
  output$downloadUI <- renderUI({
    if (controlVar$tableReady)
      div(
        textInput("write_gff_path", "Export to:", value = "~/Downloads"),
        downloadButton("downloadReport", "Save TU anno and parameter log.")
      )
  })
  
  
  output$downloadReport <- downloadHandler(
    
    filename = paste0("TU annotation report ", Sys.Date(), ".log"), 
    content = function(file) {
      
      # download TU annotation
      if (!is.null(txInput)) {
        out_file <- paste0(input$write_gff_path, 
                           "/TU_filter_anno_", 
                           gsub("(.*)\\..*$", "\\1", txInput[[1]]),
                           ".gtf")
        export.gff3(object = TU.gr, con = out_file)
      } else {
        for (i in seq_along(all.TU.expr.list)) {
          file_name <- paste0(input$write_gff_path, 
                              '/TU_filter_anno_',
                              sample.names[[i]], 
                              ".gff3")
          export.gff3(all.TU.expr.list[[i]], con = file_name)
        }
      }
      
      
      
      file.create("tmp_report.txt")
      temp_report = file.path(tempdir(), "tmp_report.txt")
      file.copy(from = "tmp_report.txt", to = temp_report, overwrite = TRUE)
      
      # Pass parameters to the report
      if (length(L.sample.list) > 0) {
        bam.table = data.frame("Names" = names(L.sample.list),
                               "Inputs" = sapply(L.sample.list, function(x) paste(x, collapse = "\n")))
      } else {
        bam.table = data.frame("Names" = NA,
                               "Inputs" = NA)
      }
      
      anno.table = data.frame(Names = c("Reference annotation",
                                        "Transcripts input"), 
                              Inputs = c(ref.path,
                                         ifelse(is.null(input$tx_file), NA, unlist(input$tx_file)) )
                              )
      
      
      params.table = data.frame(Names = c("Annotation mode",
                                          "Minimum TU overlap",
                                          "Minimum gene overlap",
                                          "TSS downstream length",
                                          "Promoter length",
                                          "TES downstream length",
                                          "Expression cutoff"),
                                
                                Inputs = c(ifelse(is.null(input$de_novo_position_features),
                                                  NA, paste(input$de_novo_position_features, sep = ';')),
                                           ovCutoff,
                                           annoCutoff,
                                           ifelse(is.null(TSS_len), NA, TSS_len),
                                           ifelse(is.null(Promoter_len), NA, Promoter_len),
                                           ifelse(is.null(TTS_len), NA, TTS_len),
                                           ifelse(is.null(input$expr_cutoff), NA, input$expr_cutoff))
                                )
      
      params = list(total_table = rbind(bam.table,
                                        anno.table,
                                        params.table) )
      
      write.table(params[[1]], file, quote = FALSE, sep = "\t", row.names = FALSE)
      
      file.remove("tmp_report.txt")
    }
  )
  
  # -------------------------------------------------------------------------------------#
  #                                                                                      #  
  #                                     Functions                                        #
  #                                                                                      #
  # -------------------------------------------------------------------------------------#
  
  ### part1 ---------------------------------------------------------------------
  # Convert bam files to the short-read genome binned list objects
  bamBinCount <<- function(bam.file, paired.end, ...){
    # Args:
    # bam.file: single file path
    # paired.end: reads single or paired end
    # outCoverage: if TRUE reads coverage on bins, otherwise count reads mid-point in bins
    bam.index = paste0(bam.file,'.bai')
    if (!file.exists(bam.index)) bam.index = Rsamtools::indexBam(files = bam.file)
    
    if (paired.end) {
      sbw = c('pos', 'qwidth','strand','rname', 'mrnm', 'mpos', 'isize')
      flag=scanBamFlag(isFirstMateRead = paired.end, isSecondaryAlignment = FALSE)
    } else {
      sbw = c('pos', 'qwidth','strand','rname')
      flag=scanBamFlag(isSecondaryAlignment = FALSE,
                       isFirstMateRead = paired.end)
    }
    param = ScanBamParam(what = sbw,flag = flag)
    srg = scanBam(bam.file, param=param, index = bam.index)
    
    all.index = rep(T,length(srg[[1]]$pos))
    
    if (paired.end) {
      chr.index = with(srg[[1]], rname == mrnm & !is.na(isize))
      chr.index[is.na(chr.index)] = FALSE
      # remove reads with large insertion
      size.index = with(srg[[1]], xor(strand == '-', isize>0) & abs(isize) < 1000 ) 
      size.index[is.na(size.index)] = FALSE
      
      all.index = all.index & size.index  & chr.index
    }
    
    srg = lapply(srg[[1]],function(x) x=x[all.index])
    
    bam.chr.list = foreach(chr = chr.names) %dopar% {
      chr.srg.plus = lapply(srg, function(x) x = x[with(srg, rname == chr & strand == '+')])
      chr.srg.minus = lapply(srg, function(x) x = x[with(srg, rname == chr & strand == '-')])
      if (!outCoverage) {
        # Get mid.points
        if (paired.end) {
          mid.pos.plus = with(chr.srg.plus, (pos + mpos) / 2 )
          mid.pos.minus = with(chr.srg.minus, (pos + mpos) / 2 )
        } else {
          mid.pos.plus = with(chr.srg.plus, pos + round(qwidth / 2))
          mid.pos.minus = with(chr.srg.minus, pos + round(qwidth / 2))
        }
        # count bins
        br = seq(0, chr.lengths[chr]+binning, by=binning)
        freq.mat = cbind(hist(mid.pos.plus, breaks=br, include.lowest = TRUE, plot = FALSE)$count,
                         hist(mid.pos.minus, breaks=br, include.lowest = TRUE, plot = FALSE)$count)
      } else {
        if (paired.end) {
          cov.plus = with(chr.srg.plus,
                          coverage(IRanges(start = pos, width=isize+qwidth),
                                   width = chr.lengths[chr] ))
          cov.minus = with(chr.srg.minus,
                           coverage(IRanges(start = pos, width=-isize+qwidth),
                                    width = chr.lengths[chr]))
        } else {
          cov.plus = with(chr.srg.plus, 
                          coverage(IRanges(start = pos, width=qwidth),
                                   width = chr.lengths[chr] ))
          cov.minus = with(chr.srg.minus, 
                           coverage(IRanges(start = pos, width=qwidth),
                                    width = chr.lengths[chr]))
        }
        binSum = function(chr.cov, binning) {
          bin = seq_len(chr.lengths[chr] %/% binning)*binning
          if (chr.lengths[chr] %% binning == 0L) bin = c(bin, unname(chr.lengths[chr]))
          diff(c(0L, cumsum(as.numeric(chr.cov))[bin]))
        }
        freq.mat = cbind(binSum(cov.plus, binning),
                         binSum(cov.minus, binning)) / binning
      }
      bam.name = strsplit(tail(unlist(strsplit(bam.file,'\\/')),1),'\\.bam') %>% 
        unlist() %>% head(1)
      colnames(freq.mat) = c(paste(bam.name, "+"), paste(bam.name, "-"))
      return(freq.mat)
    }
    names(bam.chr.list) = chr.names
    gc()
    return(bam.chr.list)
  }
  
  # Convert bam files to the chromosom coverage list objects
  bamCovRleList <<- function(bam.file, paired.end) {
    # Args:
    # bam.file: single file path
    # paired.end: reads single or paired end
    # Output strand split list of Rle coverage object
    bam.index=paste0(bam.file, '.bai')
    if (!file.exists(bam.index)) bam.index = Rsamtools::indexBam(files = bam.file)
    
    if (paired.end) {
      sbw = c('pos', 'qwidth','strand','rname', 'mrnm', 'mpos', 'isize')
      flag=scanBamFlag(isFirstMateRead = paired.end, isSecondaryAlignment=F)
    } else {
      sbw = c('pos', 'qwidth','strand','rname')
      flag = scanBamFlag(isSecondaryAlignment = FALSE,
                         isFirstMateRead = paired.end)
    }
    param = ScanBamParam(what = sbw,flag = flag)
    srg = scanBam(bam.file, param=param, index = bam.index)
    
    all.index = rep(TRUE,length(srg[[1]]$pos))
    if (paired.end) {
      chr.index = with(srg[[1]], rname == mrnm )
      chr.index[is.na(chr.index)] = FALSE
      all.index = all.index & chr.index
    }
    size.index = with(srg[[1]], xor(strand == '-', isize>0) & abs(isize)<1000 ) #remove reads with large insertion
    size.index[is.na(size.index)] = FALSE
    all.index = all.index & size.index
    
    srg = lapply(srg[[1]],function(x) x=x[all.index])
    
    bam.chr.list = foreach(chr = chr.names) %dopar% {
      chr.srg.plus = lapply(srg, function(x) x = x[with(srg, rname == chr & strand == '+')])
      chr.srg.minus = lapply(srg, function(x) x = x[with(srg, rname == chr & strand == '-')])
      
      if (paired.end) {
        cov.rle = list(with(chr.srg.plus,
                            coverage(IRanges(start = pos, width=isize+qwidth),
                                     width = chr.lengths[chr] )),
                       with(chr.srg.minus,
                            coverage(IRanges(start = pos, width=-isize+qwidth),
                                     width = chr.lengths[chr])) )
      } else {
        cov.rle = list(with(chr.srg.plus, 
                            coverage(IRanges(start = pos, width=qwidth),
                                     width = chr.lengths[chr] )),
                       with(chr.srg.minus, 
                            coverage(IRanges(start = pos, width=qwidth),
                                     width = chr.lengths[chr])) )
      }
      
      bam.name = unlist(strsplit(tail(unlist(strsplit(bam.file,'\\/')),1), '\\.bam'))[[1]] 
      names(cov.rle) = c(paste(bam.name, "+"), paste(bam.name, "-"))
      return(cov.rle)
    }
    names(bam.chr.list) = chr.names
    gc()
    return(bam.chr.list)
  }
  
  # Convert bam files to the short-reads list objects
  bam_bin_list <<- function(sample.list, FUN = bamBinCount, all.readPaired)
  {
    all.binned.counts.list = list()
    
    for (bams in seq_along(sample.list)) {
      start_time = Sys.time()
      bam.paths = sample.list[[bams]]
      paired.end = all.readPaired[[bams]]
      
      chr.lengths = seqlengths(seqinfo(BamFile(bam.paths[1])))
      bam_lvl_style = seqlevelsStyle(names(chr.lengths))
      chr.names = paste0('chr', c(1:100,'X','Y'))
      seqlevelsStyle(chr.names) = bam_lvl_style
      chr.names <<- intersect(chr.names, names(chr.lengths))
      chr.lengths <<- chr.lengths[chr.names]
      
      if (length(bam.paths) > 1) {
        binned.RNA.counts.list = list()
        for (more.bam in seq_along(bam.paths)) {
          bam2.counts.list = FUN(bam.paths[more.bam], paired.end[more.bam])
          for (chr in chr.names) {
            binned.RNA.counts.list[[chr]] = cbind(binned.RNA.counts.list[[chr]],
                                                  bam2.counts.list[[chr]])
            }
          }
        } else {binned.RNA.counts.list = FUN(bam.paths, paired.end)}
      seqlevelsStyle(names(binned.RNA.counts.list)) = "UCSC"
      
      readsNum = sum(unlist(lapply(binned.RNA.counts.list, sum)))
      if (readsNum == 0) stop(paste('Error: 0 read in sample', bams))
      
      bin.time = round(difftime(Sys.time(), start_time, units = 'mins'), 2)
      showNotification(paste("Processed ", sample.names[bams], ':'), duration = 50)
      showNotification(paste("Found", round(readsNum / 1e+6,2), "million reads"), duration = 50)
      showNotification(paste("Binning genome took", bin.time,"min"), duration = 50)
      
      sample.num = length(bam.paths)
      if (sample.num > 1) {
        binned.RNA.counts.list = lapply(binned.RNA.counts.list, function(x){
          cbind(rowSums(x[, seq(1, sample.num * 2-1, by = 2)]),
                rowSums(x[, seq(2, sample.num * 2, by = 2)]) )
        })
      }
      all.binned.counts.list = c(all.binned.counts.list, list(binned.RNA.counts.list) )
    } # the end of bam batches loop
    gc()
    names(all.binned.counts.list) = names(sample.list)
    return(all.binned.counts.list)
  }
  
  getAvgSignal_single <<- function(viterbi, obs, FUN = mean) 
  {
    viterbi = unlist(viterbi)
    mySignals = do.call("rbind", obs)
    state2pos = tapply(1:nrow(mySignals), INDEX = viterbi, function(x) mySignals[x,])
    state2val = lapply(state2pos, FUN, na.rm = TRUE)
    vit_means = do.call("rbind", state2val)
    vit_means
  }
  
  ### part2 -------------------
  # read uploaded file into ranges
  import.input.ranges <<- function(data, is.datapath)
  {
    if (is.datapath) {
      format = toupper(sapply(strsplit(data$datapath,'\\.'), tail, 1))
      import_input = import(data$datapath, format = format)
    } else {
      format = toupper(sapply(strsplit(data,'\\.'), tail, 1))
      import_input = import(data, format = format)
    }
    seqlevelsStyle(import_input) = "UCSC"
    
    import_input = import_input[seqnames(import_input)%in%paste0('chr',c(1:100,'X','Y'))]
    import_input = keepSeqlevels(import_input, unique(seqnames(import_input)), 
                                 pruning.mode="coarse")
    
    return(import_input)
  }
  
  # TU coverage Rle-----------------
  readCoverage <- function(bam.files, object.gr, targets = "full", flank.size) {
    #count reads coverage with given bamfiles
    #output a list of Rle strings of each range
    #Agrs:
    #     bam.file: file path
    #     object.gr: multiple intervals of gene or transcript
    #     targets: "full", "exon", "intron" (depends on Gene_input: total reference)
    #     flank.size: flank both ends to find precise boundary or for plotting
    object.gr = object.gr+flank.size
    object.start = start(object.gr)
    object.end = end(object.gr)
    object.strand = as.character(strand(object.gr)) == '+'
    
    if (targets != "full") {
      exon.gr = Gene_input[Gene_input$type == "exon" &
                             Gene_input$transcript_id %in%object.gr$transcript_id]
    }
    
    for (i in seq_along(bam.files)) {
      bam.index = paste0(bam.files[i],'.bai')
      if (!file.exists(bam.index)) bam.index = Rsamtools::indexBam(files = bam.files[i])
      # invisible(capture.output(is.paired.end = testPairedEndBam(bam.files[i], index = bam.index)))
      is.paired.end = testPairedEndBam(bam.files[i], index = bam.index)
      # scanBamWhat: the info that need to be extracted from a bam file.
      sbw = c('pos', 'qwidth', 'mapq', 'strand', 'rname',
              'mrnm', 'mpos', 'isize')
      sbp = ScanBamParam(what=sbw, which=object.gr,
                         flag=scanBamFlag(isUnmappedQuery = FALSE,
                                          isSecondaryAlignment = FALSE,
                                          isFirstMateRead = is.paired.end))#use first in pair reads for coverage
      
      # Scan bam file to retrieve short reads.
      sr.in.ranges = scanBam(bam.files[i], param=sbp,index = bam.index)
      
      scanBamRevOrder = function(org.gr, sbp) {
        org.grnames = with(org.gr, paste(seqnames, start, end, sep=':'))
        sbw.gr = as.data.frame(bamWhich(sbp))  # scan-bam-ed
        if ('space' %in% names(sbw.gr)) {
          sbw.grnames = with(sbw.gr, paste(space, start, end, sep=':'))
        } else if ('group_name' %in% names(sbw.gr)) {
          sbw.grnames = with(sbw.gr, paste(group_name, start, end, sep=':'))
        } else {
          stop("Cannot locate chromosome names in extracted short reads. Report
               this problem using issue tracking or discussion forum.\n")
        }
        match(org.grnames, sbw.grnames)
      }
      combind.list = function(a,b) {
        new = lapply(names(a), function(x) c(a[[x]],b[[x]]))
        names(new) = names(a)
        new
      }
      # Restore the original order.
      sr.in.ranges = sr.in.ranges[scanBamRevOrder(as.data.frame(object.gr), sbp)]
      srg.temp = list()
      for (n in seq_along(object.gr)) {
        # Filter short reads by mapping quality and insertion size.
        all.index = with(sr.in.ranges[[n]], 
                         mapq>3 & strand == as.character(strand(object.gr[n])))
        if (is.paired.end) all.index = all.index & with(sr.in.ranges[[n]], !is.na(isize))
        srg.temp = c(srg.temp, list(lapply(sr.in.ranges[[n]], 
                                           function(x) x=x[all.index])) )
      }
      names(srg.temp) = names(sr.in.ranges)
      if(!exists('srg.filtered')){
        srg.filtered = srg.temp
      }else{
        srg.filtered = lapply(names(srg.temp),function(y) 
          new=combind.list(srg.temp[[y]],srg.filtered[[y]]) )
      }
    }
    
    cov.list = list()
    for (n in seq_along(object.gr)) 
    {
      if (targets == "full") 
      {
        cov = with( srg.filtered[[n]],
                    coverage(IRanges(start = pos-object.start[n],
                                     width= ifelse(is.paired.end, abs(isize),0)+qwidth ), 
                             width = width(object.gr[n]) ))
      } else if (targets == "exon")
      {
        ov.index = with( srg.filtered[[n]], 
                         findOverlaps(IRanges(start = pos, width = qwidth), 
                                      ranges(exon.gr[exon.gr$gene_id == object.gr[n]$gene_id])) %>%
                           countQueryHits()!=0 )
        cov = with( lapply(srg.filtered[[n]], function(x) x[ov.index]), 
                    coverage(IRanges(start = pos-object.start[n],
                                     width= ifelse(is.paired.end, abs(isize),0)+qwidth ), 
                             width = width(object.gr[n]) ))
      } else if (targets == "intron")
      {
        ov.index = with( srg.filtered[[n]], 
                         findOverlaps(IRanges(start = pos, width = qwidth),
                                      ranges(exon.gr[exon.gr$gene_id == object.gr[n]$gene_id])) %>%
                           countQueryHits() == 0 )
        cov = with( lapply(srg.filtered[[n]],function(x) x[ov.index]),
                    coverage(IRanges(start = pos-object.start[n], 
                                     width= ifelse(is.paired.end, abs(isize),0)+qwidth ), 
                             width = width(object.gr[n]) ))
      }
      if (!object.strand[n]) cov=rev(cov)
      cov.list = c(cov.list, cov)
    }
    return(cov.list)
  }
  
  # Precise boundary----------------
  preciseBoundary <- function(object.gr, bam.files, k = 20, flank.size = 200, ...)
  {
    cov.list = readCoverage(bam.files, 
                            object.gr = object.gr,
                            targets = "full", 
                            flank.size = flank.size)
    
    object.gr = object.gr + flank.size
    
    foreach(n = seq_along(object.gr), .combine = c) %dopar% {
      cov = cov.list[[n]]
      temp.gr = object.gr[n]
      firstDer = splinefun(x = 1:length(cov),
                           y = sqrt(as.numeric(cov)),
                           method = 'natural')(1:length(cov), deriv = 1)
      kLDev = numeric(2 * flank.size - k)
      for (i in 1:(2 * flank.size - k) )
        kLDev[i] = mean(firstDer[i:(i + k - 1)])
      kRDev = numeric(2 * flank.size - k + 2)
      for (i in (length(firstDer) - 2 * flank.size):(length(firstDer) - k + 1) )
        kRDev[i] = mean(firstDer[i:(i + k - 1)])
      
      if (as.character(strand(temp.gr)) == '+') {
        start(temp.gr) = start(temp.gr) + which.max(kLDev) - 1
        end(temp.gr) = start(temp.gr) + which.min(kRDev) + k - 1
      } else {
        end(temp.gr) = end(temp.gr) - which.max(kLDev) + 1
        start(temp.gr) = end(temp.gr) - which.min(kRDev) - k + 1
      }
      return(temp.gr)
    }
  }
  
  # TU merging----------------------
  
  #generate gene ranges from TUs
  exon_overlap <<- function(tus, anno, exons) {
    # find exon overlapping TUs, and integrate them into one TU
    # original funtion from Michael
    # edit: merging step
    # Args:
    #   tus: TU GRanges object, raw file inpu
    #   anno: gene GRanges object
    #   exons: exon level GRanges object, ##with protein-coding genes and lincRNA (test non-coding RNA with >1 exons?)
    
    min.cutoffIdx = get_min_cutoff(tus, anno, ovCutoff=ovCutoff, annoCutoff=annoCutoff)
    no.cutoffIdx = get_min_cutoff(tus, anno, ovCutoff=0, annoCutoff=0)
    
    all_matching_tus = tus[no.cutoffIdx[[1]]]
    all_matching_tus$gene_id = anno$gene_id[no.cutoffIdx[[2]]]
    
    anno_matching_tus = tus[min.cutoffIdx[[1]]]
    anno_matching_tus$gene_id = anno$gene_id[min.cutoffIdx[[2]]]  ###give gene_id to any tu overlapping with genes, help to know related genes in downstream analysis
    
    cutoff_matching_tus = all_matching_tus[findOverlaps(query = anno_matching_tus, 
                                                        subject = all_matching_tus, 
                                                        type = "equal",ignore.strand = FALSE) %>%
                                             countSubjectHits() == 0]
    
    if (length(anno_matching_tus)>0) {
      sp = split(anno_matching_tus, anno_matching_tus$gene_id)
      message(paste(seqnames(tus)[1],': found minimum overlap', length(sp), 'of', length(anno),'genes'))
      exon_matching_tus.merged = foreach(g = names(sp), .combine="c") %dopar% {
        exons.gene_id = exons[exons$gene_id == g]
        
        tu_subset = sp[[g]]
        tu_subset$exonOV = FALSE
        
        ### extend upstream flank of first tu, in some cases the 1st TU starting sites are behind the end of 1st exons
        tu_subset_1flank = tu_subset
        tu_1.Idx = ifelse(as.character(strand(tu_subset[1])) == '+', 1, length(tu_subset))
        tu_subset_1flank[tu_1.Idx] = resize(tu_subset_1flank[tu_1.Idx],
                                            width = width(tu_subset_1flank[tu_1.Idx])+500,
                                            fix = 'end')
        # tu_subset_1flank=sort(tu_subset_1flank)
        tu_return_subset = c()
        
        mtch = findOverlaps(tu_subset_1flank, exons.gene_id,type = 'any',ignore.strand = FALSE)
        tu_subset$exonOV[queryHits(mtch)] = TRUE
        
        true.inds = which(tu_subset$exonOV == TRUE)
        if (length(true.inds) >= 1) {  ###if(length(true.inds) >= 1){ ### incase only one TU has exon overlap so all the others are unmerged, still marked as ncRNA
          min.true.ind = min(true.inds)
          max.true.ind = max(true.inds)
          
          false = c(1:length(tu_subset))[-c(min.true.ind:max.true.ind)]
          tu_return_subset_merged = tu_subset[min.true.ind]
          end(tu_return_subset_merged) = end(tu_subset[max.true.ind])
          tu_return_subset_merged$gene_id = g
          tu_return_subset_merged$transcript_id = NA
          tu_return_subset_merged$type = ifelse(as.character(exons.gene_id$gene_type)[1] == 'protein_coding','protein_coding','ncRNA')
          tu_return_subset_merged$gene_type = as.character(exons.gene_id$gene_type)[1]
          
          tu_return_subset_unmerged = tu_subset[false]
          if (length(false)>0) {
            tu_return_subset_unmerged$gene_id = g
            tu_return_subset_unmerged$transcript_id = NA
            tu_return_subset_unmerged$type = 'ncRNA'
            tu_return_subset_unmerged$gene_type = 'ncRNA'
          }
          tu_return_subset = c(tu_return_subset_merged, tu_return_subset_unmerged)
        } else {
          tu_return_subset = tu_subset[which.min(start(tu_subset))] ###some TUs don't overlap with exon and only in intron, need to be merged into one
          end(tu_return_subset) = end(tu_subset[which.max(end(tu_subset))])
          tu_return_subset$gene_id = g ###better to know which gene the tu overlaps
          tu_return_subset$transcript_id = NA
          tu_return_subset$type = 'ncRNA'
          tu_return_subset$gene_type = 'ncRNA'
        }
        return(tu_return_subset)
      }
      return.gr = c(exon_matching_tus.merged,cutoff_matching_tus)
    } else {
      return.gr = cutoff_matching_tus
    }
    if (length(return.gr)>0) {
      return(return.gr)
    } else {
      return(GRanges())
    }
  }
  
  exon_intersect <<- function(tus, anno, exons, transcripts) {
    # find exon overlapping TUs, and join exons boundary into coding TU, intersect upstream/downstream TU to usRNA/daRNA
    # original funtion from Michael
    # edit: merging step, add two classes of ncRNAs
    # Args:
    #   tus: TU GRanges object, raw file inpu
    #   anno: gene GRanges object
    #   exons: exon level GRanges object, ##with protein-coding genes and lincRNA (test non-coding RNA with >1 exons?)
    
    min.cutoffIdx = get_min_cutoff(tus, anno, ovCutoff=ovCutoff, annoCutoff=annoCutoff)
    no.cutoffIdx = get_min_cutoff(tus, anno, ovCutoff=0, annoCutoff=0)
    
    all_matching_tus = tus[no.cutoffIdx[[1]]]
    all_matching_tus$gene_id = anno$gene_id[no.cutoffIdx[[2]]]
    
    anno_matching_tus = tus[min.cutoffIdx[[1]]]
    anno_matching_tus$gene_id = anno$gene_id[min.cutoffIdx[[2]]]  ### give gene_id to any tu overlapping with genes, help to know related genes in downstream analysis
    
    # keep the tu interval couldn't reach minimum cutoff
    cutoff_matching_tus=all_matching_tus[findOverlaps(query = anno_matching_tus,
                                                      subject = all_matching_tus,
                                                      type = "equal",
                                                      ignore.strand = FALSE) %>%
                                           countSubjectHits() == 0]
    
    if (length(anno_matching_tus) > 0){
      sp = split(anno_matching_tus, anno_matching_tus$gene_id)
      message(paste(seqnames(tus)[1],': found minimum overlap', length(sp), 'of', length(anno),'genes'))
      
      exon_matching_tus.merged = foreach(g = names(sp), .combine="c",
                                         .packages = c('IRanges','GenomicRanges')) %dopar% {
        exons.gene_id = exons[exons$gene_id == g]
        transcripts.gene_id = transcripts[transcripts$gene_id == g]# all transcripts variant in the reference for gene g
        tu_subset = sp[[g]]
        tu_subset$exonOV = FALSE
        ### Previous strategy: 'extend upstream flank of first tu, in some cases the 1st TU starting sites are behind the end of 1st exons'
        ### Now new strategy: find tu overlap on TSS/1st exon, or closest 1st exon in annotated transcripts
        ### get %overlaps of tu cross transcripts, return the smallest distance between tu and transcript-100%-overlap while transcript keeps the shorter
        tu_subset.boundary = tu_subset[1];start(tu_subset.boundary) = min(start(tu_subset))
        end(tu_subset.boundary) = max(end(tu_subset))#maximum exon overlapped width
        tu_subset_1flank = tu_subset
        tu_1.Idx = ifelse(as.character(strand(tu_subset[1])) == '+', 1, length(tu_subset))
        tu_subset_1flank[tu_1.Idx] = tu_subset_1flank[tu_1.Idx]+1000
        
        mtch = findOverlaps(tu_subset_1flank, exons.gene_id, type = 'any', ignore.strand = FALSE)
        tu_subset$exonOV[queryHits(mtch)] = TRUE
        
        true.inds = which(tu_subset$exonOV == TRUE)
        if (length(true.inds) >= 1){  ### one exon overlapping can't indicate alternative TSS ### but could have execesive nc TU
          min.true.ind = min(true.inds)
          max.true.ind = max(true.inds)
          tu_return_subset_merged = tu_subset[min.true.ind]
          end(tu_return_subset_merged) = end(tu_subset[max.true.ind])
          ### get the maximum exon range, and cut excessive tu outside of exon boundary
          ### return the transcript with maximum TU overlap
          if (length(transcripts.gene_id) > 1){
            tx.idx = which.min(abs(start(transcripts.gene_id) - start(tu_subset.boundary)) + abs(end(transcripts.gene_id)-end(tu_subset.boundary)) )
            exon.boundary = transcripts.gene_id[ tx.idx ]# unmerged tu closest to exon-overlapped tu boundary
          } else {
            exon.boundary = transcripts.gene_id
          }
          tu_subset$transcript_id = NA
          tu_subset[which.min(start(tu_subset))]$transcript_id = as.character(exon.boundary$transcript_id)[1]
          mcols(exon.boundary) = mcols(tu_subset[which.min(start(tu_subset))])
          reduced.exon.boundary = IRanges::reduce(c(tu_return_subset_merged, exon.boundary))
          excessive.tu = GRanges()
          ### get nc excessive part
          # disjoin exon with reduced tu
          disjoin.exon.tu = disjoin(c(reduced.exon.boundary, exon.boundary))
          # get exon excessive regions
          excessive.tu = disjoin.exon.tu[ranges(disjoin.exon.tu) != ranges(exon.boundary)]
          excessive.tu = excessive.tu[width(excessive.tu) > binning] # discard flanking nc tu smaller than bin size, usually they are artifacts
          ### add metadata
          exon.boundary$gene_id = g
          exon.boundary$type = ifelse(as.character(exons.gene_id$gene_type)[1] == 'protein_coding', 'protein_coding', 'ncRNA')
          exon.boundary$gene_type = as.character(exons.gene_id$gene_type)[1]
          exon.boundary$exonOV = TRUE
          
          if (length(excessive.tu) > 0) {
            mcols(excessive.tu) = mcols(tu_subset[which.min(start(tu_subset))])### this will be removed later step, doesn't matter to randomly attribute
            excessive.tu$gene_id = g
            excessive.tu$transcript_id = NA
            excessive.tu$type = 'ncRNA'
            excessive.tu$gene_type = 'ncRNA'
            excessive.tu$exonOV = TRUE ### this will not be kept if tu overlaps multiple genes, see intersect.nc.tu
          }
          
          tu_return_subset = c(exon.boundary, excessive.tu)# tu_return_subset_unmerged)
        } else {
          tu_return_subset = tu_subset[which.min(start(tu_subset))] ### some TUs don't overlap with exon and only in intron, need to be merged into one
          end(tu_return_subset) = end(tu_subset[which.max(end(tu_subset))])
          
          tu_return_subset$gene_id = g ### which gene the tu overlaps
          tu_return_subset$transcript_id = NA
          tu_return_subset$type = 'ncRNA'
          tu_return_subset$gene_type = as.character(exons.gene_id$gene_type)[1]
        }
        return(tu_return_subset)
      }
      
      nc.gene.disjoin.idx = which(exon_matching_tus.merged$gene_type == 'ncRNA' & 
                                    exon_matching_tus.merged$exonOV == TRUE)
      ### reduce non-coding tus redundantly overlap multiple genes
      if (length(nc.gene.disjoin.idx) > 0)
      {
        nc.gene.disjoin.tu = exon_matching_tus.merged[nc.gene.disjoin.idx]
        intersect.nc.tu = disjoin(nc.gene.disjoin.tu)
        internal.nc.idx = findOverlaps(query = intersect.nc.tu, 
                                       subject = nc.gene.disjoin.tu) %>% countQueryHits() == 2
        ### add matadata
        mid.nc.tu = intersect.nc.tu[internal.nc.idx]
        mid.nc.tu = mid.nc.tu[width(mid.nc.tu)>binning] #### discard middle nc tu if they are too small
        if (length(mid.nc.tu) > 0)
        {
          mid.matches = findOverlaps(mid.nc.tu, tus,ignore.strand = FALSE, type = "within")
          mcols(mid.nc.tu) = suppressWarnings(mcols(tus[subjectHits(mid.matches)]))
          mid.nc.tu$exonOV = F
          mid.nc.tu$gene_type = "ncRNA"
          mid.nc.tu$gene_id = NA
          mid.nc.tu$transcript_id = NA
          mid.nc.tu$type = 'ncRNA'
          mcols(mid.nc.tu) = mcols(mid.nc.tu)[, colnames(mcols(nc.gene.disjoin.tu))]
        } else
        {
          mid.nc.tu = GRanges()
        }
        
        intersect.idx = findOverlaps(query = intersect.nc.tu, 
                                     subject = nc.gene.disjoin.tu) %>% countSubjectHits() < 2
        all.nc.gene.tu = c(nc.gene.disjoin.tu[intersect.idx], mid.nc.tu)
        exon_matching_tus.merged = c(exon_matching_tus.merged[-nc.gene.disjoin.idx], all.nc.gene.tu)
      } else
      {
        exon_matching_tus.merged = GRanges()}
      return.gr = c(exon_matching_tus.merged,cutoff_matching_tus)
    } else
    {
      return.gr = cutoff_matching_tus
    }
    if (length(return.gr) > 0)
    {
      return(return.gr)
    } else
    {
      return(GRanges())
    }
  }
  
  get_min_cutoff <<- function(tus, anno, ovCutoff, annoCutoff)
  {
    mtch = findOverlaps(tus, anno)
    ourTranscripts = queryHits(mtch) # index based on TUs
    annoTranscripts = subjectHits(mtch) # index based on GencodeAnno
    ovRanges = overlapsRanges(hits = mtch, query = IRanges::ranges(tus),
                              subject = IRanges::ranges(anno))
    
    myOverlapLens = as.numeric(width(ovRanges)) / width(tus[ourTranscripts])
    myAnnoLens = as.numeric(width(ovRanges)) / width(anno[annoTranscripts])
    if (xor(ovCutoff == 0,annoCutoff == 0))
    {
      passCutoff = which(myOverlapLens >= ovCutoff & myAnnoLens >= annoCutoff)### to remove embeded genes in other larger genes
    } else
    {
      passCutoff = which(myOverlapLens >= ovCutoff | myAnnoLens >= annoCutoff)
    }
    idxOurTranscripts = ourTranscripts[passCutoff] # index based on TUs
    idxAnno = annoTranscripts[passCutoff] # index based on GencodeAnno
    idx = list(idxOurTranscripts, idxAnno)
    return(idx)
  }
  
  set_tx_location <<- function(ncRNAs.gr, main.genes.gr)
  {
    within.sense = findOverlaps(ncRNAs.gr, main.genes.gr,  
                                type = 'within',ignore.strand = FALSE) ###add type = 'within'
    ncRNAs_antisense = ncRNAs.gr
    levels(strand(ncRNAs_antisense)) = c("-","+","*")
    within.antisense = findOverlaps(ncRNAs_antisense, main.genes.gr, 
                                    type = 'any', ignore.strand = FALSE) ###add type = 'any'
    
    if (length(queryHits(within.antisense)) > 0) 
      ncRNAs.gr$location[unique(queryHits(within.antisense))] = "antisense"
    if (length(queryHits(within.sense)) > 0) 
      ncRNAs.gr$location[unique(queryHits(within.sense))] = "gene_sense"
    if (TRUE)
    {
      ncRNAs.TSS = promoters(ncRNAs.gr, upstream=1, downstream=0)
      gene.TSS = promoters(main.genes.gr, upstream = 1, downstream = 0)
      levels(strand(ncRNAs.TSS)) = c("-","+","*")
      upstream.antisense = findOverlaps(ncRNAs.TSS, gene.TSS, 
                                        ignore.strand = FALSE, maxgap = TSS_len)
      
      plus = which(strand(gene.TSS[subjectHits(upstream.antisense)]) == "+")
      minus = which(strand(gene.TSS[subjectHits(upstream.antisense)]) == "-")
      
      plus.prompt = which(start(ncRNAs.TSS[queryHits(upstream.antisense)][plus]) <
                            start(gene.TSS[subjectHits(upstream.antisense)][plus]))
      
      minus.prompt = which(start(ncRNAs.TSS[queryHits(upstream.antisense)][minus]) >
                             start(gene.TSS[subjectHits(upstream.antisense)][minus]))
      
      prompt = c(plus[plus.prompt],minus[minus.prompt])
      
      plus.conv = which(start(ncRNAs.TSS[queryHits(upstream.antisense)][plus]) > 
                          start(gene.TSS[subjectHits(upstream.antisense)][plus]))
      
      minus.conv = which(start(ncRNAs.TSS[queryHits(upstream.antisense)][minus]) < 
                           start(gene.TSS[subjectHits(upstream.antisense)][minus]))
      
      conv = c(plus[plus.conv],minus[minus.conv])
      
      if (length(conv) > 0) 
        ncRNAs.gr$location[queryHits(upstream.antisense)[conv]] = "conRNA"
      if (length(prompt) > 0) 
        ncRNAs.gr$location[queryHits(upstream.antisense)[prompt]] = "uaRNA"
    }
    if (mergeMethod == 'exon')
    {
      downstream.regions = flank(main.genes.gr, width=TTS_len, start = FALSE) + 2 ###+2 to disjoin nc parts
      ncRNAs.TSS = promoters(ncRNAs.gr, upstream=1, downstream=0)
      downstream = findOverlaps(ncRNAs.TSS, downstream.regions)
      
      if (length(queryHits(downstream))>0)
      {
        ncRNAs.gr$location[unique(queryHits(downstream))] = "dsRNA"
      }
      # add usRNA and daRNA
      upstream.regions = promoters(main.genes.gr, upstream=Promoter_len, downstream=0) + 2 ###+2 to disjoin nc parts
      ncRNAs.TTS = flank(ncRNAs.gr,width=1,start=F)
      upstream.sense = findOverlaps(ncRNAs.TTS, upstream.regions)
      if (length(queryHits(upstream.sense))>0) 
        ncRNAs.gr$location[unique(queryHits(upstream.sense))] = "usRNA"
    }
    
    ncRNAs.TTS=flank(ncRNAs.gr,width=1,start=F)
    downstream.regions = flank(main.genes.gr, width=TTS_len, start = FALSE)
    levels(strand(ncRNAs.TTS)) = c("-","+","*")
    downstream.antisense = findOverlaps(ncRNAs.TTS, downstream.regions)
    if (length(queryHits(downstream.antisense))>0)
    {
      ncRNAs.gr$location[unique(queryHits(downstream.antisense))] = "daRNA"
    }
    
    if (T)
    {
      # add usRNA and dsRNA
      upstream.regions = promoters(main.genes.gr, upstream=Promoter_len, downstream=0)
      ncRNAs.TTS = flank(ncRNAs.gr,width=1,start=F)
      upstream.sense = findOverlaps(ncRNAs.TTS, upstream.regions)
      if (length(queryHits(upstream.sense))>0) 
        ncRNAs.gr$location[unique(queryHits(upstream.sense))] = "usRNA"
      
      ncRNAs.TSS = promoters(ncRNAs.gr, upstream=1, downstream=0)
      downstream = findOverlaps(ncRNAs.TSS, downstream.regions)
      if (length(queryHits(downstream))>0) 
        ncRNAs.gr$location[unique(queryHits(downstream))] = "dsRNA"
    }
    return(ncRNAs.gr)
  }
  
  # GRanges interval calling function with GenoStan
  STAN_anno <<- function(binned.RNA.counts.list)
  {
    sample.num = ncol(binned.RNA.counts.list[[1]]) / 2
    if (sample.num > 1)
    {
      binned.plus.merged = lapply(binned.RNA.counts.list, 
                                  function(x) x[, seq(1, sample.num * 2 - 1, by = 2)])
      binned.minus.merged = lapply(binned.RNA.counts.list,
                                   function(x) x[, seq(2, sample.num * 2, by = 2)])
    } else {
      binned.plus.merged = lapply(binned.RNA.counts.list,
                                  function(x) matrix(x[, 1], ncol = 1) )
      binned.minus.merged = lapply(binned.RNA.counts.list,
                                   function(x) matrix(x[, 2], ncol = 1) )
    }
    names(binned.plus.merged) = paste(names(binned.plus.merged), ".L", sep="")
    names(binned.minus.merged) = paste(names(binned.minus.merged), ".L", sep="")
    
    chrs = unique(do.call("rbind", strsplit(names(binned.plus.merged), "\\."))[, 1])
    pilot.mm10.regions = GRanges(seqnames = chrs,
                                 ranges = IRanges(start = 1,
                                                  end = unique(sapply(binned.plus.merged,length)*binning)),
                                 strand = "*",
                                 name = chrs)
    
    celltypes = lapply(names(binned.plus.merged), function(n) grep(n, names(binned.plus.merged)))
    names(celltypes) = names(binned.plus.merged)
    
    
    .getSizeFactors = function(obs, celltypes) 
    {
      myAvg = apply(t(sapply(celltypes,
                             function(x) apply(do.call("rbind", obs[x]), 2, sum)
                             )),
                    2, mean)
      sizeFactors = matrix(NA, nrow = length(obs), ncol = ncol(obs[[1]]))
      rownames(sizeFactors) = sapply(names(celltypes),
                                     function(x) rep(x, length(celltypes[[x]]))
                                     )
      colnames(sizeFactors) = colnames(obs[[1]])
      for (cell in 1:length(celltypes)) {
        for (currDim in 1:ncol(sizeFactors)) {
          sizeFactors[celltypes[[cell]], currDim] =
            myAvg[currDim] / sum(do.call("rbind", obs[celltypes[[cell]]])[, currDim])
        }
      }
      sizeFactors
    }
    
    sizeFactors.plus = .getSizeFactors(binned.plus.merged, celltypes)
    sizeFactors.minus = .getSizeFactors(binned.minus.merged, celltypes)
    
    # add background counts
    binned.plus = lapply(binned.plus.merged, function(x) x + 1)
    binned.minus = lapply(binned.minus.merged, function(x) x + 1)
    
    anno = GRanges()
    anno = foreach(i = seq_along(fit_method), .combine = c, .packages = c("STAN")) %dopar% {
      # watson strand
      segmentation = segment_model_fitting(binned.plus, fit_method[i],
                                           nStates, ncores, sizeFactors.plus,
                                           pilot.mm10.regions)
      expState = segmentation[2]$name
      segmentation = segmentation[which(segmentation$name == expState)]
      strand(segmentation) = "+"
      anno = c(anno, segmentation)
      # crick strand
      segmentation = segment_model_fitting(binned.minus, fit_method[i], 
                                           nStates, ncores, sizeFactors.minus,
                                           pilot.mm10.regions)
      expState = segmentation[2]$name
      segmentation = segmentation[which(segmentation$name == expState)]
      strand(segmentation) = "-"
      anno = c(anno, segmentation)
      return(anno)
    }
    
    TU_input = IRanges::reduce(anno)
    
    # fill metadata--------
    TU_input$id = 1:length(TU_input)
    TU_input$gene_id = as.character(NA)
    TU_input$exonOV = FALSE
    TU_input$type = "ncRNA"
    TU_input$gene_type = as.character(NA)
    
    # filnal output
    return(TU_input)
  }
  
  segment_model_fitting <<- function(DTASeq_data_labeled,
                                     fit_method,
                                     nStates,
                                     ncores,
                                     sizeFactors,
                                     pilot.mm10.regions)
  {
    ## Model initialization
    if(fit_method == 'Gauss')
    {
      DTASeq_data_labeled = lapply(DTASeq_data_labeled, function(x)
        apply(log(x+sqrt(x^2+1)), 2, runningMean, 2))
      hmm.data = initHMM(DTASeq_data_labeled, nStates, "IndependentGaussian", sharedCov = TRUE)
    } else if(fit_method == "NB")
    {
      hmm.data = initHMM(DTASeq_data_labeled, nStates=nStates, "NegativeBinomial",sizeFactors=sizeFactors)
    } else if(fit_method == "Poilog")
    {
      hmm.data = initHMM(DTASeq_data_labeled, nStates=nStates, "PoissonLogNormal", sizeFactors=sizeFactors)
    } else if(fit_method == "ZINB")
    {
      hmm.data = initHMM(DTASeq_data_labeled, nStates=nStates, "ZINegativeBinomial")
    }
    ## Model fitting
    hmm_fitted = fitHMM(DTASeq_data_labeled, hmm.data, maxIters=10, nCores= ncores)
    ## Calculate state path
    viterbi = getViterbi(hmm_fitted, DTASeq_data_labeled)
    ## Convert state path to GRanges object
    viterbi_granges = viterbi2GRanges(viterbi, pilot.mm10.regions, binning)
    return(viterbi_granges)
  }
  
  # call TU from bin counts --------------------------
  bin_transcribed_TU <<- function(binned.RNA.counts.list, Gene_input, bam.files, ...)
  {
    # this function wraps STAN HMM state calling and TU merging steps
    # args:
    # binned.RNA.counts.list: count list of each sample
    # Gene_input: gene reference
    
    # STAN calling --------------
    TU_input <<- sort(STAN_anno(binned.RNA.counts.list))
    # reference prep -------------
    if (!exists('gene.gr'))
    {
      gene.gr = gene_prep(Gene_input, TU_input)
    }
    current.gene.gr = gene.gr[seqnames(gene.gr) %in% chr.names]
    # merging TU -----------------------------
    TU.gr = TU_prep(Gene_input = Gene_input,
                    gene.gr = current.gene.gr,
                    TU_input = TU_input,
                    bam.files)
    return(TU.gr)
  }
  
  # extract gene granges
  gene_prep <<- function(Gene_input, TU_input)
  {
    # this function will exclude some psudeogenes for avoiding gene embeding issue
    showNotification('Preparing minimal gene transcripts...')
    seqlevelsStyle(Gene_input) = "UCSC"
    if (!'gene' %in% unique(Gene_input$type))
    {
      stop("Annotation reference does not contain any 'gene' feature.
             Please check 'type' column in the uploaded reference.")
    } else {
      gene.gr = Gene_input[Gene_input$type == 'gene' & Gene_input$gene_type %in% mergeFeatures]
      if ("gene_source" %in% colnames(mcols(gene.gr)))
      {
        if (any('havana' %in% tolower(gene.gr$gene_source))) 
          gene.gr = gene.gr[gene.gr$gene_source!='havana'] # different sources introduce duplicates
      }
      ### remove embedded genes
      if(any(grepl('^Gm',gene.gr$gene_name)))
      {
        Gm.gene.gr = gene.gr[grepl('^Gm', gene.gr$gene_name)]
        gene.gr = c(gene.gr[!grepl('^Gm', gene.gr$gene_name)],
                    Gm.gene.gr[countQueryHits(findOverlaps(Gm.gene.gr, gene.gr)) <= 1])
        gene.gr = sort(gene.gr)
      }
      gene.gr = gene.gr[countSubjectHits(findOverlaps(gene.gr, gene.gr, type = 'within')) <= 1]
      return(gene.gr)
    }
  }
  
  # merge TU ranges to genes-------------------------
  merge_annotation <<- function(Gene_input, gene.gr, TU_input, mergeMethod)
  {
    exon.gr = Gene_input[Gene_input$type == 'exon' & Gene_input$gene_id %in% gene.gr$gene_id]
    transcript.gr = Gene_input[Gene_input$type == 'transcript' & Gene_input$gene_id %in% gene.gr$gene_id]
    # split by seqnames for seperate process
    sp.gene.gr = split(gene.gr, seqnames(gene.gr))
    sp.exon.gr = split(exon.gr, seqnames(exon.gr))
    sp.transcript.gr = split(transcript.gr, seqnames(transcript.gr))
    sp.TU_input = split(TU_input, seqnames(TU_input))
    exon.TU.gr = GRanges()
    
    for (i in seqlevels(TU_input))
    {
      if (mergeMethod == 'tu')
      {
        temp.TU = exon_overlap(tus = sp.TU_input[[i]],
                              anno = sp.gene.gr[[i]],
                              exons = sp.exon.gr[[i]])
      } else if (mergeMethod == 'exon')
      {
        temp.TU = exon_intersect(tus = sp.TU_input[[i]], 
                                anno = sp.gene.gr[[i]], 
                                exons = sp.exon.gr[[i]], 
                                transcripts = sp.transcript.gr[[i]])
      } else
      {
        stop('Error: please choose a level to merge gene-TU. Stop.\n')
      }
      exon.TU.gr = c(exon.TU.gr, TU_dedup(temp.TU, sp.gene.gr[[i]]))
    }
    
    tus = sp.TU_input[[i]]
    anno = sp.gene.gr[[i]]
    exons = sp.exon.gr[[i]]
    transcripts = sp.transcript.gr[[i]]
    
    if (T)
    {
      TU_input$type = 'ncRNA'
      nc.TU.gr = TU_input[-get_min_cutoff(TU_input, gene.gr,
                                          ovCutoff = ovCutoff, 
                                          annoCutoff = annoCutoff)[[1]]]
      
      nc.TU.gr = nc.TU.gr[countQueryHits(findOverlaps(nc.TU.gr, 
                                                      IRanges::reduce(exon.TU.gr),
                                                      type='within', 
                                                      ignore.strand = F)) == 0] ###remove embedded genes
      TU.gr = sort(c(exon.TU.gr, nc.TU.gr)) 
      TU.gr$id = seq_along(TU.gr)
    }
    return(TU.gr)
  }
  
  # deduplication ---------------------------------
  TU_dedup <<- function(current.TU.gr, current.ref)
  {
    ### some duplicated TU overlap with multiple genes (often happens when short genes embeded in long genes)
    ### only keep the large gene id for unique TUs
    current.TU.plus = current.TU.gr[strand(current.TU.gr) == '+']
    current.TU.minus = current.TU.gr[strand(current.TU.gr) == '-']
    dup.ranges.plus = duplicated(ranges(current.TU.plus))
    dup.ranges.minus = duplicated(ranges(current.TU.minus))
    dup.ranges = IRanges::reduce(c(current.TU.plus[dup.ranges.plus], 
                          current.TU.minus[dup.ranges.minus]) ) # in case more than 2 duplicates
    
    if (length(dup.ranges) > 0)
    {
      all.dup.ranges = c(current.TU.plus[ranges(current.TU.plus) %in% 
                                           ranges(current.TU.plus[dup.ranges.plus])],
                         current.TU.minus[ranges(current.TU.minus) %in% 
                                            ranges(current.TU.minus[dup.ranges.minus])])
      unique.TU.gr = current.TU.gr[!current.TU.gr$gene_id %in% all.dup.ranges$gene_id]
      
      max.gene.TU.id = foreach(x = seq_along(dup.ranges), .combine = c) %dopar% {
        temp.gene = current.ref[countSubjectHits(findOverlaps(query = dup.ranges[x], subject = current.ref)) > 0]
        intct_width = width(intersect(temp.gene, dup.ranges[x]))
        if (all(intct_width == intct_width[1])) {
          temp.gene[which.max(intct_width / width(temp.gene))]$gene_id
        } else {
          temp.gene[which.max(intct_width)]$gene_id
        }
      }
      dedup.TU.gr = all.dup.ranges[all.dup.ranges$gene_id %in% max.gene.TU.id]
      dedup.TU.gr = dedup.TU.gr[!duplicated(ranges(dedup.TU.gr))]
      return(sort(c(unique.TU.gr,dedup.TU.gr)))
    } else
    {
      return(current.TU.gr)
    }
  }
  
  # annotate ncRNA -------------------------
  TU_prep <<- function (Gene_input, gene.gr, TU_input, bam.files)
  { # this function will merge and filter TU inputs by overlapping to gene intervals
    # target chromosomes, e.g. some chromosome might have 0 interval
    
    chr.names = unique(as.character(seqnames(TU_input)))
    # some library prep methods could have mate in pairs as the forward reads, check which strand has max gene coverage
    all_matched =  
      sum( width( TU_input[get_min_cutoff(TU_input, gene.gr, ovCutoff=0L, annoCutoff=0L)[[1]]] ) )
    TU_input_rev = TU_input; levels(strand(TU_input_rev)) = c("-", "+", "*")
    all_matched_rev = 
      sum( width( TU_input_rev[get_min_cutoff(TU_input_rev, gene.gr, ovCutoff=0L, annoCutoff=0L)[[1]]] ) )
    
    if (all_matched_rev > all_matched)
    {
      TU_input = TU_input_rev
      strandFlipped <<- TRUE
    }
    
    if ( is.null(txInput) )
    {
      if (mergeMethod != 'skip')
      {
        TU.gr <<- merge_annotation(Gene_input, gene.gr, TU_input, mergeMethod)
      } else { # skip TU merging
        all.cutoffIdx = get_min_cutoff(TU_input, gene.gr, ovCutoff=ovCutoff, annoCutoff=annoCutoff)
        all_matching_tus = TU_input[all.cutoffIdx[[1]]]
        all_matching_tus$gene_id = gene.gr$gene_id[all.cutoffIdx[[2]]]
        all_matching_tus$transcript_id = gene.gr$transcript_id[all.cutoffIdx[[2]]]
        all_matching_tus$type = 
          ifelse(gene.gr$gene_type[all.cutoffIdx[[2]]] == 'protein_coding',
                 'protein_coding','ncRNA')
        all_matching_tus$gene_type = gene.gr$gene_type[all.cutoffIdx[[2]]]
        # remove duplicated gene overlapped ranges
        all_matching_tus = all_matching_tus[!duplicated(all.cutoffIdx[[1]])]
        TU.gr = sort(c(all_matching_tus, TU_input[-all.cutoffIdx[[1]]])) 
      }
      
      TU.gr$location = TU.gr$type
      if (is.null(TU.gr$id)) TU.gr$id = seq_along(TU.gr)
      if (is.null(TU.gr$expr)) TU.gr$expr = 0
      mcols(TU.gr) = mcols(TU.gr)[, c('id','expr','gene_id','transcript_id','exonOV',
                                      'type','gene_type','location')]
      
      our.genes = TU.gr[TU.gr$type == 'protein_coding']
      ref.gene = gene.gr[gene.gr$gene_type == 'protein_coding']
      ref.gene = ref.gene[ countQueryHits(findOverlaps(ref.gene,TU.gr,ignore.strand=F) ) > 0]
      mcols(our.genes) = mcols(our.genes)[, c('gene_id','type')]
      mcols(ref.gene) = mcols(ref.gene)[, c('gene_id','type')]
      main.genes.gr = sort(c(our.genes, ref.gene[!ref.gene$gene_id %in% our.genes$gene_id]))
      
    } else { # uploaded TU ranges
      
      main.genes.gr = sort(gene.gr[gene.gr$gene_type == 'protein_coding'])
      TU.gr = TU_input
      matches = findOverlaps(query = TU_input, subject = main.genes.gr, ignore.strand=F)
      TU.gr$type = 'ncRNA'
      TU.gr$type[countQueryHits(matches)>0] = 'protein_coding'
      TU.gr$id = seq_along(TU.gr)
      TU.gr$expr = 0
      TU.gr$gene_id = TU.gr$transcript_id = NA
      TU.gr$gene_id[queryHits(matches)] = main.genes.gr$gene_id[subjectHits(matches)]
      TU.gr$transcript_id[queryHits(matches)] = main.genes.gr$transcript_id[subjectHits(matches)]
      TU.gr$exonOV = FALSE
      TU.gr$exonOV[countQueryHits(matches)>0] = TRUE
      TU.gr$gene_type = TU.gr$location = TU.gr$type
    } # the end of TU preprocessing
    
    # mask gene list ----------
    if ( !is.null(mask_list) )
    {
      inList = mask_list
      mask_list = read.table(inList$datapath, header = FALSE)
      rm.idx = mcols(main.genes.gr)$gene_id %in% unlist(mask_list)
      main.genes.gr = main.genes.gr[!rm.idx,]
    }
    # merge nearby nc TU ----------------
    if (!is.null(unmappable))
    {
      # reduce ncRNA seperated by unmappable regions
      ncRNAs.gr = TU.gr[which(TU.gr$type == "ncRNA")]
      # get gaps
      nc_reduce.gap = gaps(ncRNAs.gr)
      nc_reduce.gap = nc_reduce.gap[width(nc_reduce.gap) < 1500]
      unmappable.gr <<- import.input.ranges(unmappable, TRUE)
      unmappable.gr = unmappable.gr[width(unmappable.gr) < 1000000] + binning # extend ranges by the average reads insert size
      nc_reduce.gap = nc_reduce.gap[countSubjectHits(findOverlaps(query = unmappable.gr,
                                                                  subject = nc_reduce.gap, 
                                                                  type = 'any',
                                                                  ignore.strand = TRUE) ) > 0]
      
      ncRNAs.reduce.gap.Idx = findOverlaps(query = nc_reduce.gap + 1, 
                                           subject = ncRNAs.gr,
                                           type = 'any',
                                           ignore.strand = FALSE)
      ncRNAs.nonreduce.gr = ncRNAs.gr[countSubjectHits(ncRNAs.reduce.gap.Idx) == 0]
      sp.ncRNA.reduce.ls = split(ncRNAs.gr[subjectHits(ncRNAs.reduce.gap.Idx)], 
                                 queryHits(ncRNAs.reduce.gap.Idx))
      ncRNA.reduce.gr = GRanges()
      ncRNA.reduce.gr = foreach(i=1:length(sp.ncRNA.reduce.ls), .combine = c) %dopar% {
        temp.reduce.gr = sp.ncRNA.reduce.ls[[i]]
        return.reduce.gr = temp.reduce.gr[1]
        start(return.reduce.gr) = min(start(temp.reduce.gr))
        end(return.reduce.gr) = max(end(temp.reduce.gr))
        if (!is.null(return.reduce.gr$score)) 
          return.reduce.gr$score = mean(temp.reduce.gr$score)
        return(return.reduce.gr)
      }
      ncRNAs.gr = IRanges::reduce(c(ncRNA.reduce.gr, ncRNAs.nonreduce.gr))
    } else
    {
      ncRNAs.gr = TU.gr[TU.gr$type == "ncRNA"]
    }
    
    # set ncRNA TU from reduce step -----------
    if (TRUE)
    {
      if (precBound)
      {
        cat("---------------------------------- \n")
        cat("Refining non-coding TU boundary... \n")
        cat("---------------------------------- \n")
        ncRNAs.gr = preciseBoundary(ncRNAs.gr, bam.files, k=20, flank.size = binning)
      }
      # 'id','expr','gene_id','exonOV','type','gene_type','location'
      ncRNAs.gr$id = NA
      ncRNAs.gr$expr = 0
      # give old id to reduced ncRNA
      nc.ov.idx = findOverlaps(query = ncRNAs.gr, 
                               subject = TU.gr[TU.gr$type == "ncRNA"],
                               ignore.strand = F)
      ncRNAs.gr$id = base::replace(ncRNAs.gr$id, queryHits(nc.ov.idx),
                                   TU.gr[TU.gr$type == "ncRNA"]$id[subjectHits(nc.ov.idx)])
      # ncRNAs.gr$expr=base::replace(ncRNAs.gr$expr,queryHits(nc.ov.idx),TU.gr[TU.gr$type == "ncRNA"]$expr[subjectHits(nc.ov.idx)])
      ncRNAs.gr$gene_id = as.character(NA)
      ncRNAs.gr$transcript_id = as.character(NA)
      ncRNAs.gr$exonOV = FALSE
      ncRNAs.gr$type = "ncRNA"
      ncRNAs.gr$gene_type = "ncRNA"
      nc.ov.idx = findOverlaps(query = ncRNAs.gr, 
                               subject = Gene_input[Gene_input$gene_type!='protein_coding'],
                               ignore.strand = FALSE, minoverlap = binning)
      ncRNAs.gr$gene_type = base::replace(ncRNAs.gr$gene_type, 
                                          queryHits(nc.ov.idx), 
                                          Gene_input[Gene_input$gene_type != 'protein_coding']$gene_type[subjectHits(nc.ov.idx)])
      ncRNAs.gr$location = "intergenic"
      mcols(ncRNAs.gr) = mcols(ncRNAs.gr)[, c('id','expr','gene_id','transcript_id',
                                              'exonOV','type','gene_type','location')]
      coding.gr = TU.gr[TU.gr$type != "ncRNA"]
      coding.gr$location = "protein_coding"
      mcols(coding.gr) = mcols(coding.gr)[, c('id','expr','gene_id','transcript_id',
                                              'exonOV','type','gene_type','location')]
      # write location
      TU.gr = sort(c(coding.gr, set_tx_location(ncRNAs.gr, main.genes.gr) ))
    }
    return(TU.gr)
  }
  
  # count list handling---------------------------
  count_list_combine <<- function(all.binned.counts.list)
  {
    # etract bin counts on the chromosome existing in all the samples
    binned.RNA.counts.list = list()
    temp.chrs = table(unlist(lapply(all.binned.counts.list, function(x) names(x)) ) )
    temp.chrs = intersect(c(paste0('chr',1:100),'chrX','chrY'), 
                          names(temp.chrs)[temp.chrs == length(all.binned.counts.list)])
    for (chr in temp.chrs)
    {
      binned.RNA.counts.list[[chr]] = Reduce(cbind, lapply(all.binned.counts.list,
                                                           function(x) x[[chr]] ))
    }
    return(binned.RNA.counts.list)
  }
  
  count.list.split <<- function(binned.RNA.counts.list)
  {
    new.list = list()
    sample.num = ncol(binned.RNA.counts.list[[1]]) / 2
    for (i in seq(1,sample.num*2,by = 2)) 
      new.list = c(new.list, 
                   list(lapply(binned.RNA.counts.list,
                               function(x) x[,c(i,i+1)])) 
      )
    return(new.list)
  }
  
  #gene intersection------------------------------
  get_common_gene <<- function(TU.list)
  {
    # TU.genes=GRanges()
    for (i in seq_along(TU.list) )
    {
      if (i == 1)
      {
        TU.genes=sort(TU.list[[i]])
      } else
      {
        TU.genes=intersect(TU.genes, sort(TU.list[[i]]) )
      }
    }
    return(TU.genes)
  }
  
  #expr cutoff ----------------------------------
  TU_count_cal <<- function(binned.RNA.counts.list, TU.genes.gr, chr.lengths, RPK = FALSE)
  {
    # count TU gene reads from binned list
    # Args:
    # binned.RNA.counts.list
    # TU.genes.gr: sum counts in bin inside subject ranges
    # chr.lengths: chromosome length for rebiulding bins
    # RPK: logical, if returns reads per Kb (RPK)
    
    temp.sample.num = ncol(binned.RNA.counts.list[[1]]) / 2
    if(temp.sample.num == 1)
    {
      binned.RNA.counts.list=lapply(binned.RNA.counts.list, function(x) cbind(x,x))
      temp.sample.num = ncol(binned.RNA.counts.list[[1]]) / 2
    }
    
    chrs = unique(as.character(seqnames(TU.genes.gr)))
    foreach(chr=chrs, .combine = rbind) %dopar% {
      binned.temp.counts = binned.RNA.counts.list[[chr]]
      gene.temp.plus = TU.genes.gr[seqnames(TU.genes.gr) == chr & strand(TU.genes.gr) == '+']
      gene.temp.minus = TU.genes.gr[seqnames(TU.genes.gr) == chr & strand(TU.genes.gr) == '-']
      
      bin.ranges = IRanges(start = seq(1, chr.lengths[chr], by = binning), width = binning)
      plus.matches = findOverlaps(query = bin.ranges,
                                  subject = ranges(gene.temp.plus),type='any',minoverlap = binning/2)
      minus.matches = findOverlaps(query = bin.ranges,
                                   subject = ranges(gene.temp.minus),type='any',minoverlap = binning/2)
      
      chr.counts = NULL
      col_idx = seq(1,temp.sample.num*2,2)
      for (i in seq_along(gene.temp.plus))
      {
        tu.matches = queryHits(plus.matches)[subjectHits(plus.matches) == i]
        if (length(tu.matches) == 1)
        {
          chr.counts = rbind(chr.counts,
                             binned.temp.counts[tu.matches, col_idx] %>% t() )
        } else
        {
          chr.counts = rbind(chr.counts,
                             do.call(ifelse(RPK,'colMeans','colSums'), 
                                     list(binned.temp.counts[tu.matches, col_idx])) )
        }
      }
      for (i in seq_along(gene.temp.minus))
      {
        tu.matches = queryHits(minus.matches)[which(subjectHits(minus.matches) == i)]
        if (length(tu.matches) == 1)
        {
          chr.counts = rbind(chr.counts,
                             binned.temp.counts[tu.matches, col_idx+1]%>%t() )
        } else
        {
          chr.counts = rbind(chr.counts,
                             do.call(ifelse(RPK,'colMeans','colSums'),
                                     list(binned.temp.counts[tu.matches, col_idx+1])) )
        }
      }
      if (RPK) chr.counts = chr.counts * 5
      return(chr.counts)
    }
  }
  
  cal_exon_intron_count <<- function(binned.RNA.counts.list,
                                    Gene_input, 
                                    TU.genes.gr,
                                    chr.lengths,
                                    RPK = FALSE, 
                                    countExon = TRUE, ...)
  {
    # count TU reads from binned list
    # Args:
    # binned.RNA.counts.list
    # Gene_input: gene reference
    # TU.genes.gr: sum counts in bin inside subject ranges
    # chr.lengths: chromosome length for rebuilding bins
    # RPK: logical, if returns reads per Kb (RPK)
    
    temp.sample.num = ncol(binned.RNA.counts.list[[1]]) / 2
    if (temp.sample.num == 1)
    {
      binned.RNA.counts.list = lapply(binned.RNA.counts.list, function(x) cbind(x, x))
      temp.sample.num = ncol(binned.RNA.counts.list[[1]]) / 2
    }
    
    chrs = as.character(unique(seqnames(TU.genes.gr)))
    foreach(chr = chrs, .combine = rbind) %dopar%
      {
        binned.temp.counts = binned.RNA.counts.list[[chr]]
        gene.temp.plus = TU.genes.gr[seqnames(TU.genes.gr) == chr & strand(TU.genes.gr) == '+']
        gene.temp.minus = TU.genes.gr[seqnames(TU.genes.gr) == chr & strand(TU.genes.gr) == '-']
        
        exon.temp.plus = Gene_input[Gene_input$type == 'exon' &
                                      seqnames(Gene_input) == chr &
                                      strand(Gene_input) == '+']
        exon.temp.minus = Gene_input[Gene_input$type == 'exon' &
                                       seqnames(Gene_input) == chr &
                                       strand(Gene_input) == '-']
        
        exon.plus.matches = findOverlaps(query = exon.temp.plus,
                                         subject = gene.temp.plus,
                                         type='within')
        exon.minus.matches = findOverlaps(query = exon.temp.minus,
                                          subject = gene.temp.minus,
                                          type='within')
        
        bin.ranges =  IRanges(start =  seq(1, chr.lengths[chr], by=binning), width = binning)
        bin.gene.plus.matches = findOverlaps(query = bin.ranges,
                                             subject = ranges(gene.temp.plus),
                                             minoverlap = binning / 2)
        bin.gene.minus.matches = findOverlaps(query = bin.ranges,
                                              subject = ranges(gene.temp.minus),
                                              minoverlap = binning / 2)
        
        chr.counts = NULL
        col_idx = seq(1, temp.sample.num * 2, 2)
        for (i in seq_along(gene.temp.plus))
        {
          exon.id.temp = 
            IRanges::reduce(exon.temp.plus[queryHits(exon.plus.matches)[subjectHits(exon.plus.matches) == i]])
          bin.gene.matches = 
            queryHits(bin.gene.plus.matches)[subjectHits(bin.gene.plus.matches) == i]
          bin.exon.matches = 
            bin.gene.matches[countQueryHits(findOverlaps(query = bin.ranges[bin.gene.matches],
                                                         subject = ranges(exon.id.temp),
                                                         minoverlap = binning * 2/3) ) > 0]
          bin.gene.matches = bin.gene.matches[bin.gene.matches < nrow(binned.temp.counts)]
          bin.exon.matches = bin.exon.matches[bin.exon.matches < nrow(binned.temp.counts)]
          if (countExon)
          {
            if (length(bin.exon.matches) == 1)
            {
              chr.counts = rbind(chr.counts,
                                 t(binned.temp.counts[bin.exon.matches, col_idx]) )
            } else {
              chr.counts = rbind(chr.counts,
                                 do.call(ifelse(RPK, 'colMeans', 'colSums'), 
                                         list(binned.temp.counts[bin.exon.matches, col_idx])) )
            }
          } else {
            if (length(setdiff(bin.gene.matches, bin.exon.matches)) == 1)
            {
              chr.counts = 
                rbind(chr.counts,
                      t(binned.temp.counts[setdiff(bin.gene.matches, bin.exon.matches), col_idx]) )
            } else {
              chr.counts = 
                rbind(chr.counts,
                      do.call(ifelse(RPK, 'colMeans', 'colSums'),
                              list(binned.temp.counts[setdiff(bin.gene.matches, bin.exon.matches), col_idx])) )
            }
          }
        }
        
        for (i in seq_along(gene.temp.minus))
        {
          exon.id.temp = 
            IRanges::reduce(exon.temp.minus[queryHits(exon.minus.matches)[subjectHits(exon.minus.matches) == i]])
          
          bin.gene.matches = 
            queryHits(bin.gene.minus.matches)[subjectHits(bin.gene.minus.matches) == i]
          
          bin.exon.matches = 
            bin.gene.matches[findOverlaps(query = bin.ranges[bin.gene.matches],
                                          subject = ranges(exon.id.temp),
                                          minoverlap = binning * 2/3) %>% countQueryHits()>0]
          bin.gene.matches = bin.gene.matches[bin.gene.matches < nrow(binned.temp.counts)]
          bin.exon.matches = bin.exon.matches[bin.exon.matches < nrow(binned.temp.counts)]
          if (countExon)
          {
            if (length(bin.exon.matches) == 1)
            {
              chr.counts = rbind(chr.counts,
                                 binned.temp.counts[bin.exon.matches, col_idx+1] %>% t() )
            } else
            {
              chr.counts = rbind(chr.counts,
                                 do.call(ifelse(RPK,'colMeans','colSums'), 
                                         list(binned.temp.counts[bin.exon.matches, col_idx+1])) )
            }
          } else
          {
            if(length(setdiff(bin.gene.matches,bin.exon.matches)) == 1)
            {
              chr.counts = 
                rbind(chr.counts,
                      binned.temp.counts[setdiff(bin.gene.matches, bin.exon.matches), col_idx+1]%>%t() )
            } else
            {
              chr.counts = 
                rbind(chr.counts,
                      do.call(ifelse(RPK, 'colMeans', 'colSums'),
                              list(binned.temp.counts[setdiff(bin.gene.matches,bin.exon.matches), col_idx+1])) )
            }
          }
        }
        mat = chr.counts
        colnames(mat) = paste0('Sample', seq_len(temp.sample.num))
        if (RPK) mat = mat * 5
        gc()
        return(mat)
      }
  }
  
  ### part3 ------------------
  # count exon coverage 
  cutoff_expression <<- function(L.sample.list, all.TU.gr.list, all.binned.counts.list)
  {
    # inputs: bin counts and annotated TUs of all samples
    # outputs: bin counts and RPK of TUs after normalizing by exon counts across all samples
    
    .jaccardCalculation = function(minBC, all.TU.gr.list)
    {
      ncRNA.list = lapply(all.TU.gr.list, function(x) x[x$expr>minBC & x$type == 'ncRNA'])
      common.ncRNA.gr = get_common_gene(ncRNA.list)
      all.ncRNA.gr = IRanges::reduce(Reduce(c, ncRNA.list))
      jacc = length(common.ncRNA.gr) / length(all.ncRNA.gr)
      return(jacc)
    }
    
    # count all reads coverage, include spliced and non-spliced reads
    size.factors <<- sampleSizeFactors(all.TU.gr.list = all.TU.gr.list,
                                       all.binned.counts.list = all.binned.counts.list,
                                       Gene_input, L.sample.list)
    
    flat.binned.cov.list = count_list_combine(all.binned.counts.list)
    flat.binned.cov.list = lapply(flat.binned.cov.list, 
                                  function(x) sweep(x, 2, rep(size.factors, each = 2), '/') )
    size.binned.cov.list = count.list.split(flat.binned.cov.list)
    
    for (x in seq_along(all.TU.gr.list))
    {
      all.TU.gr.list[[x]] = annotatedExpression(transcriptAnno = all.TU.gr.list[[x]],
                                                binned.RNA.counts.list = size.binned.cov.list[[x]],
                                                binSize = binning)
    }
    
    if (jaccardCutoff) 
    {
      expr.array <<- seq(0, Reduce(mean,
                                   lapply(all.TU.gr.list,
                                          function(x) median(x[x$gene_type == 'ncRNA']$expr))),
                         by = 1)
      jaccardScores <<- sapply(expr.array, .jaccardCalculation, all.TU.gr.list = all.TU.gr.list)
      cutoff.BC <<- expr.array[which.max(diff(jaccardScores))]
    } else 
    {
      cutoff.BC <<- 0
    }
    # cutoff uncertain ncRNAs
    all.TU.gr.list = lapply(all.TU.gr.list, 
                            function(x) sort(c(x[x$type == 'protein_coding'],
                                               x[x$type == 'ncRNA' & x$expr > cutoff.BC])) )
    return(all.TU.gr.list)
  }
  
  sampleSizeFactors <<- function(all.TU.gr.list, all.binned.counts.list, Gene_input, L.sample.list)
  {
    .tempId = function(x) sapply(x, function(y) unlist(strsplit(y,'\\.'))[1] ) # remove gene version suffix
    common.gene.gr = get_common_gene(lapply(all.TU.gr.list, function(x) x[x$type == 'protein_coding']))
    # size factor from binned list
    Gene_exon = Gene_input[ Gene_input$type == 'exon' & Gene_input$gene_type == "protein_coding"]
    Input_ids = unique(Gene_exon$gene_id)
    
    if (exists('all.strandFlipped') && any(all.strandFlipped))
    {
      for (bam in which(all.strandFlipped))
      {
        all.binned.counts.list[[bam]] = lapply(all.binned.counts.list[[bam]],
                                        function(x) x[, rev(seq_len(ncol(x)))])
      }
    }
    gc()
    flat.binned.cov.list = count_list_combine(all.binned.counts.list)
    exon.counts = cal_exon_intron_count(binned.RNA.counts.list = flat.binned.cov.list, 
                                        Gene_input = Gene_input, 
                                        TU.genes.gr = common.gene.gr,
                                        chr.lengths,
                                        RPK = FALSE, 
                                        countExon = TRUE)
    size_factor_cal(exon.counts, RefExpr = NULL)
  }
  
  annotatedExpression <<- function(transcriptAnno, binned.RNA.counts.list, binSize = binning)
  {
    allExpr = NA
    if (ncol(binned.RNA.counts.list[[1]]) == 2)
    {
      
      allExpr = foreach(i = seq_along(transcriptAnno), .combine="c") %dopar%
        {
          sum(binned.RNA.counts.list[[as.character(seqnames(transcriptAnno[i]))]][
            round(start(transcriptAnno[i]) / binSize) : round(end(transcriptAnno[i]) / binSize),
            ifelse(strand(transcriptAnno[i]) == "+", 1, 2)
            ])
        }
      
    } else
    {
      col_idx = seq(1, ncol(binned.RNA.counts.list[[1]]), 2)
      allExpr = foreach(i=seq_along(transcriptAnno), .combine="c") %dopar%
        {
          sum(binned.RNA.counts.list[[as.character(seqnames(transcriptAnno[i]))]][
            round(start(transcriptAnno[i]) / binSize) : round(end(transcriptAnno[i]) / binSize),
            ifelse(strand(transcriptAnno[i]) == "+", col_idx, col_idx + 1)]
          ) / length(col_idx)
        }
    }
    transcriptAnno$expr = round(allExpr, 2)
    return(transcriptAnno)
  }
  
  # size.factor----------------------------------
  size_factor_cal <<- function(gene.counts, RefExpr=NULL)
  {
    # Args:
    # gene.counts: genes by row and samples by column
    # RefExpr: if return relative size factors comparing to reference counts
    # return the size factors of gene counts
    gene.counts = gene.counts[!apply(gene.counts,1,function(x) any(is.na(x)) ) & rowSums(gene.counts)>0,]
    if (!is.null(RefExpr))
    {
      gene.counts = cbind(RefExpr, gene.counts)
      # geometric mean
      geoMeans = apply(gene.counts, 1, function(x) exp(sum(log(x[x > 0]), na.rm = TRUE) / length(x)) )
      # quotients list to median (size factor)
      quoMedian = apply(sweep(gene.counts, 1, geoMeans, '/'), 2, median)
      return(quoMedian[-1] / quoMedian[1])
    } else
    {
      geoMeans = apply(gene.counts, 1, function(x) exp(sum(log(x[x > 0]), na.rm = TRUE) / length(x)) )
      return(apply(sweep(gene.counts, 1, geoMeans,'/'), 2, median))
    }
  }
})
# Run the app
shinyApp(ui = ui, server = server)
