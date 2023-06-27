#https://bioconductor.org/packages/devel/bioc/vignettes/APAlyzer/inst/doc/APAlyzer.html#analysis-of-apa-in-3utrs
#Apalyzer 
#Installing Apalyzer

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("APAlyzer")
BiocManager::install("Rsamtools")
library("APAlyzer")
library("base")
library("data.table")
library("Rsamtools")
library("tidyverse")
library("stats")
library(shiny)
library(purrr)
library(data.table)
library(base)
library(tidyverse)
library(dplyr)
library(splitstackshape)
library(shinythemes)
library(readr)
library(ggplot2)
library("repmis")
options(shiny.maxRequestSize=10000000*1024^2)

ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage(
                  "APAtizer",
                  tabPanel("DAPARS",
                           titlePanel("Use DaPars2 to analyse 3'UTR-APA from RNA-Seq data"),
                           sidebarLayout(
                             sidebarPanel(
                               fileInput("txt_files", "Select Multiple .txt Files", multiple = TRUE),
                               fileInput("txt_file2", "Select Single .txt File"),
                               actionButton("run2", "DaPars Analysis")
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Len genes",
                                          br(),
                                          textInput(inputId = "search_term", label = "Search for a gene"),
                                          downloadButton("download_datax", "Download Len Gene Data"),
                                          tableOutput("data_table_x")
                                          
                                 ),
                                 tabPanel("Short Genes",
                                          br(),
                                          textInput(inputId = "search_term2", label = "Search for a gene"),
                                          downloadButton("download_datay", "Download Short Gene Data"),
                                          tableOutput("data_table_y")
                                          
                                 ),
                                 
                               )
                             )
                           )
                           
                           ),
                  tabPanel("APA APALYZER",
                           titlePanel("Use Apalyzer to analyse 3'UTR-APA from RNA-Seq data"),
                           sidebarLayout(
                             sidebarPanel(
                               textInput("path2", "Directory path:"),
                               fileInput("txt_file3", "Select Single .txt File"),
                               actionButton("run3", "APA Analysis"),
                               br(),
                               br(),
                               br(),
                               selectInput(inputId = "output_type4", label = "Select Output Type",
                                           choices = c("NvsT_APA_UP", "NvsT_APA_DN","NvsT_APA_NC" )),
                               br(),
                               selectInput(inputId = "output_type5", label = "Select Plot Type",
                                           choices = c("IPA Volcano plot (top 40)", "IPA Volcano plot", "IPA Box"))
                               
                             ),
                             
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Number of APA events",
                                          tableOutput("data_table_2")
                                          
                                 ),
                                 
                                 tabPanel("NvsT_APA",
                                          conditionalPanel(
                                            condition = "input.output_type4 == 'NvsT_APA_UP'",
                                            br(),
                                            textInput(inputId = "search_term3", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data2", label = "Download NvsT_APA_UP"),
                                            tableOutput(outputId = "data_table_3")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type4 == 'NvsT_APA_DN'",
                                            br(),
                                            textInput(inputId = "search_term4", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data3", label = "Download NvsT_APA_DN"),
                                            tableOutput(outputId = "data_table_4")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type4 == 'NvsT_APA_NC'",
                                            br(),
                                            textInput(inputId = "search_term5", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data4", label = "Download NvsT_APA_NC"),
                                            tableOutput(outputId = "data_table_5")
                                          )
                                 ),
                                 
                                 tabPanel("Plots",
                                          conditionalPanel(
                                            condition = "input.output_type5 == 'IPA Volcano plot (top 40)'",
                                            br(),
                                            downloadButton(outputId = "download_plot1", label = "Download IPA Volcano plot (top 40)"),
                                            plotOutput(outputId = "plot1")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type5 == 'IPA Volcano plot'",
                                            br(),
                                            downloadButton(outputId = "download_plot2", label = "Download IPA Volcano plot"),
                                            plotOutput(outputId = "plot2")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type5 == 'IPA Box'",
                                            br(),
                                            downloadButton(outputId = "download_plot3", label = "Download IPA Box"),
                                            plotOutput(outputId = "plot3")
                                          )
                                 )
                               
                               )
                               
                             )
                           )
                           
                           ),
                  
                  tabPanel("IPA APALYZER",
                           titlePanel("Use Apalyzer to analyse IPA from RNA-Seq data"),
                           sidebarLayout(
                             sidebarPanel(
                               textInput("path", "Directory path:"),
                               fileInput("txt_file", "Select Single .txt File"),
                               actionButton("run", "IPA Analysis"),
                               br(),
                               br(),
                               br(),
                               selectInput(inputId = "output_type", label = "Select Output Type",
                                           choices = c("NvsT_IPA_events_UP", "NvsT_IPA_events_DN","NvsT_IPA_events_NC" )),
                               br(),
                               selectInput(inputId = "output_type2", label = "Select Output Type",
                                           choices = c("NvsT_IPA_genes_UP", "NvsT_IPA_genes_DN", "NvsT_IPA_genes_NC")),
                               br(),
                               selectInput(inputId = "output_type3", label = "Select Plot Type",
                                           choices = c("IPA Volcano plot (top 40)", "IPA Volcano plot", "IPA Box"))
                               
                             ),
                             
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Number of IPA events",
                                          tableOutput("data_table_6")
                                          
                                 ),
                                 tabPanel("NvsT_IPA_events",
                                          conditionalPanel(
                                            condition = "input.output_type == 'NvsT_IPA_events_UP'",
                                            br(),
                                            textInput(inputId = "search_term6", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data5", label = "Download NvsT_IPA_events_UP"),
                                            tableOutput(outputId = "data_table_7")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type == 'NvsT_IPA_events_DN'",
                                            br(),
                                            textInput(inputId = "search_term7", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data6", label = "Download NvsT_IPA_events_DN"),
                                            tableOutput(outputId = "data_table_8")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type == 'NvsT_IPA_events_NC'",
                                            br(),
                                            textInput(inputId = "search_term8", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data7", label = "Download NvsT_IPA_events_NC"),
                                            tableOutput(outputId = "data_table_9")
                                          )
                                          
                                 ),
                                 
                                 tabPanel("NvsT_IPA_genes",
                                          conditionalPanel(
                                            condition = "input.output_type2 == 'NvsT_IPA_genes_UP'",
                                            br(),
                                            textInput(inputId = "search_term9", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data8", label = "Download NvsT_IPA_genes_UP"),
                                            tableOutput(outputId = "data_table_10")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type2 == 'NvsT_IPA_genes_DN'",
                                            br(),
                                            textInput(inputId = "search_term10", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data9", label = "Download NvsT_IPA_genes_DN"),
                                            tableOutput(outputId = "data_table_11")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type2 == 'NvsT_IPA_genes_NC'",
                                            br(),
                                            textInput(inputId = "search_term11", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data10", label = "Download NvsT_IPA_genes_NC"),
                                            tableOutput(outputId = "data_table_12")
                                          )
                                 ),
                                 
                                 tabPanel("Plots",
                                          conditionalPanel(
                                            condition = "input.output_type3 == 'IPA Volcano plot (top 40)'",
                                            br(),
                                            downloadButton(outputId = "download_plot4", label = "Download IPA Volcano plot (top 40)"),
                                            plotOutput(outputId = "plot4")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type3 == 'IPA Volcano plot'",
                                            br(),
                                            downloadButton(outputId = "download_plot5", label = "Download IPA Volcano plot"),
                                            plotOutput(outputId = "plot5")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type3 == 'IPA Box'",
                                            br(),
                                            downloadButton(outputId = "download_plot6", label = "Download IPA Box"),
                                            plotOutput(outputId = "plot6")
                                          )
                                 )
                               )
                             )
                           )
                  )
                )
)

server <- function(input, output,session) {
  
  df_pacientes <- eventReactive(input$run,{
    file <- input$txt_file
    if (is.null(file)) {
      return(NULL)
    }
    df_pacientes <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes <- select(df_pacientes, File.Name,Case.ID,Sample.Type)
    df_pacientes$File.Name<- str_replace(df_pacientes$File.Name, "rna_seq.genomic.gdc_realn.bam", "rna_seq.genomic.gdc_realn.trim.bam")
    
    case_order <- unique(df_pacientes$Case.ID)
    df_pacientes <- df_pacientes %>% arrange(df_pacientes$Sample.Type, match(df_pacientes$Case.ID, case_order))
    
    #normal <- c("Solid Tissue Normal")
    df_pacientes$category <- ifelse(df_pacientes$Sample.Type %in% c("Solid Tissue Normal"), "Normal", "Tumor")
    df_pacientes$category = paste(df_pacientes$Case.ID, df_pacientes$category, sep="_")
    
    return(df_pacientes)
    
    
  })
  
  
  NvsT_IPA <- eventReactive(input$run, {
    datapath <- input$path
    datapath <- str_replace_all(datapath, "\\\\", "/")
    flsall <- dir(datapath,".bam")
    flsall <- paste0(datapath, '/', df_pacientes()$File.Name)
    names(flsall) <- df_pacientes()$category
    #Genomic reference
    library("repmis")
    URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
    file="hg38_REF.RData"
    source_data(paste0(URL,file,"?raw=True"))
    
    refUTRraw=refUTRraw_hg38
    dfIPAraw=dfIPA_hg38
    dfLEraw=dfLE_hg38
    PASREF=REF4PAS(refUTRraw,dfIPAraw,dfLEraw)
    UTRdbraw=PASREF$UTRdbraw
    dfIPA=PASREF$dfIPA
    dfLE=PASREF$dfLE   
    dfIPA=dfIPA
    dfLE=dfLE
    IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="forward", SeqType="ThreeMostPairEnd")
    
    x=nrow(df_pacientes())/2
    sampleTable2 = data.frame(samplename = c(names(flsall)),
                              condition = c(rep("KD",x),rep("NT",x)))
    
    NvsT_IPA=APAdiff(sampleTable2,
                     IPA_OUTraw, 
                     conKET='NT',
                     trtKEY='KD',
                     PAS='IPA',
                     CUTreads=5,
                     p_adjust_methods="fdr")
    return(NvsT_IPA)
  })
  
  Nr_IPA_events <- eventReactive(input$run,{
    Nr_IPA_events<- table(NvsT_IPA()$APAreg)  
    return(Nr_IPA_events)
  })
  # NvsT_IPA_UP
  NvsT_IPA_events_UP <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_UP <- NvsT_IPA[ which(NvsT_IPA$APAreg=='UP'),]
    
    return(NvsT_IPA_events_UP)
  })
  
  # NvsT_IPA_DN
  NvsT_IPA_events_DN <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_DN <- NvsT_IPA[ which(NvsT_IPA$APAreg=='DN'),]
    
    return(NvsT_IPA_events_DN)
  })
  
  # NvsT_IPA_NC
  NvsT_IPA_events_NC <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_NC <- NvsT_IPA[ which(NvsT_IPA$APAreg=='NC'),]
    
    return(NvsT_IPA_events_NC)
  })
  
  NvsT_IPA_genes_UP <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_UP <- NvsT_IPA_events_UP()
    NvsT_IPA_genes_UP <- distinct(NvsT_IPA_events_UP,select=c(gene_symbol))
    
    return(NvsT_IPA_genes_UP)
  })
  
  # NvsT_IPA_DN
  NvsT_IPA_genes_DN <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_DN <- NvsT_IPA_events_DN()
    NvsT_IPA_genes_DN <- distinct(NvsT_IPA_events_DN,select=c(gene_symbol))
    
    return(NvsT_IPA_genes_DN)
  })
  
  # NvsT_IPA_NC
  NvsT_IPA_genes_NC <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_NC <- NvsT_IPA_events_NC()
    NvsT_IPA_genes_NC <- distinct(NvsT_IPA_events_NC,select=c(gene_symbol))
    
    return(NvsT_IPA_genes_NC)
  })
  
  e <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    e <- APAVolcano(NvsT_IPA, PAS='IPA', Pcol = "pvalue", top=40)
    
    return(e)
  })
  f <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    f <- APAVolcano(NvsT_IPA, PAS='IPA', Pcol = "pvalue")
    
    return(f)
  })
  g <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    g <- APABox(NvsT_IPA, xlab = "APAreg", ylab = "RED", plot_title = NULL)
    
    return(g)
  })
  
  
  output$data_table_6 <- renderTable({
    Nr_IPA_events()
  })
  
  output$data_table_7 <- renderTable({
    if (input$search_term6 != "") {
      NvsT_IPA_events_UP() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term6, ignore_case = TRUE))))
    } else {
      NvsT_IPA_events_UP()
    }
  })
  
  output$data_table_8 <- renderTable({
    if (input$search_term7 != "") {
      NvsT_IPA_events_DN() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term7, ignore_case = TRUE))))
    } else {
      NvsT_IPA_events_DN()
    }
  })
  
  output$data_table_9 <- renderTable({
    if (input$search_term8 != "") {
      NvsT_IPA_events_NC() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term8, ignore_case = TRUE))))
    } else {
      NvsT_IPA_events_NC()
    }
  })
  
  output$data_table_10 <- renderTable({
    if (input$search_term9 != "") {
      NvsT_IPA_genes_UP() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term9, ignore_case = TRUE))))
    } else {
      NvsT_IPA_genes_UP()
    }
  })
  
  output$data_table_11 <- renderTable({
    if (input$search_term10 != "") {
      NvsT_IPA_genes_DN() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term10, ignore_case = TRUE))))
    } else {
      NvsT_IPA_genes_DN()
    }
  })
  
  output$data_table_12 <- renderTable({
    if (input$search_term11 != "") {
      NvsT_IPA_genes_NC() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term11, ignore_case = TRUE))))
    } else {
      NvsT_IPA_genes_NC()
    }
  })
  
  output$plot4 <- renderPlot({
    e()
  })
  output$plot5 <- renderPlot({
    f()
  })
  output$plot6 <- renderPlot({
    g()
  })
  
  
  
  output$download_data5 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_events_UP", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_events_UP(), file, row.names = FALSE)
    }
  )
  output$download_data6 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_events_DN", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_events_DN(), file, row.names = FALSE)
    }
  )
  output$download_data7 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_events_NC", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_events_NC(), file, row.names = FALSE)
    }
  )
  output$download_data8 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_genes_UP", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_genes_UP(), file, row.names = FALSE)
    }
  )
  output$download_data9 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_genes_DN", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_genes_DN(), file, row.names = FALSE)
    }
  )
  output$download_data10 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_genes_NC", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_genes_NC(), file, row.names = FALSE)
    }
  )
  output$download_plot4 <- downloadHandler(
    filename = function() {
      paste("plot4", ".png")
    },
    content = function(file) {
      ggsave(file, e(), width = 6000, height = 4000,units = c("px"),dpi = 300)
    }
  )
  output$download_plot5 <- downloadHandler(
    filename = function() {
      paste("plot5", ".png")
    },
    content = function(file) {
      ggsave(file, f(),width = 800, height = 700,units = c("px"),dpi = 300)
    }
  )
  output$download_plot6 <- downloadHandler(
    filename = function() {
      paste("plot6", ".png")
    },
    content = function(file) {
      ggsave(file, g(), dpi = 300)
      
    }
  )
  
  df_pacientes2 <- eventReactive(input$run2,{
    file <- input$txt_file2
    if (is.null(file)) {
      return(NULL)
    }
    df_pacientes2 <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes2 <- select(df_pacientes2, File.Name,Case.ID,Sample.Type)
    df_pacientes2$File.Name <- paste("WIG/", df_pacientes2$File.Name, sep = "")
    case_order <- unique(df_pacientes2$Case.ID)
    df_pacientes2 <- df_pacientes2 %>% arrange(df_pacientes2$Sample.Type, match(df_pacientes2$Case.ID, case_order))
    normal <- c("Solid Tissue Normal")
    df_pacientes2$category <- ifelse(df_pacientes2$Sample.Type %in% normal, "Normal", "Tumor")
    df_pacientes2$category = paste(df_pacientes2$Case.ID, df_pacientes2$category, sep="_")
    df_pacientes2$File.Name<- str_replace(df_pacientes2$File.Name, ".bam", "_PDUI")
    return(df_pacientes2)
    
    
  })
  
  # DPDUI
  dpdui <- eventReactive(input$run2,{
    files <- input$txt_files
    if (is.null(files)) {
      return(NULL)
    }
    
    
    df <- rbindlist(sapply(files$datapath, fread, simplify = FALSE), use.names = TRUE)
    df <- subset(df, select = -c(fit_value,Loci,Predicted_Proximal_APA))
    
    
    idx <- match(df_pacientes2()$File.Name, colnames(df))
    idx <- append(1, idx)
    
    df <- df[, ..idx]
    
    num_cols = ncol(df)
    
    colnames(df)[2:num_cols]<-df_pacientes2()$category
    
    num_cols2 = ((ncol(df)-1)/2)+1
    num_cols3 = ((ncol(df)-1)/2)+2
    num_cols4 = ncol(df)+1
    
    for (i in colnames(df)[2:num_cols2]) {
      df[[i]][is.na(df[[i]])] <- rowMeans(df[,2:num_cols2], na.rm = TRUE)[is.na(df[[i]])]
    }
    
    for (i in colnames(df)[num_cols3:num_cols]) {
      df[[i]][is.na(df[[i]])] <- rowMeans(df[,num_cols3:num_cols], na.rm = TRUE)[is.na(df[[i]])]
    }
    
    df <- data.frame(df)  
    res <- df[, grepl("_Tumor", colnames(df))] - df[, grepl("_Normal", colnames(df))]
    
    colnames(res) <- paste(colnames(df[, grepl("_Tumor", colnames(df))]),
                           colnames(df[, grepl("_Normal", colnames(df))]), sep = "-")
    
    df <-cbind(df, res)
    
    num_cols5 = ncol(df)
    
    dpdui <- df[,c(1,num_cols4:num_cols5)]
    
    return(dpdui)
  })
  
  # SHORT GENES
  short_genes <- eventReactive(input$run2,{
    dpdui <- dpdui()
    dpdui$mean <- rowMeans(dpdui[,2:length(dpdui())],na.rm=TRUE) 
    dpdui$Len_Bladder <- dpdui[,c("mean")] >= 0.2
    dpdui$Short_Bladder <- dpdui[,c("mean")] <= -0.2
    gene_len_Bladder <- dpdui[which(dpdui$Len_Bladder == 1), ]                             
    gene_short_Bladder <- dpdui[which(dpdui$Short_Bladder == 1), ]  
    
    
    return(gene_short_Bladder)
  })
  
  # LEN GENES
  len_genes <- eventReactive(input$run2,{
    dpdui <- dpdui()
    dpdui$mean <- rowMeans(dpdui[,2:length(dpdui())],na.rm=TRUE) 
    dpdui$Len_Bladder <- dpdui[,c("mean")] >= 0.2
    dpdui$Short_Bladder <- dpdui[,c("mean")] <= -0.2
    gene_len_Bladder <- dpdui[which(dpdui$Len_Bladder == 1), ]                             
    gene_short_Bladder <- dpdui[which(dpdui$Short_Bladder == 1), ]  
    
    
    return(gene_len_Bladder)
  })
  

  # Display combined data as a table
  output$data_table_x <- renderTable({
    if (input$search_term != "") {
      len_genes() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term, ignore_case = TRUE))))
    } else {
      len_genes()
    }
  })
  
  output$data_table_y <- renderTable({
    if (input$search_term2 != "") {
      short_genes() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term2, ignore_case = TRUE))))
    } else {
      short_genes()
    }
  })
  
  
  # Download combined data as a .csv file

  output$download_datax <- downloadHandler(
    filename = function() {
      paste("Len_Genes", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(len_genes(), file, row.names = FALSE)
    }
  )
  output$download_datay <- downloadHandler(
    filename = function() {
      paste("Short_Genes", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(short_genes(), file, row.names = FALSE)
    }
  )

  
  df_pacientes3 <- eventReactive(input$run3,{
    file <- input$txt_file3
    if (is.null(file)) {
      return(NULL)
    }
    df_pacientes3 <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes3 <- select(df_pacientes3, File.Name,Case.ID,Sample.Type)
    df_pacientes3$File.Name<- str_replace(df_pacientes3$File.Name, "rna_seq.genomic.gdc_realn.bam", "rna_seq.genomic.gdc_realn.trim.bam")
    
    case_order <- unique(df_pacientes3$Case.ID)
    df_pacientes3 <- df_pacientes3 %>% arrange(df_pacientes3$Sample.Type, match(df_pacientes3$Case.ID, case_order))
    
    #normal <- c("Solid Tissue Normal")
    df_pacientes3$category <- ifelse(df_pacientes3$Sample.Type %in% c("Solid Tissue Normal"), "Normal", "Tumor")
    df_pacientes3$category = paste(df_pacientes3$Case.ID, df_pacientes3$category, sep="_")
    
    return(df_pacientes3)
    
    
  })
  
  NvsT_APA <- eventReactive(input$run3,{
    datapath <- input$path2
    datapath <- str_replace_all(datapath, "\\\\", "/")
    flsall <- dir(datapath,".bam")
    flsall <- paste0(datapath, '/', df_pacientes3()$File.Name)
    names(flsall) <- df_pacientes3()$category
    #Genomic reference
    URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
    file="hg38_REF.RData"
    source_data(paste0(URL,file,"?raw=True"))
    
    refUTRraw=refUTRraw_hg38
    dfIPAraw=dfIPA_hg38
    dfLEraw=dfLE_hg38
    PASREF=REF4PAS(refUTRraw,dfIPAraw,dfLEraw)
    UTRdbraw=PASREF$UTRdbraw
    dfIPA=PASREF$dfIPA
    dfLE=PASREF$dfLE   
    
    #Analysis of APA in 3’UTRs
    refUTRraw=refUTRraw
    UTRdbraw=REF3UTR(refUTRraw)
    DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="forward")
    x=nrow(df_pacientes3())/2 
    sampleTable1 = data.frame(samplename = c(names(flsall)),
                              condition = c(rep("KD",x),rep("NT",x)))
    NvsT_APA=APAdiff(sampleTable1,DFUTRraw, 
                     conKET='NT',
                     trtKEY='KD',
                     PAS='3UTR',
                     CUTreads=5,
                     p_adjust_methods="fdr")
    return(NvsT_APA)
  })
  
  Nr_APA_events <- eventReactive(input$run3,{
    #NvsT_APA <- NvsT_APA()
    Nr_APA_events<- table(NvsT_APA()$APAreg)  
    return(Nr_APA_events)
  })
  # NvsT_APA_UP
  NvsT_APA_UP <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    NvsT_APA_UP <- NvsT_APA[ which(NvsT_APA$APAreg=='UP'),]
    
    return(NvsT_APA_UP)
  })
  
  # NvsT_APA_DN
  NvsT_APA_DN <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    NvsT_APA_DN <- NvsT_APA[ which(NvsT_APA$APAreg=='DN'),]
    
    return(NvsT_APA_DN)
  })
  
  # NvsT_APA_NC
  NvsT_APA_NC <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    NvsT_APA_NC <- NvsT_APA[ which(NvsT_APA$APAreg=='NC'),]
    
    return(NvsT_APA_NC)
  })
  
  a <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    a <- APAVolcano(NvsT_APA, PAS='3UTR', Pcol = "pvalue", top=40)
    
    return(a)
  })
  b <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    b <- APAVolcano(NvsT_APA, PAS='3UTR', Pcol = "pvalue")
    
    return(b)
  })
  d <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    d <- APABox(NvsT_APA, xlab = "APAreg", ylab = "RED", plot_title = NULL)
    
    return(d)
  })
  
  
  
  
  output$data_table_2 <- renderTable({
    Nr_APA_events()
  })
  
  output$data_table_3 <- renderTable({
    if (input$search_term3 != "") {
      NvsT_APA_UP() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term3, ignore_case = TRUE))))
    } else {
      NvsT_APA_UP()
    }
  })
  
  output$data_table_4 <- renderTable({
    if (input$search_term4 != "") {
      NvsT_APA_DN() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term4, ignore_case = TRUE))))
    } else {
      NvsT_APA_DN()
    }
  })
  
  output$data_table_5 <- renderTable({
    if (input$search_term5 != "") {
      NvsT_APA_NC() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term5, ignore_case = TRUE))))
    } else {
      NvsT_APA_NC()
    }
  })
  
  
  
  output$plot1 <- renderPlot({
    a()
  })
  output$plot2 <- renderPlot({
    b()
  })
  output$plot3 <- renderPlot({
    d()
  })
  
  
  
  
  output$download_data2 <- downloadHandler(
    filename = function() {
      paste("NvsT_APA_UP", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_APA_UP(), file, row.names = FALSE)
    }
  )
  output$download_data3 <- downloadHandler(
    filename = function() {
      paste("NvsT_APA_DN", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_APA_DN(), file, row.names = FALSE)
    }
  )
  output$download_data4 <- downloadHandler(
    filename = function() {
      paste("NvsT_APA_NC", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_APA_NC(), file, row.names = FALSE)
    }
  )
  output$download_plot1 <- downloadHandler(
    filename = function() {
      paste("plot1", ".png")
    },
    content = function(file) {
      ggsave(file, a(), width = 6000, height = 4000,units = c("px"),dpi = 300)
    }
  )
  output$download_plot2 <- downloadHandler(
    filename = function() {
      paste("plot2", ".png")
    },
    content = function(file) {
      ggsave(file, b(),width = 800, height = 700,units = c("px"),dpi = 300)
    }
  )
  output$download_plot3 <- downloadHandler(
    filename = function() {
      paste("plot2", ".png")
    },
    content = function(file) {
      ggsave(file, d(), dpi = 300)
      
    }
  )
  

  
  
  
}
# Run the app
shinyApp(ui, server)
