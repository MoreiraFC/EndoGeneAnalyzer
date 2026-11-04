
rm(list = ls())

library(shiny)
library(shinyauthr)
library(tidyverse)
library(knitr)
library(kableExtra)
library(DT)
library(dunn.test)
library(reshape2)
library(ggExtra)
library(readxl)
library(writexl)
source("r.NormFindStab5.R")

#library(kableExtra)


## cálculo de outliers
calc.outlier = function (tabela, limiar){
  cont = 1
  res = list(
    outlier = tabela$outlier,
    dados=data.frame(Gene = "", group = "", sd = 0, mean = 0)
  )
  
  for (i in unique(tabela$Genes)){
    for (j in unique(tabela$group)){
      a = filter(tabela, Genes == i, group == j)
      dp = sd(a$MeanCq, na.rm = T)
      media = mean(a$MeanCq, na.rm = T)
      out.samples = as.character(a$sample[a$MeanCq > (media+limiar*dp) | a$MeanCq < (media-limiar*dp)])
      k = tabela$Genes == i & tabela$group ==j & tabela$sample %in% out.samples
      res$outlier[k] = TRUE
      res$dados[cont,1:4] = c(i, j, dp, media)
      cont=cont+1
    }
  }
  return(res)
}

# Define UI for application that draws a histogram
ui <- fluidPage(

  tags$head(
    tags$style(
      HTML("
        .header {
          display: flex;
          justify-content: space-between;
          align-items: center;
          padding: 10px;
 #         background-color: #f5f5f5;
        }
        .header-text {
          flex: 1;
          #color: darkblue;
        }
        .header-image {
          flex: 1;
          text-align: right;
        }
      ")
    )
  ),
  
  tags$div(class = "header",
           tags$div(class = "header-text", tags$img(src = "bitmap.png", width = "350px", height = "45px", style='padding-Botton:10px'),), #h1("EndoGeneAnalyzer"),), #
           tags$div(class = "header-image",
                    tags$img(src = "logonpo2.png", width = "70px", height = "70px", style='padding-top:10px'),
                    tags$img(src = "logo_ufpa.png", width = "70px", height = "70px", style='padding-left:10px')
           )
  ),
  
  tabsetPanel(#type = "pills",
    tabPanel(style='padding-top:25px', "Data Upload", 
             fluidRow(
               column(3, style='padding-left:25px', 
                      fileInput(
                        inputId = "uploadFile", 
                        label = "Load Data",
                        multiple = F,
                        accept = NULL, width = NULL,
                        placeholder = "No file selected")),
               column(3,radioButtons('type', 'File type',
                                     choices = 
                                     c('.xls/.xlsx', '.csv/.txt'), '.xls/.xlsx')
                      ),
               column(3,radioButtons('sep', 'Separator (ignore for .xls)',
                                     c(Comma=',',
                                       Semicolon=';',
                                       Tab='\t',
                                       Space=' '), '\t')
                      ),
#             column(2, actionButton(inputId = "load", label = "Load Example Data"), style="padding-top:5px"),
             ),
             fluidRow(column(7, tableOutput("confirm.table"), style="padding-left:10px"),
                      column(5, actionButton(inputId = "confirm", label = "Confirm Data Table"), style="padding-left:30px")
                      ),
             mainPanel(tableOutput(outputId = "upload.table"))
    ),
    tabPanel("Data Summary", 
             fluidRow( column(12, plotOutput(outputId = "group.plot"))),
             fluidRow(
               column(2, style="padding-left:45px",
                 checkboxGroupInput("Sel.target", "Select Targets:", choices = NULL)
               ),
               column(2,
                      actionButton("update", "Update Target Gene(s)")
               )
             ),
             mainPanel(
               tableOutput(outputId = "data_table")
               )
             ), ## Fim do tabPanel Resumo dos dados
    tabPanel("Gene Reference Samples",
             fluidRow(column(4, style=c("padding-left:25px", "padding-top:10px"), 
                             sliderInput("out.dp", "Set outlier standard deviation:",
                                  min = 1, max = 4, value = 2, step = 0.25),
                             "Obs.: Indicates the threshold for the MeanCq of a sample to be considered an outlier (default at two standard deviations)."
                             ),

                      column(4, style="padding-top:10px",
                             fluidRow(actionButton("remove", "Remove Outliers"), actionButton("restore", "Restore Outliers")),
                             fluidRow(style="padding-top:10px",
                             "Obs.: Eliminating outliers reduces the standard deviation of the group. This may reveal additional outliers."
                             ),
                      ),
                      column(4,radioButtons('mean.outlier', 'Choose which outlier to remove',
                                            choices = 
                                              c('Only Mean', 'All Outliers'), 'Only Mean'),
                               'Obs.: "Only Mean" will show and remove only outliers in the mean of the reference genes. "All Outliers" will show and remove outliers in each gene individually.'
                      ),
                      ),
             fluidRow(column(12, style="padding-top:15px", plotOutput(outputId = "sample.plot", height = "600px"))),
             fluidRow(column(12, plotOutput(outputId = "mean.plot", height = "200px"))),
             fluidRow(column(9, style="padding-left:45px", h4("List of outliers reference samples"))),
             fluidRow(column(9, style="padding-left:45px", tableOutput(outputId = "outlier.tab")))
    ),
    tabPanel("Gene Reference Analysis",
             fluidRow(column(9, style="padding-left:25px", h4("  Gene Reference by group"))),
             fluidRow(column(12, style="padding-top:15px", HTML(
               "The first line of the table is the p-value of the Kruskall-Wallis test for each reference gene. <br/>
                      The other lines are p-values of the Wilcoxon-Mann-Whitney test between each condition/group.<br/>
                      The last column is the tests p-value of the reference genes mean. <br/>
                      Red values indicate p-values < 0.05. It is desirable that reference genes and/or the reference genes mean do not vary significantly among groups/conditions."))),
             fluidRow(
               column(12, DT::dataTableOutput("endo.p.values")),
             ),
             fluidRow(column(9, style="padding-left:25px", h4("Gene Reference Descriptive Statistics"))),
             fluidRow(column(12, style="padding-top:15px", HTML(
               "It is preferable for reference genes not to vary across different conditions or groups.  <br/>
               Additionally, it is desirable for reference genes to exhibit low variances,  <br/>
               enhancing the ability to discern small yet significant differences in the target.<br/>
               If a reference gene exhibits high standard deviation and/or has low stability and/or is signifficantly altered across group/condition, consider to remove it from the analysis."))),
             fluidRow(
               column(8, style="padding-top:15px", tableOutput("descriptive")),
               column(3,  style="padding-top:15px", checkboxGroupInput(inputId = "deg.endo", 
                                            label = "Selected Reference Genes:", 
                                            choices = NULL)),
             ),
             fluidRow(column(5, style="padding-left:25px", h4("NormFinder Analysis"), 
                             tags$a(href = "https://aacrjournals.org/cancerres/article/64/15/5245/511517/Normalization-of-Real-Time-Quantitative-Reverse", 
                                    target = "_blank", "NormFinder Article"))),
             fluidRow(column(12, style="padding-top:15px", HTML(
               "NormFinder is a method to evaluate reference gene stability using the variation of the gene expression. <br/>
              The lower the stability value, the more stable the gene is considered."))),
             fluidRow(column(9, style="padding-left:25px", h5("  NormFinder Result Ordered"))),
             fluidRow(
               column(10, tableOutput(outputId = "nf.order"))
             ),
             fluidRow(column(9, style="padding-left:25px", h5("  NormFinder Result Unordered by group"))),
             fluidRow(column(12, tableOutput(outputId = "nf.unorder"))),
             fluidRow(column(9, style="padding-left:25px; padding-bottom:50px", 
                             HTML("GroupDif - Group Differences <br/> GroupSD - Group Standard Deviation <br/> 
                                  IGroupDif - Intergroup Differences<br/> IGroupSD - Intergroup Standard Deviation" ))
                      ),
    ), ## fim do tabPanel Endogenos
    
    tabPanel("Differential Expression Analysis",
             fluidRow(column(2, selectInput(inputId = "deg.target", label = "Select a Target:", choices = NULL)),
                      column(3, checkboxGroupInput(inputId = "deg.group", label = "Select two or more groups:", choices = NULL)),
                      column(3, selectInput(inputId = "plot.type", label = "Select a plot:", 
                                            choices = c("Boxplot + Dispersion", "Boxplot", "Scatter plot", "Violin"))),
                      column(3, selectInput(inputId = "plot.color", label = "Select a color:", 
                                            choices = c("Set1", "Set2", "Set3", "Blues", "Greys", "Accent", "Purple to Green", "Red to Blue")))
             ),
             fluidRow(
                      column(6, style="padding-left:45px", plotOutput(outputId = "deg.boxplot"),
                             downloadButton(outputId = "download", label = 'Download Figure')),
                      
                      column(4, h5("Fold-Change (FC) analysis:"),
                             tableOutput(outputId = "deg.FC"),
                             fluidRow(h5("Normality test:"),
                                      tableOutput(outputId  = "shapiro.test"))
                             )
             ),
             fluidRow(column(9, style="padding-left:25px", h4("Differential Expression Analysis Results"))),
             fluidRow(column(9, style="padding-left:30px", h5("Left analysis: Pearson's t-test (2 groups/conditions) or ANOVA (3 or more groups/conditions)"))),
             fluidRow(column(9, style="padding-left:30px", h5("Right analysis: Wilcoxon-Mann-Whitney test (2 groups/conditions) or Kruskal-Wallis (3 or more groups/conditions)"))),
             fluidRow(
                      column(6, verbatimTextOutput(outputId = "param.test")),
                      column(6, verbatimTextOutput(outputId = "nparam.test")),
             ),
             fluidRow(
               column(6, style="padding-left:45px", 
                      downloadButton(outputId = "report", label = 'Download D.E.A. Report'))),
             fluidRow(column(9, style="padding-left:25px", h4("Expression Data Table"))),
             mainPanel(tableOutput("deg.table"))
            
    ),
    tabPanel("Tutorial",
             fluidRow(
               tags$img(src = "Fig_1.png", width = "900px"),
                      
             ),
             fluidRow(h4("Sample data", style="padding-left:45px")),
             fluidRow(
               column(4, 
                      fluidRow(style="padding-left:45px", h5(".xls / .xlsx")),
                      fluidRow(style="padding-left:45px", downloadButton(outputId = "xlsfile", label = 'Download .xls/.xlsx Sample File'))
                      ),
               column(4, 
                      fluidRow(style="padding-left:45px", h5(".csv / .txt")),
                      fluidRow(style="padding-left:45px", downloadButton(outputId = "csvfile", label = 'Download .csv/.txt Sample File'))
               )
             ),
             fluidRow(h2("Tutorial", style="padding-left:45px; padding-top:25px")),
             fluidRow(style="padding-left:45px", '1. Upload Mean Cq File. The first column must have the sample names, and the last column must have the group/condition name. Middle columns, genes meanCq values. After Uploading the table, confirm the data upload by clicking on the "Confirm Data Table" button.'),
             fluidRow(tags$img( style="padding-top:15px", src = "Tutorial_01.png", width = "900px"),),
             fluidRow(style="padding-left:45px; padding-top:25px", '2. Data Summary will show each gene`s mean Cq by group. Select which are the target genes (non-reference genes) by clicking on the check box, and then confirm the selection by clicking on the "Update Target Genes(s)" button.'),
             fluidRow(tags$img( style="padding-top:15px", src = "Tutorial_02.png", width = "900px"),),
             fluidRow(style="padding-left:45px; padding-top:25px", '3. Gene Reference Samples will show how the samples of each reference gene is distributed among the groups. At this step, users may consider remove outliers, if present.'),
             fluidRow(style="padding-left:75px; padding-top:10px", '3.1. The "Set outlier standard deviation" indicates the threshold for a sample to be considered an outlier (default at two standard deviations); samples away from this threshold in each group of each gene will be regarded as outlier.'),
             fluidRow(style="padding-left:75px; padding-top:10px", '3.2. The user can choose which outlier to show and remove by clicking on the radio button "Choose which outlier to remove." "Only Mean" will show and remove only outliers in the mean of the reference genes. "All Outliers" will show and remove outliers in each gene individually. This second option leans towards removing more outliers and, consequently, more samples of the analysis.'),
             fluidRow(style="padding-left:75px; padding-top:10px", '3.3. Removing outliers is an interactive process since it decreases the group standard deviation and may reveal additional outliers.'),
             fluidRow(style="padding-left:75px; padding-top:10px", '3.4. The user may restore all removed outliers by clicking on the "Restore Outliers" button.'),
             fluidRow(tags$img( style="padding-top:15px", src = "Tutorial_03.png", width = "900px"),),
             fluidRow(style="padding-left:45px; padding-top:25px", '4. The Gene Reference Analysis step aims to identify unsuitable reference genes for the analysis. If a reference gene displays a high standard deviation, low stability, or significant alterations across groups/conditions, it is advisable to consider excluding it from the analysis.'),
             fluidRow(style="padding-left:75px; padding-top:10px", '4.1. The first row of the table displays p-values from the Kruskal-Wallis test, indicating whether a gene varies among the conditions. Subsequent rows present p-values from pairwise Wilcoxon-Mann-Whitney tests between conditions. P-values less than 0.05 are highlighted in red.'),
             fluidRow(style="padding-left:75px; padding-top:10px", '4.2. The final column displays p-values for the mean expression levels of the genes. In the differential analysis, the mean of reference genes will be employed to estimate gene expression, and this mean should remain consistent across conditions.'),
             fluidRow(tags$img( style="padding-top:15px", src = "Tutorial_04.png", width = "900px"),),
             fluidRow(tags$img( style="padding-top:15px", src = "Tutorial_05.png", width = "900px"),),
             fluidRow(style="padding-left:45px; padding-top:25px", '5. Differential Expression Analysis'),
             fluidRow(style="padding-left:75px; padding-top:10px", '5.1. Select target genes and groups/conditions to test.'),
             fluidRow(style="padding-left:75px; padding-top:10px", '5.2. Gene fold-change refers to the ratio of the expression level of a gene under one condition or treatment compared to its expression level under another condition or treatment. To estimate the expression level of a gene 2^-ΔCq is used.'),
             fluidRow(style="padding-left:75px; padding-top:10px", '5.3. Normality tests are made with the Shapiro-Wilks test using ΔCq data.'),
             fluidRow(style="padding-left:75px; padding-top:10px", '5.4. Differential Analyses are performed using ΔCq data in two different ways: Parametric and a non-parametric tests. If there are only two groups/conditions, the Pearson t-test and Wilcoxon-Mann-Whitney test are performed. If there are more than two groups/conditions, ANOVA and Tukey HSD (for multiple comparisons) tests, and Kruskal-Wallis and Dunn (for multiple comparisons using the Benjamini-Hochberg adjustment method) tests are performed.'),
             fluidRow(tags$img( style="padding-top:15px; padding-bottom:50px", src = "Tutorial_06.png", width = "1100px"),)
    ),
    tabPanel("About",
             fluidRow(style="padding-top:25px",
             column(6, style="padding-top:15px", 
                    p(style="padding-left:35px", HTML("EndoGeneAnalyzer is an open-source web tool developed by researchers at the Oncology Research Center (NPO) of the Federal University of Pará for the identification of the optimal normalization gene among a set of candidates in RT-qPCR studies.<br/><br/>
Additionally, the EndoGeneAnalyzer enables its users to perform differential gene expression analysis.")),
                    fluidRow(style="padding-top:25px;padding-left:35px", h4("CITATION")),
                    p(style="padding-left:35px", HTML("Teixeira EB, Khayat AS, Assumpção PP, Casseb SM, Moreira-Nunes CA, Moreira FC. (2024) EndoGeneAnalyzer: A tool for selection and validation of reference genes. PLOS ONE 19(4): e0299993. https://doi.org/10.1371/journal.pone.0299993.")),
                    fluidRow(style="padding-top:25px;padding-left:50px", "Original Article:", tags$a("EndoGeneAnalyzer: A tool for selection and validation of reference genes", href = "https://journals.plos.org/plosone/article/citation?id=10.1371/journal.pone.0299993")),
                    fluidRow(style="padding-top:25px;padding-left:35px", h4("SUPORT")),
                    fluidRow(style="padding-top:25px;padding-left:50px", "The source code is available at:", tags$a("https://github.com/MoreiraFC/EndoGeneAnalyzer", href = "https://github.com/MoreiraFC/EndoGeneAnalyzer")),
                    fluidRow(style="padding-top:25px;padding-left:50px", "For any further questions, comments or feedback, we are happy to provide email support at", tags$a("npo.ufpa@gmail.com", href = "mailto:npo.ufpa@gmail.com")),
                    ),
              column(6,
                     fluidRow(h4(align = "center", "PARTICIPATING INSTITUTIONS AND PROGRAMS")),
                     fluidRow(align = "center",
                      tags$img( style="padding-top:15px", src = "ufpa.png", width = "75px"),
                      tags$img( style="padding-top:15px; padding-left:20px", src = "npo.png", width = "130px"),
                      tags$img( style="padding-top:15px; padding-left:20px", src = "UFC.png", width = "75px"),
                      tags$img( style="padding-top:15px; padding-left:20px", src = "npdmufc_logo.jpeg", width = "130px"),
                      
                     ),
                     fluidRow(align = "center",
                      tags$img( style="padding-top:25px", src = "PPGOCM.png", width = "170px"),
                      tags$img( style="padding-top:25px; padding-left:20px", src = "PPGBIOTEC.jpeg", width = "270px")
                     ),
                     fluidRow(style="padding-top:25px;padding-left:35px", h4("FUNDING")),
                     fluidRow(align = "center",
                     tags$img( style="padding-top:25px", src = "CNPq.png", width = "150px"),
                     tags$img( style="padding-top:25px; padding-left:20px", src = "Capes.png", width = "100px"),
                     tags$img( style="padding-top:25px; padding-left:20px", src = "FUNCAP.png", width = "130px")
                     )
                  )
             )
          )
  )        
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  ### observables 
  observe({
    updateCheckboxGroupInput(
      session = session, 
      inputId = "Sel.target", 
      choices = unique(data$tabela$Genes), 
      selected = unique(data$tabela$Genes[data$tabela$target == T])
    )
    
    updateSelectInput(
      session = session, 
      inputId = "deg.target", 
      choices = unique(data$tabela$Genes[data$tabela$target == TRUE]))
    
    updateCheckboxGroupInput(
      session = session, 
      inputId = "deg.endo", 
      choices = unique(data$tabela$Genes[data$tabela$target == FALSE]),
      selected = data$endogenos)
    
    updateCheckboxGroupInput(
      session = session, 
      inputId = "deg.group", 
      choices = unique(data$tabela$group),
      selected = unique(data$tabela$group))
    
  })
  
  
  ### Data upload
  
  file = reactive({
    req(input$uploadFile)
    
    inFile = input$uploadFile
    if (input$type == ".csv/.txt"){
      df = read.csv(inFile$datapath, header = T, sep = input$sep)
      } else {
        df = read_excel(inFile$datapath)
      }
    
    colnames(df)[c(1,ncol(df))] = c("sample", "group")
    
    return(df)
  })
  
  ### download sample files
  output$csvfile = downloadHandler(
    filename = 'sampledata.txt',
    content = function(file) { 
      df = read_excel("sampledata.xlsx")
      write.table(df, file, sep = "\t", col.names = T, row.names = FALSE, quote = F)
    })
  
  output$xlsfile = downloadHandler(
    filename = 'sampledata.xlsx',
    content = function(file) { 
      df = read_excel("sampledata.xlsx")
      write_xlsx(df, path = file)
    })
  
  
  output$upload.table = renderTable(file())
  
  output$confirm.table = renderTable({
    df = NULL
    mensagem = ''
    if (is.null(file())) return()
    #{
    #  if (!is.null(data$tabela)) 
    #    mensagem = "Example data loaded."
    #  else 
    #}
    else {
      
      if (ncol(file()) < 3) {
        mensagem = "Table has 2 comuns or fewer. If it is a .csv file, try to adjust separator.\n"
      }
      if (nrow(file()) < 3) {
        mensagem = "Table has 2 rows or fewer. Check the file format.\n"
      }
      
      if (mensagem == ''){
        flag = apply(X = file()[,2:(ncol(file())-1)], MARGIN = c(1,2), FUN = is.character)
        if (sum(flag) > 0) mensagem = "There is text at genes MeanCq columns. Check the table format or decimal separator"
      } 
      
      if (mensagem == ''){
        df = data.frame ("Number of Samples" = nrow(file()), 
                        "Number of Genes" = as.integer(ncol(file())-2),
                        "Number of Groups" = as.integer(length(unlist(unique(file()[,ncol(file())]))))
        )
      }
      else df = data.frame(Message = mensagem)
      
    }
  }, bordered = T)
  
  data = reactiveValues(tabela = NULL, original = NULL, endogenos = NULL)
  
  observeEvent(input$confirm, {
    
    data$original = gather(file(), value = "MeanCq", key = "Genes", -sample, -group)
    
    data$original$target = FALSE
    data$original$outlier = FALSE
    
    data$endogenos = unique(data$original$Genes)
    data$original$MeanCq = as.double(data$original$MeanCq)

    data$tabela = data$original
  })
  
    
 # })
  
  observeEvent(input$load, {
    
    #data$original = read.delim2("endogenos_final.txt", header = T)
    data$original = read.delim2("endogenos_completo.2.tsv", header = T)
    data$original$MeanCq = as.numeric(data$original$MeanCq)
    data$original$target = FALSE
    data$original$outlier = FALSE
    data$endogenos = unique(data$original$Genes)
    
    data$tabela = data$original
    
  })
  
  ## Data Summary
  output$group.plot = renderPlot({
    if (is.null(data$tabela)) return()
    set.seed(1)
    ggplot(data$tabela, mapping = aes(x = group, y=MeanCq )) +
      geom_boxplot(aes(fill=group)) +
      geom_jitter(aes(col=MeanCq), width = 0.3, cex = 0.5) +
      facet_wrap(~Genes)+
      ggtitle("Genes by groups")
  })
  
  observeEvent (input$update, {
    data$tabela$target = FALSE
    selected = input$Sel.target
    #data$original$target[data$original$Genes %in% selected] = TRUE
    data$tabela$target[data$tabela$Genes %in% selected] = TRUE
    #data$tabela = data$original
  })
  
  output$data_table = renderTable({
    data$tabela}, digits = 3, hover = T, striped = T) 
  
  ### Endogenous Analysis
  
  ### tabela kkw
  endo.values = reactive({
    ref.table = filter(data$tabela, Genes %in% input$deg.endo)
    
    ref.mean = ref.table %>% group_by(sample, group) %>% summarise(MeanCq = mean(MeanCq))
    ref.mean$Genes = "MeanRef"
    ref.mean$target = FALSE
    ref.mean$outlier = FALSE
    ref.table$outlier = FALSE
    #data$tabela$outlier = FALSE
    
    if (input$mean.outlier == 'All Outliers') {
      ref.table$outlier = calc.outlier(tabela = ref.table, limiar = input$out.dp)$outlier
      data$tabela$outlier = calc.outlier(tabela = data$tabela, limiar = input$out.dp)$outlier
    }
    
    ref.mean$outlier = calc.outlier(tabela = ref.mean, limiar = input$out.dp)$outlier
    sample.out = filter(ref.mean, outlier == T)$sample
    
    endo.tab = bind_rows(ref.table, ref.mean)

    k = data$tabela$sample %in% sample.out
    print(paste(sum(data$original$outlier == data$tabela$outlier), nrow(data$original)))
    data$tabela$outlier[k] = TRUE
    
    genes = unique(endo.tab$Genes)
    grupos = unique(endo.tab$group)
    dp = NULL
    meanSD = NULL
    sum.diff.mean.square = NULL
    sum.diff.sd.square = NULL
    endo.pv = NULL
    pv = NULL
    teste = NULL
    for (i in genes){
      tab = filter(endo.tab, Genes == i)
      dp[i] = sd(tab$MeanCq, na.rm = T)
      a = group_by(tab,group) %>%
        summarise(sd = sd(MeanCq, na.rm = T),
                  mean = mean(MeanCq, na.rm = T))
      a$diff.mean.gene = a$mean - mean(tab$MeanCq, na.rm = T)
      a$diff.sd.gene = a$sd - sd(tab$MeanCq, na.rm = T)
      meanSD[i] = mean(a$sd, na.rm = T)
      
      sum.diff.mean.square[i] = sum((a$mean - mean(tab$MeanCq, na.rm = T)) ^ 2, na.rm = T)
      sum.diff.sd.square[i] = sum((a$sd - sd(tab$MeanCq, na.rm = T)) ^ 2, na.rm = T)
      
      
      if (length(grupos) > 2){
        pv = NULL
        
        teste = "Kruskal-Wallis"
        kkw = kruskal.test(MeanCq ~ group, tab)
        a = dunn.test(x = tab$MeanCq, g = tab$group, kw = T, method = "BH")
        pv = as.data.frame(c(kkw$p.value, a$P.adjusted))
        row.names(pv) = c(teste, a$comparisons)
        colnames(pv) = i
        
        if (is.null(endo.pv)) {endo.pv = pv}
        else {endo.pv = bind_cols(endo.pv, pv)}
        
      } else if (length(grupos) == 2) {
        teste = "Wilcoxon-Mann-Whitney"
        wmw = wilcox.test(MeanCq ~ group, tab)
        pv[i] = wmw$p.value
        endo.pv = data.frame(t(pv))
        row.names(endo.pv) = "Wilcoxon-Mann-Whitney"
      }
    }
    endo.DT = endo.pv %>% datatable(options = list(
      dom = 'Bfrtip',
      lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
      pageLength = 15) ) %>% 
      formatRound(columns = genes, digits = 3) %>%
      formatStyle(
        columns = genes, 
        color = styleInterval(0.05, c("red", "black"))) %>%
      formatStyle(0, target = 'row', fontWeight = styleEqual(teste, 'bold'))
    
    return(list(endo.DT = endo.DT,
                endo.tab = endo.tab,
                dp = dp,
                meanSD = meanSD,
                sum.diff.sd.square = sum.diff.sd.square,
                sum.diff.mean.square = sum.diff.mean.square
    ))
  })
  
  output$endo.p.values = DT::renderDataTable({
    if (is.null(data$tabela)) return()
    endo.values()$endo.DT})
  
  ### tabela descritiva
  output$descriptive = renderTable({
    if (is.null(data$tabela)) return()
    t.data.frame(data.frame(
      Standard.Deviation=endo.values()$dp, #mean.SD=endo.values()$meanSD, 
      sum.SD.square.diff=endo.values()$sum.diff.sd.square,
      sum.mean.square.diff=endo.values()$sum.diff.mean.square
      ))
  }, striped = T, hover = T, rownames = T)
  
  
  ### tabela Norm.Finder
  nf.res = reactive({
    nf = dcast(filter(data$tabela, target == FALSE), formula = sample~Genes, value.var = "MeanCq")
    grupos = data$tabela$group
    names(grupos) = data$tabela$sample
    nf$group = grupos[nf$sample]
    
    nf = as.data.frame(t(nf))
    colnames(nf)=nf[1,]
    nf = nf[2:nrow(nf),]
    
    k = complete.cases(t(nf))
    nf = nf[,k]
    
    res = Normfinder(nf)
    return(res)
  })

  output$nf.order = renderTable({
    if (is.null(data$tabela)) return()
    #nf.res()}, rownames = T)
    nf.res()$Ordered}, striped = T, rownames = T, hover = T)
  output$nf.unorder = renderTable({
    if (is.null(data$tabela)) return()
    nf.res()$UnOrdered}, striped = T, rownames = T, hover = T)
  
  ## endogenous sample stability
  
  observeEvent (input$set.outlier.sd, {
    data$tabela$outlier = FALSE
    data$tabela$outlier = calc.outlier(data$tabela, input$out.dp)$outlier
    data$original = data$tabela
  })

  output$sample.plot = renderPlot({
    if (is.null(data$tabela)) return()
    #tabela = data$tabela
    #tab = filter(tabela, target == FALSE)
    tabela = filter(endo.values()$endo.tab, Genes != "MeanRef")
    tabela$group = factor(tabela$group, levels = sort(unique(tabela$group)))
    tabela = arrange(tabela, group, sample)
    tabela$sample = factor(tabela$sample, levels = unique(tabela$sample))
    tabela$MeanCq = as.numeric(tabela$MeanCq)
    
    ggplot(tabela) + 
      geom_point(mapping = aes(x = sample, y = MeanCq, col = outlier, shape = group), 
                 size = 2) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))+ ylim(0, 40) +
      facet_wrap(~Genes, nrow = length(unique(tabela$Genes)))
  })
  
  output$mean.plot = renderPlot({
    if (is.null(data$tabela)) return()

    tabela = filter(endo.values()$endo.tab, Genes == "MeanRef")
    tabela$group = factor(tabela$group, levels = sort(unique(tabela$group)))
    tabela = arrange(tabela, group, sample)
    tabela$sample = factor(tabela$sample, levels = unique(tabela$sample))
    tabela$MeanCq = as.numeric(tabela$MeanCq)
    
    ggplot(tabela) + 
      geom_point(mapping = aes(x = sample, y = MeanCq, col = outlier, shape = group), 
                 size = 2) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) + 
      ylim(0, 40)+ ggtitle("Mean of Selected Reference Genes")
  })
  
  observeEvent (input$remove, {
    k = data$tabela$outlier == TRUE & data$tabela$target == FALSE
    data$tabela = data$tabela[!k,]
    data$endogenos = input$deg.endo
    #    data$tabela = filter(data$tabela, outlier == FALSE)
    # data$tabela = data$original
  })
  
  observeEvent (input$restore, {
    data$endogenos = input$deg.endo
    target = unique(data$tabela$Genes[data$tabela$target == T])
    data$tabela = data$original
    data$tabela$target[data$tabela$Genes %in% target] = TRUE
  })
  
  output$outlier.tab = renderTable({
    if (is.null(data$tabela)) return()
    filter(endo.values()$endo.tab, outlier == TRUE, target == FALSE) %>% arrange(Genes)}, striped = T, hover = T, digits = 3)
  
  ### Differential Analysis
  de.tab = reactive({
    tab = filter(data$tabela, 
                 Genes %in% c(input$deg.target, input$deg.endo),
                 group %in% input$deg.group)
    tab.d = dcast(tab, sample~Genes, value.var = "MeanCq")
    row.names(tab.d) = tab.d$sample
    tab.d = tab.d[,2:ncol(tab.d)]
    
    #tab.d = tab.d[complete.cases(tab.d),]
    
    if (length(input$deg.endo) > 1){
      tab.d$mean = apply(tab.d[, input$deg.endo], 1, mean, na.rm=T)
    } else { tab.d$mean = tab.d[, input$deg.endo]}
    tab.d$deltaCq = tab.d[, input$deg.target] - tab.d$mean
    tab.d$twoDeltaCq = 2^(-tab.d$deltaCq)
    grupo = tab$group
    names(grupo) = tab$sample
    tab.d$group = grupo[row.names(tab.d)]
    return(tab.d)
  })
  
  
  output$param.test = renderPrint({
    if (input$deg.target == "") {
      print("Select a target gene in Data Summary Tab")
    } else if (length(input$deg.group) == 1) {
      print("Select two or more groups")
    } else if (length(input$deg.group) == 2) {
      t.test(formula = deltaCq ~ group, data = de.tab())
    } else if (length(input$deg.group) > 2) {
      model = aov(formula = deltaCq ~ group, data = de.tab())
      print("ANOVA Results:")
      print(summary(model))
      print(TukeyHSD(model))
    }
  })
  
  output$nparam.test = renderPrint({
    if (input$deg.target == "") {
      print("Select a target gene in Data Summary Tab")
    } else if (length(input$deg.group) == 1) {
      print("Select two or more groups")
    } else if (length(input$deg.group) == 2) {
      wilcox.test (formula = deltaCq ~ group, data = de.tab())
    } else if (length(input$deg.group) > 2) {
      dunn.test(x = de.tab()$deltaCq, g = de.tab()$group, kw = T, method = "BH")
    }
  })
  
  output$shapiro.test = renderTable({
    if (is.null(data$tabela)) return()
    df = data.frame()
    if (input$deg.target == "") {
      df = data.frame(message = "Select a target gene in Data Summary Tab")
    } else{
      
      tab = filter (de.tab(), group %in% input$deg.group)
      df = as.data.frame(tab %>% group_by(group) %>% summarise(count=n()))
      df$shapiro.p.value = 0
      for (i in 1:nrow(df)){
        a = filter(tab, group == df$group[i])
        df$shapiro.p.value[i] = shapiro.test(a$deltaCq)$p.value
        
      }
    }
    
    return(df)
    
  }, striped = T, hover = T)
  

  plot.deg.boxplot = function() {
    if (is.null(data$tabela)) return()
    
    cores = c("Set1", "Set2", "Set3", "Blues", "Greys", "Accent", "PRGn", "RdBu")
    names(cores) = c("Set1", "Set2", "Set3", "Blues", "Greys", "Accent", "Purple to Green", "Red to Blue")
    
    
    if (input$deg.target == "") {
      ggplot(data$tabela, mapping = aes(x = Genes, y=MeanCq)) +
        geom_boxplot(aes(fill=Genes)) +
        geom_jitter(aes(col=MeanCq), width = 0.3, cex = 1) +
        scale_fill_brewer(palette = cores[input$plot.color])
      
    } else {
      
      if (input$plot.type == "Boxplot") {
        
        ggplot(de.tab(), mapping = aes(x = group, y=-deltaCq)) +
          geom_boxplot(aes(fill=group))+
          scale_fill_brewer(palette = cores[input$plot.color])
        
      } else if(input$plot.type == "Scatter plot") {
        set.seed(1)
        ggplot(de.tab(), mapping = aes(x = group, y=-deltaCq)) +
          geom_jitter(aes(col=group), width = 0.3, cex = 1.5)+
          scale_fill_brewer(palette = cores[input$plot.color])
        
      } else if(input$plot.type == "Violin") {
        
        ggplot(de.tab(), mapping = aes(x = group, y=-deltaCq)) +
          geom_violin (aes(fill=group))+
          scale_fill_brewer(palette = cores[input$plot.color])
        
      } else {
        set.seed(1)
        ggplot(de.tab(), mapping = aes(x = group, y=-deltaCq)) +
          geom_boxplot(aes(fill=group)) +
          geom_jitter(aes(col=-deltaCq), width = 0.3, cex = 1)+
          scale_fill_brewer(palette = cores[input$plot.color])
      }
    }
  }
  
  output$deg.boxplot = renderPlot (plot.deg.boxplot())
  
  output$download = downloadHandler(
    filename = function() { paste(input$dataset, 'Figure.boxplot.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot.deg.boxplot(), device = "pdf")
    })
  
  output$deg.FC = renderTable({
    if (is.null(data$tabela)) return()
    df = data.frame()
    if (input$deg.target == "") {
      df = data.frame(message = "Select a target gene in Data Summary Tab")
    } else if (length(input$deg.group) == 1) {
      df = data.frame(message = "Select two or more groups")
    } else if (length(input$deg.group) > 1) {
      
      grupos = input$deg.group
      cont = 0
      for (i in 1:(length(grupos)-1)){
        for (j in (i+1):length(grupos)){
          meanA = mean(filter(de.tab(), group == grupos[i])$twoDeltaCq, na.rm = T)
          meanB = mean(filter(de.tab(), group == grupos[j])$twoDeltaCq, na.rm = T)
          cont=cont+1
          df[cont,1:2] = c(meanA/meanB, meanB/meanA)
          row.names(df)[cont] = paste(grupos[i], grupos[j], sep="-")
        }
      }
      colnames(df) = c("FC", "1/FC")
    }
    return(df)
  }, rownames = T, striped = T)
  
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "DEA_report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "DEA_report.Rmd")
      file.copy("DEA_report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(n = de.tab())
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  output$deg.table = renderTable({
    if (is.null(data$tabela)) return()
    if (input$deg.target == "") return()
    a = de.tab()
    a = rename(a, 
               !!paste("\u0394", "Cq", sep = "") := deltaCq,
               !!paste(c("2", "^", "-\u0394", "Cq"), collapse = "") := twoDeltaCq)
    
    return(a)
  }, striped = T, hover = T, rownames = T, digits = 4)
  
  #  })
}

# Run the application 
shinyApp(ui = ui, server = server)
