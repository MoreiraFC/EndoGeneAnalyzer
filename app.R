
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
      dp = sd(a$MeanCt, na.rm = T)
      media = mean(a$MeanCt, na.rm = T)
      out.samples = as.character(a$sample[a$MeanCt > (media+limiar*dp) | a$MeanCt < (media-limiar*dp)])
      k = tabela$Genes == i & tabela$group ==j & tabela$sample %in% out.samples
      res$outlier[k] = TRUE
      res$dados[cont,1:4] = c(i, j, dp, media)
      cont=cont+1
    }
  }
  return(res)
}

#user_base <- tibble::tibble(
#  user = c("npo", "fabiano", "alunoAndre"),
#  password = sapply(c("GCnp0!", "sn00py", "endogeno"), sodium::password_store),
#  permissions = c("admin", "standard", "standard"),
#  name = c("User One", "User Two", "User Three")
#)

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

#  headerPanel(windowTitle = "EndoGeneAnalyzer", 
#              title = div(
#                class = "pull-right", tags$img(height = 70, width = 70, src = "logonpo2.png", style='padding-top:10px'),
#                tags$img(height = 70, width = 70, src = "logo_ufpa.png", style='padding-left:10px'))
#              ),
#    titlePanel(
#      title = div(
#      class = "pull-right", tags$img(height = 70, width = 70, src = "logonpo2.png", style='padding-top:10px'),
#      tags$img(height = 70, width = 70, src = "logo_ufpa.png", style='padding-left:10px'),
      
      #class = "pull-right", shinyauthr::logoutUI(id = "logout")
 #   ), "Nome da Ferramenta"
#  ),
  
  # login section
#  shinyauthr::loginUI(id = "login", title = "", user_title = "Usuário", pass_title = "Senha"),
  
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
             column(2, actionButton(inputId = "load", label = "Load Example Data"), style="padding-top:5px"),
             ),
             fluidRow(column(5, tableOutput("confirm.table"), style="padding-left:10px"),
                      column(3, actionButton(inputId = "confirm", label = "Confirm Data Table"), style="padding-left:10px")
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
                                  min = 1, max = 4, value = 2, step = 0.25)),
                      #column(2, style="padding-top:10px", actionButton("set.outlier.sd", "Update outlier SD")),
                      column(2, style="padding-top:10px",
                             fluidRow(actionButton("remove", "Remove Outliers")),
                             fluidRow(
                             "Obs.: Removing outliers decreases the group standard deviation. 
                                           This action may make visible other outliers"
                             ),
                      ),
                      column(3,radioButtons('mean.outlier', 'Choose which outlier to remove',
                                            choices = 
                                              c('Only Mean', 'All Outliers'), 'Only Mean')
                      ),
                      ),
             fluidRow(column(12, style="padding-top:15px", plotOutput(outputId = "sample.plot", height = "600px"))),
             fluidRow(column(12, plotOutput(outputId = "mean.plot", height = "200px"))),
             fluidRow(column(9, style="padding-left:45px", h4("List of outliers reference samples"))),
             fluidRow(column(9, style="padding-left:45px", tableOutput(outputId = "outlier.tab")))
    ),
    tabPanel("Gene Reference Analysis",
             fluidRow(column(9, style="padding-left:25px", h4("  Gene Reference by group"))),
             fluidRow(
               column(12, DT::dataTableOutput("endo.p.values"))
             ),
             fluidRow(column(9, style="padding-left:25px", h4("Gene Reference Descriptive Statistics"))),
             fluidRow(
               column(8, tableOutput("descriptive")),
               column(3, checkboxGroupInput(inputId = "deg.endo", 
                                            label = "Select one or more Genes:", 
                                            choices = NULL)),
             ),
             fluidRow(column(9, style="padding-left:25px",h4("Normfinder Analysis"))),
             fluidRow(column(9, style="padding-left:25px", h5("  Normfinder Result Ordered"))),
             fluidRow(
               column(10, tableOutput(outputId = "nf.order"))
             ),
             fluidRow(column(9, style="padding-left:25px", h5("  Normfinder Result Unordered by group"))),
             fluidRow(column(12, tableOutput(outputId = "nf.unorder")))
    ), ## fim do tabPanel Endogenos
    
    tabPanel("Differential Analysis",
             fluidRow(column(2, selectInput(inputId = "deg.target", label = "Select a Target:", choices = NULL)),
                      column(3, checkboxGroupInput(inputId = "deg.group", label = "Select two or more groups:", choices = NULL)),
             ),
             fluidRow(
                      column(6, style="padding-left:45px", plotOutput(outputId = "deg.boxplot"),
                             downloadButton(outputId = "download", label = 'Download Figure')),
                      
                      column(4, style="padding-left:45px",
                             tableOutput(outputId = "deg.FC"),
                             fluidRow(h5("Normality test:"),
                                      tableOutput(outputId  = "shapiro.test"))
                             )
             ),
             fluidRow(h1("")),
             fluidRow(
                      column(6, verbatimTextOutput(outputId = "param.test")),
                      column(6, verbatimTextOutput(outputId = "nparam.test")),
             ),
             mainPanel(tableOutput("deg.table"))
            
    )
  )        
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
#  credentials <- shinyauthr::loginServer(
#    id = "login",
#    data = user_base,
#    user_col = user,
#    pwd_col = password,
#    sodium_hashed = TRUE,
#    log_out = reactive(logout_init())
#  )
  
  # Logout to hide
#  logout_init <- shinyauthr::logoutServer(
#    id = "logout",
#    active = reactive(credentials()$user_auth)
#  )
  
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
  
  output$upload.table = renderTable(file())
  
  output$confirm.table = renderTable({
    if (is.null(file())) return()
    else {
      df = NULL
      mensagem = ''
      if (ncol(file()) < 3) {
        mensagem = "Table has 2 comuns or fewer. If it is a .csv file, try to adjust separator.\n"
      }
      if (nrow(file()) < 3) {
        mensagem = "Table has 2 rows or fewer. Check the file format.\n"
      }
      
      if (mensagem == ''){
        flag = apply(X = file()[,2:(ncol(file())-1)], MARGIN = c(1,2), FUN = is.character)
        if (sum(flag) > 0) mensagem = "There is text at genes MeanCt columns. Check the table format or decimal separator"
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
    
    data$original = gather(file(), value = "MeanCt", key = "Genes", -sample, -group)
    
    data$original$target = FALSE
    data$original$outlier = FALSE
    
    data$endogenos = unique(data$original$Genes)
    data$original$MeanCt = as.double(data$original$MeanCt)

    data$tabela = data$original
  })
  
    
 # })
  
  observeEvent(input$load, {
    
    #data$original = read.delim2("endogenos_final.txt", header = T)
    data$original = read.delim2("endogenos_completo.2.txt", header = T)
    data$original$MeanCt = as.numeric(data$original$MeanCt)
    data$original$target = FALSE
    data$original$outlier = FALSE
    data$endogenos = unique(data$original$Genes)
    
    #data.table$MeanCt = as.double(data.table$MeanCt)
#    data.table$outlier = calc.outlier(data.table, limiar = input$out.dp)$outlier

    data$tabela = data$original
    })
  
  ## Data Summary
  output$group.plot = renderPlot({
    if (is.null(data$tabela)) return()
    set.seed(1)
    ggplot(data$tabela, mapping = aes(x = group, y=MeanCt )) +
      geom_boxplot(aes(fill=group)) +
      geom_jitter(aes(col=MeanCt), width = 0.3, cex = 0.5) +
      facet_wrap(~Genes)+
      ggtitle("Genes by groups")
  })
  
  observeEvent (input$update, {
    data$tabela$target = FALSE
    selected = input$Sel.target
    data$original$target[data$original$Genes %in% selected] = TRUE
    data$tabela = data$original
  })
  
  output$data_table = renderTable({
    data$tabela}, digits = 3, hover = T, striped = T) 
  
  ### Endogenous Analysis
  
  ### tabela kkw
  endo.values = reactive({
    b = filter(data$tabela, Genes %in% input$deg.endo)
    
    a = b %>% group_by(sample, group) %>% summarise(MeanCt = mean(MeanCt))
    a$Genes = "MeanRef"
    a$target = FALSE
    a$outlier = FALSE
    b$outlier = FALSE
    #data$tabela$outlier = FALSE
    
    if (input$mean.outlier == 'All Outliers') {
      b$outlier = calc.outlier(tabela = b, limiar = input$out.dp)$outlier
      data$tabela$outlier = calc.outlier(tabela = data$tabela, limiar = input$out.dp)$outlier
    }
    
    a$outlier = calc.outlier(tabela = a, limiar = input$out.dp)$outlier
    med.out = filter(a, outlier == T)$sample
    
    endo.tab = bind_rows(b, a)

    k = data$tabela$sample %in% med.out
    print(paste(sum(data$original$outlier == data$tabela$outlier), nrow(data$original)))
    data$original$outlier[k] = TRUE
    
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
      dp[i] = sd(tab$MeanCt)
      a = group_by(tab,group) %>%
        summarise(sd = sd(MeanCt),
                  mean = mean(MeanCt))
      a$diff.mean.gene = a$mean - mean(tab$MeanCt)
      a$diff.sd.gene = a$sd - sd(tab$MeanCt)
      meanSD[i] = mean(a$sd)
      
      sum.diff.mean.square[i] = sum((a$mean - mean(tab$MeanCt)) ^ 2)
      sum.diff.sd.square[i] = sum((a$sd - sd(tab$MeanCt)) ^ 2)
      
      
      if (length(grupos) > 2){
        pv = NULL
        
        teste = "Kruskal-Wallis"
        kkw = kruskal.test(MeanCt ~ group, tab)
        a = dunn.test(x = tab$MeanCt, g = tab$group, kw = T, method = "BH")
        pv = as.data.frame(c(kkw$p.value, a$P.adjusted))
        row.names(pv) = c(teste, a$comparisons)
        colnames(pv) = i
        
        if (is.null(endo.pv)) {endo.pv = pv}
        else {endo.pv = bind_cols(endo.pv, pv)}
        
      } else if (length(grupos) == 2) {
        teste = "Wilcoxon-Mann-Whitney"
        wmw = wilcox.test(MeanCt ~ group, tab)
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
    nf = dcast(filter(data$tabela, target == FALSE), formula = sample~Genes, value.var = "MeanCt")
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
    tabela$MeanCt = as.numeric(tabela$MeanCt)
    
    ggplot(tabela) + 
      geom_point(mapping = aes(x = sample, y = MeanCt, col = outlier, shape = group), 
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
    tabela$MeanCt = as.numeric(tabela$MeanCt)
    
    ggplot(tabela) + 
      geom_point(mapping = aes(x = sample, y = MeanCt, col = outlier, shape = group), 
                 size = 2) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) + 
      ylim(0, 40)+ ggtitle("Mean of Selected Gene Reference")
  })
  
  observeEvent (input$remove, {
    k = data$original$outlier == TRUE & data$original$target == FALSE
    data$original = data$original[!k,]
    data$endogenos = input$deg.endo
    #    data$tabela = filter(data$tabela, outlier == FALSE)
    data$tabela = data$original
  })
  
  output$outlier.tab = renderTable({
    if (is.null(data$tabela)) return()
    filter(data$tabela, outlier == TRUE, target == FALSE) %>% arrange(Genes)}, striped = T, hover = T, digits = 3)
  
  ### Differential Analysis
  de.tab = reactive({
    tab = filter(data$tabela, 
                 Genes %in% c(input$deg.target, input$deg.endo),
                 group %in% input$deg.group)
    tab.d = dcast(tab, sample~Genes, value.var = "MeanCt")
    row.names(tab.d) = tab.d$sample
    tab.d = tab.d[,2:ncol(tab.d)]
    
    #tab.d = tab.d[complete.cases(tab.d),]
    
    if (length(input$deg.endo) > 1){
      tab.d$mean = apply(tab.d[, input$deg.endo], 1, mean, na.rm=T)
    } else { tab.d$mean = tab.d[, input$deg.endo]}
    tab.d$deltaCt = tab.d[, input$deg.target] - tab.d$mean
    tab.d$twoDeltaCt = 2^(-tab.d$deltaCt)
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
      t.test(formula = deltaCt ~ group, data = de.tab())
    } else if (length(input$deg.group) > 2) {
      model = aov(formula = deltaCt ~ group, data = de.tab())
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
      wilcox.test (formula = deltaCt ~ group, data = de.tab())
    } else if (length(input$deg.group) > 2) {
      dunn.test(x = de.tab()$deltaCt, g = de.tab()$group, kw = T, method = "BH")
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
        df$shapiro.p.value[i] = shapiro.test(a$deltaCt)$p.value
        
      }
    }
    
    return(df)
    
  }, striped = T, hover = T)
  
  
  
  
  
  plot.deg.boxplot = function() {
    if (is.null(data$tabela)) return()
    
    if (input$deg.target == "") {
      ggplot(data$tabela, mapping = aes(x = Genes, y=MeanCt)) +
        geom_boxplot(aes(fill=Genes)) +
        geom_jitter(aes(col=MeanCt), width = 0.3, cex = 1) 
      } else {
    
      set.seed(1)
      ggplot(de.tab(), mapping = aes(x = group, y=-deltaCt)) +
        geom_boxplot(aes(fill=group)) +
        geom_jitter(aes(col=-deltaCt), width = 0.3, cex = 1)
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
          meanA = mean(filter(de.tab(), group == grupos[i])$twoDeltaCt, na.rm = T)
          meanB = mean(filter(de.tab(), group == grupos[j])$twoDeltaCt, na.rm = T)
          cont=cont+1
          df[cont,1:2] = c(meanA/meanB, meanB/meanA)
          row.names(df)[cont] = paste(grupos[i], grupos[j], sep="-")
        }
      }
      colnames(df) = c("FC", "1/FC")
    }
    return(df)
  }, rownames = T, striped = T)
  
  output$deg.table = renderTable({
    if (is.null(data$tabela)) return()
    if (input$deg.target == "") return()
    a = de.tab()
    a = rename(a, 
               !!paste("\u0394", "Ct", sep = "") := deltaCt,
               !!paste(c("2", "^", "-\u0394", "Ct"), collapse = "") := twoDeltaCt)
    
    return(a)
  }, striped = T, hover = T, rownames = T, digits = 4)
  
  #  })
}

# Run the application 
shinyApp(ui = ui, server = server)
