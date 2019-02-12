#####################
# CHLAMYNET WEB APP #
#####################

# Copyright (C) 2018  Francisco J. Romero-Campero
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public
# License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Authors: Francisco J. Romero-Campero
# 
# Contact: Francisco J. Romero-Campero - fran@us.es 
# Date: February 2019

## Load libraries
library(shiny)
library(gplots)
library(ggplot2)
library(igraph)
library(topGO)

## Data loading

## Load network
network.data <- read.table(file="data/chlamynet_data.tsv",header = TRUE,as.is=TRUE,fill=TRUE)
head(network.data)

## Adjust y positions
network.data$y <- - network.data$y

## Extract gene ids
genes <- sort(network.data$name)

## Extract Pfam structure
pfam <- unique(network.data$pfam)

## Read network in gml format
chlamynet <- read.graph(file = "data/chlamynet.gml",format = "gml")
nodes <- network.data$name

## Reading in data frame with annotation
annotation.data <- read.table(file="data/Creinhardtii_223_annotation_info.txt",sep="\t")

## Load gene expression levels in FPKM
fpkm.data.read <- read.table(file="data/fpkm_data.txt",header=TRUE)

## Setup for GO enrichment based on PFAM annotation
## Load Chlamydomonas reinhardtii gene to GO annotation
geneID2GO <- readMappings(file = "data/Chlamy_PFAM_based_GO_annotation.map")

## Set the background gene set to the whole gene set in Chlamy 
## We create a vector with all elements 1 and name it with the gene 
## identifiers 
gene.names <- attributes(geneID2GO)[[1]]
gene.background.pfam <- rep(1,length(gene.names))
names(gene.background.pfam) <- gene.names

## The function used to selec genes is defined according to the assignment of 0 to the genes
## of interest
cre.gene.selec <- function(gene.list)
{
  return(gene.list == 0)
}

## The function used to select genes is defined according to the assignment of 0 to the genes
## of interest
ath.gene.selec <- function(gene.list)
{
  return(gene.list == 0)
}

ontology.type <- "BP"

## Load A. thaliana orthology
chlamy.athaliana <- read.table(file="data/Creinhardtii_athaliana.txt")
chlamy.names <- as.vector(chlamy.athaliana[[1]])
athaliana.names <- as.vector(chlamy.athaliana[[2]])
names(athaliana.names) <- chlamy.names

## Loading the probe names and ATG identifiers and their correspondence
probe.atg <- read.table(file="data/probe_names_atg.txt")
probe.names <- as.vector(probe.atg[[1]])
atg.names <- as.vector(probe.atg[[2]])
names(probe.names) <- atg.names

###################
## CHLAMYNET APP ##
###################

ui <- fluidPage(
  ## Application title
  titlePanel(tags$b("ChlamyNET, a ", tags$i("Chlamydomonas reinhardtii"), " Gene Co-expression Network ")),
  
  # Extra line breaks
  tags$br(),
  tags$br(),
  
  ## Introduction to the tool 
  fluidRow(
    column(2,
           tags$img(src="ChlamyNet_Logo.jpg",alt="ChlamyNET_logo", width=160,height=160,align="middle")
    ),
    
    column(3,
           tags$ul(style="list-style-type:none",
                   tags$li(tags$a(href="http://viridiplantae.ibvf.csic.es/ChlamyNet/index.html", 
                                  tags$b("Home"))),

                   tags$li(tags$a(href="http://viridiplantae.ibvf.csic.es/ChlamyNet/clusters.html",
                                  tags$b("ChlamyNet Exploration"))),
                   
                   tags$li(tags$a(href="http://viridiplantae.ibvf.csic.es/ChlamyNet/clusters.html",
                     tags$b("Clusters GO Term Enrichment"))),
                   
                   tags$li(tags$a(href="http://viridiplantae.ibvf.csic.es/ChlamyNet/transcription_factors.html",
                     tags$b("Transcription Factors and Regulators"))),
                   
                   tags$li(tags$a(href="http://viridiplantae.ibvf.csic.es/ChlamyNet/tutorial.html",
                     tags$b("Tutorial"))),
                   
                   tags$li(tags$a(href="http://viridiplantae.ibvf.csic.es/ChlamyNet/case_studies.html",
                     tags$b("Case Studies"))),
                   
                   tags$li(tags$a(href="http://viridiplantae.ibvf.csic.es/ChlamyNet/links.html",
                     tags$b("Related Links")))
           )),

    column(7,
           p(tags$b("ChlamyNET 2.0 "), "is a web-based tool developed using ", 
             tags$b(tags$a(href="https://shiny.rstudio.com/","shiny,")), " an R package for the development 
     of web apps from R code. The aim of ChlamyNET is to facilitate the studies over the 
     Chlamydomonas transcriptome. Please, take a look at the ", 
     tags$b(tags$a(href="http://viridiplantae.ibvf.csic.es/ChlamyNet/tutorial.html","tutorial")),  " or ", 
     tags$b(tags$a(href="https://www.youtube.com/watch?v=QuySUPid-rg","watch this video tutorial.")), 
     "A user can search for a set of genes of interest using the ", tags$b("Gene Selection panel."), " The analysis 
     of individual genes or sets of genes including their neighbouring genes can be performed using 
     the ", tags$b("Gene Expression"), " and ", tags$b("GO Enrichment Analysis"), " panels. Please, take a look at the ",
     tags$b(tags$a(href="http://viridiplantae.ibvf.csic.es/ChlamyNet/case_studies.html","case studies")),
     " for some examples on how to use ChlamyNET 2.0.") 
     )
  ),
  
  ## Extra line breaks for separation
  tags$br(),
  tags$br(),
  
  ## Sidebar layout is chosen. Panels for analysis will be in the sidebar and in the
  ## main panel a visualization of the network and results of the analysis (tables and graphs)
  ## will be displayed. 
  sidebarLayout(
    sidebarPanel(
  
      ## ChlamyNET contains FOUR different panels for analysis
      
      ## PANEL 1: GENE SELECTION ##
      wellPanel(
        
        ## Panel title
        tags$h3(tags$b("Gene Selection:")),
        
        ## Radio buttons to chose between two modes of gene selection
        ## using gene ID (CreXX.gXXXXXX) or PFAM ID (PFXXXXX)
        radioButtons(inputId = "gene_selection_mode",
                    label = "Mode", 
                     choices = c("Gene ID", "Protein Domain PFAM ID","Gene List"),
                     selected = "Gene ID"),
        
        ## Dynamic panel for selecting genes based on their ID
        conditionalPanel(condition = "input.gene_selection_mode == 'Gene ID'",
                         ## Select gene ID
                         selectizeInput(inputId = "selected.gene",
                                        label = "Gene ID",
                                        choices = genes,
                                        selected = "Cre08.g370400",
                                        multiple = TRUE),
                         ## Select distance of genes to consider
                         selectInput(inputId = "distance", 
                                     label = "Select Co-expressed Genes at Distance:", 
                                     choices = 0:3,
                                     selected = 0,
                                     multiple = FALSE),
                         ## Button to trigger selections based on gene ID
                         actionButton(inputId = "button_gene_id",label="Select Genes")),
        
        conditionalPanel(condition = "input.gene_selection_mode == 'Gene List'",
                         ## Enter gene list
                         textAreaInput(inputId = "gene.list", label= "Set of genes", width="90%", 
                                       height = "200px",placeholder = "Insert set of genes",
                                       value = "Cre01.g011150
Cre04.g224600
Cre07.g332250
Cre14.g620850"),
                         ## Select distance of genes to consider
                         selectInput(inputId = "distance_gene_list", 
                                     label = "Select Co-expressed Genes at Distance:", 
                                     choices = 0:3,
                                     selected = 0,
                                     multiple = FALSE),
                         ## Button to trigger selections based on gene ID
                         actionButton(inputId = "button_gene_list_id",label="Select Genes")),
        
        
        ## Dynamic panel for selecting genes based on their PFAM domain IDs
        conditionalPanel(condition = "input.gene_selection_mode == 'Protein Domain PFAM ID'",
                         ## Select PFAM ID
                         selectizeInput(inputId = "selected_pfam",
                                        label = "PFAM Protein Domain ID",
                                        choices = pfam,
                                        selected = "PF00643,PF06203",#"PF00010",
                                        multiple = TRUE),
                         ## Action button to trigger identification of genes with the selected
                         ## PFAM ID
                         actionButton(inputId = "go_pfam",label="Search Genes")),
        conditionalPanel(condition = "output.number_pfam_genes > 0 && input.gene_selection_mode == 'Protein Domain PFAM ID'",
                         tags$br(), # extra line break for separation
                         ## Check box to select genes with the given PFAM ID
                         checkboxGroupInput(inputId = "pfam_genes",
                                            label = "Select Genes:"),
                         ## Select distance 
                         selectInput(inputId = "distance_pfam", 
                                     label = "Select Co-expressed Genes at Distance:", 
                                     choices = 0:3,
                                     selected = 0,
                                     multiple = FALSE),
                         ## Action button to trigger gene selection based on PFAM ID
                         actionButton(inputId = "select_neighbours_pfam",label="Select Genes")
                           
                        )
      ), #end panel for gene selection

      
      ## PANEL 2: GENE EXPRESSION ANALYSIS ##
      wellPanel(
        
        # Title, explanatory text and extra line breaks for separation
        tags$h3(tags$b("Gene Expression Analysis:")),
        tags$b("Generate a line graph of the expression levels of the selected genes:"),
        tags$br(),
        tags$br(),
        
        # Action button for line graph
        actionButton(inputId = "button_linegraph",label="Line Graph"),
        
        # Extra line breaks for separation and explanatory text
        tags$br(),
        tags$br(),
        tags$b("Generate a heatmap of the correlation between the expression levels of the selected genes:"),
        tags$br(),
        tags$br(),
        
        # Action button for heatmap
        actionButton(inputId = "button_heatmap",label="Heatmap")
        
      ), #end panel for expression analysis
      
      ## PANEL 3: GO ENRICHMENT ANALYSIS ##
      wellPanel(
        
        # Title
        tags$h3(tags$b("GO Enrichment Analysis:")),
        
        # GO enrichment mode
        radioButtons(inputId = "go_mode",
                     label = "GO Enrichment",
                     choices = c("Based on Arabidopsis orthology","Based on PFAM Annotation")),
      
        # Action button to trigger GO enrichment analysis
        actionButton(inputId = "go_enrichment",label = "GO Enrichment")
        
      ), #end panel for GO enrichment

      width = 4
    ),
    
    # In the main panel the network, output graphs and tables are displayed. 
    mainPanel(
      plotOutput("networkPlot", height = "800px"), #click="plot_click",
      uiOutput(outputId = "message"),
      uiOutput(outputId = "download"),
      dataTableOutput(outputId = "output_table"),
      plotOutput("output_graph"),

      width = 8
    )
  ),
  
  
  
  tableOutput("gene_info_table"),
  tableOutput("go_table"),
  plotOutput("heatmapPlot"),
  plotOutput("linePlot"),
  plotOutput("go_enrichment_plot")
  
)

## Server function for ChlamyNET
server <- function(input, output, session) {

  ## Initial/default visualization of ChlamyNET network
  output$networkPlot <- renderPlot({
    ggplot(network.data, aes(x,y)) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks.y = element_blank()) + 
      geom_point(color=network.data$color,size=1)
  })

  ## Reactive responding to gene selection by ID when clicking button_gene_id
  selected_gene_id <- eventReactive(input$button_gene_id,{
      target.gene <- input$selected.gene
      
      selection <- c()
      for(i in 1:length(target.gene))
      {
        selection <- c(selection,
                       nodes[as.vector(ego(graph = chlamynet,nodes = target.gene[i],order=input$distance)[[1]])])
      }
      selection <- unique(selection)
      #print(selection)
      #subset(network.data, name %in% nodes[as.vector(ego(graph = chlamynet,nodes = target.gene,order=input$distance)[[1]])])
      subset(network.data, name %in% selection)
  })

  ## Reactive responding to gene list selection by ID when clicking button_gene_list_id
  selected_gene_list_id <- eventReactive(input$button_gene_list_id,{
    target.gene.list <- as.vector(unlist(
      strsplit(input$gene.list, split="\n",
               fixed = TRUE)[1]))
      
    selection <- c()
    genes.not.chlamynet <- setdiff(x = target.gene.list, y = network.data$name)
    target.gene.list <- intersect(target.gene.list,network.data$name)
    
    for(i in 1:length(target.gene.list))
    {
      selection <- c(selection,
                     nodes[as.vector(ego(graph = chlamynet,nodes = target.gene.list[i],order=input$distance_gene_list)[[1]])])
    }
    selection <- unique(selection)
    #print(paste("distance:",input$distance_gene_list))
    #print(selection)
    #subset(network.data, name %in% nodes[as.vector(ego(graph = chlamynet,nodes = target.gene,order=input$distance)[[1]])])
    subset(network.data, name %in% selection)
  })
  
  
  ## Gene selection by PFAM
  selected_gene_pfam <- eventReactive(input$go_pfam,{
    ## Extract x and y position for the selected gene
    pfam_selection <- subset(network.data, pfam %in% input$selected_pfam)
    pfam_genes_selected <<- pfam_selection$name
    #print("selection:")
    #print(pfam_genes_selected)
    pfam_selection
  })
  
  ## Reactive to determine selected genes with a given PFAM ID
  my_pfam_genes <- reactive({
    pfam_selection <- subset(network.data, pfam %in% input$selected_pfam)
    pfam_genes_selected <<- pfam_selection$name
    return(pfam_genes_selected)
  })
  
  ## Observe to update the list of genes with a given PFAM ID
  observe({
    updateCheckboxGroupInput(session, "pfam_genes",
                      choices = my_pfam_genes()
    )
  })
  
  ## Determination of the number of genes with a given PFAM ID to be
  ## used in a dynamic conditional panel
  output$number_pfam_genes <- eventReactive(input$go_pfam,{
    ## Extract x and y position for the selected gene
    pfam_selection <- subset(network.data, pfam == input$selected_pfam)
    pfam_genes_selected <<- pfam_selection$name
    #print("selection:")
    ##print(pfam_genes_selected)
    length(pfam_genes_selected)
  })
  
  ## Necessary for dinamical panel
  outputOptions(output, "number_pfam_genes", suspendWhenHidden = FALSE)

  ## Visualization of selected genes by ID
  observeEvent(input$button_gene_id, {
        #print("aquí")
        #print(selected_gene_id())
        gene.list <<- selected_gene_id()$name
        output$networkPlot <- renderPlot({
          ggplot(network.data, aes(x,y)) + 
            theme(panel.background = element_blank(), 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks.y = element_blank()) + 
            geom_point(color=network.data$color,size=1) + 
            geom_point(data=selected_gene_id(), aes(x,y), size=3, fill=selected_gene_id()$color,colour="black",pch=21) 
        })
        
        ## Generate annotation table
        genes.annotation.data <- subset(annotation.data, V1 %in% gene.list)
        
        colnames(genes.annotation.data) <- c("Gene id", "PFAM", "PANTHER", "KOG", "EC", "K", "Arabidopsis", "Common name", "Description")
        genes.annotation.data <- genes.annotation.data[c("Gene id", "Description", "Common name", "Arabidopsis", "PFAM", "PANTHER", "KOG", "EC", "K")]
        
        output$message <- renderUI({
          tagList(
            tags$p(tags$b("The table below presents the annotation of the selected genes.
                          This information can be downloaded by clicking on the following button:")),
            tags$br())
        })
        
        output$download <- renderUI({
          tagList(
            downloadButton(outputId= "downloadData", "Get Selected Genes Annotation"),
            tags$br(),
            tags$br()
          )
        })
        
        output$downloadData <- downloadHandler(
          filename = function() 
          {
              return('selected_genes_annotation.txt')
          },
          content = function(file) {
            write.table(as.data.frame(genes.annotation.data), file=file, sep="\t",quote=FALSE )
        })
        
        
        ## Construct data frame with links
        
        genes.ids <- as.vector(genes.annotation.data[["Gene id"]])
        genes.pfam <- as.vector(genes.annotation.data[["PFAM"]])
        genes.kog <- as.vector(genes.annotation.data[["KOG"]])
        genes.ec <- as.vector(genes.annotation.data[["EC"]])
        genes.k <- as.vector(genes.annotation.data[["K"]])
        genes.pthr <- as.vector(genes.annotation.data[["PANTHER"]])
        
        genes.with.links <- vector(mode="character",length=length(genes.ids))
        genes.pfam.with.links <- vector(mode="character",length=length(genes.pfam))
        genes.kog.with.links <- vector(mode="character",length=length(genes.kog))
        genes.ec.with.links <-  vector(mode="character",length=length(genes.ec))
        genes.k.with.links <- vector(mode="character",length=length(genes.ec))
        genes.pthr.with.links <- vector(mode="character",length=length(genes.pthr))
        
        for(i in 1:length(genes.ids))
        {
          current.gene <- genes.ids[i]

          current.phytozome.link <- paste0("https://phytozome.jgi.doe.gov/pz/portal.html#!results?search=0&crown=1&star=1&method=4614&searchText=",
                                           current.gene,
                                           "&offset=0")

          phytozome.href <- paste(c("<a href=\"",
                                    current.phytozome.link,
                                    "\" target=\"_blank\">Phytozome</a>"),
                                  collapse="")

          current.circadianet.link <- paste0("http://viridiplantae.ibvf.csic.es/circadiaNet/genes/cre/",
                                             current.gene,
                                             ".html")


          circadianet.href <- paste(c("<a href=\"",
                                      current.circadianet.link,
                                      "\" target=\"_blank\">CircadiaNET</a>"),
                                    collapse="")

          new.gene.id.element <- paste(c(current.gene,
                                         paste(phytozome.href,
                                         circadianet.href,sep=",")),collapse = " ")
          
          current.pfam <- genes.pfam[i]

          if(current.pfam != "")
          {
            pfams.ids <- strsplit(current.pfam,split=",")[[1]]
            pfam.href <- vector(mode="character",length=length(pfams.ids))
            
            for(j in 1:length(pfams.ids))
            {
              pfam.link <- paste0("https://pfam.xfam.org/family/",pfams.ids[j])
              
              pfam.href[j] <- paste(c("<a href=\"",
                                        pfam.link,
                                        "\" target=\"_blank\">",
                                        pfams.ids[j], "</a>"),
                                      collapse="")
            }
            pfams.hrefs <- paste(pfam.href,collapse=",")
          } else
          {
            pfams.hrefs <- ""
          }
          
          
          ## Create kog links
          current.kog <- genes.kog[i]
          
          if(current.kog != "")
          {
              kog.link <- paste0("http://eggnogdb.embl.de/#/app/results?seqid=Q6CPW9&target_nogs=",
                                 current.kog)
              kog.href <- paste(c("<a href=\"",
                                      kog.link,
                                      "\" target=\"_blank\">",
                                      current.kog, "</a>"),
                                    collapse="")
            } else
          {
            kog.href <- ""
          }
          
          ## Create ec link
          current.ec <- genes.ec[i]
          
          if(current.ec != "")
          {
            ec.link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ec:",
                               current.ec)
            ec.href <- paste(c("<a href=\"",
                                ec.link,
                                "\" target=\"_blank\">",
                                current.ec, "</a>"),
                              collapse="")
          } else
          {
            ec.href <- ""
          }
          
          ## Create k link
          current.k <- genes.k[i]
          
          if(current.k != "")
          {
            k.link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ko:",
                              current.k)
            k.href <- paste(c("<a href=\"",
                               k.link,
                               "\" target=\"_blank\">",
                               current.k, "</a>"),
                             collapse="")
          } else
          {
            k.href <- ""
          }
          
          ## Create pthr link
          current.pthr <- genes.pthr[i]
          
          if(current.pthr != "")
          {
            pthr.link <- paste0("http://www.pantherdb.org/panther/family.do?clsAccession=",
                             current.pthr)
            pthr.href <- paste(c("<a href=\"",
                              pthr.link,
                              "\" target=\"_blank\">",
                              current.pthr, "</a>"),
                            collapse="")
          } else
          {
            pthr.href <- ""
          }
          
          genes.with.links[i] <- new.gene.id.element
          genes.pfam.with.links[i] <- pfams.hrefs
          genes.kog.with.links[i] <- kog.href
          genes.ec.with.links[i] <- ec.href
          genes.k.with.links[i] <- k.href
          genes.pthr.with.links[i] <- pthr.href
        } 
        
        genes.annotation.data.with.links <- 
          data.frame(genes.with.links,
                     genes.annotation.data[,2:4],
                     genes.pfam.with.links,
                     genes.pthr.with.links,
                     genes.kog.with.links,
                     genes.ec.with.links,
                     genes.k.with.links)

        colnames(genes.annotation.data.with.links) <- colnames(genes.annotation.data)
        
        output$output_table <- renderDataTable({
            as.data.frame(genes.annotation.data.with.links)
          },escape=FALSE)
        
    })

## ---------------------------------------
  
  ## Visualization of selected genes by ID
  observeEvent(input$button_gene_list_id, {
    #print("aquí con lista")
    #print(selected_gene_id())
    gene.list <<- selected_gene_list_id()$name
    output$networkPlot <- renderPlot({
      ggplot(network.data, aes(x,y)) + 
        theme(panel.background = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks.y = element_blank()) + 
        geom_point(color=network.data$color,size=1) + 
        geom_point(data=selected_gene_list_id(), aes(x,y), size=3, fill=selected_gene_list_id()$color,colour="black",pch=21) 
    })
    
    ## Generate annotation table
    genes.annotation.data <- subset(annotation.data, V1 %in% gene.list)
    
    colnames(genes.annotation.data) <- c("Gene id", "PFAM", "PANTHER", "KOG", "EC", "K", "Arabidopsis", "Common name", "Description")
    genes.annotation.data <- genes.annotation.data[c("Gene id", "Description", "Common name", "Arabidopsis", "PFAM", "PANTHER", "KOG", "EC", "K")]

    target.gene.list <- as.vector(unlist(
      strsplit(input$gene.list, split="\n",
               fixed = TRUE)[1]))
    
    genes.not.chlamynet <- setdiff(x = target.gene.list, y = network.data$name)
    
        
    if(length(genes.not.chlamynet) > 0)
    {
      msg <- paste(c("The following genes are not in ChlamyNET: ", genes.not.chlamynet, 
                     ". Check if the IDs are correct and correspond to Chlamydomonas reinhardtii annotation v5. Alternatively, these genes may not be differentially expressed in the conditions integrated in ChlamyNET."),collapse=" ")
    } else
    {
      msg <- ""
    }
    
    output$message <- renderUI({
      tagList(
        tags$p(tags$b(msg)), 
        tags$p(tags$b("The table below presents the annotation of the selected genes.
                      This information can be downloaded by clicking on the following button:")),
        tags$br())
    })
    
    output$download <- renderUI({
      tagList(
        downloadButton(outputId= "downloadData", "Get Selected Genes Annotation"),
        tags$br(),
        tags$br()
      )
    })
    
    output$downloadData <- downloadHandler(
      filename = function() 
      {
        return('selected_genes_annotation.txt')
      },
      content = function(file) {
        write.table(as.data.frame(genes.annotation.data), file=file, sep="\t",quote=FALSE )
      })
    
    
    ## Construct data frame with links
    
    genes.ids <- as.vector(genes.annotation.data[["Gene id"]])
    genes.pfam <- as.vector(genes.annotation.data[["PFAM"]])
    genes.kog <- as.vector(genes.annotation.data[["KOG"]])
    genes.ec <- as.vector(genes.annotation.data[["EC"]])
    genes.k <- as.vector(genes.annotation.data[["K"]])
    genes.pthr <- as.vector(genes.annotation.data[["PANTHER"]])
    
    genes.with.links <- vector(mode="character",length=length(genes.ids))
    genes.pfam.with.links <- vector(mode="character",length=length(genes.pfam))
    genes.kog.with.links <- vector(mode="character",length=length(genes.kog))
    genes.ec.with.links <-  vector(mode="character",length=length(genes.ec))
    genes.k.with.links <- vector(mode="character",length=length(genes.ec))
    genes.pthr.with.links <- vector(mode="character",length=length(genes.pthr))
    
    for(i in 1:length(genes.ids))
    {
      current.gene <- genes.ids[i]
      
      current.phytozome.link <- paste0("https://phytozome.jgi.doe.gov/pz/portal.html#!results?search=0&crown=1&star=1&method=4614&searchText=",
                                       current.gene,
                                       "&offset=0")
      
      phytozome.href <- paste(c("<a href=\"",
                                current.phytozome.link,
                                "\" target=\"_blank\">Phytozome</a>"),
                              collapse="")
      
      current.circadianet.link <- paste0("http://viridiplantae.ibvf.csic.es/circadiaNet/genes/cre/",
                                         current.gene,
                                         ".html")
      
      
      circadianet.href <- paste(c("<a href=\"",
                                  current.circadianet.link,
                                  "\" target=\"_blank\">CircadiaNET</a>"),
                                collapse="")
      
      new.gene.id.element <- paste(c(current.gene,
                                     paste(phytozome.href,
                                           circadianet.href,sep=",")),collapse = " ")
      
      current.pfam <- genes.pfam[i]
      
      if(current.pfam != "")
      {
        pfams.ids <- strsplit(current.pfam,split=",")[[1]]
        pfam.href <- vector(mode="character",length=length(pfams.ids))
        
        for(j in 1:length(pfams.ids))
        {
          pfam.link <- paste0("https://pfam.xfam.org/family/",pfams.ids[j])
          
          pfam.href[j] <- paste(c("<a href=\"",
                                  pfam.link,
                                  "\" target=\"_blank\">",
                                  pfams.ids[j], "</a>"),
                                collapse="")
        }
        pfams.hrefs <- paste(pfam.href,collapse=",")
      } else
      {
        pfams.hrefs <- ""
      }
      
      
      ## Create kog links
      current.kog <- genes.kog[i]
      
      if(current.kog != "")
      {
        kog.link <- paste0("http://eggnogdb.embl.de/#/app/results?seqid=Q6CPW9&target_nogs=",
                           current.kog)
        kog.href <- paste(c("<a href=\"",
                            kog.link,
                            "\" target=\"_blank\">",
                            current.kog, "</a>"),
                          collapse="")
      } else
      {
        kog.href <- ""
      }
      
      ## Create ec link
      current.ec <- genes.ec[i]
      
      if(current.ec != "")
      {
        ec.link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ec:",
                          current.ec)
        ec.href <- paste(c("<a href=\"",
                           ec.link,
                           "\" target=\"_blank\">",
                           current.ec, "</a>"),
                         collapse="")
      } else
      {
        ec.href <- ""
      }
      
      ## Create k link
      current.k <- genes.k[i]
      
      if(current.k != "")
      {
        k.link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ko:",
                         current.k)
        k.href <- paste(c("<a href=\"",
                          k.link,
                          "\" target=\"_blank\">",
                          current.k, "</a>"),
                        collapse="")
      } else
      {
        k.href <- ""
      }
      
      ## Create pthr link
      current.pthr <- genes.pthr[i]
      
      if(current.pthr != "")
      {
        pthr.link <- paste0("http://www.pantherdb.org/panther/family.do?clsAccession=",
                            current.pthr)
        pthr.href <- paste(c("<a href=\"",
                             pthr.link,
                             "\" target=\"_blank\">",
                             current.pthr, "</a>"),
                           collapse="")
      } else
      {
        pthr.href <- ""
      }
      
      genes.with.links[i] <- new.gene.id.element
      genes.pfam.with.links[i] <- pfams.hrefs
      genes.kog.with.links[i] <- kog.href
      genes.ec.with.links[i] <- ec.href
      genes.k.with.links[i] <- k.href
      genes.pthr.with.links[i] <- pthr.href
    } 
    
    genes.annotation.data.with.links <- 
      data.frame(genes.with.links,
                 genes.annotation.data[,2:4],
                 genes.pfam.with.links,
                 genes.pthr.with.links,
                 genes.kog.with.links,
                 genes.ec.with.links,
                 genes.k.with.links)
    
    colnames(genes.annotation.data.with.links) <- colnames(genes.annotation.data)
    
    output$output_table <- renderDataTable({
      as.data.frame(genes.annotation.data.with.links)
    },escape=FALSE)
    
  })
  
  
  
  
  
## ---------------------------------------
  
  
  
  
  ## Visualization of identified genes with the given PFAM domains
  observeEvent(input$go_pfam, {
      #print("pfam")
      #print(selected_gene_pfam())
      gene.list <<- selected_gene_pfam()$name
      output$networkPlot <- renderPlot({
        ggplot(network.data, aes(x,y)) + 
          theme(panel.background = element_blank(), 
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks.y = element_blank()) + 
          geom_point(color=network.data$color,size=1) + 
          geom_point(data=selected_gene_pfam(), aes(x,y), size=3, fill=selected_gene_pfam()$color,colour="black",pch=21) 
      })
    })
  
    ## Visualization of selected genes by ID
    observeEvent(eventExpr = input$select_neighbours_pfam, {
      
      ## Determine the neighbours of the identified genes with the given PFAM
      found.pfam.genes <- input$pfam_genes
      neighbours.pfam <- c()
      for(i in 1:length(found.pfam.genes))
      {
        neighbours.pfam <- rbind(neighbours.pfam,
                                 subset(network.data, name %in% nodes[as.vector(ego(graph = chlamynet,nodes = found.pfam.genes[i] ,order=input$distance_pfam)[[1]])]))
      }
      
      ## Set gene list for further analysis
      gene.list <<- neighbours.pfam$name
      
      ## Network visualization
      output$networkPlot <- renderPlot({
        ggplot(network.data, aes(x,y)) + 
          theme(panel.background = element_blank(), 
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks.y = element_blank()) + 
          geom_point(color=network.data$color,size=1) + 
          geom_point(data=neighbours.pfam, aes(x,y), size=3, fill=neighbours.pfam$color,colour="black",pch=21) 
      })
      
      ## Generate annotation table
      genes.annotation.data <- subset(annotation.data, V1 %in% gene.list)
      
      colnames(genes.annotation.data) <- c("Gene id", "PFAM", "PANTHER", "KOG", "EC", "K", "Arabidopsis", "Common name", "Description")
      genes.annotation.data <- genes.annotation.data[c("Gene id", "Description", "Common name", "Arabidopsis", "PFAM", "PANTHER", "KOG", "EC", "K")]
      
      ## Construct data frame with links
      
      genes.ids <- as.vector(genes.annotation.data[["Gene id"]])
      genes.pfam <- as.vector(genes.annotation.data[["PFAM"]])
      genes.kog <- as.vector(genes.annotation.data[["KOG"]])
      genes.ec <- as.vector(genes.annotation.data[["EC"]])
      genes.k <- as.vector(genes.annotation.data[["K"]])
      genes.pthr <- as.vector(genes.annotation.data[["PANTHER"]])
      
      genes.with.links <- vector(mode="character",length=length(genes.ids))
      genes.pfam.with.links <- vector(mode="character",length=length(genes.pfam))
      genes.kog.with.links <- vector(mode="character",length=length(genes.kog))
      genes.ec.with.links <-  vector(mode="character",length=length(genes.ec))
      genes.k.with.links <- vector(mode="character",length=length(genes.ec))
      genes.pthr.with.links <- vector(mode="character",length=length(genes.pthr))
      
      for(i in 1:length(genes.ids))
      {
        current.gene <- genes.ids[i]
        
        current.phytozome.link <- paste0("https://phytozome.jgi.doe.gov/pz/portal.html#!results?search=0&crown=1&star=1&method=4614&searchText=",
                                         current.gene,
                                         "&offset=0")
        
        phytozome.href <- paste(c("<a href=\"",
                                  current.phytozome.link,
                                  "\" target=\"_blank\">Phytozome</a>"),
                                collapse="")
        
        current.circadianet.link <- paste0("http://viridiplantae.ibvf.csic.es/circadiaNet/genes/cre/",
                                           current.gene,
                                           ".html")
        
        
        circadianet.href <- paste(c("<a href=\"",
                                    current.circadianet.link,
                                    "\" target=\"_blank\">CircadiaNET</a>"),
                                  collapse="")
        
        new.gene.id.element <- paste(c(current.gene,
                                       paste(phytozome.href,
                                             circadianet.href,sep=",")),collapse = " ")
        
        current.pfam <- genes.pfam[i]
        
        if(current.pfam != "")
        {
          pfams.ids <- strsplit(current.pfam,split=",")[[1]]
          pfam.href <- vector(mode="character",length=length(pfams.ids))
          
          for(j in 1:length(pfams.ids))
          {
            pfam.link <- paste0("https://pfam.xfam.org/family/",pfams.ids[j])
            
            pfam.href[j] <- paste(c("<a href=\"",
                                    pfam.link,
                                    "\" target=\"_blank\">",
                                    pfams.ids[j], "</a>"),
                                  collapse="")
          }
          pfams.hrefs <- paste(pfam.href,collapse=",")
        } else
        {
          pfams.hrefs <- ""
        }
        
        
        ## Create kog links
        current.kog <- genes.kog[i]
        
        if(current.kog != "")
        {
          kog.link <- paste0("http://eggnogdb.embl.de/#/app/results?seqid=Q6CPW9&target_nogs=",
                             current.kog)
          kog.href <- paste(c("<a href=\"",
                              kog.link,
                              "\" target=\"_blank\">",
                              current.kog, "</a>"),
                            collapse="")
        } else
        {
          kog.href <- ""
        }
        
        ## Create ec link
        current.ec <- genes.ec[i]
        
        if(current.ec != "")
        {
          ec.link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ec:",
                            current.ec)
          ec.href <- paste(c("<a href=\"",
                             ec.link,
                             "\" target=\"_blank\">",
                             current.ec, "</a>"),
                           collapse="")
        } else
        {
          ec.href <- ""
        }
        
        ## Create k link
        current.k <- genes.k[i]
        
        if(current.k != "")
        {
          k.link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ko:",
                           current.k)
          k.href <- paste(c("<a href=\"",
                            k.link,
                            "\" target=\"_blank\">",
                            current.k, "</a>"),
                          collapse="")
        } else
        {
          k.href <- ""
        }
        
        ## Create pthr link
        current.pthr <- genes.pthr[i]
        
        if(current.pthr != "")
        {
          pthr.link <- paste0("http://www.pantherdb.org/panther/family.do?clsAccession=",
                              current.pthr)
          pthr.href <- paste(c("<a href=\"",
                               pthr.link,
                               "\" target=\"_blank\">",
                               current.pthr, "</a>"),
                             collapse="")
        } else
        {
          pthr.href <- ""
        }
        
        genes.with.links[i] <- new.gene.id.element
        genes.pfam.with.links[i] <- pfams.hrefs
        genes.kog.with.links[i] <- kog.href
        genes.ec.with.links[i] <- ec.href
        genes.k.with.links[i] <- k.href
        genes.pthr.with.links[i] <- pthr.href
      } 
      
      genes.annotation.data.with.links <- 
        data.frame(genes.with.links,
                   genes.annotation.data[,2:4],
                   genes.pfam.with.links,
                   genes.pthr.with.links,
                   genes.kog.with.links,
                   genes.ec.with.links,
                   genes.k.with.links)
      
      colnames(genes.annotation.data.with.links) <- colnames(genes.annotation.data)
      

      output$message <- renderUI({
        tagList(
          tags$p(tags$b("The table below presents the annotation of the selected genes.
                          This information can be downloaded by clicking on the following button:")),
          tags$br())
      })

      output$download <- renderUI({
        tagList(
          downloadButton(outputId= "downloadData", "Get Selected Genes Annotation"),
          tags$br(),
          tags$br()
        )
      })
      
      output$downloadData<- downloadHandler(
        filename = 'selected_genes_annotation.txt',
        content = function(file) {
          write.table(as.data.frame(genes.annotation.data), file=filename, sep="\t",quote=FALSE )
        })
      
      output$output_table <- renderDataTable({
        as.data.frame(genes.annotation.data.with.links)
      },escape=FALSE)
      
    })
  
    
    ## Generate Heatmap
    observeEvent(eventExpr = input$button_heatmap, {
      #print("Generate heatmap")
      #gene.list <- selected.gene.pos()$name
      #print("Selected genes:")
      #print(gene.list)
      ## Creating a matrix that stores the expression level of the genes of interest
      number.of.genes <- length(gene.list)
      
      genes.data <- matrix(nrow=number.of.genes,ncol=20)
      
      for(i in 1:number.of.genes)
      {
        ## We analyse 20 different conditions
        for(j in 1:20)
        {
          genes.data[i,j] <- fpkm.data.read[fpkm.data.read[["gene_id"]] == gene.list[i],j+1]
        }
      }
      
      ## Naming the rows
      rownames(genes.data) <- gene.list
      
      ## Compute the gene correlation (note that is necessary to traspose the matrix)
      gene.correlation <- cor(t(genes.data))
      
      output$message <- renderUI({
        tagList(
          tags$p(tags$b("The heatmap and table below presents the correlation between the expression
                         profiles of the selected genes. This information can be downloaded by 
                         clicking on the corresponding buttons:")),
          tags$br())
      })
      
      output$download <- renderUI({
        tagList(
          downloadButton(outputId= "downloadCorrelation", "Get Selected Genes Correlation"),
          tags$br(),
          tags$br()
        )
      })
      
      ## Heatmap figure
      output$output_graph <- renderPlot({
        heatmap.2(gene.correlation,dendrogram="row",labCol=c(""),density.info="none",trace="none",col=greenred(100),margins = c(10,10),cexRow=1)
      },width = 10*96,height = 10*96)
      
      ## Heatmap table
      gene.correlation.table <- data.frame(colnames(gene.correlation),gene.correlation)
      colnames(gene.correlation.table) <- c("Gene ID",colnames(gene.correlation))
      
      output$output_table <- renderDataTable({
        gene.correlation.table
      },escape=TRUE)
      
    })
    
    ## Generate line graph 
    observeEvent(eventExpr = input$button_linegraph,{
      #print("line graph")
      #gene.list <- selected_gene_pfam()$name
      #print("Selected pfam genes")
      #print(gene.list)
      ## Create a matrix that stores the expression level of the genes of interest
      number.of.genes <- length(gene.list)
      
      genes.data <- matrix(nrow=number.of.genes,ncol=20)
      
      for(i in 1:number.of.genes)
      {
        ## We analyse 20 different conditions
        for(j in 1:20)
        {
          genes.data[i,j] <- fpkm.data.read[fpkm.data.read[["gene_id"]] == gene.list[i],j+1]
        }
      }
      
      rownames(genes.data) <- gene.list
      table.colnames <- colnames(fpkm.data.read)
      colnames(genes.data) <- table.colnames[2:21]
      
      genes.data.table <- data.frame(gene.list,genes.data)
      colnames(genes.data.table) <- c("Gene ID", colnames(genes.data))
      
      output$message <- renderUI({
        tagList(
          tags$p(tags$b("The line graph and table below presents the expression profiles 
                         of the selected genes. This information can be downloaded by 
                         clicking on the corresponding buttons:")),
          tags$br())
      })
      
      output$download <- renderUI({
        tagList(
          downloadButton(outputId= "downloadProfiles", "Get Selected Genes Profiles"),
          tags$br(),
          tags$br()
        )
      })
      
      output$output_table <- renderDataTable({
        genes.data.table
      },escape=TRUE)
      
      ## We use as many different colors as input genes.
      line.colors <- rainbow(number.of.genes)
      
      output$output_graph <- renderPlot({
        ## First line with all the details of the graph
        plot(genes.data[1,],type="l",col=line.colors[1],
             ylim=c(min(genes.data) - 0.1, max(genes.data) + 0.1),
             xlab="",ylab="",lwd=3,axes=FALSE) #,tick=FALSE,labels=FALSE
        
        ## The rest of the lines to be added to the graph started in the previous line
        if(nrow(genes.data) > 1)
        {
          for(k in 2:nrow(genes.data))
          {
            lines(genes.data[k,],type="l",col=line.colors[k],lwd=3)
          }
        }
        
        ## Add info on the axis
        axis(2,font=2)
        mtext(side=2,"Expression Level (FPKM)",cex=1.5,line=3)
        
        ## Vertical lines to divide the plot into different conditions
        lines(c(4,4),c(min(genes.data) - 0.1, max(genes.data) + 0.1),lty=2)
        lines(c(7,7),c(min(genes.data) - 0.1, max(genes.data) + 0.1),lty=2)
        lines(c(10,10),c(min(genes.data) - 0.1, max(genes.data) + 0.1),lty=2)
        lines(c(12,12),c(min(genes.data) - 0.1, max(genes.data) + 0.1),lty=2)
        lines(c(16,16),c(min(genes.data) - 0.1, max(genes.data) + 0.1),lty=2)
        lines(c(18,18),c(min(genes.data) - 0.1, max(genes.data) + 0.1),lty=2)
        
        ## Add studied conditions on the x axis
        axis(1,font=2,tick=FALSE,at=2.3,labels="Copper",las=2)
        axis(1,font=2,tick=FALSE,at=5.6,labels="Iron",las=2)
        axis(1,font=2,tick=FALSE,at=8.1,labels="Oxidative",las=2)
        axis(1,font=2,tick=FALSE,at=8.9,labels="Stress",las=2)
        axis(1,font=2,tick=FALSE,at=11,labels="Nitrogen",las=2)
        axis(1,font=2,tick=FALSE,at=14,labels="Sulfur",las=2)
        axis(1,font=2,tick=FALSE,at=17,labels="sor1",las=2)
        axis(1,font=2,tick=FALSE,at=19,labels="Minerals",las=2)
        
        ## Add legend
        legend("topright",legend=gene.list,col=line.colors,lwd=3)
      },height = 5*96) 
    })
    
    
    ## GO enrichment
    observeEvent(input$go_enrichment,{
      
      #print("go:")
      #print(input$go_mode)
      ## Selected genes
      cre.gene.names <- gene.list
      
      if(input$go_mode == "Based on PFAM Annotation")
      {
        #print("pfam")
        ## Set the background gene set to the whole gene set in Chlamy 
        ## We create a vector with all elements 1 and name it with the gene 
        ## identifiers 
        gene.names <- attributes(geneID2GO)[[1]]
        gene.background <- rep(1,length(gene.names))
        names(gene.background) <- gene.names
        
        ## In order to distinguish between background genes and genes of interest the value 0 is 
        ## associated to our genes of interest
        gene.background[cre.gene.names] <- 0
        
        ## Construct the topGOdata object
        sampleGOdata <- new("topGOdata",
                            description = "Chlamydomonas session", ontology = ontology.type,
                            allGenes = gene.background, geneSel = cre.gene.selec,
                            nodeSize = 10,
                            annot = annFUN.gene2GO, gene2GO = geneID2GO)
        
        resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
        
        ## Table containing the 20 most significant GO terms to be
        ## save as a text file. To avoid an error when no significant
        ## GO term is found a tryCatch sentence is used.
        
        ms <- 1
        number.top.nodes <- 20
        
        while(ms == 1 | ms == 2)
        {
          
          ms <- tryCatch(
            {
              allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                                 ranksOf = "classicFisher", topNodes = number.top.nodes, 
                                 numChar=100)
              ms <- 3
            },warning = function(war) {
              
              # warning handler picks up where error was generated
              print(paste("MY_WARNING:  ",war))
              print(number.top.nodes)
              return(1)
              
            }, error = function(err) {
              
              # error handler picks up where error was generated
              print(paste("MY_ERROR:  ",err))
              return(2)
              
            }, finally = {
              print(3)
              number.top.nodes <- number.top.nodes - 1
            }
          )
          
        }
        
        ## Generate for each GO term the Cre genes associated with it
        sig.go.terms <- allRes[["GO.ID"]]
        
        go.to.gene <- read.table(file="data/Chlamy_GO_to_gene_annotation.txt")
        go.term.names <- as.vector(go.to.gene[[1]])
        genes.go.terms <- as.vector(go.to.gene[[2]])
        names(genes.go.terms) <- go.term.names
        
        genes.with.GO <- vector(length=length(sig.go.terms))
        
        for(i in 1:length(sig.go.terms))
        {
          cre.go.terms <- strsplit(genes.go.terms[sig.go.terms[i]],split=",")[[1]]
          cre.go.names <- intersect(cre.go.terms,cre.gene.names)
          
          chlamy.go.names.string <- paste(cre.go.names,collapse=" ")
          genes.with.GO[i] <- chlamy.go.names.string
        }
        
        allRes[["Chlamy.GO"]] <- genes.with.GO
        
        ## Extracting only the relevant columns
        relevant.res <- allRes[,c("GO.ID","Term","classicFisher","Chlamy.GO")]
        colnames(relevant.res) <- c("GO term", "Description", "p value", "Genes")
        
        relevant.res <- relevant.res[(relevant.res[["Genes"]] != ""),]

        
          
        relevant.res.2 <- relevant.res
        
        ## Associating links to the genes
        genes.with.link <- vector(length=nrow(relevant.res))
        for(i in 1:nrow(relevant.res))
        {
         chlamy.genes <- strsplit(relevant.res[i,4],split=" ")[[1]]
         gene.href <- vector(length=length(chlamy.genes))
         for(j in 1:length(chlamy.genes))
           {
           # gene.url <- paste("https://phytozome.jgi.doe.gov/pz/portal.html#!results?search=0&crown=1&star=1&method=4614&searchText=",chlamy.genes[j],sep="")
           gene.href[j] <- paste(c("<a href=\"",gene.url,"\" target=\"_blank\">",chlamy.genes[j],"</a>"),collapse="")
         }
        
         gene.href.string <- paste(gene.href,collapse=" ")
         genes.with.link[i] <- gene.href.string
        }
        relevant.res[["Genes"]] <- genes.with.link



        ## Associating links to the GO terms
        go.links <- vector(length=nrow(relevant.res))
        for(i in 1:nrow(relevant.res))
        {
           go.term.url <- paste("http://amigo.geneontology.org/amigo/term/",relevant.res[i,1],sep="")
           href.go.term <- paste(c("<a href=\"",go.term.url,"\" target=\"_blank\">",relevant.res[i,1],"<a>"),collapse="")
           go.links[i] <- href.go.term
        }
        
        relevant.res[["GO term"]] <- go.links
        
        output$message <- renderUI({
          tagList(
            tags$p(tags$b("The table and graph below presents the enriched GO terms in the
                          selected genes. Both can be downloaded by 
                          clicking on the corresponding buttons:")),
            tags$br())
        })
        
        output$download <- renderUI({
          tagList(
            downloadButton(outputId= "downloadGOTable", "Get GO Enrichment Data"),
            tags$br(),
            tags$br()
          )
        })
        
        output$output_table <- renderDataTable({
          relevant.res
        },escape=FALSE)
        
        ## Plotting the DAG with the first 20 significant GO terms
        output$output_graph <- renderPlot({
          showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all')
        },height = 1000)
        
        
        
        table.data <- relevant.res
        
        
        
      } else if(input$go_mode == "Based on Arabidopsis orthology")
      {
        #print("Arabidopsis orthology")
        ## Set the background gene set to the whole gene set in A. thaliana 
        ## We create a vector with all elements 1 and name it with the gene 
        ## identifiers 
        gene.background <- rep(1,length(probe.names))
        names(gene.background) <- probe.names
        
        cre.gene.names <- gene.list
        
        ## Extract A. thaliana orthologs
        athaliana.gene.names <- athaliana.names[cre.gene.names]
        athaliana.gene.names <- athaliana.gene.names[athaliana.gene.names != "*"]
        
        ## Convert to probe names
        ath.gene.names <- probe.names[athaliana.gene.names]
        
        ## In order to distinguish between background genes and genes of interest the value 0 is 
        ## associated to our genes of interest
        gene.background[ath.gene.names] <- 0
        
        #print("apply topGO")
        
        ## Construct the topGOdata object
        sampleGOdata <- new("topGOdata",
                            description = "Arabidopsis session", ontology = ontology.type,
                            allGenes = gene.background, geneSel = ath.gene.selec,
                            nodeSize = 10,
                            annot = annFUN.db, affyLib = "ath1121501.db")
        
        
        resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
        
        ## Table containing the 100 most significant GO terms to be
        ## save as a text file. To avoid an error when no significant
        ## GO term is found a tryCatch sentence is used.
        
        ms <- 1
        number.top.nodes <- 100
        
        while(ms == 1 | ms == 2)
        {
          
          ms <- tryCatch(
            {
              allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                                 ranksOf = "classicFisher", topNodes = number.top.nodes, numChar=100)
              ms <- 3
            },warning = function(war) {
              
              # warning handler picks up where error was generated
              print(paste("MY_WARNING:  ",war))
              print(number.top.nodes)
              return(1)
              
            }, error = function(err) {
              
              # error handler picks up where error was generated
              print(paste("MY_ERROR:  ",err))
              return(2)
              
            }, finally = {
              print(3)
              number.top.nodes <- number.top.nodes - 1
            }
          )
          
        }
        
        ## Generating for each GO term the Cre genes associated with it
        sig.go.terms <- allRes[["GO.ID"]]
        genes.with.GO <- vector(length=length(sig.go.terms))
        genes.with.GO.no.links <- vector(length=length(sig.go.terms))
        
        for(i in 1:length(sig.go.terms))
        {
          ath.go.terms <- unique(org.At.tairGO2ALLTAIRS[[sig.go.terms[i]]])
          ath.go.names <- intersect(ath.go.terms,athaliana.gene.names)
          chlamy.go.names <- vector()
          for(j in 1:length(ath.go.names))
          {
            chlamy.go.names <- c(chlamy.go.names, chlamy.names[ath.go.names[j] == athaliana.names])
            chlamy.go.names <- unique(chlamy.go.names)
          }
          
          href.chlamy.go.names <- vector(length=length(chlamy.go.names))
          for(k in 1:length(chlamy.go.names))
          {
            url.chlamy.gene <- paste("https://phytozome.jgi.doe.gov/pz/portal.html#!results?search=0&crown=1&star=1&method=4614&searchText=",chlamy.go.names[k],sep="")
            href.chlamy.go.names[k] <- paste(c("<a href=\"",url.chlamy.gene,"\" target=\"_blank\">",chlamy.go.names[k],"</a>"),collapse="")
          }
          
          chlamy.go.names.string <- paste(href.chlamy.go.names,collapse=" ")
          genes.with.GO[i] <- chlamy.go.names.string
          genes.with.GO.no.links[i] <- chlamy.go.names
        }
        
        allRes[["Chlamy.GO"]] <- genes.with.GO
        
        go.links <- vector(length=length(sig.go.terms))
        
        # Adding link to the GO id column
        for(i in 1:length(sig.go.terms))
        {
          go.url <- paste("http://amigo.geneontology.org/amigo/term/",sig.go.terms[i],sep="")
          go.href <- paste(c("<a href=\"",go.url,"\" target=\"_blank\">",sig.go.terms[i],"<a>"),collapse="")
          go.links[i] <- go.href
        }
        
        allRes[["GO.ID"]] <- go.links
        
        ## Extracting only the relevant columns
        relevant.res <- allRes[,c("GO.ID","Term","classicFisher","Chlamy.GO")]
        colnames(relevant.res) <- c("GO term", "Description", "p value", "Genes")
        
        table.data <- relevant.res
        
        
        output$message <- renderUI({
          tagList(
            tags$p(tags$b("The table and graph below presents the enriched GO terms in the
                          selected genes. Both can be downloaded by 
                          clicking on the corresponding buttons:")),
            tags$br())
        })
        
        output$download <- renderUI({
          tagList(
            downloadButton(outputId= "downloadGOTable", "Get GO Enrichment Data"),
            tags$br(),
            tags$br()
          )
        })
        
        output$output_table <- renderDataTable({
          relevant.res
        },escape=FALSE)
        
        ## Plotting the DAG with the first 20 significant GO terms
        output$output_graph <- renderPlot({
          showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all')
        },height = 1000)
        
        
        
        
        
        # ## Plotting the DAG with the first 20 significant GO terms
        # output$go_enrichment_plot <- renderPlot({
        #   showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = 20, useInfo = 'all')
        # },height = 2000)
      }
      
    })
  
}



# Run the application 
shinyApp(ui = ui, server = server)

