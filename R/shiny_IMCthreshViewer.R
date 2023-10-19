
#' @importFrom shinydashboard dashboardPage
#' @importFrom shinydashboard dashboardHeader
#' @importFrom shinydashboard dashboardSidebar
#' @importFrom shinydashboard sidebarMenu
#' @importFrom shinydashboard menuItem
#' @importFrom shiny icon
#' @importFrom shinydashboard dashboardBody
#' @importFrom shinyjs extendShinyjs
#' @importFrom shinyjs useShinyjs
#' @importFrom shinydashboard tabItems
#' @importFrom shinydashboard tabItem
#' @importFrom shiny fluidPage
#' @importFrom shiny fluidRow
#' @importFrom shiny column
#' @importFrom shiny br
#' @importFrom shiny h1
#' @importFrom shinydashboard box
#' @importFrom shiny fileInput
#' @importFrom shiny uiOutput


#' @importFrom shinyjs addClass
#' @importFrom shiny reactiveValues
#' @importFrom shiny reactive
#' @importFrom shiny renderUI
#' @importFrom shiny observeEvent
#' @importFrom shiny selectInput
#' @importFrom flowCore read.FCS
#' @importFrom flowCore featureNames
#' @importFrom EBImage readImage
#' @importFrom stats quantile
#' @importFrom shiny actionButton
#' @importFrom shiny checkboxInput
#' @importFrom shiny debounce
#' @importFrom shiny renderPlot
#' @importFrom shiny renderPrint
#' @importFrom shiny plotOutput
#' @importFrom shiny renderTable
#' @importFrom shiny tableOutput
#' @importFrom shiny stopApp
#' @importFrom grDevices dev.set
#' @importFrom grDevices dev.off
#' @importFrom grDevices dev.next
#' @importFrom shiny shinyApp
#' @importFrom EBImage renderDisplay
#' @importFrom dplyr %>%

NULL

#'
#' shiny_IMCthreshViewer
#'
#' App
#'
#'
#'
#' @return app
#' @encoding UTF-8
#' @export
#'
#'
#' @examples
#' #shiny_IMCthreshViewer()




shiny_IMCthreshViewer <- function() {

  jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"

options(shiny.maxRequestSize = 2000*1024^2)



ui = dashboardPage(


  dashboardHeader(title = "Shiny - Threshold"),
  dashboardSidebar(width = 350,
                   sidebarMenu(id = "tabs",
                               menuItem(
                                 text = "IMCthreshViewer",
                                 tabName = "Threshold",
                                 icon = icon("eye",verify_fa = FALSE, lib= "font-awesome"))
                   )),
  dashboardBody(shinyjs::useShinyjs(),shinyjs::extendShinyjs(text = jscode, functions = "collapse"),
                tabItems(


                  # Tab1 : Conversion csv to fcs --------------------------------------------

                  tabItem("Threshold",
                          fluidPage(
                            fluidRow(
                              column(6, align="center", offset = 3,
                                     br(),br(),br(),
                                     h1("IMC Thresholds Viewer", align="center"),
                                     br(), br()
                              )
                            )
                          ),


                          fluidPage(
                            fluidRow(
                              column(12, offset = 2,
                                     box(
                                       id = "box_files_norm",
                                       title = "Input Files",
                                       status = "danger",
                                       solidHeader = TRUE,
                                       collapsible = TRUE,
                                       width = 8,
                                       fileInput("fcs_files", "FCS File:", accept = ".fcs", multiple = F),
                                       fileInput("tiff_files", "TIFF Files :", accept = ".tiff", multiple = T),
                                       uiOutput("show_marqueurs"),
                                       uiOutput("valid")
                                     )
                              )
                            )
                          ),


                          fluidPage(
                            fluidRow(
                              column(1, align="center",
                                     uiOutput("selec_cols")
                              ),
                              column(1, align="center",
                                     br(),
                                     uiOutput("selec_mark")
                              ),
                              column(2, align="center",
                                     br(),
                                     uiOutput("selec_thresh_img")
                              ),
                              column(8, align="center",
                                     br(), br(),
                                     uiOutput("box_plot"),
                                     uiOutput("box_dotplot")
                              ),
                              column(8, align="center",
                                     uiOutput("montre_glob")
                              )
                            )
                          )


                  )
                )
  ),
  title = "IMCthreshViewer"
)
server = function(input, output, session) {


  shinyjs::addClass(selector = "body", class = "sidebar-collapse")

  # Variables globales utilisées------------------------------------------------------

  global <- reactiveValues(
    data = NULL, #flowframe
    markers = NULL,
    namefile_tiffs_sans_mask = NULL,
    path_tiffs_sans_mask = NULL,
    mask = NULL, # mask tiff
    showMaskSegm = FALSE, # remplace input$show_mask_segmentation
    interactivityOK = NULL, 
    h = NULL, #height images
    w = NULL, #width images
    bits = NULL, #nb of bits
    res_mark = NULL,
    res_col = NULL,
    res_thresh_im = NULL,
    Xpos_mask = NULL,
    Ypos_mask= NULL,
    nrow = NULL,
    new_idx = NULL,
    old_idx = NULL,
    mark_and_threshfcs = list(),
    listemarq = NULL, # pour etre dependant du bouton validation sinon bug (modif dans l'input$listemarq)
    idx_thresh = list() #idx des thresh si strat gating
  )

  # _____1 ------------------------------------------------------------------


  # ----------------------------------------------------------
  # ----------- CHARGEMENT FCS
  # ----------------------------------------------------------
  observeEvent(input$fcs_files, {

    global$interactivityOK = NULL
    global$Xpos_mask = NULL
    global$Ypos_mask = NULL
    global$mark_and_threshfcs = list()

    # Récupération des fichiers
    cheminfcs <- input$fcs_files
    # Récupération des data
    global$data <- flowCore::read.FCS(cheminfcs$datapath, truncate_max_range = F)
    # Volet sélection des marqueurs issus du fcs
    markerstmp <- as.vector(featureNames(global$data))
    # Adjust to steinbock
    steintest <- length(which(markerstmp =="centroid-0"))
    if(steintest == 1){
      tochange <- as.data.frame(global$data@exprs)
      coltochange <- colnames(tochange)
      coltochange[which(coltochange == "centroid-0")] <- "Y_position"
      coltochange[which(coltochange == "centroid-1")] <- "X_position"
      coltochange[which(coltochange == "DNA1")] <- "DNA1_Ir191"
      colnames(tochange) <- coltochange
      global$data = NULL
      global$data <- .newff(matrice = tochange)
      markerstmp <- as.vector(featureNames(global$data))
    }
    global$markers <- markerstmp[-c(which(markerstmp == "Y_position"), which(markerstmp == "X_position"))]
    output$show_marqueurs <- renderUI(selectInput(width = "400px", "listemarq", "List of fcs markers :", choices=global$markers, selected = global$markers[[1]], multiple=TRUE))
  })


  # ----------------------------------------------------------
  # ----------- CHARGEMENT TIFFS
  # ----------------------------------------------------------
  observeEvent(input$tiff_files, {

    global$interactivityOK = NULL
    global$new_idx = NULL
    global$old_idx = NULL

    # Récupération des fichiers
    chemintiff <- input$tiff_files

    if(length(grep(".*_mask.tiff", chemintiff$name)) != 1){

      message("--------------------------------------")
      message("NO MASK.")
      message("--------------------------------------")
      global$namefile_tiffs_sans_mask <- chemintiff$name
      global$path_tiffs_sans_mask <- chemintiff$datapath
      tiff_alea <- chemintiff$datapath[[1]]
      tmp = EBImage::readImage(tiff_alea)
      tmp1 = tmp * 2^16-1
      global$w = dim(tmp1)[1] # image width
      global$h = dim(tmp1)[2] # image height

    }else{

      # Récupération du masque tiff
      whichonemasktiff <- grep(".*_mask.tiff", chemintiff$name)
      global$namefile_tiffs_sans_mask <- chemintiff$name[-whichonemasktiff]
      global$path_tiffs_sans_mask <- chemintiff$datapath[-whichonemasktiff]
      tiff_mask <- chemintiff$datapath[[whichonemasktiff]]
      mask_tmp = EBImage::readImage(tiff_mask)
      global$mask = mask_tmp * 2^16-1
      global$w = dim(global$mask)[1] # image width
      global$h = dim(global$mask)[2] # image height

      # Check out the number of bits
      other_tiff <- chemintiff$datapath[[whichonemasktiff+1]]
      if(length(other_tiff)==0){ other_tiff <- chemintiff$datapath[[whichonemasktiff-1]] }
      imgtestbits_tmp <- readImage(other_tiff)
      imgtestbits_tmp2 <- imgtestbits_tmp* 2^16-1
      imgtestbits <- exp(quantile(log(imgtestbits_tmp2), probs = 0.98, na.rm=T)[[1]])
      if(imgtestbits > 10000){ #On est en mode 32bits
        global$bits = 32
        message("32-bit tiffs.")
      }else{
        global$bits = 16
        message("16-bit tiffs.")
      }


    }

    # Bouton validation ( /!\ s'affiche dès le fichier tiff, donc possibilite de bug si validation sans fcs)
    output$valid <- renderUI(actionButton(inputId = "validation", label = "OK", icon("angle-double-right"),
                                          style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))

  })



  # ----------------------------------------------------------
  # ----------- AFFICHAGE DES BOX
  # ----------------------------------------------------------
  observeEvent(input$validation, {

    global$listemarq = NULL # pour etre dependant du bouton validation sinon bug (modif dans l'input$listemarq)
    global$idx_thresh = list()
    global$listemarq <- input$listemarq

    # Affichage des choix de couleurs :
    output$selec_cols <- .renderColorChoiceImage(noms_fichiers_sans_masques = global$namefile_tiffs_sans_mask)
    # Affichage des choix de marqueurs :
    output$selec_mark <- .renderMarkerChoiceImage(noms_des_fichiers_sans_masques = global$namefile_tiffs_sans_mask)
    # Affichage des choix thresh img :
    output$selec_thresh_img <- .renderThreshChoiceImage(noms_files_sans_mask = global$namefile_tiffs_sans_mask,
                                                       chemins_files_sans_mask = global$path_tiffs_sans_mask,
                                                       thresh_image = global$res_thresh_im,
                                                       imgbits = global$bits)

    if(!is.null(global$mask)){
      output$mask_segmentation_version_tiff <- renderUI({
        shiny::checkboxInput("show_mask_segmentation", label="Segmentation Mask", value = FALSE, width = "100px")})
    }else{ # Si pas de masque dans le tiff
      output$mask_segmentation_version_tiff <- renderUI({""})
      global$showMaskSegm = FALSE
    }

    global$interactivityOK <- "OK"

    if(length(global$listemarq) == 1){# Si un seul marqueur de sélectionné

      # Si premier passage de ce marqueur :
      if(is.null(global$mark_and_threshfcs[[global$listemarq]])){
        global$mark_and_threshfcs[[global$listemarq]] <- 0
        val = 0
      }else{
        val <- global$mark_and_threshfcs[[global$listemarq]]
      }

      # Affichage Box pour gestion des thresholds image/masque + affichage nombre de cellules resultant dans le masque
      output$box_plot <- .renderForOnlyOneMarqueur(dat = global$data, mark = global$listemarq, bas = 0, gauche = 0, val_thresh_save = val)


    }else{# Si plusieurs marqueurs sélectionnés

      # Si premier passage de ces marqueurs, récupération des valeurs des thresholds enregistrées :
      vals = list()
      for(cpt in 1:length(global$listemarq)){
        if(is.null(global$mark_and_threshfcs[[global$listemarq[[cpt]]]])){
          global$mark_and_threshfcs[[global$listemarq[[cpt]]]] <- 0
          vals[[cpt]] = 0
        }else{
          vals[[cpt]] <- global$mark_and_threshfcs[[global$listemarq[[cpt]]]]
        }
      }

      output$box_plot <- .renderForMultipleMarqueurs(dat = global$data, marks = global$listemarq, bas = 0, gauche = 0, vals_thresh_save = vals)

    }

  })



  observeEvent(input$show_mask_segmentation,{ global$showMaskSegm = input$show_mask_segmentation })


  # ----------------------------------------------------------
  # ----------- RECUPERATION DES MARQUEURS COCHES ET DES COULEURS
  # ----------------------------------------------------------


  res_mark_reactive <- reactive({
    if(!is.null(global$interactivityOK)){
      lapply(
        X = 1:length(global$namefile_tiffs_sans_mask),
        FUN = function(i){
          observeEvent(input[[paste0("mark", i)]], {
            global$res_mark <- sapply(1:length(global$namefile_tiffs_sans_mask), function(i) {paste0("", input[[paste0('mark', i)]])})

            if(global$res_mark[[i]] == TRUE){ # Coche marqueur
              if(is.null(global$new_idx)){ # un seul select
                global$new_idx = i # init
                global$old_idx = i # init
              }else{ # Plusieurs select
                global$new_idx = i # Dernier select
                global$old_idx = append(global$old_idx,i) #On stock les infos des selects
              }
            }else{ # decoche un marqueur
              if(!is.null(global$old_idx) && !is.null(global$new_idx)){ #si les var on ete init
                if(i == global$new_idx){ # Si le marq decocher correspond au dernier select
                  if(length(global$old_idx)>1){ # Si il y a d'autres marqueurs de cochés
                    global$old_idx = global$old_idx[-(length(global$old_idx))] #on enleve le dernier decocher
                    global$new_idx <- global$old_idx[[length(global$old_idx)]] #le dernier select devient celui a afficher
                  }else{ # Decoche du dernier marqueur select
                    global$new_idx = NULL # reinit
                    global$old_idx = NULL # reinit
                  }
                }else{ # Si on decoche un autre marqueur
                  global$old_idx = global$old_idx[-(which(global$old_idx==i))]
                }
              }
            }
          })
        }
      )
    }
  })

  res_col_reactive <- reactive({
    if(!is.null(global$interactivityOK)){
      lapply(
        X = 1:length(global$namefile_tiffs_sans_mask),
        FUN = function(i){
          observeEvent(input[[paste0("col", i)]], {
            global$res_col <- sapply(1:length(global$namefile_tiffs_sans_mask), function(i) {paste0("", input[[paste0('col', i)]])})
          })
        }
      )
    }
  })

  res_thresh_reactive <- reactive({
    if(!is.null(global$interactivityOK)){
      lapply(
        X = 1:length(global$namefile_tiffs_sans_mask),
        FUN = function(i){
          observeEvent(input[[paste0("thresh_im", i)]], {
            global$res_thresh_im <- sapply(1:length(global$namefile_tiffs_sans_mask), function(i) {paste0("", input[[paste0('thresh_im', i)]])})
          })
        }
      )
    }
  }) %>% debounce(2000) # Pour éviter bug de boucles, délai de 2sec



  observeEvent(input$Thresh, { # Si un seul marqueur /!\

    ff <- global$data
    idx <-which(ff@exprs[, global$listemarq] >= input$Thresh)

    if(length(idx) == 1){ # S'il ne reste qu'une cellule
      global$Xpos_mask <- ff@exprs[,"X_position"][which(ff@exprs[,global$listemarq] == max(ff@exprs[,global$listemarq]))]
      global$Ypos_mask <- ff@exprs[,"Y_position"][which(ff@exprs[,global$listemarq] == max(ff@exprs[,global$listemarq]))]
      global$nrow <- 1
    }else{
      ff@exprs <- ff@exprs[idx,]
      global$Xpos_mask <- ff@exprs[,"X_position"]
      global$Ypos_mask <- ff@exprs[,"Y_position"]
      global$nrow <- nrow(ff)
    }

    global$mark_and_threshfcs[[global$listemarq]] <- input$Thresh

  })


  ################################ WATCHOUTTTTTTTTTTTTTTTTTTTTT
  res_thresh2_reactive <- reactive({
    if(!is.null(global$interactivityOK)){
      lapply(
        X = 1:length(global$listemarq),
        FUN = function(i){
          observeEvent(input[[paste0("Thresh", i)]], {

            ff <- global$data
            global$idx_thresh[[i]] <-which(ff@exprs[, global$listemarq[[i]]] >= input[[paste0("Thresh", i)]])

            global$mark_and_threshfcs[[global$listemarq[[i]]]] <- input[[paste0("Thresh", i)]]

            finalf <- global$data
            finalf@exprs <- finalf@exprs[Reduce(intersect, global$idx_thresh),]
            global$Xpos_mask <- finalf@exprs[,"X_position"]
            global$Ypos_mask <- finalf@exprs[,"Y_position"]
            global$nrow <- nrow(finalf)

          })
        }
      )

    }
  }) %>% debounce(2000)






  currentImage <- reactive({
    if(!is.null(global$interactivityOK)){
      if(length(which(global$res_mark == TRUE))>0){
        res <- .getImage(resu_mark =global$res_mark, pathtiffs_sansmask = global$path_tiffs_sans_mask, resu_col = global$res_col,
                        mask = global$mask, segmentation = global$showMaskSegm, nbbits = global$bits,
                        lesthresh = global$res_thresh_im)
        .plotmask(img = res, showmaskfcs = input$show_mask, ff= global$data, h = global$h, w = global$w, gauche = 0,
                  bas = 0, Xpos_mask =global$Xpos_mask, Ypos_mask = global$Ypos_mask)
      }else{""}
    }else{""}
  })


  # Affichage image
  observeEvent(global$res_mark,{
    if(length(which(global$res_mark == TRUE))>0){
      output$plotting <- renderDisplay({currentImage()})
    }else{ # Quand rien n'est sélectionné
      output$plotting <- NULL
    }
  })

  currentdotplot <- reactive({
    if(!is.null(global$interactivityOK)){
      if(length(global$listemarq) == 1){
        .plotdot(data_ff = global$data, chans = c(global$listemarq,global$new_idx), thresh = input$Thresh,
                 tiffname_sans_mask = global$namefile_tiffs_sans_mask, actual_mark_fcs = global$listemarq, markers = global$markers)
      }else{
        .plotdot(data_ff = global$data, chans = c(global$listemarq[[1]],global$new_idx), thresh = input$Thresh,
                 tiffname_sans_mask = global$namefile_tiffs_sans_mask, actual_mark_fcs = global$listemarq[[1]], markers = global$markers)
      }
    }
  })

  output$plottingDot <- renderPlot({currentdotplot()})



  # ----------------------------------------------------------
  # ----------- AFFICHAGE PLOTS MASQUE DEP
  # ----------------------------------------------------------
  observeEvent(input$show_mask, {

    res_mark_reactive()
    res_col_reactive()
    res_thresh_reactive()
    res_thresh2_reactive()

    # Affichage du nombre de cellules restant avec le threshold appliqué :
    output$numbercells <- renderPrint({paste0(global$nrow, " Cells.", sep="")})

    output$box_dotplot <- renderUI(
      box(
        id = "box_dot",
        solidHeader = F,
        collapsible = F,
        width = 12,
        plotOutput("plottingDot"),
        actionButton(label = "Show threshold values", inputId = "Show_Glob", icon("angle-double-right"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
      )
    )

  })

  # ----------------------------------------------------------
  # ----------- AFFICHAGE DES THRESHOLDS
  # ----------------------------------------------------------
  observeEvent(input$Show_Glob,{

    test <- .gettableaufinal(markers = sapply(global$namefile_tiffs_sans_mask, function(x){gsub(pattern = "\\.tiff$", "",x)}),
                              minsImg = global$res_thresh_im[1,],
                              maxImg = global$res_thresh_im[2,],
                              maxautoImg = sapply(global$path_tiffs_sans_mask, function(x){round(.calcthresh(x, bits = global$bits)[[1]],2)}),
                              listMinMask = global$mark_and_threshfcs)

    output$montre_glob <- renderUI({ tableOutput("table") })
    output$table <- renderTable(test)

  })



  # _____2 ------------------------------------------------------------------

  # Fonction pour que le process s'arrete lorsque l'on ferme la fenetre
  session$onSessionEnded(function() {
    dev.set(dev.next())
    dev.off()
    gc()
    rm(list=ls())
    gc()
    stopApp()
  })


}
shinyApp(ui, server)
}
