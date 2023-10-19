#' @importFrom EBImage readImage
#' @importFrom stats quantile


.calcthresh <- function(path_tiff, bits){
  if(bits == 16){
    tiff = readImage(path_tiff)
    img <- tiff * 2^16-1
    thresh_max <- exp(quantile(log(img), probs = 0.98, na.rm=T)[[1]]) #98% de la somme du signal = thresh de base
  }else{
    img = readImage(path_tiff)
    thresh_max <- exp(quantile(log(img), probs = 0.993, na.rm=T)[[1]])
  }
  valmaxsignal <- round(max(img),1)
  fin <- c(thresh_max, valmaxsignal)
  return(fin)
}




.gettableaufinal <- function(markers, minsImg, maxImg, maxautoImg, listMinMask){

  bilantmp <- data.frame(markers, minsImg, maxImg, maxautoImg)
  dftmp <- as.data.frame(listMinMask)
  dftmp[2,] <- colnames(dftmp)
  df <- t(dftmp)
  colnames(df) <- c("threshminmask", "markers")
  bilan <- merge(bilantmp, df, all=T)
  colnames(bilan) <- c("Markers", "Thresh Image min", "Thresh Image max", "Thresh Image max auto", "Thresh Min Mask")

  return(bilan)
}




#' @importFrom graphics text

.gettext <- function(datff, channels, thresh){
  tot <- nrow(datff)
  popavecthresh <- length(which(datff@exprs[,channels[[1]]]>=thresh))
  prop1 <- round((popavecthresh/tot)*100, 2)
  propf <- paste(prop1, "% \n of Total Cells", sep = "")
  xx <- (range(datff@exprs[,channels[[1]]])[[2]]/10)*7
  yy <- (range(datff@exprs[,channels[[2]]])[[2]]/10)*8
  return(text(xx, yy, propf))
}




#' @importFrom EBImage abind
#' @importFrom EBImage paintObjects
#' @importFrom EBImage getFrames
#' @importFrom EBImage colorMode


.getImage <- function(resu_mark, pathtiffs_sansmask, resu_col, mask, segmentation, lesthresh, nbbits){
  if(length(resu_mark)>1){
    which_one <- which(resu_mark == TRUE)
    tiffs_sel <- pathtiffs_sansmask[c(which_one)]
    cols_sel <- resu_col[c(which_one)]

    tiffs_to_plot = list()
    for(i in 1:length(tiffs_sel)){
      tiffs_to_plot[[i]] <- .prepareImage(path_tiff=tiffs_sel[[i]], col =cols_sel[[i]],
                                         thresh_min = lesthresh[,which_one[[i]]][[1]], thresh_max = lesthresh[,which_one[[i]]][[2]],
                                         bitsimg = nbbits)
    }

    imgdat <- NULL
    imgdat <- EBImage::abind(tiffs_to_plot[1:length(tiffs_to_plot)], along =4)
    EBImage::colorMode(imgdat) <- "Color"
    imgdat2 <- Reduce("+", getFrames(imgdat, type="render"))

    if(segmentation == TRUE){
      segmented <- EBImage::paintObjects(mask, imgdat2, col=c('white'), opac = c(0.6,1))
      img_resize <- segmented
    }else{
      img_resize <- imgdat2
    }
    return(img_resize)
  }
}




#' @importFrom Biostrings pairwiseAlignment


.plotdot <- function(data_ff, chans, thresh, tiffname_sans_mask, actual_mark_fcs, markers){

  # On retrouve le marqueur d'ADN
  dnachan <- tryCatch(markers[[grep("dna",markers, ignore.case = T)[[1]]]],  error = function(e) "PASDADN")
  marktoplot <- dnachan

  if(length(chans) == 2){ # Au moins un marqueur de coché
    mark_tmp <- tools::file_path_sans_ext(tiffname_sans_mask[[as.numeric(chans[[2]])]]) # Recup du marqueur coché (tiff file)
    # On retrouve le marqueur dans le fcs par similarité (permet d'avoir des noms diff entre fichiers et param fcs):
    idxmark <- which.max(pairwiseAlignment(rep(mark_tmp, length(markers)), markers, type = "overlap")@score)
    mark <- markers[[idxmark]]
    if(mark != actual_mark_fcs){ marktoplot <- mark }
  }

  if(marktoplot == "PASDADN"){
    return("")
  }else{
    .baseplotdens(ff = data_ff, canaux = c(chans[[1]], marktoplot), fcsthresh = thresh)
  }

}



#' @importFrom grDevices densCols
#' @importFrom grDevices adjustcolor
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot
#' @importFrom graphics abline

.baseplotdens <- function(ff, canaux, fcsthresh){
  df <- as.data.frame(ff@exprs)
  colPalette <- colorRampPalette(colors = c("blue", "turquoise", "green", "yellow", "orange", "red"))
  col <- suppressWarnings(grDevices::densCols(df[,canaux], colramp = colPalette))
  col <- grDevices::adjustcolor(col)
  graphics::plot(df[,canaux], col = col, pch = 19)
  abline(v = fcsthresh, col = "red", lty = 2, lwd = 3)
  .gettext(datff = ff, channels = canaux, thresh = fcsthresh)
}



#' @importFrom EBImage readImage
#' @importFrom EBImage paintObjects
#' @importFrom EBImage channel
#' @importFrom EBImage display
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom graphics par


.plotmask <- function(img, showmaskfcs, ff, h, w, gauche, bas, Xpos_mask, Ypos_mask, col_mask = "red" ){
  if(showmaskfcs){
    haut = -bas
    droite = -gauche
    tmp <- tempfile()
    png(tmp, width = w, height = h)
    par(bg = NA, mar = rep(0,4))
    plot(x = Xpos_mask, y = Ypos_mask,
         pch = 16, col = "red", cex = 0.5,
         axes=F, ann=F,
         xaxs = "i", yaxs = "i",
         ylim = c(h+bas,0+haut),
         xlim=c(0+gauche,w+droite))
    dev.off()
    imgfcs <- readImage(tmp, 'png')
    grayimage<-channel(imgfcs,"gray")
    imgf  = paintObjects(grayimage, tgt = img, col = col_mask)
    EBImage::display(imgf)
  }else{
    EBImage::display(img)
  }
}



#' @importFrom EBImage readImage
#' @importFrom EBImage normalize
#' @importFrom EBImage colormap
#' @importFrom grDevices colorRampPalette
#' @importFrom stats quantile


.prepareImage <- function(path_tiff, thresh_min = 0, thresh_max = NULL, col = NULL, bitsimg){

  tiff = EBImage::readImage(path_tiff)
  if(bitsimg == 16){
    img = tiff * 2^16-1
    if(is.null(thresh_max)){
      thresh_max <- exp(quantile(log(img), probs = 0.98, na.rm=T)[[1]])
    }
  }else{
    img = tiff
    if(is.null(thresh_max)){
      thresh_max <- exp(quantile(log(img), probs = 0.993, na.rm=T)[[1]])
    }
  }
  res_img_tmp <- EBImage::normalize(img, inputRange=c(thresh_min,thresh_max))
  if(!is.null(col)){
    jet.colors = colorRampPalette(c("black", col))
  }else{
    jet.colors = colorRampPalette(c("black", "white"))
  }
  res_img = colormap(res_img_tmp, jet.colors(256))

  return(res_img)
}




#' @importFrom colourpicker colourInput
#' @importFrom shiny renderUI
#' @importFrom shiny tagList
#' @importFrom shiny tags
#' @importFrom shiny div


.renderColorChoiceImage <- function(noms_fichiers_sans_masques){
  return(
    renderUI({
      lapply(1:length(noms_fichiers_sans_masques), function(x){
        tagList(tags$style(type = 'text/css', '#big_slider2 {height: 110px; width: 60%;}'),
                div(id = 'big_slider2', colourInput(inputId = paste("col", x, sep=""), label = "", value = "white", palette = "limited")
                ))
      })
    })
  )
}





#' @importFrom shinydashboard box
#' @importFrom EBImage displayOutput
#' @importFrom shiny splitLayout
#' @importFrom shiny sliderInput
#' @importFrom shiny uiOutput
#' @importFrom shiny verbatimTextOutput
#' @importFrom shiny renderUI
#' @importFrom shiny h2


.renderForMultipleMarqueurs <- function(dat, marks, bas, gauche, vals_thresh_save){
  express <- dat@exprs
  threshs_max_fcs <- lapply(1:length(marks), function(y){ round(max(express[,y]),1) })
  return(
    renderUI(
      box(
        id = "box_options",
        solidHeader = F,
        collapsible = F,
        width = 12,
        EBImage::displayOutput("plotting", width = "100%", height = "700px"),
        shiny::splitLayout(
          shiny::checkboxInput("show_mask", label="Centroids from FCS", value = TRUE, width = "100px"),
          lapply(1:length(marks), function(x){
            shiny::sliderInput(inputId = paste("Thresh", x, sep=""), min = 0, label = h2(paste0("Intensity threshold ", marks[[x]], sep="")),
                               max = threshs_max_fcs[[x]], value = vals_thresh_save[[x]], step = 0.01, width = "400px")
          }),
          shiny::uiOutput("mask_segmentation_version_tiff")),
        shiny::verbatimTextOutput("numbercells")
      )
    )
  )
}



#' @importFrom shinydashboard box
#' @importFrom EBImage displayOutput
#' @importFrom shiny splitLayout
#' @importFrom shiny sliderInput
#' @importFrom shiny uiOutput
#' @importFrom shiny verbatimTextOutput
#' @importFrom shiny renderUI
#' @importFrom shiny h4


.renderForOnlyOneMarqueur <- function(dat, mark, bas, gauche, val_thresh_save){
  express <- dat@exprs
  thresh_max_fcs = round(max(express[,mark]),1)
  return(
    renderUI(
      box(
        id = "box_options",
        solidHeader = F,
        collapsible = F,
        width = 12,
        EBImage::displayOutput("plotting", width = "100%", height = "700px"),
        shiny::splitLayout(
          shiny::checkboxInput("show_mask", label="Centroids from FCS", value = TRUE, width = "100px"),
          shiny::sliderInput("Thresh", min = 0, label = h4("Intensity threshold"), max = thresh_max_fcs, value = val_thresh_save, step = 0.01, width = "400px"),
          shiny::uiOutput("mask_segmentation_version_tiff")),
        shiny::verbatimTextOutput("numbercells")
      )
    )
  )
}



#' @importFrom shiny renderUI
#' @importFrom shiny tagList
#' @importFrom shiny checkboxInput
#' @importFrom shiny div



.renderMarkerChoiceImage <- function(noms_des_fichiers_sans_masques){
  return(
    renderUI({
      lapply(1:length(noms_des_fichiers_sans_masques), function(x){
        tagList(tags$style(type = 'text/css', '#big_slider1 {height: 100px; width: 100%;}'),
                div(id = 'big_slider1',
                    checkboxInput(inputId = paste("mark", x, sep=""), label = basename(noms_des_fichiers_sans_masques[[x]]), value = FALSE)
                ))
      })
    })
  )
}




#' @importFrom shiny renderUI
#' @importFrom shiny tagList
#' @importFrom shiny sliderInput
#' @importFrom shiny div


.renderThreshChoiceImage <- function(noms_files_sans_mask, chemins_files_sans_mask, thresh_image, imgbits){
  return(
    renderUI({
      lapply(1:length(noms_files_sans_mask), function(x){
        tocalc <- .calcthresh(path_tiff = chemins_files_sans_mask[[x]], bits = imgbits)
        min <- 0
        max <- tocalc[[1]]
        if(!is.null(thresh_image)){
          if(thresh_image[1,][[x]] != 0){
            min <- thresh_image[1,][[x]]
          }
          if(thresh_image[2,][[x]] != tocalc[[1]]){
            max <- thresh_image[2,][[x]]
          }
        }
        tagList(tags$style(type = 'text/css', '#big_slider3 {height: 110px;}'),
                div(id = 'big_slider3',
                    sliderInput(inputId = paste("thresh_im", x, sep=""), label = paste("Thresh max auto : ", round(tocalc[[1]], 2)),
                                min = 0, max = tocalc[[2]], value = c(min, max), dragRange = FALSE, step=0.1)
                ))
      })
    })
  )
}


#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom methods new

.newff <- function(matrice){

  metadata <- data.frame(name=colnames(matrice),desc=colnames(matrice))
  metadata$minRange <- apply(matrice,2,min)
  metadata$maxRange <- apply(matrice,2,max)
  Descr <- apply(matrice, 2, function(x){ paste0(round(min(x), 1), ",", round(max(x), 0) +1, sep="")})
  Descr_list <- as.list(Descr)
  names(Descr_list) <- sapply(c(1:length(Descr_list)), function(x) {paste0("$P",x,"D")})
  ff <- new("flowFrame",exprs=as.matrix(matrice), parameters=Biobase::AnnotatedDataFrame(metadata), description = Descr_list)

  return(ff)
}
