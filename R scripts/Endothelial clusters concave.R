
library(concaveman)

endo_concavity <- getOption("endo.concavity", 2)
extract_concave_coords <- function(hull) {
  if (is.list(hull) && !is.null(hull$polygon)) {
    hull <- hull$polygon
  }
  hull <- as.matrix(hull)
  if (!is.matrix(hull) || ncol(hull) < 2) {
    stop("need at least two columns of coordinates")
  }
  if (nrow(hull) >= 3) {
    x <- hull[, 1]
    y <- hull[, 2]
    x_next <- c(x[-1], x[1])
    y_next <- c(y[-1], y[1])
    signed_area <- sum(x * y_next - x_next * y) / 2
    if (!is.na(signed_area) && signed_area < 0) {
      hull <- hull[rev(seq_len(nrow(hull))), , drop = FALSE]
    }
  }
  hull
}

polyAreaList <- list()   

for (ROIImage in unique(datan$ROI)){
  print(ROIImage)
  
  # Extract ROI data
  ROIsubset <- datan[ROI==ROIImage & Phenotype == "Endo",..PointDatacols]
  
  if(nrow(ROIsubset)>15){
    # For each cell type i
    X <- ppp(x = ROIsubset[,im], 
             y = ROIsubset[,jm],
             window = owin(c(ROIsubset[,min(im)],
                             ROIsubset[,max(im)]),
                           c(ROIsubset[,min(jm)],
                             ROIsubset[,max(jm)])),
             # marks = ROIsubset[,..PointDatalocation],
             unitname = "micrometer")
    
    Xconnect <- connected(X,20)
    
    Endothresholds <- names(table(Xconnect$marks)[table(Xconnect$marks)>10])
    
    datan[ROI==ROIImage & Phenotype == "Endo", Endocluster := as.character(Xconnect$marks)]
    #convex hull computation
    if (length(Endothresholds)) {
      # areas for this ROI
      roiAreas <- lapply(Endothresholds, function(cl) {
        pts <- X[Xconnect$marks == cl]
        # need ≥3 points for a polygon
        coords <- cbind(pts$x, pts$y)
        uniqueCoords <- unique(coords)
        # need ≥3 unique points for a polygon
        if (nrow(uniqueCoords) >= 3) {
          hullCoords <- extract_concave_coords(
            concaveman::concaveman(uniqueCoords, concavity = endo_concavity)
          )
          hullWindow <- owin(poly = list(x = hullCoords[, 1],
                                         y = hullCoords[, 2]))
          areaHull <- area.owin(hullWindow)           # μm²
          
        } else {
          areaHull <- NA_real_
        }
        data.table(ROI = ROIImage,
                   Endocluster = as.character(cl),
                   Endoarea = areaHull)
      })
      polyAreaList[[ROIImage]] <- rbindlist(roiAreas)
    }
    #datan[, Endoclusterviable := NULL]
    
    datan[ROI==ROIImage & Phenotype == "Endo" & Endocluster %in% Endothresholds, Endoclusterviable := T]
    
  }
}

EndopolyAreas <- rbindlist(polyAreaList, use.names = TRUE, fill = TRUE)
datan <- EndopolyAreas[datan, on = .(ROI, Endocluster)]

###vascularization model

#datan[!is.na(TLSarea), .N, by = .(ROI, TLSarea)][N > 1][order(ROI, -N)]
datan[, TLS_true_count := sum(TLSclusterviable, na.rm = TRUE), by = ROI]
datan[, ROI_Endoarea_convex := sum(unique(TLSarea), na.rm=TRUE), by = ROI]
datan[, ROI_Endoarea_concave := sum(unique(Endoarea), na.rm=TRUE), by = ROI]
datan[, patient_endoarea_convex := sum(unique(TLSarea), na.rm=TRUE), by=case]
datan[, patient_endoarea_concave := sum(unique(Endoarea), na.rm=TRUE), by=case]