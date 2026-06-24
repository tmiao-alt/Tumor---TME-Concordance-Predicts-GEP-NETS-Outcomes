
required_packages <- c("data.table", "spatstat.geom")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Please install required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(data.table)
library(spatstat.geom)


datan <- as.data.table(get("datan", inherits = TRUE))

required_columns <- c("ROI", "jm", "im", "ci", "Phenotype_refined")
missing_columns <- setdiff(required_columns, names(datan))
if (length(missing_columns) > 0) {
  stop(
    "`datan` is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

censor <- function(data, percentile = 0.99, cap = NA) {
  data <- as.matrix(data)
  if (is.na(cap)) {
    highlim <- quantile(data, percentile, na.rm = TRUE)
  } else {
    highlim <- cap
  }
  lowlim <- quantile(data, 1 - percentile, na.rm = TRUE)
  data[data < lowlim] <- lowlim
  data[data > highlim] <- highlim
  as.data.table(data)
}

progresstrack <- function(position, list, history) {
  progress <- match(position, list)
  pct <- round(progress / length(list) * 100, 1)
  if (pct >= history + 1 || progress == length(list)) {
    message(pct, "% complete: ", position)
    history <- floor(pct)
  }
  history
}

if (!exists("correctedcelltypes", inherits = TRUE)) {
  correctedcelltypes <- sort(unique(as.character(datan$Phenotype_refined)))
  correctedcelltypes <- correctedcelltypes[!is.na(correctedcelltypes) & correctedcelltypes != "Unknown"]
}

correctedcelltypes <- as.character(correctedcelltypes)
correctedcelltypes <- correctedcelltypes[
  !is.na(correctedcelltypes) &
    correctedcelltypes != "Unknown" &
    correctedcelltypes %in% unique(as.character(datan$Phenotype_refined))
]

if (length(correctedcelltypes) == 0) {
  stop("No usable cell types found in `correctedcelltypes` or `datan$Phenotype_refined`.", call. = FALSE)
}

# Generate the _pos columns in data.
for (type in 1:length(correctedcelltypes)){
  datan[, paste0(correctedcelltypes[type],"_pos") := 0]
  datan[Phenotype_refined==correctedcelltypes[type], paste0(correctedcelltypes[type],"_pos") := 1]
}
# All core clustering distances
cap <- 75

# Enter target cells (major corrected cell types plus additional _pos markers)

targetpositives <- c(correctedcelltypes)

PointDatacols <- c("ROI", "jm", "im", "ci", paste0(targetpositives,"_pos"))
PointDatalocation <- c("jm","im","ci")

ROIsourcedata <- copy(datan[,..PointDatacols])
NNtarget <- 5
count <- 0
ROIs <- unique(datan$ROI)
for (ROIImage in ROIs)
{
  count <- progresstrack(ROIImage, ROIs, count)
  # Extract ROI data
  ROIsubset <- ROIsourcedata[ROI==ROIImage]
  
  # Make dummy cells
  ROIdummy <- copy(ROIsubset[1:(NNtarget+1)])
  ROIdummy[,ci := 1000000:(1000000+NNtarget)]
  ROIdummy[, which(colnames(ROIdummy) %like% "_pos") := 1]
  ROIdummy[, which(colnames(ROIdummy) %in% c("im","jm")) := 1000000:(1000000+NNtarget)]
  
  ROIsubset <- rbind(ROIsubset, ROIdummy)
  
  # For each cell type i
  X <- ppp(x = ROIsubset[,jm],
           y = ROIsubset[,im],
           window = owin(c(ROIsubset[,min(jm)],
                           ROIsubset[,max(jm)]),
                         c(ROIsubset[,min(im)],
                           ROIsubset[,max(im)])),
           # marks = ROIsubset[,..PointDatalocation],
           unitname = "micrometer")
  # Prepare dataset for cell subset
  
  NNsubi <- ROIsubset[,..PointDatalocation]
  
  Xcentroids <- vector("list",length(targetpositives))
  Ycentroids <- vector("list",length(targetpositives))
  NNcentroiddistances <- vector("list", length(targetpositives))
  
  # Compare to each cell type j
  for(j in targetpositives){
    # Generate points
    index_j <- ROIsubset[[paste0(j,"_pos")]]==1
    
    # Y <- ppp(x = ROIsubset[index_j,jm],
    #          y = ROIsubset[index_j,im],
    #          window = owin(c(ROIsubset[,min(jm)],
    #                          ROIsubset[,max(jm)]),
    #                        c(ROIsubset[,min(im)],
    #                          ROIsubset[,max(im)])),
    #          # marks = ROIsubset[index_j,..PointDatalocation],
    #          unitname = "micrometer")
    Y <- X[which(index_j==T)]
    
    # Find spatial nearest neighbors of i in j
    NNdataj <- as.data.table(nncross(X,Y,k=1:NNtarget,
                                     iX=as.integer(ROIsubset$ci),
                                     iY=as.integer(ROIsubset[index_j]$ci),
                                     #is.sorted.X=T, is.sorted.Y=T, sortby="y"
    ))
    
    # Find cells by object number instead of index
    NNdataj[, which(names(NNdataj) %like% "which") := lapply(.SD, function(x) ROIsubset[index_j][x,ci]), .SDcols = names(NNdataj) %like% "which"]
    
    # Name the columns for cell type j
    colnames(NNdataj) <- paste0(j, ".", colnames(NNdataj))
    
    # Store data
    NNsubi <- cbind(NNsubi, NNdataj)
    
    coltargets <- names(NNdataj)[names(NNdataj) %like% ".which"]
    NNminiX <- NNdataj[,..coltargets]
    NNminiY <- NNdataj[,..coltargets]
    NNminiX[,(coltargets) := lapply(.SD, function(x) ROIsubset[match(x, ci),jm]), .SDcols=coltargets]
    NNminiY[,(coltargets) := lapply(.SD, function(x) ROIsubset[match(x, ci),im]), .SDcols=coltargets]
    
    # Calculate X,Y centroid for each cell type
    Xcentroids[[which(targetpositives==j)]] <- NNminiX[,rowMeans(NNminiX, na.rm=T)]
    
    Ycentroids[[which(targetpositives==j)]] <- NNminiY[,rowMeans(NNminiY, na.rm=T)]
    
    # Calculate distance from each cell to each centroid for each cell type
    NNcentroiddistances[[which(targetpositives==j)]] <- sqrt((Xcentroids[[which(targetpositives==j)]]-NNsubi$jm)^2 + (Ycentroids[[which(targetpositives==j)]]-NNsubi$im)^2)
  }

  # Store new distance data for centroids
  newnncentroidcols <- apply(expand.grid(c("NN.centroid.X.","NN.centroid.Y.","NN.centroiddist."), targetpositives), 1, paste, collapse="")
  
  idx <- order(c(seq_along(Xcentroids), seq_along(Ycentroids), seq_along(NNcentroiddistances))) 
  
  NNsubi[, (newnncentroidcols) := (c(Xcentroids, Ycentroids, NNcentroiddistances))[idx]]
  
  # Substitute spatial calculations back in to ROI subset data
  spatialcolumns <- names(NNsubi)[(names(NNsubi) %like% "dist") | (names(NNsubi) %like% "NN")]
  spatialdata <- NNsubi[,..spatialcolumns]
  spatialdata <- spatialdata[1:(nrow(spatialdata)-NNtarget-1)]
  
  spatialdata <- censor(spatialdata, cap=cap)/cap
  
  for (target in targetpositives){
    regexin <- paste0("^",target,"[.]")
    distancecols <- colnames(spatialdata)[grep(regexin, colnames(spatialdata))]
    
    spatialdata[, paste0(target, ".density") := (NNtarget-rowSums(.SD))/NNtarget, .SDcols=distancecols]
  }
  
  
  rmcols <- c(grep("centroid\\.", colnames(spatialdata), value = T),grep("\\.dist(?!.1)",colnames(spatialdata), value = T, perl=T))
  
  spatialdata[, (rmcols) := NULL]
  
  datan[ROI==ROIImage, (colnames(spatialdata)) := spatialdata]
}


corespatialmarkers <- c("CD4.density","CD8.density","Treg.density","B.density","Mac.density","Myel.density",
                        "Epithelial.density", "Cancer.density","Fibro.density","Endo.density") 

set.seed(1236)

spatialclusters <- kmeans(datan[,..corespatialmarkers], 10) 
datan[, spatialcluster := paste0("TME_spatial_", spatialclusters$cluster)]
spcluster <- heatmapwrap(datan, corespatialmarkers,"spatialcluster",cluster_cols = FALSE, cluster_rows = TRUE, display_numbers = F, fontsize=14)
ggsave(filename = file.path("C:/Users/tmiao/Desktop/NEN/R/Figures/", "spcluster_refined.tiff"),
       plot = spcluster, width = 8, height = 6, dpi = 300)