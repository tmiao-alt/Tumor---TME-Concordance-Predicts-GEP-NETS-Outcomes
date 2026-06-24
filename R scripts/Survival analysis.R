biomarkergroups <- list(
  list(c(quote(origin == "Pancreas")), "Pan_origin",  "datan",  "case"),
  list(c(quote(origin == "Ileum")),    "SB_origin",  "datan",  "case"),
  list(c(quote(origin == "Pancreas" & tissue == "L")),        "P_Lmet",     "datan",  "case"),
  list(c(quote(origin == "Ileum" & tissue == "L")),        "SB_Lmet",     "datan",  "case"),
  list(c(quote(origin == "Ileum" & tissue == "T")),        "SB_T",     "datan",  "case"),
  list(c(quote(origin == "Pancreas" & tissue == "T")),        "P_T",     "datan",  "case"),
  list(c(quote(origin == "Ileum" & tissue == "LN")),       "SB_LNmet",    "datan",  "case"),
  list(c(quote(origin == "Pancreas" & tissue == "LN")),       "P_LNmet",    "datan",  "case"),
  list(c(quote(tissue == "T")),        "tumor",    "datan",  "case"),
  list(c(quote(tissue == "L")),        "Lmet",    "datan",  "case"),
  list(c(quote(tissue == "LN")),        "LNmet",    "datan",  "case")
)


for(condition in biomarkergroups){
  survivalcondition <- condition[[2]]
  
  survivalvars <- c("sex","origin","tissue","OS","PFS","status")
  
  glmvars <- c("sex","origin","tissue","status")
  
  cellselectindex <- 1:nrow(datan)
  
  groupvariable <- condition[[4]]
  #length(unique(paste0(origin, get(groupvariable))))/length(unique(get(groupvariable)))
  if("metadata" %in% condition[[3]]){
    metaselectsubset <- metadata[eval(condition[[1]][[match("metadata", condition[[3]])]]), ROI]
    #condition[[1]][[1]], list(list("a","b","c"),list(1,2,3))
    cellselectindex <- intersect(cellselectindex, datan[ROI %in% metaselectsubset, which=T])  
  }  
  
  if ("Primarytumor" %in% condition[[3]]) {
    pt_cases <- Primarytumor[eval(condition[[1]][[match("Primarytumor", condition[[3]])]]), ROI]
    datselectindex <- datan[ROI %in% pt_cases, which = TRUE]
    cellselectindex <- intersect(cellselectindex, datselectindex)
  }
  # 
  # if ("N_removed" %in% condition[[3]]) {
  # Nrm_cases <- N_removed[eval(condition[[1]][[match("N_removed", condition[[3]])]]), ROI]
  # datselectindex <- datan[ROI %in% Nrm_cases, which = TRUE]
  # cellselectindex <- intersect(cellselectindex, datselectindex)
  # }
  # 
  if("datan" %in% condition[[3]]){
    datselectindex <- datan[eval(condition[[1]][[match("datan", condition[[3]])]]), which=T]
    cellselectindex <- intersect(cellselectindex, datselectindex)
  }
  
  subdatacols <- c(
    "ROI",
    "Object",
    "spatialcluster",
    "case",
    "Endoarea",
    "Endocluster",
    "patient_endoarea_concave",
    "ROI_Endoarea_concave",
    "patient_endoarea_convex",
    "EndothelialAreaPercent",
    "patient_endoarea_percent",
    "CD8.density",
    "CD4.density",
    "Mac.density",
    "M2mac.density",
    "Fibro.density",
    "Treg.density",
    "B.density",
    "Epithelial.density",
    "Cancer.density",
    "Endo.density", "IL6", "IL10", "PDL1", "TGFB", "Ki67"
  )
  subdata <- datan[cellselectindex, names(datan)%like%"pos"|names(datan)%in%subdatacols, with=F]
  interestedcell <- c("Fibro", "Mac", "CD4", "CD8", "Treg", "B", "Endo")
  
  # Compute mean density per cell type per ROI
  density_cols   <- paste0(interestedcell, ".density")   # e.g. "Fibro.density"
  mean_col_names <- paste0(interestedcell, "_meandensity")
  
  # Step 1: Filter to Cancer cells, compute mean density per ROI for each cell type
  cancer_means <- subdata[Cancer_pos == 1,
                          lapply(.SD, mean, na.rm = TRUE),
                          by = ROI,
                          .SDcols = density_cols]
  setnames(cancer_means, density_cols, mean_col_names)
  
  subdata <- cancer_means[subdata, on = "ROI"]
  subdata_archive <- if (!exists("subdata_archive")) list() else subdata_archive
  subdata_archive[[survivalcondition]] <- subdata
 
  survivaldat <- metadata[ROI %in% unique(subdata$ROI),
                          lapply(.SD, function(x) unique(na.omit(x))[1]),
                          .SDcols = c(union(survivalvars, glmvars)),
                          by = c(groupvariable)]
  survivaldat <- merge(survivaldat, subdata[, .(Cells = .N), by = groupvariable], by = groupvariable)
  
  # merge correlation results to survivaldat
  # ── Per-condition, per-patient correlation ─────────────────────────────────
  corr_condition_cols <- character(0)
  
  corr_cols_for_survival <- c("CD8.density", "Mac.density", "CD4.density", "Fibro.density",
                              "Treg.density",
                              "B.density",
                              "PDL1", "TGFB", "IL6", "IL10", "Ki67")
  corr_targets           <- c("PDL1", "TGFB", "IL6", "IL10", "Ki67")
  
  cancer_subdata <- subdata[Cancer_pos == 1]
  
  missing_corr_cols <- setdiff(corr_cols_for_survival, names(cancer_subdata))
  
  if (nrow(cancer_subdata) >= 3L && length(missing_corr_cols) == 0L) {
    
    correlation_condition <- tryCatch(
      correlate_dt2(
        cancer_subdata,
        cols       = corr_cols_for_survival,
        by         = groupvariable,        # "case" → one correlation per patient
        method     = "spearman",
        mode       = "target",
        target     = corr_targets,
        slope_type = "ols",
        p_adjust   = "BH"
      ),
      error = function(e) {
        warning("Correlation failed for [", survivalcondition, "]: ", e$message)
        NULL
      }
    )
    
    if (!is.null(correlation_condition) && nrow(correlation_condition) > 0L) {
      
      correlation_condition[, pair_id := paste(var_x, var_y, sep = "_")]
      correlation_condition[, corr_strata := fifelse(
        p_adj <= 0.05 & cor > 0, "pos_sig",
        fifelse(p_adj <= 0.05 & cor < 0, "neg_sig", "not_sig")
      )]
      
      dup_check <- correlation_condition[, .N, by = c(groupvariable, "pair_id")]
      if (any(dup_check$N > 1)) {
        stop("Duplicate groupvariable + pair_id combinations found in correlation_condition")
      }
      corr_wide <- dcast(
        correlation_condition,
        as.formula(paste(groupvariable, "~ pair_id")),
        value.var     = "corr_strata"
      )
      
      survivaldat       <- merge(survivaldat, corr_wide, by = groupvariable, all.x = TRUE)
      corr_condition_cols <- setdiff(names(corr_wide), groupvariable)
      
      
    }
    
  } else {
    if (length(missing_corr_cols) > 0L)
      warning("[", survivalcondition, "] Skipping correlation — missing columns: ",
              paste(missing_corr_cols, collapse = ", "))
    else
      warning("[", survivalcondition, "] Skipping correlation — fewer than 3 cancer cells after filtering.")
  }
  
  if (length(corr_condition_cols)) {
    survivaldat[, (corr_condition_cols) := lapply(.SD, function(x) {
      vals <- as.character(x)
      #vals[is.na(vals)] <- "not_sig"
      factor(vals, levels = c("neg_sig", "not_sig", "pos_sig"))
    }), .SDcols = corr_condition_cols]
  }
  
  cellsubtypes <- grep("_pos$", colnames(datan), value = T)
  cellsubtypes <- sub("_pos", "", cellsubtypes)
  
  for (i in cellsubtypes){
    cell_counts <- subdata[get(paste0(i, "_pos")) == 1, .N, by = groupvariable]
    survivaldat[, paste0(i, "_Cells") := 0L]
    survivaldat[cell_counts, on = groupvariable, (paste0(i, "_Cells")) := i.N]
    survivaldat[, paste0(i, "_Percent") := get(paste0(i, "_Cells"))/get("Cells")]
  }
  
  for (i in unique(subdata$spatialcluster)){
    cell_counts <- subdata[spatialcluster == i, .N, by = groupvariable]
    survivaldat[, paste0(i, "_Cells") := 0L]
    survivaldat[cell_counts, on = groupvariable, (paste0(i, "_Cells")) := i.N]
    survivaldat[, paste0(i, "_Percent") := get(paste0(i, "_Cells"))/get("Cells")]
  }
  
 
  density_features <- c("CD8.density", "CD4.density", "Mac.density", "Fibro.density", "Epithelial.density","Cancer.density", "Endo.density")
  if (all(density_features %in% names(subdata))) {

    spatial7_counts <- subdata[
      Cancer.density > 0.45 & Fibro.density > 0.25 & (Mac.density + CD8.density > 0.45) & CD4.density < 0.1,
      .N,
      by = groupvariable
    ]

    spatial4_counts <- subdata[
      Cancer.density >0.7 & Fibro.density>0.6 & Endo.density > 0.5,
      .N,
      by = groupvariable
    ]
    survivaldat[, spatialbiomarker7_Cells := 0L]
    survivaldat[spatial7_counts, on = groupvariable, spatialbiomarker7_Cells := i.N]
    survivaldat[, spatialbiomarker7_Percent := spatialbiomarker7_Cells/Cells]
    
    
    survivaldat[, spatialbiomarker4_Cells := 0L]
    survivaldat[spatial4_counts, on = groupvariable, spatialbiomarker4_Cells := i.N]
    survivaldat[, spatialbiomarker4_Percent := spatialbiomarker4_Cells/Cells]
  }
  
  
  # survivaldat <- merge(survivaldat,
  #     subdata[, .(patient_endoarea_concave = patient_endoarea_concave[1L]),by = groupvariable],
  #     by= groupvariable, all.x=TRUE)
  survivaldat <- merge(survivaldat,
                       subdata[, .(patient_endoarea_percent = patient_endoarea_percent[1L]),by = groupvariable],
                       by= groupvariable, all.x=TRUE)
  
  #add mean density of interested cells around cancer cells to survivaldat
  meandensity_patient <- unique(subdata[, c("ROI", groupvariable, mean_col_names), with = FALSE])[
    , lapply(.SD, mean, na.rm = TRUE), by = groupvariable, .SDcols = mean_col_names]
  survivaldat <- merge(survivaldat, meandensity_patient, by = groupvariable, all.x = TRUE)
  
  # ── Cancer cell immunomodulatory expression ──────────────────────────────
  cancer_immuno_markers   <- c("IL6", "IL10", "PDL1", "TGFB", "Ki67")
  cancer_immuno_col_names <- paste0("Cancer_", cancer_immuno_markers, "_mean")
  
  if (all(cancer_immuno_markers %in% names(subdata))) {
    
    # ROI-level mean for cancer cells only
    cancer_immuno_roi <- subdata[
      Cancer_pos == 1,
      lapply(.SD, mean, na.rm = TRUE),
      by = c("ROI", groupvariable),
      .SDcols = cancer_immuno_markers
    ]
    setnames(cancer_immuno_roi, cancer_immuno_markers, cancer_immuno_col_names)
    
    # aggregate ROI means to patient level
    cancer_immuno_patient <- cancer_immuno_roi[
      , lapply(.SD, mean, na.rm = TRUE),
      by = groupvariable,
      .SDcols = cancer_immuno_col_names
    ]
    
    survivaldat <- merge(survivaldat, cancer_immuno_patient, by = groupvariable, all.x = TRUE)
    
  } else {
    cancer_immuno_col_names <- character(0)   # graceful fallback
    warning("One or more of cytokines not found in subdata; skipping cancer immuno stratification.")
  }

  
  survivaldat[is.nan(survivaldat)] <- 0
  cellcols    <- grep("Cells$",   names(survivaldat), value = TRUE)
  percentcols <- grep("Percent$", names(survivaldat), value = TRUE)
  areacols    <- "patient_endoarea_percent"   
  
  survivaldat[, (cellcols) := NULL]
  if (length(percentcols)) survivaldat[, (percentcols) := round(.SD, 4), .SDcols = percentcols]
  
  # threshold BOTH Percent and Area features
  contvars <- c(percentcols, areacols, mean_col_names, cancer_immuno_col_names)                                
  for (i in contvars) {
    classif <- as.integer(survivaldat[[i]] > median(survivaldat[[i]], na.rm = TRUE))
    survivaldat[, paste0(i, "_Median_Thres") := classif]
    classif <- as.integer(survivaldat[[i]] > mean(survivaldat[[i]], na.rm = TRUE))
    survivaldat[, paste0(i, "_Mean_Thres") := classif]
  }

  conditions <- c(intersect(survivalvars, c("EBV","MHCI","MHCII","SEX","STAGEGR","PDL1CopyNum","PDL1Amp","EarlyRelapse","LateRelapse","Refractory")),grep("Thres", colnames(survivaldat), value=T))
  conditions <- conditions[(survivaldat[, lapply(.SD, function(x) length(unique(na.omit(x)))==1), .SDcols=conditions]!=1)]
  conditions <- unique(c(conditions, corr_condition_cols))
  
  codes <- list(c("OS","status"), c("PFS", "status"))
  
  codes <- codes[which(!sapply(codes, function(x) !(all(x %in% survivalvars))))]
  
  
  non_corr_conditions <- setdiff(conditions, corr_condition_cols)
  if (length(non_corr_conditions)) {
    survivaldat[, (non_corr_conditions) := lapply(.SD, function(x) factor(x)), .SDcols = non_corr_conditions]
  }
  
  
  glmconditions <- union(setdiff(grep("Percent",colnames(survivaldat), value=T), grep("Thres", colnames(survivaldat), value=T)), "patient_endoarea_percent" )
  #setdiff(x,y) find the elements which are in the first Object but not in the second Object.
  #union(x,y)takes two objects like Vectors, dataframes, etc. as arguments and results in a third object with the combination of the data of both the objects.
  glmconditions <- glmconditions[(survivaldat[, lapply(.SD, function(x) length(unique(na.omit(x)))==1), .SDcols=glmconditions]!=1)]
  
  coxresults <- data.table(index=1:(length(conditions)*length(codes)),
                           groups=vector("character",length(conditions)*length(codes)),
                           conditions=vector("character",length(conditions)*length(codes)),
                           codes=vector("character",length(conditions)*length(codes)),
                           pvalues=vector("numeric",length(conditions)*length(codes)))
  glmresults <- data.table(index=1:(length(glmconditions)*length(glmvars)),
                           groups=vector("character",length(glmconditions)*length(glmvars)),
                           conditions=vector("character",length(glmconditions)*length(glmvars)),
                           codes=vector("character",length(glmconditions)*length(glmvars)),
                           pvalues=vector("numeric",length(glmconditions)*length(glmvars)),
                           coef=vector("numeric",length(glmconditions)*length(glmvars)))
  survformulas <- list()
  glmobjects <- list()
  
  counter <- 1
  for (i in conditions){
    for (j in codes){
      timevar <- j[1]
      statusvar <- j[2]
      coxresults[counter, c("groups","conditions","codes","pvalues") := data.table(survivalcondition, i, j[1],NA)]
      if(!(all(is.na(survivaldat[!is.na(survivaldat[[i]])][[statusvar]])))){
        survformula <- as.formula(paste0("Surv(", j[1], ", ", j[2], ") ~ ", i))
        fit10 <- survfit(survformula, data = survivaldat[, c(j[1], j[2], i, groupvariable), with=F], id=survivaldat[[groupvariable]])
        fit10$call$formula <- survformula
        coxout <- coxph(survformula, data = survivaldat[, c(j[1], j[2], i), with=F], id=survivaldat[[groupvariable]])
        coef_est <- summary(coxout)$coefficients[ , "coef"]
        hr_est   <- exp(coef_est)
        coxresults[counter, `:=`(
          coef = list(coef_est),
          HR = list(hr_est),
          ratio = max(fit10$n)/sum(fit10$n)
        )]
        if(!is.na(summary(coxout)$waldtest[3])){
          coxresults[counter, c("pvalues") := list(summary(coxout)$waldtest[3])] 
          attr(survformula, ".Environment") <- NULL
          survformulas[[counter]] <- survformula
          survP <- ggsurvplot(
            fit10,                     # survfit object with calculated statistics.
            data = survivaldat,             # data used to fit survival curves.
            risk.table = TRUE,       # show risk table.
            pval = TRUE,             # show p-value of log-rank test.
            pval.size = 14,
            conf.int = FALSE,         # show confidence intervals for
            # point estimates of survival curves.
            xlab = "Months",   # customize X axis label.
            ylab = j,
            title = survivalcondition,
            break.time.by = 10,     # break X axis in time intervals by 500.
            ggtheme = theme_bw(base_size = 25), # customize plot and risk table with a theme.
            risk.table.y.text.col = T,# colour risk table text annotations.
            risk.table.height = 0.25, # the height of the risk table
            risk.table.y.text = FALSE,# show bars instead of names in text annotations
            # in legend of risk table.
            risk.table.fontsize = 7,
            ncensor.plot = FALSE,      # plot the number of censored subjects at time t
            ncensor.plot.height = 0.1,
            conf.int.style = "step",  # customize style of confidence intervals
            # surv.median.line = "hv",  # add the median survival pointer.
            tables.theme = theme_cleantable(),
            palette="nejm"
          )
          outfn <- sprintf("%s_%s_%s.pdf",
                           survivalcondition,
                           gsub("\\W+", "", i),
                           timevar)
          pdf(file.path(survPlocation, outfn),
              width = 10,
              height = 10)
          print(survP)
          dev.off()
        } else {
          coxresults[counter, c("pvalues") := list(NA)]
          survformulas[[counter]] <- NA
        }
      } else {
        coxresults[counter, c("pvalues") := list(NA)]
        survformulas[[counter]] <- NA
      }
      
      counter <- counter+1
      
    }
  }
  
  counter <- 1
  
  for (i in glmconditions){
    for (j in glmvars){
      glmresults[counter, c("groups","conditions","codes","pvalues","coef") := data.table(survivalcondition, i, j[1],NA,NA)]
      if(length(unique(survivaldat[[j]][!is.na(survivaldat[[i]])]))!=1){
        glmdat <- survivaldat[, c(i,j,groupvariable), with=F]
        glmdat <- glmdat[apply(glmdat[, c(i), with=F], 1, function(x) all(is.finite(x)))]
        glmout <- glm.cluster(as.formula(paste0("`",i,"` ~ `",j,"`")), data = glmdat, cluster=glmdat[[groupvariable]])
        if(!any(is.na(summary(glmout$glm_res)$coefficients[-1,4]))){
          glmresults[counter, c("pvalues","coef") := list(min(summary(glmout$glm_res)$coefficients[-1,4]),summary(glmout$glm_res)$coefficients[2,1])]
          if(any(summary(glmout$glm_res)$coefficients[-1,4]<0.05)){# | counter %in% interestingsamples){
            glmout$glm_res <- strip_glm_env(glmout$glm_res)
            glmobjects[[counter]] <- glmout
            # p <- eval(parse(text=paste0("ggplot(survivaldat, aes(x=factor(",
            #                             j,"), y=",i,"))+geom_boxplot()+labs(x='",j,"', y='",i,"')")))
          } else {
            glmobjects[[counter]] <- NA
          }
        }
      } else {
        glmresults[counter, c("pvalues","coef") := list(NA, NA)]
        glmobjects[[counter]] <- NA
      }
      
      counter <- counter+1
      
    }
  }
  
  coxresults[, adjp2 := p.adjust(pvalues,method="BH")]
  coxresults[, adjp := p.adjust(pvalues)]
  glmresults[, adjp := p.adjust(pvalues)]
  
  survivalbiomarkerexportlist <- list(coxresults=coxresults,survformulas=survformulas, glmresults=glmresults, glmobjects=glmobjects, survivaldat=survivaldat)
  survivalbiomarkerexportlist <- setNames(survivalbiomarkerexportlist, paste0(survivalcondition, names(survivalbiomarkerexportlist)))
  
  filetime <- format(Sys.time(), "%b%e_%H%M")
  
  save(list=names(survivalbiomarkerexportlist), envir = list2env(survivalbiomarkerexportlist), file=paste0(biomarkersavelocation, paste(survivalcondition, filetime, "survivalbiomarkerobjects.RData", sep="_")))

}