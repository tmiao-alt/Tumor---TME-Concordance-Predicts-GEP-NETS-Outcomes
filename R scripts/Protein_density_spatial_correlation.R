# -- helpers -----------------------------------------------------------------
.select_cols <- function(dt, cols) {
  if (length(cols) == 2L && all(cols %in% names(dt)) && match(cols[1], names(dt)) <= match(cols[2], names(dt))) {
    i1 <- match(cols[1], names(dt)); i2 <- match(cols[2], names(dt))
    return(names(dt)[i1:i2])
  }
  
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("Columns not found: ", paste(miss, collapse=", "))
  cols
}

.numeric_only <- function(dt, cols) {
  num <- cols[sapply(cols, function(x) is.numeric(dt[[x]]))]
  dropped <- setdiff(cols, num)
  if (length(dropped)) warning("Dropping non-numeric columns: ", paste(dropped, collapse=", "))
  num
}

# core worker for a single x/y vector; 
.cor_worker <- function(x, y, method, slope_type) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3L) {
    return(list(n=n, cor=NA_real_, p_value=NA_real_, conf_low=NA_real_, conf_high=NA_real_,
                K=NA_real_, b=NA_real_, r_squared=NA_real_, trend="none"))
  }
  x <- x[ok]; y <- y[ok]
  
  # correlation & test
  ct <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE, conf.level = 0.95))
  r  <- unname(ct$estimate); p <- unname(ct$p.value)
  ci <- if (!is.null(ct$conf.int)) unname(ct$conf.int) else c(NA_real_, NA_real_)
  
  # slope options
  if (slope_type == "rank") {
    xr <- rank(x, ties.method = "average"); yr <- rank(y, ties.method = "average")
    fit <- lm(yr ~ xr)
  } else { # "ols"
    fit <- lm(y ~ x)
  }
  co <- coef(fit)
  K  <- unname(co[2]); b <- unname(co[1]); r2 <- summary(fit)$r.squared
  
  list(n=n, cor=r, p_value=p, conf_low=ci[1], conf_high=ci[2],
       K=K, b=b, r_squared=r2, trend = ifelse(is.na(r), "none", ifelse(r>0,"increasing", ifelse(r<0,"decreasing","none"))))
}

# -- main --------------------------------------------------------------------
correlate_dt2 <- function(
    dt,
    cols,                               # vector of names or c(start, end)
    method = c("pearson","spearman"),
    by = character(0),                   # e.g., c("ROI") or c("case","origin")
    mode = c("pairwise","target"),
    target = NULL,                       # string or vector when mode="target"
    na_method = c("pairwise","complete"),
    slope_type = c("ols","rank"),        # "rank" = slope on ranks
    p_adjust = c("none","BH","bonferroni","holm","BY","hochberg","hommel")
) {
  stopifnot(is.data.table(dt))
  method    <- match.arg(method)
  mode      <- match.arg(mode)
  na_method <- match.arg(na_method)
  slope_type<- match.arg(slope_type)
  p_adjust  <- match.arg(p_adjust)
  
  use_arg <- if (na_method == "pairwise") "pairwise.complete.obs" else "complete.obs"
  
  sel  <- .select_cols(dt, cols)
  sel  <- .numeric_only(dt, sel)
  if (length(sel) < 2L) stop("Need at least two numeric columns after filtering.")
  
  # build (x,y) pairs
  if (mode == "target") {
    if (is.null(target)) stop("Provide target when mode='target'.")
    targets <- intersect(if (length(target)==2L && all(target %in% names(dt))) .select_cols(dt, target) else target, sel)
    if (!length(targets)) stop("No valid targets among selected columns.")
    xs <- setdiff(sel, targets)
    combs <- CJ(var_x = xs, var_y = targets, unique = TRUE)
  } else {
    cm <- utils::combn(sel, 2); combs <- data.table(var_x = cm[1,], var_y = cm[2,])
  }
  
  run_group <- function(gdt) {
    combs[, {
      w <- .cor_worker(gdt[[var_x]], gdt[[var_y]], method = method, slope_type = slope_type)
      c(list(var_x=var_x, var_y=var_y), w)
    }, by = seq_len(nrow(combs))][, seq_len := NULL][]
  }
  
  if (length(by)) {
    miss_by <- setdiff(by, names(dt))
    if (length(miss_by)) stop("Grouping columns not found: ", paste(miss_by, collapse=", "))
    res <- dt[, run_group(.SD), by = by]
  } else {
    res <- run_group(dt)
  }
  
  if (length(by)) {
    data.table::setorderv(res, c(by, "var_x", "var_y"))
  } else {
    data.table::setorder(res, "var_x", "var_y")
  }
  # adjust p-values within each group (if grouped), else globally
  if (length(by)) {
    res[, p_adj := p.adjust(p_value, method = p_adjust), by = by]
  } else {
    res[, p_adj := p.adjust(p_value, method = p_adjust)]
  }
  
  # metadata
  res[, `:=`(method = method, slope_type = slope_type, use = use_arg, adjust = p_adjust)]
  setcolorder(res, intersect(c(by, "method","use","adjust","slope_type","var_x","var_y",
                               "n","cor","p_value","p_adj","conf_low","conf_high","K","b","r_squared","trend"), names(res)))
  if (length(by)) {
    data.table::setorderv(res, c(by, "var_x", "var_y"))
  } else {
    data.table::setorderv(res, c("var_x", "var_y"))
  }
  res[]
}