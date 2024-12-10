PredictorResponseBivarLevels <- function(pred.resp.df, Z = NULL, qs = c(0.25, 0.5, 0.75), both_pairs = TRUE, z.names = NULL) {
  var.pairs <- dplyr::distinct(dplyr::select_(pred.resp.df, ~variable1, ~variable2))
  if (both_pairs) {
    var.pairs.rev <- dplyr::tibble(
      variable1 = var.pairs$variable2,
      
      variable2 = var.pairs$variable1
    )
    var.pairs <- rbind(var.pairs, var.pairs.rev)
  }
  
  if (is.null(z.names)) {
    z.names <- colnames(Z)
    if (is.null(z.names)) {
      z.names <- paste0("z", 1:ncol(Z))
      colnames(Z) <- z.names
    }
  }
  
  df <- data.frame()
  for (i in 1:nrow(var.pairs)) {
    var1 <- as.character(unlist(var.pairs[i, "variable1"]))
    var2 <- as.character(unlist(var.pairs[i, "variable2"]))
    preds <- pred.resp.df[pred.resp.df$variable1 == var1 & pred.resp.df$variable2 == var2, ]
    if (nrow(preds) == 0) {
      preds <- pred.resp.df[pred.resp.df$variable1 == var2 & pred.resp.df$variable2 == var1, ]
      preds.rev <- dplyr::tibble(
        variable1 = preds$variable2,
        variable2 = preds$variable1,
        z1 = preds$z2,
        z2 = preds$z1,
        est = preds$est,
        se = preds$se
      )
      preds <- preds.rev
      preds <- dplyr::arrange_(preds, ~z2, ~z1)
    }
    
    ngrid <- sqrt(nrow(preds))
    preds.plot <- preds$est
    se.plot <- preds$se
    
    hgrid <- matrix(preds.plot, ngrid, ngrid)
    se.grid <- matrix(se.plot, ngrid, ngrid)
    z1 <- preds$z1[1:ngrid]
    z2 <- preds$z2[seq(1, by = ngrid, length.out = ngrid)]
    
    quants <- quantile(Z[, var2], qs)
    
    ## relation of z1 with outcome at different levels of z2
    se.grid.sub <- hgrid.sub <- matrix(NA, ngrid, length(qs))
    for (k in seq_along(quants)) {
      sub.sel <- which.min(abs(z2 - quants[k]))
      hgrid.sub[, k] <- hgrid[, sub.sel]
      se.grid.sub[, k] <- se.grid[, sub.sel]
    }
    colnames(hgrid.sub) <- colnames(se.grid.sub) <- paste0("q", seq_along(qs))
    hgrid.df <- tidyr::gather(data.frame(hgrid.sub), quantile, 'est', convert = TRUE)
    se.grid.df <- tidyr::gather(data.frame(se.grid.sub), quantile, 'se')
    
    df.curr <- data.frame(variable1 = var1, variable2 = var2, z1 = z1, quantile = factor(hgrid.df$quantile, labels = qs), est = hgrid.df$est, se = se.grid.df$se, stringsAsFactors = FALSE)
    df <- rbind(df, df.curr)
  }
  df <- dplyr::tbl_df(df) %>%
    dplyr::arrange_(~variable1, ~variable2)
  df
}
