basicForPlots <- function (dat, types = c(), align = TRUE, binsize = "10K", cohort, bac, sample) 
{
  if ((length(types) == 3) && (align == TRUE)) {
    mat = t(matrix(c(3, 3, 1, 2), nrow = 2))
    layout(mat)
    align = FALSE
  }
  forbin = sprintf("for%s", binsize)
  limfor = quantile(dat[[forbin]], 0.999)
  if (is.element(1, types)) {
    if (align) {
      par(mfcol = c(1, 1))
    }
    par(mar = c(4, 4, 3, 3))
    dfor = density(dat[[forbin]], to = limfor)
    plot(dfor, col = 0, main = sprintf("%s Bins: Read Density",binsize))
    lines(dfor, col = 1, lw = 2)
  }
  if (is.element(2, types)) {
    if (align) {
      par(mfcol = c(1, 1))
    }
    par(mar = c(4, 4, 3, 3))
    plot(density(abs(dat$forw[, 2]), from = quantile(abs(dat$forw[,2]), 0.001),to = quantile(abs(dat$forw[,2]), 0.999)),
         col = 0, main = "Fragment Length Density")
    lines(density(dat$forw[, 2], from = 1, to = quantile(dat$forw[,2], 0.999)), 
          col = 1, lw = 2)
  }
  if (is.element(3, types)) {
    if (align) {
      par(mfcol = c(2, 1))
    }
    par(mar = c(2, 2, 2, 2))
    plot(dat[[forbin]], xlab = "Genome Location", 
         main = sprintf("%s | %s | %s || Read counts in %s bins",cohort,bac,sample,binsize), 
         cex = 0.3, pch = 20,col = rgb(0, 0, 0, 0.3), 
         #ylim = c(0, limfor), 
         ylim = c(0, 500))
    outli = which(dat[[forbin]] > limfor)
    points(outli, rep(limfor, length(outli)), col = 4, cex = 0.6)
  }
}


compDataForRefPlots <- function (chr, dat, types, binsize = "10K", align = TRUE) {
  if ((length(types) >= 3) && (align == TRUE)) {
    mat = matrix(c(1, 2, 3, 1, 2, 3, 1, 2, 3, 0, 4, 5, 0, 
                   4, 5), nrow = 3)
    layout(mat)
    align = FALSE
  }
  forbin = sprintf("for%s", binsize)
  mapbin = sprintf("map%s", binsize)
  gcbin = sprintf("gc%s", binsize)
  limboth = quantile(dat[[forbin]], 0.997) + 20
  outli = which(dat[[forbin]] > limboth)
  if (is.element(1, types)) {
    if (align) {
      par(mfcol = c(3, 1))
    }
    par(mar = c(2, 2, 2, 2))
    plot(dat[[forbin]] + dat[[revbin]], xlab = "Genome Location", 
         main = sprintf("Read counts in %s bins", binsize), 
         cex = 0.3, pch = 20, ylim = c(0, limboth))
    points(outli, rep(limboth, length(outli)), col = 4, cex = 0.6)
    plot(1 - chr[[mapbin]], main = sprintf("Mappability %s",binsize), cex = 0.3, pch = 20)
    plot(chr[[gcbin]], main = sprintf("GC %s Bins", binsize),cex = 0.3, pch = 20)
  }
  if (is.element(2, types)) {
    if (align) {
      par(mfcol = c(1, 1))
    }
    par(mar = c(4, 4, 3, 3))
    maxsizebin = min(c(length(dat[[forbin]]), 30000))
    sampbin = setdiff(sample(length(dat[[forbin]]), maxsizebin), 
                      outli)
    plot(1 - chr[[mapbin]][sampbin], (dat[[forbin]] + dat[[revbin]])[sampbin], 
         pch = 20, cex = 1, col = rgb(0, 0, 0, 0.2), main = "Map vs. Counts bins", 
         ylab = "Counts", xlab = "% Mappable", ylim = c(0, 
                                                        limboth))
    abline(0, median((dat[[forbin]] + dat[[revbin]])[chr[[mapbin]][sampbin] > 
                                                       0.99]), col = 4, lw = 2)
  }
  if (is.element(3, types)) {
    fowardUnmap = mean(chr$isrep[dat$forw[, 1]])
    reverseUnmap = mean(chr$isrep[dat$reve[, 1]])
    Unmap = mean(chr$isrep)
    cat(sprintf("Ratio between unmappable bases and unmappable reads rate %g (Forward Strand)\n", 
                Unmap/fowardUnmap))
    cat(sprintf("Ratio between unmappable bases and unmappable reads rate %g (Reverse Strand)\n", 
                Unmap/reverseUnmap))
  }
  if (is.element(4, types)) {
    if (align) {
      par(mfcol = c(1, 1))
    }
    par(mar = c(4, 4, 3, 3))
    map_cutoff = 0.1
    maxsizeb = min(sum(chr[[mapbin]] < map_cutoff), 8000)
    sampbinb = setdiff(sample(which(chr[[mapbin]] < map_cutoff), 
                              maxsizeb), outli)
    plot(chr[[gcbin]][sampbinb], (dat[[forbin]] + dat[[revbin]])[sampbinb] + 
           rnorm(length(sampbinb), 0, 0.3), pch = 20, cex = 0.6, 
         col = rgb(0, 0, 0, 0.3), main = "GC effect", ylab = "Counts", 
         xlab = "GC", ylim = c(0, limboth))
    lines(smooth.spline(chr[[gcbin]][sampbinb], (dat[[forbin]] + 
                                                   dat[[revbin]])[sampbinb], spar = 1), col = 4, lw = 2)
  }
}
###################
makeGCL2 <- function (gcline, locs_leng, minlens = 5, maxlens = 400,
		sampline, margin = 1, max_frag_for_loc = 20) 
{
	maxfraglenC = 700
		stopifnot(maxlens < maxfraglenC)
		nlocs = numeric(maxfraglenC^2)
		nfrag = numeric(maxfraglenC^2)
		tmpgc = gcline
		tmpgc[is.na(gcline)] = FALSE
		sampline_vec = logical(length(gcline))
		sampline_vec[sampline] = TRUE
		sampline_vec[is.na(gcline)] = FALSE
		res = .C("makeGCL2", nlocs = as.integer(nlocs), 
				nfrag = as.integer(nfrag), as.integer(tmpgc), as.integer(length(tmpgc)), 
				as.integer(locs_leng[, 1]), as.integer(nrow(locs_leng)), 
				as.integer(sampline_vec), as.integer(minlens), as.integer(maxlens), 
				as.integer(margin), as.integer(max_frag_for_loc))
		fragmat = matrix(res$nfrag, maxfraglenC, maxfraglenC)
		locsmat = matrix(res$nlocs, maxfraglenC, maxfraglenC)
		rownames(locsmat) = 0:(maxfraglenC - 1)
		rownames(fragmat) = 0:(maxfraglenC - 1)
		colnames(locsmat) = 0:(maxfraglenC - 1)
		colnames(fragmat) = 0:(maxfraglenC - 1)
		return(list(frag = fragmat, locs = locsmat))
}
###################
c.plotCondMean <- function (cMean, ci = FALSE, newplot = TRUE, normRange = 0, meanLine = FALSE, 
		lt = 1,cohort,bac,sample, ...) 
{
	if (cMean$norm) {
		ci_corr = cMean$tot_mean
			dat_mean = 1
			ylabel = "Normalized Rate"
	}
	else {
		ci_corr = 1
			dat_mean = cMean$tot_mean
			ylabel = "Rate"
	}
	xvals = as.numeric(names(cMean$tab_means))
	if (normRange > 0) {
		xvals = xvals/normRange
	}
	if (newplot) {
#		plot(xvals, cMean$tab_means, col = 0, main = sprintf("%s | %s | %s \nTV = %g",cohort,bac,sample,cMean[[3]]), 
		plot(xvals, cMean$tab_means, col = 0,
				ylim = c(0, max(cMean$tab_means)), xlab = "GC", ylab = ylabel)
	}
	lines(xvals, cMean$tab_means, lwd=2, ...)
	if (meanLine) {
		abline(dat_mean, 0, lt = 2)
	}
	if (ci) {
		tot_mean = cMean$tot_mean
			segments(xvals, cMean$tab_means - (sqrt(cMean$tab_means)/sqrt(cMean$cts *ci_corr)), 
					xvals, cMean$tab_means + (sqrt(cMean$tab_means)/sqrt(cMean$cts * ci_corr)),col=2)
	}
	lines(xvals, cMean$tab_means, lwd=2, ...)
	title(main=sprintf("%s | %s | %s \nTV = %g",cohort,bac,sample,cMean$TV))
}

