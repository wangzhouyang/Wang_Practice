

fastcount = function(xvar,yvar) {
  naline= is.na(xvar)
  naline[is.na(yvar)] = TRUE
  xvar[naline] = 0
  yvar[naline] = 0
  useline = !naline;
  # Table must be initialized for -1's
  tablex = numeric(max(xvar)+1)
  tabley = numeric(max(xvar)+1)
  stopifnot(length(xvar) == length(yvar))

  res = .C("fastcount",PACKAGE="GCcorrect",tablex = as.integer(tablex),tabley = as.integer(tabley),
    as.integer(xvar), as.integer(yvar),as.integer(useline),as.integer(length(xvar)))
  xuse = which(res$tablex>0)
  xnames = xuse - 1
  resb = rbind(res$tablex[xuse], res$tabley[xuse])  
  colnames(resb) = xnames
  return(resb)
}


fastcountR = function(xvar,yvar) {
  xnames= sort(unique(xvar))
  tab = matrix(0,nrow = 2, ncol =length(xnames))
  for (x in 1:length(xnames)){
    tmp = xvar==xnames[x]
    tab[1,x] = sum(tmp,na.rm = TRUE)    
    tab[2,x] = sum(yvar[tmp],na.rm= TRUE)
  }
  colnames(tab) = xnames
  return(tab)
}


compTV = function(frag,locs){
  grandmean = sum(frag)/sum(locs+1)
  TV = 0.5*sum((abs(frag/(locs+1) - grandmean) * (locs+1))/ sum(locs+1) )
  return(TV)
}

compSTV = function(frag,locs){
  grandmean = sum(frag)/sum(locs+1)
  TV = 0.5*sum((abs(frag/(locs+1) - grandmean) * (locs+1))/ sum(locs+1) )
  return(TV/grandmean)
}

####
#
# Estimate conditional mean
#
####
getCondMean = function(xvar,yvar,cutoff = -1,jump = 1,norm = FALSE){
  # We keep bin-names integer, but should shift xvar by (jump-1)/2
  # To account for bin centers
  center_shift = (jump-1)/2
  if (jump>1){
    xvar  = floor(xvar/jump)*jump  
  }
  if (cutoff>0){
    yvar[yvar>cutoff] = cutoff
  }
  tab = fastcount(xvar,yvar)

  cts = tab[1,]
  nreads = tab[2,]
  names(cts) = as.numeric(names(cts))+center_shift
  names(nreads) = as.numeric(names(cts))
  tab_means = nreads/cts
  tot_mean = sum(nreads)/sum(cts)
  TV = compTV(nreads,cts)

  if (norm) {
    tab_means = tab_means / tot_mean
    TV = TV / tot_mean
  }
  return(list(tab_means= tab_means,nreads = nreads,cts = cts,TV = TV,tot_mean = tot_mean,norm = norm))
}

sumCondMeans = function(cma, cmb) {
  if(cma$norm){
    readsa = cma$tab_means * cma$cts * cma$tot_mean
  } else {
    readsa = cma$tab_means * cma$cts
  }
  if(cmb$norm){
    readsb = cmb$tab_means * cmb$cts * cmb$tot_mean
  } else {
    readsb = cmb$tab_means * cmb$cts
  }
	# Make sure that they are aligned
  allnames = names(c(cma$cts,cmb$cts))
  allbins = sapply(allnames,as.numeric)
  bin_ord = order(allbins)
  
  sortcts = c(cma$cts,cmb$cts)[bin_ord]
  sortreads = c(readsa,readsb)[bin_ord]
  sortnames = allnames[bin_ord]

  rm_double = logical(length(bin_ord))
  for (i in 2:length(sortnames)){
    if (sortnames[i] == sortnames[i-1]) {
      rm_double[i] = TRUE
      sortcts[i-1] = sortcts[i-1] + sortcts[i]
      sortreads[i-1] = sortreads[i-1]+ sortreads[i]
    }
  }
  new_reads = sortreads[!rm_double]
  new_cts = sortcts[!rm_double]
  tab_means = new_reads/new_cts
  tot_mean = sum(new_reads)/sum(new_cts)
  TV = compTV(new_reads,new_cts)

  if (cma$norm) {
    tab_means = tab_means / tot_mean
    TV = TV / tot_mean
  }
  return(list(tab_means= tab_means,cts = new_cts,TV = TV,tot_mean = tot_mean,norm = cma$norm))
}

getTVScore = function(xvar,yvar,cutoff = -1,jump = 1,norm=FALSE){
  cm = getCondMean(xvar,yvar,cutoff ,jump,norm)
  return(cm$TV)
}

  
####
#
# Estimate table (conditional mean for each fragment length)
#
####

# Generate an intepolated version of condMean
smoothCondMean = function(tab_means, olocations, nlocations, span=0.2, zero_pad = c(), eps = 0.001,plot = FALSE){
  origpts = c(olocations,zero_pad)
  rates = c(tab_means,rep(0,length(zero_pad)))
  interp = loess(rates~origpts,span = span)
  
  predvec= predict(interp,nlocations)
  # Make sure everything positive, and threshold at max*epsilons.
  predvec[predvec< max(predvec,na.rm=TRUE)*eps] = 0
  if (plot){
    plot(nlocations,predvec,type = 'l')
    title(sprintf('Span: %f',span))
  }
  names(predvec) = nlocations
  return(predvec)
}

shift = function(sline,offset,fill = NA) {
  if (offset>0) {
    return(c(sline[(1 + offset):length(sline)],rep(fill,offset)))
  }
  else {
    return(c(rep(fill,abs(offset)),sline[1:(length(sline)+offset)])) 
  }
}


#makeFullTable = function(isgcline,locs_leng , minlens,maxlens, sampline,margin, max_frag_for_loc) {
#  maxfraglenC = 700 # Hardcoded into C file
#  stopifnot(maxlens<maxfraglenC)
#  nlocs = numeric(maxfraglenC^2)
#  nfrag = numeric(maxfraglenC^2)

#  tmpgc = isgcline # Remove NA's
#  tmpgc[is.na(tmpgc)] = FALSE
#  res = .C("makeFullTable",PACKAGE = "GCcorrect",nlocs = as.integer(nlocs),nfrag = as.integer(nfrag), 
#		   gcline = as.integer(tmpgc),gclen = as.integer(length(tmpgc)),locs=as.integer(locs_leng[,1]),
#		   lengs=as.integer(locs_leng[,2]),loclenglen = as.integer(nrow(locs_leng)), 
#		   sampline = as.integer(sampline),minlens = as.integer(minlens),maxlens = as.integer(maxlens), 
#		   margin = as.integer(margin),max_frag = as.integer(max_frag_for_loc))
#  
#  return (list(frag = matrix(res$nfrag,maxfraglenC,maxfraglenC),locs = matrix(res$nlocs,maxfraglenC,maxfraglenC)))
#}



#getFullTablePreds = function(isgcline, predtable, predbeg,predlen , minlens, maxlens,  margin,strand="F") {
#  maxfraglenC = 700
#  stopifnot(length(predtable)== maxfraglenC^2)

#  predictions = numeric(predlen)

#  tmpgc = isgcline # Remove NA's
#  tmpgc[is.na(tmpgc)] = FALSE
#  gclen = length(tmpgc)
  # Missing strand based information... as.integer(strand=="F")
#  res = .C("getFullTablePreds", PACKAGE = "GCcorrect", preds = as.double(predictions),as.integer(predbeg),as.integer(predlen), 
#		       as.double(predtable), as.integer(minlens), as.integer(maxlens),
#		       as.integer(tmpgc), as.integer(gclen), as.integer(margin ) )
#  return(res)
#  
#}
    
makeGCLens = function(gcline,locs_leng , minlens=5,maxlens=400, sampline,margin = 1, max_frag_for_loc=20) {
  # This function estimates conditional means of multiple GC windows of different lengths starting at the same position
  # It is essentially like makeFullTable but does not distinguish according to fragment length.

  # Note - margin=1 discards 1 bp at each side of the window. margin = 0 does not discard anything...
  # In the notion of the paper:
  # For model W_{a,l} we have a=margin and l=length. 
  
  maxfraglenC = 700 # Hardcoded into C file
  stopifnot(maxlens<maxfraglenC)
  nlocs = numeric(maxfraglenC^2)
  nfrag = numeric(maxfraglenC^2)

  tmpgc = gcline # Remove NA's
  tmpgc[is.na(gcline)] = FALSE
  sampline_vec= logical(length(gcline))
  sampline_vec[sampline] = TRUE
  sampline_vec[is.na(gcline)] = FALSE # Perhaps make a larger window around NAs of GC
  res = .C("makeGCLens", PACKAGE = "GCcorrect",nlocs = as.integer(nlocs),nfrag = as.integer(nfrag), 
		   as.integer(tmpgc),as.integer(length(tmpgc)),as.integer(locs_leng[,1]),
		   as.integer(nrow(locs_leng)), 
		   as.integer(sampline_vec),as.integer(minlens),as.integer(maxlens), 
		   as.integer(margin),as.integer(max_frag_for_loc))
  fragmat = matrix(res$nfrag,maxfraglenC,maxfraglenC)
  locsmat = matrix(res$nlocs,maxfraglenC,maxfraglenC)
  rownames(locsmat) = 0:(maxfraglenC-1)
  rownames(fragmat) = 0:(maxfraglenC-1)
  colnames(locsmat) = 0:(maxfraglenC-1)
  colnames(fragmat) = 0:(maxfraglenC-1)
  return (list(frag = fragmat,locs = locsmat))
}

scoreGCLens = function(begdata, maxlens, minlens =  5,margin = 1,scale = TRUE) { 
  tvs = numeric(maxlens)
  for (a in minlens:maxlens) {
    if(scale){
      tvs[a] = compSTV(begdata$frag[a,1:(a+1)],begdata$locs[a,1:(a+1)])
      ylab = "TV Score"
    }
    else{
      tvs[a] = compTV(begdata$frag[a,1:(a+1)],begdata$locs[a,1:(a+1)])
      ylab = "TV distance"
    }
  }
  return(tvs)
}


#genTabAndPredict = function(chr, dat, samp, strand, minlen, maxlen, margin= 3, smooth = 0.7) {
#  require("fields")
#  est_samp = logical(length(chr$isrep))
#  est_samp[samp$singleLocSamp] = TRUE
#  est_samp = est_samp & (!chr$isrep)
#  if(strand=="F") {
#    tab = makeFullTable(chr$isgc,dat$forw , min=minlen,max = maxlen, sampline = est_samp,margin = margin, max_frag_for_loc=5)
#  } else {
#    # Still need to fix this for reverse strand
#    tab = makeFullTable(chr$isgc,dat$reve , min=minlen,max = maxlen, sampline = est_samp,margin = margin, max_frag_for_loc=5)
#  }
#  baseline_reg = 10
#  ratesmooth =  image.smooth(t(as.matrix(tab$frag))/(t(as.matrix(tab$locs))+baseline_reg),theta=smooth)
#  preds = getFullTablePreds (chr$isgc, 0 ,length(chr$isgc), as.numeric(t(ratesmooth$z)), min = minlen, max = maxlen,  margin = margin, strand= "F")
#  preds$preds[is.na(chr$isgc)] = 0;
#  preds$preds[chr$isrep] = 0;
#  preds$preds[is.na(chr$isrep)] = 0;
#  return(list(tab = tab, preds = preds$preds))
#}

condMeanAndPredict = function(chr, samp, strand = "F", winstart, winend,readlen, smooth = 0.2,jump=1,cutoff=3,gcfilename = "gcline") {
  winsize = winend-winstart
  gc = prepareGCChrom(chr,winsize,filename = gcfilename)
  if(strand=="F") {
    cMean = getCondMean(gc[samp$singleLocSamp+winstart],samp$forsamped,jump=jump,cutoff=cutoff)
  } else if (strand == "R") {
    cMean = getCondMean(gc[samp$singleLocSamp-winend+readlen],samp$revsamped,jump=jump,cutoff=cutoff)
  }

  res = predictLine(cMean,gc,winsize,winstart,chr$isrep,strand = strand,readlen=readlen, smooth=smooth)
  return(res)
}

predictLine = function(cMean,gcline,winsize,winstart,isrep,strand = "F",readlen ,smooth=0.2) {
  winend = winsize+winstart
  
  res = list()
  res$unsmooth = cMean

  oldpts = as.numeric(names(cMean$tab_means))/winsize
  newpts = (0:winsize)/winsize
  res$predvec = smoothCondMean(cMean$tab_means,oldpts,newpts,span=smooth,
    zero_pad = c(-0.1,-0.05,1,1.05))  

  # Leaner version of predvec:
  tmp_predvec = res$predvec
  names(tmp_predvec) = c()
  
  # Get Predictions

  tmppreds = tmp_predvec[gcline+1]

  if (strand=="F"){    # The gc of location winstart+1 should be first
    preds = shift(tmppreds,winstart,fill=0)
  }
  if(strand == "R") {
    preds = shift(tmppreds,readlen-winend,fill=0)
  }
  preds[is.na(preds)] = 0;
  preds[isrep] = 0;
  preds[is.na(isrep)] = 0;

  res$preds = preds
  return(res)
}


unEvenBin = function(predline,binstart,binend){
  binline = numeric(length(binstart))
  for (i in 1:length(binstart)) {
    binline[i] = sum(predline[binstart[i]:binend[i]])
  }
  return(binline)
}

# add plot method 
predLoess = function(dat, chr, binsize = "10K", workrange, predict = FALSE , map_cutoff = 0.1, span = 0.3, plot=FALSE,direct = "F"){
  
  forbin = sprintf('for%s',binsize)
  revbin = sprintf('rev%s',binsize)
  mapbin = sprintf('map%s',binsize)
  gcbin = sprintf('gc%s',binsize)
  
  workrange = workrange[chr[[mapbin]][workrange] < map_cutoff]
  if (direct=="F"){
    adjcounts = dat[[forbin]][workrange]/(1-chr[[mapbin]][workrange])
  } else if(direct == "R"){
    adjcounts = dat[[revbin]][workrange]/(1-chr[[mapbin]][workrange])

  } else if(direct == "B"){
    adjcounts = (dat[[forbin]][workrange]+dat[[revbin]][workrange])/(1-chr[[mapbin]][workrange])
  }
  loess = loess(adjcounts~chr[[gcbin]][workrange],span=span)
  if (plot){
    limboth = quantile(adjcounts,0.997)+ 20
    plot(chr[[gcbin]][workrange],adjcounts+rnorm(length(workrange),0,0.3),pch=20,cex= 0.6, col =rgb(0,0,0,0.3), main = 'GC vs. Counts - Bins',ylab = 'Counts', xlab = 'GC',ylim= c(0,limboth) )
  }
  res = list()
  res$loess = loess
  pts = seq(0.2,0.8,0.01)
  res$predvec = pmax(predict(loess,pts),0)
  res$predvec[is.na(res$predvec)] = 0
  names(res$predvec) = pts

  if (predict){
    preds = pmax(predict(loess,chr[[gcbin]]),0)
    preds[is.na(preds)] = 0
    adjpreds = preds*(1-chr[[mapbin]])
    res$preds = adjpreds
  }
  return(res)
}

uncorrLoess = function(dat, chr, binsize = "10K", workrange,  map_cutoff = 0.1, cnt_cutoff = 0,span = 0.3, plot=FALSE,direct = "F"){
  
  forbin = sprintf('for%s',binsize)
  revbin = sprintf('rev%s',binsize)
  mapbin = sprintf('map%s',binsize)
  gcbin = sprintf('gc%s',binsize)

  if (direct=="F"){
    counts = dat[[forbin]]
  } else if(direct == "R"){
    counts = dat[[revbin]]
  } else if(direct == "B"){
    counts = dat[[forbin]]+dat[[revbin]]
  }

  if (cnt_cutoff) {
    workrange = workrange[(counts[workrange] < cnt_cutoff) & (chr[[mapbin]][workrange] < map_cutoff)]
  } else {
    workrange = workrange[(chr[[mapbin]][workrange] < map_cutoff)]
  }

  counts = counts[workrange]

  loesst = loess(counts~chr[[gcbin]][workrange],span=span)

  if (plot){
    limboth = quantile(counts,0.997)+ 20
    plot(chr[[gcbin]][workrange],counts+rnorm(length(workrange),0,0.3),pch=20,cex= 0.6, col =rgb(0,0,0,0.3), main = 'GC vs. Counts - Bins',ylab = 'Counts', xlab = 'GC',ylim= c(0,limboth) )
  }
  res = list()
  res$loess = loesst
  pts = seq(0.2,0.8,0.01)
  res$predvec = pmax(predict(loesst,pts),0)
  res$predvec[is.na(res$predvec)] = 0
  names(res$predvec) = pts

  return(res)
}
