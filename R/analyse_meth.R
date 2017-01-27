#' A Function Analyses CNV data
#'
#' This function analyses methylome data.
#' @param gene A vector describing the gene (line of a bed file).
#' @param cnv_data A matrix of CN values.
#' @param cnv_platform A data frame describing probe positions.
#' @param cols A color vectors indexed by by samples names.
#' @param PLOT A boolean defining if graphical output must be dispayed on the graphical output.
#' @param up_str   An integer specifying up stream size (in bp).
#' @param dwn_str  An integer specifying down stream size (in bp).
#' @importFrom graphics legend
#' @importFrom graphics axis
#' @export
analyse_cnv = function(gene, cnv_data, cnv_platform, cols, PLOT=FALSE, up_str=5000, dwn_str=5000) {
  # params
  sample_names = colnames(cnv_data)
  if (missing(cols)) {
    cols=c()
    cols[sample_names] = 1
  }

  # get gene properties
  chr =              gene[[1]]
  beg =   as.numeric(gene[[2]])
  end =   as.numeric(gene[[3]])
  gene_name =        gene[[4]]
  id_probe =         gene[[5]]
  strand =           gene[[6]]

  # get cnv infos
  if (strand == "-") {
    off_set_beg = dwn_str
    off_set_end = up_str
  } else {
    off_set_beg = up_str
    off_set_end = dwn_str
  }  
  probe_before = rev(rownames(cnv_platform)[!is.na(cnv_platform$chrom)  & cnv_platform$chrom == chr & cnv_platform$loc < beg-off_set_beg])[1]
  probe_after =      rownames(cnv_platform)[!is.na(cnv_platform$chrom)  & cnv_platform$chrom == chr & cnv_platform$loc > end+off_set_end] [1]
  if (is.na(probe_after)) {
    probe_after =    rev(rownames(cnv_platform)[!is.na(cnv_platform$chrom)  & cnv_platform$chrom == chr]) [1]
  }
  idx = which(rownames(cnv_platform) == probe_before):which(rownames(cnv_platform) == probe_after)
  probes_cnv = rownames(cnv_platform)[idx]
  pos_cnv = cnv_platform[probes_cnv, ]$loc
  probes_cnv_close = rownames(cnv_platform)[idx]

  cnv_val = apply(cnv_data[probes_cnv_close, sample_names],2, function(c) {
    if (length(unique(c)) == 1) {
      return(unique(c)[1])
    } else {
      return(NA)
    }
  })

  if (PLOT) { 
    max_cnv = 8
    y_base = 0
    plot(0, 0, col=0, xaxt="n",
      xlim=range(c(pos_cnv)), ylim=c(0, max_cnv)                          , 
      ylab="CN",                                                          , 
      xlab=paste("chr", gene[[1]], ":", gene[[2]], "-", gene[[3]], sep=""),
      main=paste(gene_name, sep="")
    )
    axis(1, at=gene[2:3])
    polygon(c(beg,end,end,beg), c(y_base-0.1, y_base-0.1, y_base+0.1, y_base+0.1), col=1)
    if (strand == "-") {
      arrows(end, y_base, end, y_base+0.3, length=0)
      arrows(end, y_base+0.3, end - (end-beg)/3, y_base+0.3, length=0.03)
    } else {
      arrows(beg, y_base, beg, y_base+0.3, length=0)
      arrows(beg, y_base+0.3, beg + (end-beg)/3, y_base+0.3, length=0.03)
    }  
    points(pos_cnv , rep(y_base      , length(pos_cnv )), pch=16, col=adjustcolor(2, alpha.f=0.6))
    matplot(pos_cnv, cnv_data[probes_cnv,sample_names], type="l", lty=1, col=adjustcolor(cols[sample_names], alpha.f=0.5), lwd=3, add=TRUE)
    
  }
  return(cnv_val)
}


#' A Function Analyses CNV and Transcriptome data
#'
#' This function analyses methylome data.
#' @param gene A vector describing the gene (line of a bed file).
#' @param trscr_res A vector of expression values.
#' @param cnv_res A vector of CN values.
#' @param meth_idx A vector of sample for which methylom analysis is available.
#' @param cols A color vectors indexed by by samples names.
#' @param idxes A list of index defining sample groups. 
#' @param ctrl_idx A vector of sample used to define threshold 
#' @param PLOT A boolean defining if graphical output must be dispayed on the graphical output.
#' @importFrom graphics legend
#' @export
analyse_trscr_cnv = function(gene, trscr_res, cnv_res, meth_idx, ctrl_idx, cols, idxes, PLOT=TRUE) {
  # get gene properties
  chr =              gene[[1]]
  beg =   as.numeric(gene[[2]])
  end =   as.numeric(gene[[3]])
  gene_name =        gene[[4]]
  id_probe =         gene[[5]]
  strand =           gene[[6]]
  
  trscr_res_orig = trscr_res
  trscr_res = trscr_res[names(cnv_res)]

  # thresholds
  if (missing(cols) | missing(idxes)) {
    trscr_thres_hi = trscr_thres_lo = max(trscr_res_orig[ctrl_idx])
    cnv_thres_hi = 3.5
    cnv_thres_lo = 2.5
    trscr_lo_cnv_lo_idx = names(trscr_res)[!is.na(cnv_res) & trscr_res < trscr_thres_lo & cnv_res < cnv_thres_lo]
    trscr_hi_cnv_lo_idx = names(trscr_res)[!is.na(cnv_res) & trscr_res >= trscr_thres_hi & cnv_res < cnv_thres_lo]
    trscr_hi_cnv_hi_idx = names(trscr_res)[!is.na(cnv_res) & trscr_res > trscr_thres_hi & cnv_res >= cnv_thres_hi]
    idxes = list(trscr_hi_cnv_hi=trscr_hi_cnv_hi_idx,
                 trscr_hi_cnv_lo=trscr_hi_cnv_lo_idx,
                 trscr_lo_cnv_lo=trscr_lo_cnv_lo_idx)
  } else {
    idxes = idxes
  }
  
  # cols
  sample_names = names(cnv_res)
  if (missing(cols)) {
    cols = rep(NA, length(cnv_res))
    names(cols) = names(cnv_res)
    cols[!is.na(cnv_res)] = "grey"  
    cols[idxes$trscr_hi_cnv_hi] = "red"
    cols[idxes$trscr_hi_cnv_lo] = "green"
    cols[idxes$trscr_lo_cnv_lo] = "blue"
  } else {
    cols = cols[sample_names]
  }

  # cols

  idx = names(cnv_res)
  idx2 = meth_idx
  idx2 = idx2[idx2 %in% idx]
  cex = pch = c()
  pch[idx] = cex[idx] = 1
  pch[idx2] = 16
  cex[idx2] = 3

  if (PLOT) {
    plot(cnv_res[idx], trscr_res[idx], pch=pch, cex=cex, col=adjustcolor(cols[idx], alpha.f=0.3), 
      main=paste(gene_name, "@", id_probe, sep=""), xlab="CN", ylab="log2(expr)")
    # abline(h=trscr_thres_hi, v=c(cnv_thres_hi, cnv_thres_lo), col="grey", lty=2)
    legend("bottomright", legend=paste(names(idxes), " (", sapply(idxes, length), ")", sep=""), col=c("red", "green", "blue"), pch=1)
  }
  return(list(idx=idxes, cols=cols))
}

#' A Function Analyses Methylome data
#'
#' This function analyses methylome data.
#' @param gene A vector describing the gene (line of a bed file).
#' @param meth_data A matrix of beta values.
#' @param meth_platform A data frame describing CpG positions.
#' @param cols A color vectors indexed by by samples names.
#' @param PLOT A boolean defining if graphical output must be dispayed on the graphical output.
#' @param JUST_PROBE_POS A boolean defining if function returns only probes positions.
#' @param up_str   An integer specifying up stream size (in bp).
#' @param dwn_str  An integer specifying down stream size (in bp).
#' @param win_size An integer specifying slidding window size (in bp).
#' @param wig_size An integer specifying wiggle size (in bp).
#' @return A matrix of convolved probes signal around the TSS of the selected gene.
#' @importFrom grDevices adjustcolor
#' @importFrom graphics arrows
#' @importFrom graphics matplot
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @examples
#' cols = as.numeric(sunexp_design$sex) * 2
#' names(cols) = rownames(sunexp_design)
#' gene = genes[1,]
#' res1 = analyse_meth(gene, sunexp_data, sunexp_platform, cols, PLOT=TRUE)
#' legend("topright", col=as.numeric(unique(sunexp_design$sex)) * 2, 
#'        legend=unique(sunexp_design$sex), lty=1)
#' @export
analyse_meth = function(gene, meth_data, meth_platform, cols, PLOT=TRUE, up_str=5000, dwn_str=5000, win_size=1500, wig_size=50, JUST_PROBE_POS=FALSE) {
  # print(gene)
  sample_names = colnames(meth_data)
  if (missing(cols)) {
    cols=c()
    cols[sample_names] = 1
  } else {
    cols = cols[sample_names]
  }
  # get gene properties
  chr =            gene[[1]]
  strand =         gene[[6]]
  gene_name =      gene[[4]]
  beg = as.numeric(gene[[2]])
  end = as.numeric(gene[[3]])

  # get meth infos
  if (strand == "-") {
  # if (TRUE) {
    off_set_beg = dwn_str
    off_set_end = up_str
    tss = end
    reg_probes_names = rownames(meth_platform)[
      !is.na(meth_platform$MAPINFO) & !is.na(meth_platform$CHR) &
      meth_platform$CHR == chr &
      meth_platform$MAPINFO > end-dwn_str &
      meth_platform$MAPINFO <= end+up_str
    ]
  } else {
    off_set_beg = up_str
    off_set_end = dwn_str
    tss = beg
    reg_probes_names = rownames(meth_platform)[
      !is.na(meth_platform$MAPINFO) & !is.na(meth_platform$CHR) &
      meth_platform$CHR == chr &
      meth_platform$MAPINFO >= beg-up_str &
      meth_platform$MAPINFO < beg+dwn_str
    ]
  }
  if (JUST_PROBE_POS) {
    return(reg_probes_names)
  }
  reg_probes_pos = meth_platform[reg_probes_names, "MAPINFO"]
  tss_shift = tss %% wig_size

  if (length(reg_probes_pos) == 0) {
    return(NULL)
  }

  # get independant regions form pos
  foo = reg_probes_pos[-1] - reg_probes_pos[-length(reg_probes_names)]
  ends_idx = which(foo > win_size)
  starts_idx = ends_idx + 1
  starts_idx = c(1, starts_idx)
  ends_idx = c(ends_idx, length(reg_probes_names))
  reg_infos = data.frame(starts_idx, ends_idx)
  reg_infos = data.frame(starts_idx, ends_idx, beg=reg_probes_pos[starts_idx], end=reg_probes_pos[ends_idx])
  reg_infos$len =  reg_infos$end - reg_infos$beg
  reg_infos$nb_probes =  reg_infos$ends_idx - reg_infos$starts_idx + 1

  # let's do it
  wig_data_sparse = lapply(1:nrow(reg_infos), function(i) {
    tmp_reg = reg_infos[i,]
    tmp_starts_idx = tmp_reg$starts_idx
    tmp_ends_idx = tmp_reg$ends_idx
    tmp_probe_idx = tmp_starts_idx:tmp_ends_idx
    tmp_probe_pos = reg_probes_pos[tmp_probe_idx]
    tmp_d = meth_data[reg_probes_names[tmp_probe_idx],]
    if (length(tmp_probe_idx) > 1) {
      tmp_d = t(tmp_d)
    } else {
      tmp_d = t(t(tmp_d))
    }
    # bases for convotutions
    if (strand == "-") {
      tmp_pos =  (reg_probes_pos[tmp_starts_idx]-win_size/2+1):(reg_probes_pos[tmp_ends_idx]+win_size/2)
      wig_starts = (floor(min((tmp_pos - tss_shift) / wig_size)):floor(max((tmp_pos - tss_shift) / wig_size))) * wig_size + tss_shift + 1
    } else {
      tmp_pos =  (reg_probes_pos[tmp_starts_idx]-win_size/2):(reg_probes_pos[tmp_ends_idx]+win_size/2-1)
      wig_starts = (floor(min((tmp_pos - tss_shift) / wig_size)):floor(max((tmp_pos - tss_shift) / wig_size))) * wig_size + tss_shift
    }
    # convolution matrix for sliding window stuffs
    win_mat = sapply(tmp_pos, function(p){
      if (strand == "-") {
        vec = p - tmp_probe_pos > -750 & p - tmp_probe_pos <= 750
      } else {
        vec = p - tmp_probe_pos >= -750 & p - tmp_probe_pos < 750
      }
      ret = (vec / sum(vec))
      return(ret)
    })
    # convolution matrix for wiggled stuffs
    wig_mat = sapply(wig_starts, function(w){
      vec = tmp_pos - w >= 0 & tmp_pos - w < wig_size
      ret = (vec / sum(vec))
      return(ret)
    })
    # convolved data
    mat_prod = tmp_d %*% win_mat %*% wig_mat
    ret = cbind(wig_starts, t(mat_prod))
    return(ret)
  })
  wig_data_sparse = do.call(rbind, wig_data_sparse)
  rownames(wig_data_sparse) = paste("pos_", wig_data_sparse[,1], sep="")

  # from sparse matrix to full matrix
  if (strand == "-") {
    wig_coords = ((tss - tss_shift - off_set_beg) / wig_size):((tss - tss_shift + off_set_beg)/wig_size) * 50 + tss_shift + 1
  } else {
    wig_coords = ((tss - tss_shift - off_set_beg) / wig_size):((tss - tss_shift + off_set_beg)/wig_size) * 50 + tss_shift
  }
  wig_data_full = matrix(NA, nrow=length(wig_coords), ncol=length(sample_names))
  rownames(wig_data_full) = paste("pos_", wig_coords, sep="")
  colnames(wig_data_full) = sample_names
  wig_sparse_idx = wig_data_sparse[,1][wig_data_sparse[,1] %in% wig_coords]
  wig_data_full[paste("pos_", wig_sparse_idx, sep=""),sample_names] = wig_data_sparse[paste("pos_", wig_sparse_idx, sep=""),sample_names]

  # wig_data_full = NULL

  # plot
  if (PLOT) { 
    # display parameters
    y_base = -0.1
    # base
    plot(0, 0, col=0, xlim=c(tss-up_str, tss+up_str), ylim=c(-0.5,1), 
                      ylab="beta", 
                      xlab=paste("chr", unique(meth_platform[reg_probes_names, ]$CHR), sep=""), 
                      main=gene_name)
    # gene
    polygon(c(beg,end,end,beg), c(y_base-0.01, y_base-0.01, y_base+0.01, y_base+0.01), col=1)
    arr_len = (end-beg)/3
    arr_len = 500
    if (strand == "-") {
      arrows(end, y_base, end, y_base+0.05, length=0)
      arrows(end, y_base+0.05, end - arr_len, y_base+0.05, length=0.03)
    } else {
      arrows(beg, y_base, beg, y_base+0.05, length=0)
      arrows(beg, y_base+0.05, beg + arr_len, y_base+0.05, length=0.03)
    }  

    # CpG Island
    sub_meth_platform = meth_platform[reg_probes_names,]
    cpg_islands = get_cpg_islands(sub_meth_platform)
    # print(paste("cpg_islands:", cpg_islands))
    pos_cpg_islands = cpg_islands$cen[chr == cpg_islands$chrom]
    # points(pos_cpg_islands, rep(y_base -0.4, length(pos_cpg_islands)), pch=16, col=adjustcolor(2, alpha.f=0.3))
    apply(cpg_islands, 1, function(cpg_island) {
      polygon(c(cpg_island[["beg"]],cpg_island[["end"]],cpg_island[["end"]],cpg_island[["beg"]]), c(y_base-0.41, y_base-0.41, y_base-0.39, y_base-0.39), col=2)
    })
    # CpG probes
    points(meth_platform[reg_probes_names, "MAPINFO"], rep(y_base - 0.2, length(reg_probes_names)), pch=16, col=adjustcolor(4, alpha.f=0.3))
    # raw probes signal
    pos_x = meth_platform[reg_probes_names,]$MAPINFO
    meth_data = meth_data[reg_probes_names, ]
    if (length(reg_probes_names) == 1) {
      meth_data = t(meth_data)
    }
    apply(meth_data, 2, function(c) {
      points(pos_x, c, col=adjustcolor("grey",0.3))     
    })
    # convolved signal
    matplot(wig_coords, wig_data_full, type="l", lty=1, col=adjustcolor(cols[sample_names], alpha.f=0.5), lwd=3, add=TRUE)
    legend("bottomright", legend=paste("(", table(cols), ")", sep=""), col=names(table(cols)), lty=1)
  }
    
  if (strand == "-") {
    wig_data_full = wig_data_full[nrow(wig_data_full):1, ]
  }
  
  return(wig_data_full)
}

get_cpg_islands = function(meth_platform) {
  cpg_islands = sort(unique(as.character(meth_platform$UCSC_CpG_Islands_Name)))
  cpg_island = cpg_islands[1]
  foo = apply(t(cpg_islands), 2, function(cpg_island) {
    # print(cpg_island)
    s = strsplit(cpg_island, ":")[[1]]
    chrom = substr(s[1],4,1000)
    len = as.numeric(strsplit(s[2], "-")[[1]])
    # print(s)
    # print(len)
    beg = len[1]
    end = len[2]
    return(list(cpg_island=cpg_island, chrom=chrom, beg=beg, end=end))
  })
  foo = do.call(rbind,foo)
  cpg_islands = data.frame(lapply(data.frame(foo, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
  cpg_islands$cen = (cpg_islands$end + cpg_islands$beg) / 2 
  cpg_islands$len =  cpg_islands$end - cpg_islands$beg 
  return(cpg_islands)
}








# layout(matrix(1:32,4))
# bar = apply(stable_coords[1:32,],1,analyse_meth, TRUE)
# baz = apply(tsps_coords,1,analyse_meth)
# gene = tsps_coords[1,]; analyse_meth(gene)
