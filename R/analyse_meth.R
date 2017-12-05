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
#' @param method  A function used to extract CN score on each region
#' @importFrom graphics legend
#' @importFrom graphics axis
#' @export
analyse_cnv = function(gene, cnv_data, cnv_platform, cols, PLOT=FALSE, up_str=5000, dwn_str=5000, method=NULL, FULL=FALSE) {
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
  if (is.na(probe_before)) {
    probe_before =    rownames(cnv_platform)[!is.na(cnv_platform$chrom)  & cnv_platform$chrom == chr][1]
  }
  idx = which(rownames(cnv_platform) == probe_before):which(rownames(cnv_platform) == probe_after)
  if (length(idx) < 2) {
    ret = rep(NA, length(sample_names))
    names(ret) = sample_names
    return(ret)
  }
  probes_cnv = rownames(cnv_platform)[idx]
  pos_cnv = cnv_platform[probes_cnv, ]$loc
  probes_cnv_close = rownames(cnv_platform)[idx]

  if (FULL) {
    cnv_val = apply(cnv_data[probes_cnv_close, sample_names],2, function(c) {
      mean = mean(c, na.rm=TRUE)
      sd = sd(c, na.rm=TRUE)
      len = length(c)
      nb_na = sum(is.na(c))      
      return(list(mean, sd, len, nb_na))
    })    
    cnv_val = do.call(rbind, cnv_val)
  } else {
    if (!is.null(method)) {
      cnv_val = apply(cnv_data[probes_cnv_close, sample_names],2, function(c) {
        return(method(c))
      })    
    } else {
      cnv_val = apply(cnv_data[probes_cnv_close, sample_names],2, function(c) {
        if (length(unique(c)) == 1) {
          return(unique(c)[1])
        } else {
          return(NA)
        }
      })    
    }    
  }

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
#' @param PLOT A boolean defining if graphical output must be dispayed on the graphical output.
#' @param JUST_PROBE_POS A boolean defining if function returns only probes positions.
#' @param GAUSSIAN_MIXTURE A boolean specifying if gaussian mixture model is use to convolve signal.
#' @param up_str   An integer specifying up stream size (in bp).
#' @param dwn_str  An integer specifying down stream size (in bp).
#' @param win_size An integer specifying slidding window size (in bp).
#' @param wig_size An integer specifying wiggle size (in bp).
#' @param probe_idx A vector specifying probes associated to the gene.
#' @param mat_mult A function to multiply matrices.
#' @param mat_mult A function to multiply matrices.
#' @param pf_pos_colname string matching the name of the column in the platform that contain the chromosome on which we find a probes.
#' @param pf_pos_colname string matching the name of the column in the platform that contain the position information of probes.
#' @return A matrix of convolved probes signal around the TSS of the selected gene.
#' @importFrom grDevices adjustcolor
#' @importFrom graphics arrows
#' @importFrom graphics matplot
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom dmprocr get_probe_names
#' @examples
#' cols = as.numeric(sunexp_design$sex) * 2
#' names(cols) = rownames(sunexp_design)
#' gene = genes[1,]
#' res1 = analyse_meth(gene, sunexp_data, sunexp_platform, cols, PLOT=TRUE)
#' legend("topright", col=as.numeric(unique(sunexp_design$sex)) * 2, 
#'        legend=unique(sunexp_design$sex), lty=1)
#' @export
analyse_meth = function(
  gene                                ,
  meth_data                           , 
  meth_platform                       , 
  probe_idx                           , 
  pf_chr_colname="Chromosome"         ,
  pf_pos_colname="Start"              ,
  PLOT=TRUE                           , 
  up_str=5000                         , 
  dwn_str=5000                        , 
  win_size=1500                       ,
  wig_size=50                         ,
  JUST_PROBE_POS=FALSE                , 
  GAUSSIAN_MIXTURE=TRUE               ,
  mat_mult = function(A,B) {A %*% B} 
 ) {

  meth_platform = meth_platform[,c(pf_chr_colname,pf_pos_colname)]

  if (substr(meth_platform[1, pf_chr_colname], 1, 3) != "chr") {
    meth_platform[,pf_chr_colname] = paste0("chr",meth_platform[,pf_chr_colname])
  }
  if (substr(gene[[1]], 1, 3) != "chr") {
    gene[[1]] = paste0("chr",gene[[1]])
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
  } else {
    off_set_beg = up_str
    off_set_end = dwn_str
    tss = beg
  }

  ## Compute probes associated with the gene 
  if (missing(probe_idx)) {
    probe_idx = get_probe_names(gene=gene, meth_platform=meth_platform, pf_chr_colname=pf_chr_colname, pf_pos_colname=pf_pos_colname, up_str=up_str, dwn_str=dwn_str)
  }

  if (length(probe_idx) == 0) {
    warning(paste0("No probes for gene ", gene[[4]],"(",gene[[5]],")."))
    return(NULL)
  }

  if (JUST_PROBE_POS) {
    return(list(probe_idx=probe_idx))
  }

  # print(gene[[4]])
  sample_names = colnames(meth_data)


  gene_probe_pos = meth_platform[probe_idx,pf_pos_colname]
  tss_shift = tss %% wig_size

  # NOW it's two step sparse computing of methylome signal: i) mask, ii) convolution.

  # i) mask
  # get independant regions form pos (mask definition)
  foo = gene_probe_pos[-1] - gene_probe_pos[-length(probe_idx)]
  ends_idx = which(foo > win_size)
  starts_idx = ends_idx + 1
  starts_idx = c(1, starts_idx)
  ends_idx = c(ends_idx, length(probe_idx))
  reg_infos = data.frame(starts_idx, ends_idx)
  reg_infos = data.frame(starts_idx, ends_idx, beg=gene_probe_pos[starts_idx], end=gene_probe_pos[ends_idx])
  reg_infos$len =  reg_infos$end - reg_infos$beg
  reg_infos$nb_probes =  reg_infos$ends_idx - reg_infos$starts_idx + 1
  # print(paste(nrow(reg_infos), "regions found."))
  # Go!
  preproc_matrices = lapply(1:nrow(reg_infos), function(i) {
    tmp_reg = reg_infos[i,]
    tmp_starts_idx = tmp_reg$starts_idx
    tmp_ends_idx = tmp_reg$ends_idx
    tmp_probe_idx = tmp_starts_idx:tmp_ends_idx
    tmp_probe_pos = gene_probe_pos[tmp_probe_idx]
    tmp_d = meth_data[probe_idx[tmp_probe_idx],]
    if (length(tmp_probe_idx) > 1) {
      tmp_d = t(tmp_d)
    } else {
      tmp_d = t(t(tmp_d))
    }
    # bases for convotutions
    if (strand == "-") {
      tmp_pos =  (gene_probe_pos[tmp_starts_idx]-win_size/2+1):(gene_probe_pos[tmp_ends_idx]+win_size/2)
      wig_starts = (floor(min((tmp_pos - tss_shift) / wig_size)):floor(max((tmp_pos - tss_shift) / wig_size))) * wig_size + tss_shift + 1
    } else {
      tmp_pos =  (gene_probe_pos[tmp_starts_idx]-win_size/2):(gene_probe_pos[tmp_ends_idx]+win_size/2-1)
      wig_starts = (floor(min((tmp_pos - tss_shift) / wig_size)):floor(max((tmp_pos - tss_shift) / wig_size))) * wig_size + tss_shift
    }
    # convolution matrix for sliding window stuffs
    if (GAUSSIAN_MIXTURE) {
      win_mat = t(sapply(tmp_probe_pos, function(p){
        dnorm(tmp_pos, p, win_size/8)
      }))
    } else {
      win_mat = sapply(tmp_pos, function(p){
        if (strand == "-") {
          vec = p - tmp_probe_pos > -win_size/2 & p - tmp_probe_pos <= win_size/2
        } else {
          vec = p - tmp_probe_pos >= -win_size/2 & p - tmp_probe_pos < win_size/2
        }
        ret = vec
        return(ret)
      })
    }
    # density of probes
    if (length(tmp_probe_idx) == 1) {
      reg_probe_density = win_mat
    } else {
      reg_probe_density = apply(win_mat, 2, sum)
    }      
    # normalize win_mat
    if (length(tmp_probe_idx) == 1) {
      win_mat = t(rep(1, length(tmp_pos)))
    } else {
      win_mat = apply(win_mat, 2, function(c) {c / sum(c)})      
    }      
    # convolution matrix for wiggled stuffs
    wig_mat = sapply(wig_starts, function(w){
      vec = tmp_pos - w >= 0 & tmp_pos - w < wig_size
      if (sum(vec) != 0) {
        ret = (vec / sum(vec))
      } else {
        ret = vec
      }
      return(ret)
    })
    l = list(
      tmp_starts_idx=tmp_starts_idx, 
      tmp_ends_idx=tmp_ends_idx, 
      tmp_pos=tmp_pos, 
      reg_probe_density=reg_probe_density, 

      wig_starts=wig_starts, 
      tmp_d=tmp_d, 
      win_mat=win_mat, 
      wig_mat=wig_mat
    )
    return(l)
  })

  # let's do it 2 (time comsuming)
  wig_data_sparse = lapply(1:nrow(reg_infos), function(i) {
    # print(i)
    l = preproc_matrices[[i]]
    preproc_matrices
    wig_starts = l$wig_starts
    tmp_d      = l$tmp_d
    win_mat    = l$win_mat
    wig_mat    = l$wig_mat
    # convolved data
    mat_prod = mat_mult(mat_mult(tmp_d, win_mat), wig_mat)
    ret = cbind(wig_starts, t(mat_prod))
    return(ret)
  })
  wig_data_sparse = do.call(rbind, wig_data_sparse)
  rownames(wig_data_sparse) = paste("pos_", wig_data_sparse[,1], sep="")

  # let's do it 2
  wig_probe_density_split = unlist(lapply(preproc_matrices, function(l) {
    reg_probe_density = l$reg_probe_density
    wig_mat    = l$wig_mat
    # convolved data
    ret = mat_mult(reg_probe_density, wig_mat)
    return(ret)
  }))

  wig_starts_split = unlist(lapply(preproc_matrices, function(l) {
    l$wig_starts
  }))


  # from splited densities to full density
  b_coords = (tss - off_set_beg):(tss + off_set_beg)
  probe_density = rep(0, length(b_coords))
  idx = unlist(sapply(preproc_matrices, function(i) {i$tmp_pos})) - min(b_coords) + 1
  density_val = unlist(sapply(preproc_matrices, function(i) {i$reg_probe_density}))
  density_val = density_val[idx>0 & idx<=length(probe_density)]
  idx = idx[idx>0 & idx<=length(probe_density)]
  probe_density[idx] = density_val


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


  # from wig_probe_density_split to wig_probe_density
  names(wig_probe_density_split) = names(wig_starts_split) = paste0("pos_", wig_starts_split)
  wig_probe_density = rep(0, length(wig_coords))
  names(wig_probe_density) = paste0("pos_", wig_coords)
  wig_probe_density_split = wig_probe_density_split[names(wig_probe_density_split) %in% names(wig_probe_density)]
  wig_starts_split = wig_starts_split[names(wig_starts_split) %in% names(wig_probe_density)]
  wig_probe_density[names(wig_probe_density_split)] = wig_probe_density_split


  if (strand == "-") {
    wig_data_full = wig_data_full[nrow(wig_data_full):1, ]
    wig_coords = rev(wig_coords)
    wig_probe_density = rev(wig_probe_density)
  }

  return(list(
    gene            = gene                     , 
    meth_data       = meth_data[probe_idx,]    , 
    meth_platform   = meth_platform[probe_idx,], 
    probe_idx       = probe_idx                , 
    pf_chr_colname  = pf_chr_colname           , 
    pf_pos_colname  = pf_pos_colname           , 
    up_str          = up_str                   , 
    dwn_str         = dwn_str                  , 
    sample_names    = sample_names             ,
    tss             = tss                      ,
    wig_probe_density = wig_probe_density,
    wig_starts_split = wig_starts_split,
    wig_probe_density_split = wig_probe_density_split,
    wig_coords =wig_coords,
    wig_probe_density=wig_probe_density, 
    # preproc_matrices=preproc_matrices,
    wig_data_full=wig_data_full))
}


#' A Function that Plots Methylome Analysis Output
#'
#' This function plots methylome analysis output
#' @param meth_analyse A list, corresponding to analyse_meth function output
#' @param cols A color vectors indexed by by samples names.
#' @param main A character string specifying the title of the plot.
#' @param alpha.f=0.1 a numeric specifying transparency of convolved signal.
#' @importFrom grDevices adjustcolor
#' @importFrom graphics arrows
#' @importFrom graphics matplot
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @export
plot_meth = function(
  meth_analyse,
  cols         ,
  main         ,
  alpha.f=0.1                       
 ) {

  gene            = meth_analyse$gene
  meth_data       = meth_analyse$meth_data      
  meth_platform   = meth_analyse$meth_platform  
  probe_idx       = meth_analyse$probe_idx      
  pf_chr_colname  = meth_analyse$pf_chr_colname 
  pf_pos_colname  = meth_analyse$pf_pos_colname 
  up_str          = meth_analyse$up_str         
  dwn_str         = meth_analyse$dwn_str         
  sample_names    = meth_analyse$sample_names     
  wig_starts_split = meth_analyse$wig_starts_split    
  tss    = meth_analyse$tss  
  wig_coords = meth_analyse$wig_coords
  wig_probe_density = meth_analyse$wig_probe_density       
  wig_probe_density_split = meth_analyse$wig_probe_density_split
  wig_data_full = meth_analyse$wig_data_full

  # check chr notation
  if (substr(meth_platform[1, pf_chr_colname], 1, 3) != "chr") {
   meth_platform[,pf_chr_colname] = paste0("chr",meth_platform[,pf_chr_colname])
  }
  if (substr(gene[[1]], 1, 3) != "chr") {
   gene[[1]] = paste0("chr",gene[[1]])
  }

  # get gene properties
  chr =            gene[[1]]
  strand =         gene[[6]]
  gene_name =      gene[[4]]
  beg = as.numeric(gene[[2]])
  end = as.numeric(gene[[3]])



  if (missing(main)) {
    main=""
  } else {
    main = paste0(main, " ")
  }
  if (missing(cols)) {
    cols=c()
    cols[sample_names] = 1
  } else {
    cols = cols[sample_names]
  }

  # display parameters
  y_base = -0.1
  # base
  plot(0, 0, col=0, xlim=c(tss-max(up_str, dwn_str), tss+max(up_str, dwn_str)), ylim=c(-0.5,1), 
                    ylab="beta", 
                    xlab=unique(meth_platform[probe_idx, pf_chr_colname]), 
                    main=paste0(main, gene_name))
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

  # # CpG Island
  # sub_meth_platform = meth_platform[probe_idx,]
  # cpg_islands = methmybeachup:::get_cpg_islands(sub_meth_platform)
  # # print(paste("cpg_islands:", cpg_islands))
  # pos_cpg_islands = cpg_islands$cen[chr == cpg_islands$chrom]
  # # points(pos_cpg_islands, rep(y_base -0.4, length(pos_cpg_islands)), pch=16, col=adjustcolor(2, alpha.f=0.3))
  # apply(cpg_islands, 1, function(cpg_island) { polygon(c(cpg_island[["beg"]], cpg_island[["end"]], cpg_island[["end"]], cpg_island[["beg"]]), c(y_base-0.41, y_base-0.41, y_base-0.39, y_base-0.39), col=2)})
  # # CpG probes
  # points(meth_platform[probe_idx,pf_pos_colname], rep(y_base - 0.2, length(probe_idx)), pch=16, col=adjustcolor(4, alpha.f=0.3))

  # raw probes signal
  pos_x = meth_platform[probe_idx,pf_pos_colname]
  meth_data = meth_data[probe_idx, ]
  if (length(probe_idx) == 1) {
    meth_data = t(meth_data)
  }
  apply(meth_data, 2, function(c) {
    points(pos_x, c, col=adjustcolor("grey",0.3))     
  })

  # density of probes
  # lines(b_coords, (probe_density/max(probe_density)/2) - 0.5, type="l", col=4)
  lines(wig_coords, (wig_probe_density/max(wig_probe_density)/2) - 0.5, type="l", col=2)
  # lines(wig_starts_split, (wig_probe_density_split/max(wig_probe_density_split)/2) - 0.5, type="l", col=2)


  # convolved signal
  matplot(wig_coords, wig_data_full, type="l", lty=1, col=adjustcolor(cols[sample_names], alpha.f=0.1), lwd=3, add=TRUE)
  legend("bottomright", legend=paste("(", table(cols), ")", sep=""), col=names(table(cols)), lty=1)

}


get_cpg_islands = function(meth_platform, cgi_colname="UCSC_CpG_Islands_Name") {
  if (cgi_colname == "UCSC_CpG_Islands_Name") {
    cpg_islands = sort(unique(as.character(meth_platform[[cgi_colname]])))
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
  } else if (cgi_colname == "CGI_Coordinate"){
    stop("get_cpg_islands. Not yet!")    
  }
  return(cpg_islands)
}








# layout(matrix(1:32,4))
# bar = apply(stable_coords[1:32,],1,analyse_meth, TRUE)
# baz = apply(tsps_coords,1,analyse_meth)
# gene = tsps_coords[1,]; analyse_meth(gene)
