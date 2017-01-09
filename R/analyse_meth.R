#' A Function Analyses Methylome data
#'
#' This function analyses methylome data.
#' @param gene A vector describing the gene (line of a bed file).
#' @param sunexp_data A matrix of beta values.
#' @param sunexp_platform A data frame describing CpG positions.
#' @param cols A color vectors indexed by by samples names.
#' @param PLOT A boolean defining if graphical output must be dispayed on the graphical output.
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
#' legend("topright", col=as.numeric(unique(sunexp_design$sex)) * 2, legend=unique(sunexp_design$sex), lty=1)
#' @export
analyse_meth = function(gene, sunexp_data, sunexp_platform, cols, PLOT=TRUE, up_str=5000, dwn_str=5000, win_size=1500, wig_size=50) {
  # print(gene)
  sample_names = colnames(sunexp_data)
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
    reg_probes_names = rownames(sunexp_platform)[
      !is.na(sunexp_platform$MAPINFO) & !is.na(sunexp_platform$CHR) &
      sunexp_platform$CHR == chr &
      sunexp_platform$MAPINFO > end-dwn_str &
      sunexp_platform$MAPINFO <= end+up_str
    ]
  } else {
    off_set_beg = up_str
    off_set_end = dwn_str
    tss = beg
    reg_probes_names = rownames(sunexp_platform)[
      !is.na(sunexp_platform$MAPINFO) & !is.na(sunexp_platform$CHR) &
      sunexp_platform$CHR == chr &
      sunexp_platform$MAPINFO >= beg-up_str &
      sunexp_platform$MAPINFO < beg+dwn_str
    ]
  }
  reg_probes_pos = sunexp_platform[reg_probes_names, "MAPINFO"]
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
    tmp_d = sunexp_data[reg_probes_names[tmp_probe_idx],]
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
                      xlab=paste("chr", unique(sunexp_platform[reg_probes_names, ]$CHR), sep=""), 
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
    sub_sunexp_platform = sunexp_platform[reg_probes_names,]
    cpg_islands = get_cpg_islands(sub_sunexp_platform)
    # print(paste("cpg_islands:", cpg_islands))
    pos_cpg_islands = cpg_islands$cen[chr == cpg_islands$chrom]
    # points(pos_cpg_islands, rep(y_base -0.4, length(pos_cpg_islands)), pch=16, col=adjustcolor(2, alpha.f=0.3))
    apply(cpg_islands, 1, function(cpg_island) {
      polygon(c(cpg_island[["beg"]],cpg_island[["end"]],cpg_island[["end"]],cpg_island[["beg"]]), c(y_base-0.41, y_base-0.41, y_base-0.39, y_base-0.39), col=2)
    })
    # CpG probes
    points(sunexp_platform[reg_probes_names, "MAPINFO"], rep(y_base - 0.2, length(reg_probes_names)), pch=16, col=adjustcolor(4, alpha.f=0.3))
    # raw probes signal
    pos_x = sunexp_platform[reg_probes_names,]$MAPINFO
    meth_data = sunexp_data[reg_probes_names, ]
    if (length(reg_probes_names) == 1) {
      meth_data = t(meth_data)
    }
    apply(meth_data, 2, function(c) {
      points(pos_x, c, col=adjustcolor("grey",0.3))     
    })
    # convolved signal
    matplot(wig_coords, wig_data_full, type="l", lty=1, col=cols[sample_names], add=TRUE)
  }
    
  if (strand == "-") {
    wig_data_full = wig_data_full[nrow(wig_data_full):1, ]
  }
  
  return(wig_data_full)
}

get_cpg_islands = function(sunexp_platform) {
  cpg_islands = sort(unique(as.character(sunexp_platform$UCSC_CpG_Islands_Name)))
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
