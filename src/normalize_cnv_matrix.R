

load("dataMASTER.RData")
load("gencode19_gns_lite.RData")
GENCODE = gencode19_gns_lite


DB = dataMASTER

library(circlize)
library(GenomicRanges)
chr_df = read.chromInfo()$df
chr_df = chr_df[chr_df$chr %in% paste0("chr", c(1:22, "X", "Y")), ]
chr_df[, 1] = gsub("chr", "", chr_df[, 1])
chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
library(EnrichedHeatmap)
chr_window = makeWindows(chr_gr, w = 1e6)


average_in_window = function(window, gr, v, method = "weighted", empty_v = NA) {

    if(missing(v)) v = rep(1, length(gr))
    if(is.null(v)) v = rep(1, length(gr))
    if(is.atomic(v) && is.vector(v)) v = cbind(v)

    v = as.matrix(v)
    if(is.character(v) && ncol(v) > 1) {
        stop("`v` can only be a character vector.")
    }

    if(length(empty_v) == 1) {
        empty_v = rep(empty_v, ncol(v))
    }

    u = matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))

    mtch = as.matrix(findOverlaps(window, gr))
    intersect = pintersect(window[mtch[,1]], gr[mtch[,2]])
    w = width(intersect)
    v = v[mtch[,2], , drop = FALSE]
    n = nrow(v)

    ind_list = split(seq_len(n), mtch[, 1])
    window_index = as.numeric(names(ind_list))
    window_w = width(window)

    if(is.character(v)) {
        for(i in seq_along(ind_list)) {
            ind = ind_list[[i]]
            if(is.function(method)) {
                u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
            } else {
                tb = tapply(w[ind], v[ind], sum)
                u[window_index[i], ] = names(tb[which.max(tb)])
            }
        }
    } else {
        if(method == "w0") {
            gr2 = reduce(gr, min.gapwidth = 0)
            mtch2 = as.matrix(findOverlaps(window, gr2))
            intersect2 = pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])

            width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
            ind = unique(mtch2[, 1])
            width_setdiff = width(window[ind]) - width_intersect

            w2 = width(window[ind])

            for(i in seq_along(ind_list)) {
                ind = ind_list[[i]]
                x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
                u[window_index[i], ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
            }

        } else if(method == "absolute") {
            for(i in seq_along(ind_list)) {
                u[window_index[i], ] = colMeans(v[ind_list[[i]], , drop = FALSE])
            }
            
        } else if(method == "weighted") {
            for(i in seq_along(ind_list)) {
                ind = ind_list[[i]]
                u[window_index[i], ] = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
            }
        } else {
            if(is.function(method)) {
                for(i in seq_along(ind_list)) {
                    ind = ind_list[[i]]
                    u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
                }
            } else {
                stop("wrong method.")
            }
        }
    }

    return(u)
}

CNV_list = as(DB[["cnv"]], "GRangesList")

lt = lapply(CNV_list, function(gr) {
	average_in_window(chr_window, gr, log2(gr$CovRatio), empty_v = 0)
})


m = do.call(cbind, lt)
m = t(m)
m[is.infinite(m)] = 0
rownames(m) = names(CNV_list)

hc = hclust(dist(m))

m = m[hc$order, ]

# Heatmap(m, cluster_rows = FALSE, cluster_columns = FALSE, 
# 	column_split = factor(as.vector(seqnames(chr_window)), levels = c(1:22, "X", "Y")),
# 	column_gap = unit(0, "points"), border = TRUE, column_title_gp = gpar(fontsize = 10))

saveRDS(list(matrix = m, chr_window = chr_window), file = "normalized_CNV.rds")



CNV_list = as(DB[["cnv_germline"]], "GRangesList")

lt = lapply(CNV_list, function(gr) {
    average_in_window(chr_window, gr, gr$CN, empty_v = 2)
})


m = do.call(cbind, lt)
m = t(m)
m[is.infinite(m)] = 0
rownames(m) = names(CNV_list)

hc = hclust(dist(m))

m = m[hc$order, ]

# Heatmap(m, cluster_rows = FALSE, cluster_columns = FALSE, 
#   column_split = factor(as.vector(seqnames(chr_window)), levels = c(1:22, "X", "Y")),
#   column_gap = unit(0, "points"), border = TRUE, column_title_gp = gpar(fontsize = 10))

saveRDS(list(matrix = m, chr_window = chr_window), file = "normalized_CNV_germline.rds")

