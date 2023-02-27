
tabItem = function(..., .list = NULL) {
	lt = list(...)
	if(!is.null(.list)) {
		lt = c(lt, .list)
	}
	do.call(shinydashboard::tabItem, lt)
}


tabItems = function(..., .list = NULL) {
	lt = list(...)
	if(!is.null(.list)) {
		lt = c(lt, .list)
	}
	do.call(shinydashboard::tabItems, lt)
}

menuItem = function(..., .list = NULL) {
	lt = list(...)
	if(!is.null(.list)) {
		lt = c(lt, .list)
	}
	do.call(shinydashboard::menuItem, lt)
}

parse_samples = function(s, experiment_id = NULL) {
	s = strsplit(s, "\\s+")[[1]]

	l = META$PatientID %in% s
	if(!any(l)) {
		l = rownames(META) %in% s
	}
	samples = rownames(META)[l]
	if(!is.null(experiment_id)) {
		obj = DB[[experiment_id]]
		if(inherits(obj, "RaggedExperiment")) {
			samples = intersect(samples, names(obj@assays))
		} else if(inherits(obj, "RangedSummarizedExperiment")) {
			samples = intersect(samples, colnames(obj))
		} else {
			samples = intersect(samples, colnames(obj))
		}
	}
	samples
}

parse_gr_str = function(s) {
	s = gsub("chr", "", s)
	s = gsub(",", "", s)
	s = strsplit(s, ":|-")[[1]]
	if(length(s) == 1) {
		list(chr = s)
	} else if(length(s) == 2) {
		list(chr = s[1], start = as.numeric(s[2]))
	} else if(length(s) == 3) {
		list(chr = s[1], start = as.numeric(s[2]), end = as.numeric(s[3]))
	} else {
		list()
	}
}


make_rainfall = function(gr, layout = "rect") {

	if(length(gr) > 5000) {
		gr = gr[sample(1:length(gr), 5000)]
	}
	seqlevelsStyle(gr) = "UCSC"
	names(gr) = NULL
	gr = as.data.frame(gr)
	if(all(gr[, 3] == gr[, 2])) gr[, 3] = gr[, 3] + 1
		
	circos.clear()
	if(layout == "rect") {
		gr_density = circlize::genomicDensity(gr, window.size = 2e6, chr.len = read.chromInfo(species = GENOME)$chr.len)

		gtrellis_layout(n_track = 4, nrow = 4, compact = TRUE,
		    track_axis = c(FALSE, TRUE, TRUE, FALSE), 
		    track_height = unit.c(2*grobHeight(textGrob("chr1")), 
		                          unit(1, "null"), 
		                          unit(0.5, "null"), 
		                          unit(3, "mm")), 
		    track_ylim = c(0, 1, 0, 8, c(0, max(gr_density[[4]])), 0, 1),
		    track_ylab = c("", "log10(inter_dist)", "density", ""))

		# track for chromosome names
		add_track(panel_fun = function(gr) {
		    # the use of `get_cell_meta_data()` will be introduced later
		    chr = get_cell_meta_data("name")  
		    grid.rect(gp = gpar(fill = "#EEEEEE"))
		    grid.text(chr)
		})

		# track for rainfall plots
		gr_rainfall = circlize::rainfallTransform(gr)
		gr_rainfall$dist[is.na(gr_rainfall$dist)] = 8e7
		add_points_track(gr_rainfall, log10(gr_rainfall[[4]]),
		  pch = 16, size = unit(2, "pt"), gp = gpar(col = "red"))

		# track for genomic density
		
		add_lines_track(gr_density, gr_density[[4]], area = TRUE, 
		  gp = gpar(fill = "pink"))

		# track for ideogram
		cytoband_df = circlize::read.cytoband(species = "hg19")$df
		add_track(cytoband_df, panel_fun = function(gr) {
		    cytoband_chr = gr
		    grid.rect(cytoband_chr[[2]], unit(0, "npc"),
		              width = cytoband_chr[[3]] - cytoband_chr[[2]], height = unit(1, "npc"),
		              default.units = "native", hjust = 0, vjust = 0,
		              gp = gpar(fill = circlize::cytoband.col(cytoband_chr[[5]])))
		    grid.rect(min(cytoband_chr[[2]]), unit(0, "npc"),
		              width = max(cytoband_chr[[3]]) - min(cytoband_chr[[2]]), height = unit(1, "npc"),
		              default.units = "native", hjust = 0, vjust = 0,
		              gp = gpar(fill = "transparent"))
		})
	} else if(layout == "circular") {
		circos.initializeWithIdeogram(species = GENOME)
		circos.genomicRainfall(gr, col = "red", pch = 16, cex = 0.5)
		circos.genomicDensity(gr, col = "pink", border = "black")
	}
}


make_oncoprint = function(genes, samples, experiments) {

	if("rna" %in% experiments) {
		use_rna = TRUE
	} else {
		use_rna = FALSE
	}

	experiments = setdiff(experiments, "rna")

	ml = vector("list", length(experiments))
	names(ml) = experiments
	ml = lapply(ml, function(x) {
		m = matrix(FALSE, nrow = length(genes), ncol = length(samples))
		rownames(m) = genes
		colnames(m) = samples
		m
	})

	for(e in experiments) {
		rl = DB[[e]]@assays
		for(sn in samples) {
			if(!is.null(rl[[sn]])) {
				ind = genes %in% rl[[sn]]$Gene
				ml[[e]][ind, sn] = TRUE
			}
		}
	}

	ht = oncoPrint(ml, name = "oncoprint", alter_fun = ALTER_FUN, col=ALTER_COL, show_column_names=FALSE, remove_empty_columns = TRUE)

	if(use_rna) {
		expr = assay(DB[["rna"]], "TPM")

		gn = structure(names = GENCODE$gene_id, GENCODE$gene_name)

		rownames(expr) = gn[rownames(expr)]

		expr = expr[genes, intersect(samples, colnames(expr))]

		ht = ht + Heatmap(expr, name = "TPM", show_column_names = FALSE)
	}

	ht
}


ALTER_COL = c(indel="goldenrod1", indel_germline="goldenrod4",
	        snv="darkorchid1", snv_germline="darkorchid4",
	        fusion="limegreen", DUP="salmon", AMP="red", DEL="deepskyblue", 
	        HDEL="blue", LOH="grey40")

ALTER_FUN = list(
  background = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#CCCCCC", col=NA)),
  snv = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=ALTER_COL["snv"], col=NA)),
  snv_germline = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.75, gp=gpar(fill=ALTER_COL["snv_germline"], col=NA)),
  indel = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.6, gp=gpar(fill=ALTER_COL["indel"], col=NA)),
  indel_germline = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.45, gp=gpar(fill=ALTER_COL["indel_germline"],
                                           col=NA)),
  fusion = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.45, gp=gpar(fill=ALTER_COL["fusion"], col=NA)),

  DUP = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=ALTER_COL["dup"], col=NA)),
  AMP = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.8, gp=gpar(fill=ALTER_COL["amp"], col=NA)),
  DEL = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.7, gp=gpar(fill=ALTER_COL["del"], col=NA)),
  HDEL = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.6, gp=gpar(fill=ALTER_COL["hdel"], col=NA)),
  LOH = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.3, gp=gpar(fill=ALTER_COL["LOH"], col=NA))
)

plotOutput = function(...) {
	x = shiny::plotOutput(...)
	x$children[[1]] = p("generating plot or waiting for update...", class="waiting_msg")
	x
}


tableOutput = function(...) {
	x = shiny::tableOutput(...)
	x$children[[1]] = p("generating table or waiting for update...", class="waiting_msg")
	x
}

