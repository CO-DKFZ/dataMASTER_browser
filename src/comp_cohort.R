
COMPONENT$cohort = list()

COMPONENT$cohort$ui = div(
	h3("OncoPrint"),
	div(
		column(3,
			textAreaInput("cohort_genes", label = "Genes of interest", value = "TP53\nAPC\nRB1\n", rows = 6)
		),
		column(3,
			textAreaInput("cohort_patients", label = "Patients of interest", value = "BD26XP\nPP8L7S\n8NFHC7\nRX17A3\nFY5ZQW\nBA95UQ\nFKPP\n", rows = 6),
			actionButton("cohort_submit", "Submit")
		),
		div(style = "clear:both")
	),
	br(),
	tabsetPanel(
		id = qq("cohort_tabset"),
		type = "tabs",
		tabPanel("Mutations",
			uiOutput("cohort_mutations_oncoprint_ui"),
	   		style = "padding:16px 0px"
		),
		tabPanel("CNV",
			uiOutput("cohort_cnv_oncoprint_ui"),
			style = "padding:16px 0px"
		),
		tabPanel("Expression",
			uiOutput("cohort_expression_ui"),
			style = "padding:16px 0px"
		)
	)
)

COMPONENT$cohort$server = function(input, output, session) {

	observeEvent(input$cohort_submit, {

		genes = strsplit(input$cohort_genes, "\\s")[[1]]
		samples = parse_samples(input$cohort_patients)	

		if(length(genes) == 0) {
			output$cohort_mutations_oncoprint_ui = renderUI({
				p("No gene is found.")
			})
			return(NULL)
		}
		if(length(samples) == 0) {
			output$cohort_mutations_oncoprint_ui = renderUI({
				p("No sample is found.")
			})
			return(NULL)
		}

		output$cohort_mutations_oncoprint_ui = renderUI({
			div(
				div(
					column(6,
						selectInput("cohort_mutations_oncoprint_select", label = "Mutation types", 
							choices = structure(c("snv", "snv_germline", "indel", "indel_germline", "fusion"), 
								names = c("snv", "snv_germline", "indel", "indel_germline", "fusion")),
							selected = c("snv", "snv_germline", "indel", "indel_germline", "fusion"), 
							multiple = TRUE, width = 500),
					),
					column(2,
						actionButton(qq("cohort_mutations_oncoprint_submit"), "Make oncoPrint", style="margin-top:24px")
					),
					div(style = "clear:both;"),
					style = "width:900px;margin-left:-15px;"
				),
				originalHeatmapOutput(heatmap_id = "cohort_mutations_oncoprint", response = c("click"), width = 600, height = 600),
	   			HeatmapInfoOutput(heatmap_id = "cohort_mutations_oncoprint", output_ui = htmlOutput("cohort_mutations_oncoprint_output"), width = "90%")
	   		)
		})

		output$cohort_cnv_oncoprint_ui = renderUI({
			div(
				div(
					column(6,
						selectInput("cohort_cnv_oncoprint_select", label = "CNV types", 
							choices = structure(c("DUP", "AMP", "DEL", "HDEL", "CNN"), 
								names = c("DUP", "AMP", "DEL", "HDEL", "CNN")),
							selected = c("DUP", "AMP", "DEL", "HDEL", "CNN"), 
							multiple = TRUE, width = 500),
					),
					column(2,
						actionButton(qq("cohort_cnv_oncoprint_submit"), "Make oncoPrint", style="margin-top:24px")
					),
					div(style = "clear:both;"),
					style = "width:900px;margin-left:-15px;"
				),
				originalHeatmapOutput(heatmap_id = "cohort_cnv_oncoprint", response = c("click"), width = 600, height = 600),
	   			HeatmapInfoOutput(heatmap_id = "cohort_cnv_oncoprint", output_ui = htmlOutput("cohort_cnv_oncoprint_output"), width = "90%")
	   		)
		})

		output$cohort_expression_ui = renderUI({
			div(
				selectInput("cohort_expression_type", label = "Type", choices = c("TPM" = "TPM", "FPKM" = "FPKM"), selected = "TPM"),
				plotOutput("cohort_expression_plot")
			)
		})

		output[["cohort_expression_plot"]] = renderPlot({
			data_type = input[["cohort_expression_type"]]

			showNotification(qq("make heatmap."), duration = 4, type = "message")

			rg = rowRanges(DB[["rna"]])
			m = assays(DB[["rna"]])[[data_type]]

			lg = rg$Gene %in% genes
			m2 = m[rg$Gene %in% genes, colnames(m) %in% samples, drop = FALSE]
			m2 = log2(m2+1)
			Heatmap(m2, row_labels = rg$Gene[lg], name = qq("log2(@{data_type}+1)"))
		})
	})

	observeEvent(input$cohort_mutations_oncoprint_submit,  {
		genes = strsplit(input$cohort_genes, "\\s")[[1]]
		samples = parse_samples(input$cohort_patients)

		# oncoprint for mutations
		experiments = input$cohort_mutations_oncoprint_select
		if(length(experiments)) {
			generate_interactive_oncoprint(input, output, session, genes, samples, experiments)
		}

	})

	observeEvent(input$cohort_cnv_oncoprint_submit,  {
		genes = strsplit(input$cohort_genes, "\\s")[[1]]
		samples = parse_samples(input$cohort_patients)
	
		cnv_types = input$cohort_cnv_oncoprint_select
		if(length(cnv_types) && length(samples)) {
			generate_interactive_oncoprint_cnv(input, output, session, genes, samples, cnv_types)
		}

	})
}


generate_interactive_oncoprint = function(input, output, session, genes, samples, experiments) {
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
				if(e == "fusion") {
					ind = sapply(genes, function(x) any(grepl(paste0("\\b", x, "\\b"), rl[[sn]]$Gene)))
				} else {
					ind = genes %in% rl[[sn]]$Gene
				}
				ml[[e]][ind, sn] = TRUE
			}
		}
	}

	ht = oncoPrint(ml, name = "oncoprint", alter_fun = ALTER_FUN, col=ALTER_COL, 
		show_column_names = TRUE, remove_empty_columns = TRUE)

	click_action = function(df, output) {
		if(is.null(df)) {
			return(NULL)
		}
		gene = NULL
		if(df$heatmap == "oncoprint") {
			if(inherits(ht, "Heatmap")) {
				gene = rownames(ht@matrix)[df$row_index]
				sample = colnames(ht@matrix)[df$column_index]
				value = ht@matrix[df$row_index, df$column_index]
			} else {
				if(df$heatmap == "oncoprint") {
					ht2 = ht@ht_list[[1]]
					gene = rownames(ht2@matrix)[df$row_index]
					sample = colnames(ht2@matrix)[df$column_index]
					value = ht2@matrix[df$row_index, df$column_index]
				}
			}
		}
		if(!is.null(gene)) {
			output$cohort_mutations_oncoprint_output = renderUI({
				if(value != "") {
					experiments = strsplit(value, ";")[[1]]
					grl = lapply(experiments, function(e) {
						gr = DB[[e]]@assays[[sample]]
						if(e == "fusion") {
							gr[grepl(paste0("\\b", gene, "\\b"), gr$Gene)]
						} else {
							gr[gr$Gene == gene]
						}
					})
					names(grl) = experiments
					grl = lapply(grl, as.data.frame)

					tbl = list()
					for(i in seq_along(experiments)) {
						tb = grl[[i]]
						tb = t(as.matrix(tb))
						rownames(tb) = paste0("<b>", rownames(tb), "</b>")
						tbl[[i]] = box(
							title = paste0(experiments[i], ": ", sample),
							width = 4,
							HTML(knitr::kable(tb, format = "html", escape = FALSE)),
							style = "overflow-x: scroll; border:#CCC 1px solid; border-radius: 4px"
						)
					}

					do.call(fluidRow, tbl)
				} else {
					p("No alteration.")
				}
			})	
		}
	}

	makeInteractiveComplexHeatmap(input, output, session, ht, 
		heatmap_id = "cohort_mutations_oncoprint",
		click_action = click_action)
}

# this function has a large duplication in the previous one. It can be optimized in the future
generate_interactive_oncoprint_cnv = function(input, output, session, genes, samples, cnv_types) {
	full_samples = samples
	full_genes = genes

	samples = intersect(names(DB[["cnv"]]@assays), full_samples)
	genes = intersect(GENCODE$gene_name, full_genes)
	cnv_lt = DB[["cnv"]]@assays[samples]
	cnv_gr = unlist(cnv_lt)
	cnv_gr = cnv_gr[cnv_gr$Type %in% cnv_types]
	gene_gr = GENCODE[GENCODE$gene_name %in% genes]

	ov = findOverlaps(cnv_gr, gene_gr)

	ovdf = data.frame(gene = gene_gr$gene_name[subjectHits(ov)], sample = names(cnv_gr)[queryHits(ov)], cnv_type = cnv_gr$Type[queryHits(ov)])
	ovdf = unique(ovdf)

	ml = list()
	for(t in cnv_types) {
		l = ovdf[, 3] == t
		if(any(l)) {
			ml[[t]] = matrix(0, nrow = length(full_genes), ncol = length(full_samples), dimnames = list(full_genes, full_samples))
			tbb = xtabs(~ gene + sample, data = ovdf[l, , drop = FALSE])
			ml[[t]][rownames(tbb), colnames(tbb)] = tbb
		}
	}

	ht = oncoPrint(ml, name = "oncoprint", alter_fun = ALTER_FUN, col=ALTER_COL, 
		show_column_names = TRUE, remove_empty_columns = TRUE)

	click_action = function(df, output) {
		if(is.null(df)) {
			return(NULL)
		}
		gene = NULL
		if(df$heatmap == "oncoprint") {
			if(inherits(ht, "Heatmap")) {
				gene = rownames(ht@matrix)[df$row_index]
				sample = colnames(ht@matrix)[df$column_index]
				value = ht@matrix[df$row_index, df$column_index]
			} else {
				if(df$heatmap == "oncoprint") {
					ht2 = ht@ht_list[[1]]
					gene = rownames(ht2@matrix)[df$row_index]
					sample = colnames(ht2@matrix)[df$column_index]
					value = ht2@matrix[df$row_index, df$column_index]
				}
			}
		}
		if(!is.null(gene)) {
			output$cohort_cnv_oncoprint_output = renderUI({
				if(value != "") {
					cnv_gr = DB[["cnv"]]@assays[[sample]]
					cnv_gr = cnv_gr[cnv_gr$Type %in% cnv_types]
					g_gr = GENCODE[GENCODE$gene_name == gene]

					ov = findOverlaps(cnv_gr, g_gr)

					tb = as.data.frame(cnv_gr[unique(queryHits(ov))])
					tb = t(as.matrix(tb))
					rownames(tb) = paste0("<b>", rownames(tb), "</b>")
					fluidRow(
						box(title = qq("CNV: @{gene} in @{sample}"),
							width = 4,
							HTML(knitr::kable(tb, format = "html", escape = FALSE)),
							style = "overflow-x: scroll; border:#CCC 1px solid; border-radius: 4px"
						)
					)
					
				} else {
					p("No CNV")
				}
			})	
		}
	}

	makeInteractiveComplexHeatmap(input, output, session, ht, 
		heatmap_id = "cohort_cnv_oncoprint",
		click_action = click_action)
}
