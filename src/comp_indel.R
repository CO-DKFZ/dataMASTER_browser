
experiment_id = "indel"
experiment_name = "Indel"

COMPONENT[[experiment_id]] = list()

df = as.data.frame(DB[["indel"]])

indel_function_type = unique(as.vector(df$Function))
indel_function_type = indel_function_type[!is.na(indel_function_type)]
names(indel_function_type) = indel_function_type


nr = elementNROWS(DB[[experiment_id]]@assays)

COMPONENT[[experiment_id]]$ui = div(
	id = qq("@{experiment_id}_workspace"),
	h3(qq("Experiment: @{experiment_name}")),
	div(
		column(7,
			selectInput(qq("@{experiment_id}_config_function_type"), label = "Function type", 
				choices = indel_function_type, selected = setdiff(indel_function_type, c("unknown", "synonymous")), 
				multiple = TRUE, width = 600),
			tags$hr(),
			uiOutput(qq("@{experiment_id}_sample_filter_output")),
		),
		column(5,
			textAreaInput(qq("@{experiment_id}_config_sample_list"), label = "Patient/Sample list", 
				rows = 10, width = 400),
			tags$hr(),
			textOutput(qq("@{experiment_id}_config_nsample_selected")),
			actionButton(qq("@{experiment_id}_config_submit"), label = "Submit")
		),
		div(style = "clear:both;"),
		style = "width:1000px"
	),
	br(),
	tabsetPanel(
		id = qq("@{experiment_id}_tabset"),
		type = "tabs",
		tabPanel(qq("Summary"),
			uiOutput(qq("@{experiment_id}_summary_ui")),
			style = "padding:16px 0px"
		),
		tabPanel("Genome-wide distribution",
			uiOutput(qq("@{experiment_id}_genome_wide_ui")),
			style = "padding:16px 0px"
		),
		tabPanel("Top most mutated genes",
			uiOutput(qq("@{experiment_id}_top_genes_ui")),
			style = "padding:16px 0px"
		),
		tabPanel("Top most mutated samples",
			uiOutput(qq("@{experiment_id}_top_samples_ui")),
			style = "padding:16px 0px"
		),
		tabPanel(qq("@{experiment_name} table"),
			uiOutput(qq("@{experiment_id}_table_ui")),
			style = "padding:16px 0px"
		)
	)
)

COMPONENT[[experiment_id]]$server = local({

	experiment_id = experiment_id
	experiment_name = experiment_name
	
	function(input, output, session) {

	obj = EXPERIMENTS[[experiment_id]]
	df = as.data.frame(obj)
	
	CD = CD_LIST[[experiment_id]]

	observeEvent(input[[qq("@{experiment_id}_config_function_type")]], {

		if(length(input[[qq("@{experiment_id}_config_function_type")]]) == 0) {
			updateSelectInput(session, qq("@{experiment_id}_config_function_type"), selected = "nonsynonymous")
			function_type = "nonsynonymous"
		} else {
			function_type = input[[qq("@{experiment_id}_config_function_type")]]
		}

		df = df[df$Function %in% function_type, , drop = FALSE]
		l = rownames(CD) %in% unique(df$group_name)
		
		updateTextAreaInput(session, qq("@{experiment_id}_config_sample_list"),
			value = paste(rownames(CD)[l], collapse = "\n"))
		output[[qq("@{experiment_id}_config_nsample_selected")]] = renderText({
			paste0(sum(l), " samples selected")
		})
		output[[qq("@{experiment_id}_sample_filter_output")]] = renderUI(NULL)

		l2 = sapply(CD_SELECTED, function(x) length(unique(META[[x]][l])) > 1)
		CD_SELECTED = CD_SELECTED[l2]
		output[[qq("@{experiment_id}_sample_filter_output")]] = renderUI({
			div(
				column(5,
					selectInput(qq("@{experiment_id}_sample_filter_field"), label = "Filter samples by", 
								choices = structure(CD_SELECTED, names = CD_SELECTED), 
								selected = NULL, multiple = TRUE)
				),
				column(7,
					uiOutput(qq("@{experiment_id}_sample_filter_output_filed"))
				),
				div(style = "clear:both;"),
				style = "border:1px solid #CCC;padding:8px 0px;border-radius:5px;"
			)
		})

		output[[qq("@{experiment_id}_sample_filter_output_filed")]] = renderUI(NULL)
	})

	observeEvent(input[[qq("@{experiment_id}_sample_filter_field")]], {
		fields = input[[qq("@{experiment_id}_sample_filter_field")]]
		if(length(fields) == 0) {
			output[[qq("@{experiment_id}_sample_filter_output_filed")]] = renderUI(NULL)
		} else {
			ui_list = list()
			for(nm in fields) {
				if(is.numeric(CD[[nm]])) {
					rg = range(CD[[nm]], na.rm = TRUE)
					if(length(input[[qq("@{experiment_id}_sample_filter_@{nm}")]])) {
						value = input[[qq("@{experiment_id}_sample_filter_@{nm}")]]
					} else {
						value = rg
					}
					ui_list[[nm]] = sliderInput(qq("@{experiment_id}_sample_filter_@{nm}"), 
						label = nm, 
						min = rg[1], max = rg[2], 
						value = value)
				} else {
					le = unique(CD[[nm]])
					selected = input[[qq("@{experiment_id}_sample_filter_@{nm}")]]
					ui_list[[nm]] = selectInput(qq("@{experiment_id}_sample_filter_@{nm}"),
						label = nm,
						choices = structure(le, names = le), selected = selected, multiple = TRUE)
				}
			}
			ui_list = c(ui_list, list(actionButton(qq("@{experiment_id}_sample_filter_field_submit"), label = "Apply sample filter")))

			names(ui_list) = NULL
			output[[qq("@{experiment_id}_sample_filter_output_filed")]] = renderUI({
				do.call(div, ui_list)
			})
		}
	}, ignoreNULL = FALSE)

	observeEvent(input[[qq("@{experiment_id}_sample_filter_field_submit")]], {

		if(length(input[[qq("@{experiment_id}_config_function_type")]]) == 0) {
			function_type = "nonsynonymous"
		} else {
			function_type = input[[qq("@{experiment_id}_config_function_type")]]
		}

		df = df[df$Function %in% function_type, , drop = FALSE]
		l = rownames(CD) %in% unique(df$group_name)
		
		for(nm in CD_SELECTED) {
			if(length(input[[qq("@{experiment_id}_sample_filter_@{nm}")]]) > 0) {
				if(is.numeric(CD[[nm]])) {
					rg = input[[qq("@{experiment_id}_sample_filter_@{nm}")]]
					l2 = CD[[nm]] >= rg[1] & CD[[nm]] <= rg[2]
					l2[is.na(l2)] = FALSE
					l = l & l2
				} else {
					v = input[[qq("@{experiment_id}_sample_filter_@{nm}")]]
					if(length(v)) {
						l2 = CD[[nm]] %in% v
						l2[is.na(l2)] = FALSE
						l = l & l2
					}
				}
			}
		}
		updateTextAreaInput(session, qq("@{experiment_id}_config_sample_list"),
			value = paste(rownames(CD)[l], collapse = "\n"))
		output[[qq("@{experiment_id}_config_nsample_selected")]] = renderText({
			paste0(sum(l), " samples selected")
		})
	})

	observeEvent(input[[qq("@{experiment_id}_config_sample_list")]], {
		samples = parse_samples(input[[qq("@{experiment_id}_config_sample_list")]], experiment_id)

		output[[qq("@{experiment_id}_config_nsample_selected")]] = renderText({
			paste0(length(samples), " samples selected")
		})
	})

	observeEvent(input[[qq("@{experiment_id}_config_submit")]], {
		
		updateTabItems(session, "sidebar_menu", experiment_id)
		updateTextInput(session, qq("@{experiment_id}_config_sample_list"), value = input[[qq("@{experiment_id}_config_sample_list")]])

		if(length(input[[qq("@{experiment_id}_config_function_type")]]) == 0) {
			updateSelectInput(session, qq("@{experiment_id}_config_function_type"), selected = "nonsynonymous")
			function_type = "nonsynonymous"
		} else {
			function_type = input[[qq("@{experiment_id}_config_function_type")]]
		}

		samples = parse_samples(input[[qq("@{experiment_id}_config_sample_list")]], experiment_id)
		df = df[df$Function %in% function_type, , drop = FALSE]

		if(length(samples) == 0) {
			output[[qq("@{experiment_id}_summary_ui")]] = renderUI({
				p("No sample is found.")
			})
			return(NULL)
		}

		output[[qq("@{experiment_id}_summary_ui")]] = renderUI({
			div(
				box(title = qq("Number of @{experiment_name}s per sample"),
					plotOutput(qq("@{experiment_id}_summary_n_snv"))
				),
				box(title = "Number of mutated genes per sample",
					plotOutput(qq("@{experiment_id}_summary_n_gene"))
				),
				box(title = qq("Number of @{experiment_name}s per chromosome"),
					plotOutput(qq("@{experiment_id}_summary_n_snv_by_chr")),
					width = 12
				),
				box(title = "Genomic locations",
					plotOutput(qq("@{experiment_id}_summary_genomic_locations"))
				),
				box(title = "Function",
					plotOutput(qq("@{experiment_id}_summary_function"))
				),
				box(title = "Length of the deletion/insertion",
					plotOutput(qq("@{experiment_id}_summary_width")),
					width = 12
				),
				div(style="clear:both;")
			)
		})

		output[[qq("@{experiment_id}_summary_n_snv")]] = renderPlot({
			showNotification(qq("make overview plot for @{experiment_id}."), duration = 4, type = "message")
			n = tapply(df$group_name, df$group_name, length)
			n = sort(n)
			par(mar = c(4, 4, 1, 1))
			plot(n, pch = 16, cex = 0.5, log = "y", ylab = qq("# @{experiment_name}s (in log10 scale)"), xlab = "Ordered samples")
			
			ind = names(n) %in% samples
			points(which(ind), n[ind], pch = 16, cex = 1, col = "red")
			legend("topleft", pch = 16, legend = c("All samples", "Selected samples"), col = c("black", "red"), border = NA)
		})

		output[[qq("@{experiment_id}_summary_n_gene")]] = renderPlot({
			
			n = tapply(df$Gene, df$group_name, function(x) length(unique(x)))
			n = sort(n)
			par(mar = c(4, 4, 1, 1))
			plot(n, pch = 16, cex = 0.5, log = "y", ylab = "# mutated genes (in log10 scale)", xlab = "Ordered samples")
			
			ind = names(n) %in% samples
			points(which(ind), n[ind], pch = 16, cex = 1, col = "red")
			legend("topleft", pch = 16, legend = c("All samples", "Selected samples"), col = c("black", "red"), border = NA)
		})

		output[[qq("@{experiment_id}_summary_n_snv_by_chr")]] = renderPlot({
			df = df[df$group_name %in% samples, , drop = FALSE]

			par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
			x = table(df$seqnames)
			barplot(x, ylab = qq("# @{experiment_name}s"))

			chr_len = getChromInfoFromUCSC(GENOME)
			chr_len = structure(chr_len[, 2], names = gsub("chr", "", chr_len[, 1]))
			chr_len = chr_len[ names(x) ]

			fc = (x/sum(x)) / (chr_len/sum(chr_len))
			fc = structure(as.vector(fc), names = names(fc))
			par(mar = c(4, 6, 1, 1))
            plot(log2(fc), type = "h", ylab = TeX("log2(fold enrichment): log2 ( \\frac{(n_indel/n_all_indel}{(chr_len/genome_len)} )"), axes = FALSE)			
            axis(side = 1, at = seq_along(fc), labels = names(fc))
			axis(side = 2)
			graphics::box()
		})

		output[[qq("@{experiment_id}_summary_genomic_locations")]] = renderPlot({
			df = df[df$group_name %in% samples, , drop = FALSE]
			par(mar = c(4, 4, 1, 1))
			pie(table(df$Location))
		})

		output[[qq("@{experiment_id}_summary_function")]] = renderPlot({
			df = df[df$group_name %in% samples, , drop = FALSE]
			par(mar = c(4, 4, 1, 1))
			pie(table(df$Function))
		})

		output[[qq("@{experiment_id}_summary_width")]] = renderPlot({
			df = df[df$group_name %in% samples, , drop = FALSE]
			par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))

			plot(table(nchar(df$Ref)), xlab = "Number of bases in deletion", ylab = "Count")
			plot(table(nchar(df$Alt)), xlab = "Number of bases in insertion", ylab = "Count")
		})

		output[[qq("@{experiment_id}_genome_wide_ui")]] = renderUI({
			div(
				radioButtons(qq("@{experiment_id}_genome_wide_radio"), "Visualization layout", 
					choices = c("Rectangular" = "rect", "Circular" = "circular"), selected = "rect",
					inline = TRUE),
				box(title = qq("Rainfall plot of all @{experiment_name}s in selected samples"),
					plotOutput(qq("@{experiment_id}_rainfall"), height = 700),
					width = 10
				),
				div(style = "clear:both;")
			)
		})

		output[[qq("@{experiment_id}_rainfall")]] = renderPlot({
			showNotification(qq("make genome-wide plot for @{experiment_id}."), duration = 4, type = "message")
			df = df[df$group_name %in% samples, , drop = FALSE]
			make_rainfall(GRanges(seqnames = df$seqnames, ranges = IRanges(df$start, df$end)), layout = input[[qq("@{experiment_id}_genome_wide_radio")]])
		})

		output[[qq("@{experiment_id}_top_genes_ui")]] = renderUI({
			div(
				selectInput(qq("@{experiment_id}_top_genes_n_top"), label = "Number of top genes", choices = c("10" = 10, "20" = 20, "50" = 50, "100" = 100), selected = 20),
				box(title = "Top genes",
					tableOutput(qq("@{experiment_id}_top_genes")),
					width = 3
				),
				box(title = qq("OncoPrint of the top genes (in @{length(samples)} samples)"),
					plotOutput(qq("@{experiment_id}_top_genes_oncoprint")),
					width = 8
				),
				tags$script(HTML(qq('
	$("#@{experiment_id}_top_genes_n_top").change(function() {
		var n = parseInt($(this).val());
		$("#@{experiment_id}_top_genes_oncoprint").height(n*20);
	});'))),
				div(style = "clear:both;")
			)
		})

		observeEvent(input[[qq("@{experiment_id}_top_genes_n_top")]], {

			showNotification(qq("generate table and plot for top genes."), duration = 4, type = "message")
		
			df = df[df$group_name %in% samples, , drop = FALSE]
			tb = df[, c("group_name", "Gene")]
			xn = tapply(tb$group_name, tb$Gene, length) # number of SNVs per gene
			tb = unique(tb)

			x = tapply(tb$group_name, tb$Gene, length) # number of samples per gene
			p = x/length(samples)

			k = as.numeric(input[[qq("@{experiment_id}_top_genes_n_top")]])
			ind = order(p, xn[names(p)], decreasing = TRUE)[1:min(length(p), k)]

			output[[qq("@{experiment_id}_top_genes")]] = renderTable({
				dt = data.frame(Gene = names(p)[ind], Fraction = paste0(round(p[ind]*100, 1), "%"), "#Samples" = x[ind], check.names = FALSE)
				dt[, 1] = qq("<a href=\"#\" onclick=\"link_to_other_tab($(this), 'region', {region_search:'@{dt[,1]}'}, 'region_search_submit');false\">@{dt[,1]}</a>", collapse = FALSE)
				dt
			}, sanitize.text.function = function(x) x)

			tb = tb[tb$Gene %in% names(p)[ind], , drop = FALSE]
			mm = xtabs(~ Gene + group_name, tb)

			if(ncol(mm) > 2000) {
				mm = mm[, sample(ncol(mm), 2000), drop = FALSE]
			}
			output[[qq("@{experiment_id}_top_genes_oncoprint")]] = renderPlot({
				ht = oncoPrint(list(indel = mm), name = experiment_name,
					alter_fun = list(indel = alter_graphic("rect", width = 1, height = 0.9, fill = "red")),
					show_column_names = FALSE)
				draw(ht)
			})
		})

		

		output[[qq("@{experiment_id}_top_samples_ui")]] = renderUI({
			div(
				selectInput(qq("@{experiment_id}_top_samples_n_top"), label = "Number of top samples", choices = c("10" = 10, "20" = 20, "50" = 50, "100" = 100), selected = 20),
				box(title = "Top samples",
					tableOutput(qq("@{experiment_id}_top_samples")),
					width = 4
				),
				box(title = qq("OncoPrint (samples are on rows, genes are on columns)"),
					plotOutput(qq("@{experiment_id}_top_samples_oncoprint")),
					width = 8
				),
				tags$script(HTML(qq('
	$("#@{experiment_id}_top_samples_n_top").change(function() {
		var n = parseInt($(this).val());
		$("#@{experiment_id}_top_samples_oncoprint").height(n*20);
	});'))),
				div(style = "clear:both;")
			)
		})

		
		observeEvent(input[[qq("@{experiment_id}_top_samples_n_top")]], {

			showNotification(qq("generate table and plot for top samples."), duration = 4, type = "message")

			df = df[df$group_name %in% samples, , drop = FALSE]
			tb = df[, c("group_name", "Gene")]
			xn = tapply(tb$group_name, tb$group_name, length) # number of SNVs per sample
			
			tb = unique(tb)
			x = tapply(tb$Gene, tb$group_name, length) # number of genes per sample
			
			k = input[[qq("@{experiment_id}_top_samples_n_top")]]
			ind = order(xn, x, decreasing = TRUE)[1:min(length(xn), k)]
			
			output[[qq("@{experiment_id}_top_samples")]] = renderTable({
				dt = data.frame(Sample = names(x)[ind], "#Indels" = xn[ind], "#Mutated genes" = x[ind], check.names = FALSE)
				# dt[, 1] = qq("<a href=\"#\" onclick=\"link_to_other_tab($(this), 'gene', {gene_search:'@{dt[,1]}'}, 'gene_search_submit');false\">@{dt[,1]}</a>", collapse = FALSE)
				dt
			}, sanitize.text.function = function(x) x)

			tb = tb[tb$group_name %in% names(x)[ind], , drop = FALSE]
			mm = xtabs(~ group_name + Gene, tb)

			if(ncol(mm) > 2000) {
				mm = mm[, sample(ncol(mm), 2000), drop = FALSE]
			}
			output[[qq("@{experiment_id}_top_samples_oncoprint")]] = renderPlot({
				ht = oncoPrint(list(indel = mm), name = experiment_name,
					alter_fun = list(indel = alter_graphic("rect", width = 1, height = 0.9, fill = "red")),
					show_column_names = FALSE)
				draw(ht)
			})
		})


		output[[qq("@{experiment_id}_table_ui")]] = renderUI({
			div(
				div(
					column(8, textInput(qq("@{experiment_id}_table_gr_range"), "Genomic ranges", placeholder = "1:1000-2000")),
					column(2, actionButton(qq("@{experiment_id}_table_gr_range_submit"), "Query", style="margin-top:24px")),
					div(style = "clear:both;"),
					style = "width:500px;margin-left:-15px;"
				),
				div(
					DTOutput(qq("@{experiment_id}_table")),
					style = "margin:10px 15px; padding: 10px; background-color:white;"
				)
			)
		})

		output[[qq("@{experiment_id}_table")]] = renderDT({
			df = df[df$group_name %in% samples, , drop = FALSE]
			df = df[, -1]
			colnames(df)[1:2] = c("Sample", "Chr")
			df[, "AF_T"] = round(df[, "AF_T"], 2)
			df[, "AF_C"] = round(df[, "AF_C"], 2)
			rownames(df) = NULL
			df
		}, server = TRUE, options = list(scrollX = TRUE))

		observeEvent(input[[qq("@{experiment_id}_table_gr_range_submit")]], {
			df = df[df$group_name %in% samples, , drop = FALSE]
			df = df[, -1]
			colnames(df)[1:2] = c("Sample", "Chr")
			df[, "AF_T"] = round(df[, "AF_T"], 2)
			df[, "AF_C"] = round(df[, "AF_C"], 2)
			rownames(df) = NULL
			df

			lt = parse_gr_str(input[[qq("@{experiment_id}_table_gr_range")]])

			if(length(lt) == 0) {
				output[[qq("@{experiment_id}_table")]] = renderDT(NULL)
			} else {

				if(length(lt) == 1) {
					df2 = df[df$Chr == lt$chr, , drop = FALSE]
				} else if(length(lt) == 2) {
					df2 = df[df$Chr == lt$chr & df$start == lt$start, ,drop = FALSE]
				} else  if(length(lt) == 3) {
					df2 = df[df$Chr == lt$chr & df$start >= lt$start & df$end <= lt$end, ,drop = FALSE]
				}
				output[[qq("@{experiment_id}_table")]] = renderDT({
					df2
				}, server = TRUE, options = list(scrollX = TRUE))
			}
		})

	})
}

})
