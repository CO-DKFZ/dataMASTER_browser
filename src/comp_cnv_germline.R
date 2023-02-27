

experiment_id = "cnv_germline"
experiment_name = "cnv_germline"

COMPONENT[[experiment_id]] = list()

df = as.data.frame(DB[[experiment_id]])

cnv_type = unique(as.vector(df$Type))
cnv_type = cnv_type[!is.na(cnv_type)]
names(cnv_type) = cnv_type


nr = elementNROWS(DB[[experiment_id]]@assays)

COMPONENT[[experiment_id]]$ui = div(
	id = qq("@{experiment_id}_workspace"),
	h3(qq("Experiment: @{experiment_name}")),
	div(
		column(7,
			selectInput(qq("@{experiment_id}_config_function_type"), label = "Type", 
				choices = cnv_type, selected = cnv_type, 
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
		tabPanel("CNV plot",
			uiOutput(qq("@{experiment_id}_cnv_plot_ui")),
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
			updateSelectInput(session, qq("@{experiment_id}_config_function_type"), selected = c("DUP", "DEL"))
			function_type = c("DUP", "DEL")
		} else {
			function_type = input[[qq("@{experiment_id}_config_function_type")]]
		}

		df = df[df$Type %in% function_type, , drop = FALSE]
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

		updateTabItems(session, "sidebar_menu", experiment_id)
		updateTextInput(session, qq("@{experiment_id}_config_sample_list"), value = input[[qq("@{experiment_id}_config_sample_list")]])

		if(length(input[[qq("@{experiment_id}_config_function_type")]]) == 0) {
			function_type = "nonsynonymous"
		} else {
			function_type = input[[qq("@{experiment_id}_config_function_type")]]
		}

		df = df[df$Type %in% function_type, , drop = FALSE]
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
		
		if(length(input[[qq("@{experiment_id}_config_function_type")]]) == 0) {
			updateSelectInput(session, qq("@{experiment_id}_config_function_type"), selected = "nonsynonymous")
			function_type = "nonsynonymous"
		} else {
			function_type = input[[qq("@{experiment_id}_config_function_type")]]
		}

		samples = parse_samples(input[[qq("@{experiment_id}_config_sample_list")]], experiment_id)
		df = df[df$Type %in% function_type, , drop = FALSE]

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
				box(title = qq("Width distribution of @{experiment_name}s"),
					plotOutput(qq("@{experiment_id}_summary_width_dist"))
				),
				
				box(title = qq("Number of @{experiment_name}s per chromosome"),
					plotOutput(qq("@{experiment_id}_summary_n_snv_by_chr")),
					width = 12
				),
				box(title = "Confidence",
					plotOutput(qq("@{experiment_id}_summary_Confidence"))
				),
				box(title = "Function",
					plotOutput(qq("@{experiment_id}_summary_function"))
				),
				
				box(title = "Distribution of CN",
					plotOutput(qq("@{experiment_id}_summary_CN"))
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

		output[[qq("@{experiment_id}_summary_width_dist")]] = renderPlot({
			df = df[df$group_name %in% samples, , drop = FALSE]
			par(mar = c(4, 4, 1, 1))
			plot(density(log10(df$width)), xlab = "log10 width", ylab = "Density")
		})

		output[[qq("@{experiment_id}_summary_n_snv_by_chr")]] = renderPlot({
			df = df[df$group_name %in% samples, , drop = FALSE]

			par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
			x = table(as.vector(df$seqnames))
			barplot(x, ylab = qq("# @{experiment_name}s"))

			chr_len = getChromInfoFromUCSC(GENOME)
			chr_len = structure(chr_len[, 2], names = gsub("chr", "", chr_len[, 1]))
			chr_len = chr_len[ names(x) ]

			fc = (x/sum(x)) / (chr_len/sum(chr_len))
			fc = structure(as.vector(fc), names = names(fc))
			par(mar = c(4, 6, 1, 1))
            plot(log2(fc), type = "h", ylab = TeX("log2(fold enrichment): log2 ( \\frac{(n_sv/n_all_sv}{(chr_len/genome_len)} )"), axes = FALSE)			
            axis(side = 1, at = seq_along(fc), labels = names(fc))
			axis(side = 2)
			graphics::box()
		})

		output[[qq("@{experiment_id}_summary_Confidence")]] = renderPlot({
			df = df[df$group_name %in% samples, , drop = FALSE]
			par(mar = c(4, 4, 1, 1))
			pie(table(df$Confidence))
		})

		output[[qq("@{experiment_id}_summary_function")]] = renderPlot({
			df = df[df$group_name %in% samples, , drop = FALSE]
			par(mar = c(4, 4, 1, 1))
			pie(table(df$Type))
		})

		output[[qq("@{experiment_id}_summary_CN")]] = renderPlot({
			df = df[df$group_name %in% samples, , drop = FALSE]
			par(mar = c(4, 4, 1, 1))
			plot(table(df$CN), xlab = "CN", ylab = "Count", main = "")
		})

		output[[qq("@{experiment_id}_genome_wide_ui")]] = renderUI({
			div(
				div(
					column(9, selectInput(qq("@{experiment_id}_heatmap_chr"), label = "Select chromosomes", 
						choices = structure(CHROMOSOME, names = CHROMOSOME), selected = CHROMOSOME, 
						multiple = TRUE, width = "100%")),
					column(3, actionButton(qq("@{experiment_id}_heatmap_chr_submit"), "Make heatmap", style="margin-top:24px")),
					div(style = "clear:both;"),
					style = "width:800px"
				),
				box(title = qq("Genomic heatmap of all @{experiment_name}s in selected samples"),
					plotOutput(qq("@{experiment_id}_heatmap"), height = 800),
					width = 12
				),
				div(style = "clear:both;")
			)
		})

		output[[qq("@{experiment_id}_heatmap")]] = renderPlot({
			showNotification(qq("make cnv heatmap plot for @{experiment_id}."), duration = 4, type = "message")
			df = df[df$group_name %in% samples, , drop = FALSE]
			m = CNV_GERMLINE_NORMALIZED$matrix
			chr_window = CNV_GERMLINE_NORMALIZED$chr_window

			all_chr = as.vector(seqnames(chr_window))
			l_chr = all_chr %in% CHROMOSOME

			m = m[rownames(m) %in% samples, l_chr, drop = FALSE]

			if(nrow(m) > 1000) {
				m = m[sort(sample(nrow(m), 1000)), , drop = FALSE]
			}

			ht = Heatmap(m, name = "CN", col = colorRamp2(c(0, 2, 4, 6), c("blue", "white", "pink", "red")),
				cluster_rows = FALSE, cluster_columns = FALSE, 
				column_split = factor(all_chr[l_chr], levels = intersect(CHROMOSOME, all_chr[l_chr])),
				column_gap = unit(0, "points"), border = TRUE, column_title_gp = gpar(fontsize = 10),
				show_row_names = FALSE)
			draw(ht)

		})

		observeEvent(input[[qq("@{experiment_id}_heatmap_chr_submit")]], {
			showNotification(qq("make cnv heatmap plot for @{experiment_id}."), duration = 4, type = "message")
			df = df[df$group_name %in% samples, , drop = FALSE]
			m = CNV_GERMLINE_NORMALIZED$matrix
			chr_window = CNV_GERMLINE_NORMALIZED$chr_window

			all_chr = as.vector(seqnames(chr_window))
			l_chr = all_chr %in% input[[qq("@{experiment_id}_heatmap_chr")]]

			m = m[rownames(m) %in% samples, l_chr, drop = FALSE]

			if(nrow(m) > 1000) {
				m = m[sort(sample(nrow(m), 1000)), , drop = FALSE]
			}
			output[[qq("@{experiment_id}_heatmap")]] = renderPlot({
				ht = Heatmap(m, name = "CN", col = colorRamp2(c(0, 2, 4, 6), c("blue", "white", "pink", "red")),
					cluster_rows = FALSE, cluster_columns = FALSE, 
					column_split = factor(all_chr[l_chr], levels = intersect(CHROMOSOME, all_chr[l_chr])),
					column_gap = unit(0, "points"), border = TRUE, column_title_gp = gpar(fontsize = 10),
					show_row_names = FALSE)
				draw(ht)
			})
		})

		output[[qq("@{experiment_id}_cnv_plot_ui")]] = renderUI({
			div(
				div(
					column(8, textInput(qq("@{experiment_id}_cnv_plot_sample"), "Select a sample", placeholder = samples[1])),
					column(2, actionButton(qq("@{experiment_id}_cnv_plot_sample_submit"), "Make plot", style="margin-top:24px")),
					div(style = "clear:both;"),
					style = "width:500px;margin-left:-15px;"
				),
				div(
					plotOutput(qq("@{experiment_id}_cnv_plot"), height = 500),
					style = "margin:10px 0px; padding: 10px; background-color:white;"
				),
				tags$script(HTML(qq('
	$("#@{experiment_id}_cnv_plot_sample").typeahead({
	    source: @{jsonlite::toJSON(samples)},
	    lookup: function (event) {
	        var that = this, items;
	        if (that.ajax) {
	            that.ajaxer();
	        } else {
	            that.query = that.$element.val();

	            if (!that.query) {
	                return that.shown ? that.hide() : that;
	            }

	            items = that.grepper(that.source);

	            if (!items) {
	                return that.shown ? that.hide() : that;
	            }
	            return that.render(items.slice(0, that.options.items)).show();
	        }
	    }
	});
				')))
			)
		})

		observeEvent(input[[qq("@{experiment_id}_cnv_plot_sample_submit")]], {
			sample = input[[qq("@{experiment_id}_cnv_plot_sample")]]

			gr = DB[[experiment_id]]@assays[[sample]]
			if(is.null(gr)) {
				output[[qq("@{experiment_id}_cnv_plot")]] = renderPlot({
					plot(NULL, xlim = c(0, 1), ylim = c(0, 1), type = "n", axes = FALSE, ann = FASE)
					text(0.5, 0.5, "No sample is found.", cex = 3)
				})
				return(NULL)
			}
			gr = gr[width(gr) > 100]

			seqlevelsStyle(gr) = "UCSC"

			output[[qq("@{experiment_id}_cnv_plot")]] = renderPlot({
				
				max = max(gr$CN)*1.1

				gtrellis_layout(n_track = 1, nrow = 4, compact = TRUE, species = GENOME, 
				    track_ylim = c(0, max),
				    track_ylab = "CN", title = sample,
				    add_name_track = TRUE, add_ideogram_track = TRUE)

				add_segments_track(gr, gr$CN)
				add_points_track(gr, gr$CN, track = current_track())
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
			rownames(df) = NULL
			df
		}, server = TRUE, options = list(scrollX = TRUE))

		observeEvent(input[[qq("@{experiment_id}_table_gr_range_submit")]], {
			df = df[df$group_name %in% samples, , drop = FALSE]
			df = df[, -1]
			colnames(df)[1:2] = c("Sample", "Chr")
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
