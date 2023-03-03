

experiment_id = "rna"
experiment_name = "RNA"

COMPONENT[[experiment_id]] = list()

COMPONENT[[experiment_id]]$ui = div(
	id = qq("@{experiment_id}_workspace"),
	h3(qq("Experiment: @{experiment_name}")),
	div(
		column(7,
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
		),
		column(5,
			radioButtons(qq("@{experiment_id}_config_expr_type"), label = "Date type", choices = c("TPM" = "TPM", "FPKM" = "FPKM"), selected = "TPM", inline = TRUE),
			textAreaInput(qq("@{experiment_id}_config_sample_list"), label = "Patient/Sample list", 
				value = paste(rownames(CD_LIST[[experiment_id]]), collapse = "\n"),
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
		tabPanel(qq("Selected genes"),
			uiOutput(qq("@{experiment_id}_selected_genes_ui")),
			style = "padding:16px 0px"
		),
		tabPanel(qq("Dimension reduction"),
			uiOutput(qq("@{experiment_id}_dimention_reduction_ui")),
			style = "padding:16px 0px"
		)
	)
)

COMPONENT[[experiment_id]]$server = local({

	experiment_id = experiment_id
	experiment_name = experiment_name

	l_pc = GENCODE$gene_type == "protein_coding"
	
	function(input, output, session) {

	obj = EXPERIMENTS[[experiment_id]]
	
	CD = CD_LIST[[experiment_id]]

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

		l = rep(TRUE, nrow(CD))
		
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

	cache = new.env()
	cache$dimension_reduction_pos = list()

	observeEvent(input[[qq("@{experiment_id}_config_submit")]], {
		data_type = input[[qq("@{experiment_id}_config_expr_type")]]
		samples = parse_samples(input[[qq("@{experiment_id}_config_sample_list")]], experiment_id)

		cache$dimension_reduction_pos = list()

		if(length(samples) == 0) {
			output[[qq("@{experiment_id}_summary_ui")]] = renderUI({
				p("No sample is found.")
			})
			return(NULL)
		}
		
		output[[qq("@{experiment_id}_summary_ui")]] = renderUI({
			div(
				box(title = qq("Global distribution"),
					plotOutput(qq("@{experiment_id}_summary_global_distribution"), height = 600),
					width = 12
				),
				div(style="clear:both;")
			)
		})

		output[[qq("@{experiment_id}_summary_global_distribution")]] = renderPlot({

			showNotification(qq("Global distributino of all samples. This may take long."), duration = 8, type = "message")

			mat = assay(DB[["rna"]], data_type)
			mat = mat[l_pc, samples, drop = FALSE]
			if(nrow(mat) > 2000) {
				mat = mat[sample(nrow(mat), 2000), , drop = FALSE]
			}
			if(ncol(mat) > 1000) {
				mat = mat[, sample(ncol(mat), 1000), drop = FALSE]
			}
			m2 =log2(mat+1)
			densityHeatmap(m2, ylab = qq("log2(@{data_type} + 1)"), show_column_names = FALSE,
				column_order = order(colMeans(m2)))
		})

		output[[qq("@{experiment_id}_selected_genes_ui")]] = renderUI({
			div(
				box(selectInput(qq("@{experiment_id}_selected_genes_which"), label = "Select genes", 
						choices = c("highest 1000 expressed genes" = 1, "median 1000 expressed genes" = 2, "lowest 1000 expressed genes" = 3, "top 1000 most variable genes" = 4),
						selected = 4, width = 300),
					plotOutput(qq("@{experiment_id}_selected_genes"), height = 600),
					width = 12
				),
				div(style="clear:both;")
			)
		})

		output[[qq("@{experiment_id}_selected_genes")]] = renderPlot({

			showNotification(qq("Heatmap of selected genes. This may take long."), duration = 8, type = "message")

			mat = assay(DB[["rna"]], data_type)
			mat = mat[l_pc, samples, drop = FALSE]
			mat = log2(mat + 1)
			e = rowMeans(mat)

			which = input[[qq("@{experiment_id}_selected_genes_which")]]

			if(which == 1) {
				ind = order(-e)[1:1000]
			} else if(which == 3) {
				ind = order(e)[1:1000]
			} else if(which == 3) {
				nr = nrow(mat)
				ind = order(e)[seq(round(nr/2)-500, round(nr/2)+500)]
			} else {
				ind = order(-rowSds(mat))[1:1000]
			}

			m2 = mat[ind, , drop = FALSE]
			if(ncol(m2) > 1000) {
				m2 = m2[, sample(ncol(m2), 1000), drop = FALSE]
			}
			Heatmap(m2, name = qq("log2(@{data_type}+1)"), show_row_names = FALSE, show_column_names = FALSE,
				show_row_dend = FALSE, show_column_dend = FALSE)

		})

		output[[qq("@{experiment_id}_dimention_reduction_ui")]] = renderUI({
			div(
				column(3,
					selectInput(qq("@{experiment_id}_dimension_reduction_method"), label = "Method", 
						choices = c("UMAP" = "UMAP", "t-SNE" = "t-SNE", "PCA" = "PCA"),
						selected = "UMAP", width = 300),
					actionButton(qq("@{experiment_id}_dimension_reduction_rerun"), label = "Force rerun")
				),
				column(3,
					selectInput(qq("@{experiment_id}_dimension_reduction_color_by"), label = "Color by meta column",
						choices = structure(c("none", CD_SELECTED), names = c("none", CD_SELECTED)),
						selected = "none", width = 300),


					HTML(qq('<div class="form-group" style="width:300px;">
  <label class="control-label" id="@{experiment_id}_dimension_reduction_color_by_expr-label" for="@{experiment_id}_dimension_reduction_color_by_expr">Or color by gene expression</label>
  <input id="@{experiment_id}_dimension_reduction_color_by_expr" type="text" class="form-control" value=""/>
</div>')),
					tags$script(HTML(qq('
	$("#@{experiment_id}_dimension_reduction_color_by").change(function() {
		$("#@{experiment_id}_dimension_reduction_color_by_expr").val("");
		Shiny.setInputValue("@{experiment_id}_dimension_reduction_color_by_expr_g", "");
	});
	$("#@{experiment_id}_dimension_reduction_color_by_expr").focusout(function() {
		Shiny.setInputValue("@{experiment_id}_dimension_reduction_color_by_expr_g", $("#@{experiment_id}_dimension_reduction_color_by_expr").val());
	});
				')))
				),
				div(style="clear:both;"),
				plotOutput(qq("@{experiment_id}_dimention_reduction"), height = 800, width = 800)
			)
		})

		cache$dimension_reduction_submit = -1

		observe({

			mat = assay(DB[["rna"]], data_type)
			mat = mat[l_pc, samples, drop = FALSE]
			mat = log2(mat + 1)
			ind = order(-rowSds(mat))[1:1000]
			m2 = mat[ind, , drop = FALSE]

			method = input[[qq("@{experiment_id}_dimension_reduction_method")]]
			color_by = input[[qq("@{experiment_id}_dimension_reduction_color_by")]]; if(is.null(color_by)) color_by = "none"
			gene = input[[qq("@{experiment_id}_dimension_reduction_color_by_expr_g")]]

			if(!(identical(gene, NULL) || identical(gene, ""))) {
				if(!gene %in% GENCODE$gene_name) {
					showNotification(qq("Cannot find gene '@{gene}'."), duration = 10, type = "error")
					return(NULL)
				}
			}

			col_v = NULL
			if(!(identical(gene, NULL) || identical(gene, ""))) {
				col_v = log2(mat[GENCODE$gene_id[GENCODE$gene_name == gene], ] + 1)
				col_v = as.vector(scale(col_v))
				color_by = "by_expr"
			} else {
				if(color_by != "none") {
					col_v = CD[samples, color_by]
				}
			}

			output[[qq("@{experiment_id}_dimention_reduction")]] = renderPlot({
				if(cache$dimension_reduction_submit != input[[qq("@{experiment_id}_dimension_reduction_rerun")]] || is.null(cache$dimension_reduction_pos[[method]])) {
					
					showNotification(qq("Dimensioin reduction. This may take long."), duration = ifelse(ncol(m2) > 500, 10, 5), type = "message")

					cache$dimension_reduction_submit = input[[qq("@{experiment_id}_dimension_reduction_rerun")]]
					if(color_by == "none") {
						pos = dimension_reduction(m2, method = method,
							scale = TRUE)
						cache$dimension_reduction_pos[[method]] = pos
					} else {
						col = ComplexHeatmap:::default_col(col_v, TRUE)
						if(is.function(col)) {
							col2 = col(col_v)
						} else {
							col2 = col[col_v]
						}
						col2[is.na(col2)] = "#EEEEEE"
						pos = dimension_reduction(m2, method = method,
							scale = TRUE, col = col2)
						cache$dimension_reduction_pos[[method]] = pos
					}
				} else {
					pos = cache$dimension_reduction_pos[[method]]
					if(color_by == "none") {
						plot(pos[, 1], pos[, 2], pch = 16, cex = 1, xlab = paste0(method, " 1"), ylab = paste0(method, " 2"), 
							main = qq("@{method} on a matrix with 1000 rows, rows are scaled, with 10 PCs"))
					} else {
						col = ComplexHeatmap:::default_col(col_v, TRUE)
						if(is.function(col)) {
							col2 = col(col_v)
						} else {
							col2 = col[col_v]
						}
						col2[is.na(col2)] = "#EEEEEE"
						plot(pos[, 1], pos[, 2], pch = 16, cex = 1, xlab = paste0(method, " 1"), ylab = paste0(method, " 2"), 
							main = qq("@{method} on a matrix with 1000 rows, rows are scaled, with 10 PCs"), col = col2)
					}
				}
			})
		})

	})
}

})
