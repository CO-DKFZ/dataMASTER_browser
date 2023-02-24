
### overview
CD = colData(DB)
CD_SELECTED = c("ProjectID", "TumorID", "ControlID", "DNASeq", "RNASample",
	            "Ploidy", "TumorCellContent", "ValidClinicalEval",
	            "SevereDegradation", "TissueMaterial", "LOH.HRD", "LST", "TAI",
	            "MSI")
CD_SELECTED = intersect(CD_SELECTED, colnames(CD))
box_list = lapply(CD_SELECTED, function(x) {
	if(is.numeric(CD[[x]])) {
		box(
			title= x,
			status = "primary",
			plotOutput(qq("overview_@{x}_plot")),
			width = 3,
			height = 600
		)
	} else if(is.character(CD[[x]]) || is.factor(CD[[x]]) || is.logical(CD[[x]])) {
		box(
			title = x,
			status = "primary",
			DTOutput(qq("overview_@{x}_table")),
			width = 3,
			height = 600
		)
	}
})


ui = div(
	h3("Overview of the MASTER"),
	do.call(fluidRow, box_list)
)

server = function(input, output, session) {
	
	for(x in CD_SELECTED) {
		if(is.numeric(CD[[x]])) {
			local({
				x = x
				output[[qq("overview_@{x}_plot")]] = renderPlot({
					par(mar = c(4, 4, 1, 1))
					plot(sort(CD[[x]]), pch = 16, cex = 0.5, xlab = "Samples", ylab = x)
				})
			})
		} else {
			local({
				x = x
				output[[qq("overview_@{x}_table")]] = renderDT({
					tb = table(CD[[x]])
					tb = sort(tb, decreasing = TRUE)
					datatable(as.data.frame(tb), options = list(dom = 'tipr'), rownames = FALSE)
				})
			})
		}
	}
}

COMPONENT$overview = list(ui = ui, server = server)


##### experiment

ui = div(
	plotOutput("experiment_overview")
)

server = function(input, output, session) {
	output$experiment_overview = renderPlot({
		upsetSamples(DB)
	})
}

COMPONENT$experiment = list(ui = ui, server = server)


###### generate snv, indel, sv, cnv, fusion

experiment_pane_generator = function(experiment_id, experiment_name) {

	has_gene_annotation = !is.null(DB[[experiment_id]]@assays[[1]]$Gene)

	ui = div(
		h3(experiment_name),
		fluidRow(
			box(title = "Configuration",
				width = 10,
				status = "primary",
				sliderInput(qq("@{experiment_id}_config_nsample"), label = qq("Number of @{experiment_id}"), min = 1, max = max(elementNROWS(DB[[experiment_id]]@assays)), value = c(1, max(elementNROWS(DB[[experiment_id]]@assays)))),
				textAreaInput(qq("@{experiment_id}_config_text"), label = "Patient/sample list"),
				textOutput(qq("@{experiment_id}_config_nsample_selected")),
				actionButton(qq("@{experiment_id}_config_submit"), label = "Submit")
			),
			box(title = qq("Number of @{experiment_name}"),
				width = 4,
				status = "primary",
				plotOutput(qq("@{experiment_id}_overview"))
			),
			box(title = "Rainfall plot",
				width = 8,
				status = "primary",
				plotOutput(qq("@{experiment_id}_rainfall"))
			),
			if(has_gene_annotation) box(title = "Top 20 genes",
				width = 4,
				status = "primary",
				tableOutput(qq("@{experiment_id}_top_genes")),
				actionButton(qq("@{experiment_id}_top_genes_oncoprint"), label = "Make oncoprint")
			) else NULL,
			box(title = "Top 20 samples",
				width = 4,
				status = "primary",
				tableOutput(qq("@{experiment_id}_top_samples"))
			),
			if(has_gene_annotation) box(width = 12,
				status = "primary",
				plotOutput(qq("@{experiment_id}_oncoprint"))
			) else NULL
		)
	)

	server = function(input, output, session) {

		n = elementNROWS(DB[[experiment_id]]@assays)

		observeEvent(input[[qq("@{experiment_id}_config_nsample")]], {
				
			config_nsample = input[[qq("@{experiment_id}_config_nsample")]]

			l = n >= config_nsample[1] & n <= config_nsample[2]
			s = names(DB[[experiment_id]]@assays)[l]

			updateTextAreaInput(session, qq("@{experiment_id}_config_text"), value = paste(s, collapse = "\n"))
			output[[qq("@{experiment_id}_config_nsample_selected")]] = renderText({
				paste0(sum(l), " samples selected")
			})

		})

		output[[qq("@{experiment_id}_overview")]] = renderPlot({
			showNotification(qq("make overview plot for @{experiment_id}."), duration = 4, type = "message")

			par(mar = c(4, 4, 1, 1))
			plot(sort(n), pch = 16, cex = 0.5, col = "red", log = "y", ylab = qq("#@{experiment_name}"), xlab = "Samples")
		})

		observeEvent(input[[qq("@{experiment_id}_config_submit")]], {
				
			samples = parse_samples(input[[qq("@{experiment_id}_config_text")]])

			output[[qq("@{experiment_id}_overview")]] = renderPlot({
				showNotification(qq("make overview plot for @{experiment_id}."), duration = 4, type = "message")
				
				n = elementNROWS(DB[[experiment_id]]@assays)
				names(n) = names(DB[[experiment_id]]@assays)
				n = sort(n)
				par(mar = c(4, 4, 1, 1))
				plot(n, pch = 16, cex = 0.5, log = "y", ylab = "# Indel", xlab = "Samples")
				
				ind = names(n) %in% samples
				points(which(ind), n[ind], pch = 16, cex = 1, col = "red")
				legend("topleft", pch = 16, legend = c("All samples", "Selected samples"), col = c("black", "red"))
			})

			output[[qq("@{experiment_id}_rainfall")]] = renderPlot({
				showNotification(qq("make rainfall plot for @{experiment_id}."), duration = 4, type = "message")

				gr = unlist(DB[[experiment_id]]@assays[samples])
				make_rainfall(gr)
			})

		})

		output[[qq("@{experiment_id}_rainfall")]] = renderPlot({
			showNotification(qq("make rainfall plot for @{experiment_id}."), duration = 4, type = "message")

			gr = unlist(DB[[experiment_id]]@assays)
			make_rainfall(gr)
		})

		if(has_gene_annotation) {
			output[[qq("@{experiment_id}_top_genes")]] = renderTable({
				x = DB[[experiment_id]]
				pa = x@assays@partitioning
				samples = rep(names(pa), times = width(pa))
				genes = rowData(x)$Gene

				tb = table(unlist(tapply(genes, samples, unique)))
				n_sample = dim(x)[2]
				tb2 = tb/n_sample

				ind = order(-tb)[1:20]
				tb2 = tb2[ind]
				df = data.frame(gene = names(tb2), fraction = paste0(round(as.vector(tb2)*100, 1), "%"), frequency = as.vector(tb)[ind])
				df[, 1] = qq("<a href=\"#\" onclick=\"link_to_other_tab($(this), 'gene', {gene_search:'@{df[,1]}'}, 'gene_search_submit');false\">@{df[,1]}</a>", collapse = FALSE)
				df
			}, sanitize.text.function = function(x) x)
		}

		output[[qq("@{experiment_id}_top_samples")]] = renderTable({
			x = DB[[experiment_id]]
			pa = x@assays@partitioning
			samples = rep(names(pa), times = width(pa))
			
			tb = table(samples)
			n_sample = dim(x)[2]
			tb2 = tb/n_sample

			ind = order(-tb)[1:20]
			tb2 = tb2[ind]
			df = data.frame(sample = names(tb2), fraction = paste0(round(as.vector(tb2)*100, 1), "%"), frequency = as.vector(tb)[ind])
			df
		}, sanitize.text.function = function(x) x)
		

		if(has_gene_annotation) {
			observeEvent(input[[qq("@{experiment_id}_top_genes_oncoprint")]], {
				showNotification(qq("make oncoPrint top 20 genes in @{experiment_id}."), duration = 4, type = "message")

				output[[qq("@{experiment_id}_oncoprint")]] = renderPlot({

					x = DB[[experiment_id]]
					pa = x@assays@partitioning
					samples = rep(names(pa), times = width(pa))
					genes = rowData(x)$Gene

					tb = table(unlist(tapply(genes, samples, unique)))
					
					ind = order(-tb)[1:20]
					tb = tb[ind]

					genes = names(tb)
					samples = rownames(colData(x))
					make_oncoprint(genes, samples, experiment_id)
				})
			})
		}

	}

	list(ui = ui, server = server)
}

# COMPONENT$snv = experiment_pane_generator("snv", "SNV")
# COMPONENT$snv_germline = experiment_pane_generator("snv_germline", "SNV Germline")
# COMPONENT$indel = experiment_pane_generator("indel", "Indel")
# COMPONENT$indel_germline = experiment_pane_generator("indel_germline", "Indel Germline")
COMPONENT$sv = experiment_pane_generator("sv", "SV")
COMPONENT$sv_germline = experiment_pane_generator("sv_germline", "SV Germline")
COMPONENT$cnv = experiment_pane_generator("cnv", "CNV")
COMPONENT$cnv_germline = experiment_pane_generator("cnv_germline", "CNV Germline")
COMPONENT$fusion = experiment_pane_generator("fusion", "Gene Fusion")


##### rna

ui = div(
	h3("rna"),
	plotOutput("rna_distribution")
)

server = function(input, output, session) {
	output[["rna_distribution"]] = renderPlot({
		mat = assay(DB[["rna"]], "TPM")
		mat = mat[sample(nrow(mat), 1000), sample(1:ncol(mat), 500)]
		densityHeatmap(log2(mat + 1), ylab = "log2(TPM + 1)", show_column_names = FALSE,
			column_order = order(colMeans(mat)))
	})
}

COMPONENT$rna = list(ui = ui, server = server)




####### mutcat_SBS

ui = div(
	plotOutput("mutcat_SBS_plot", height = 800)
)

server = function(input, output, session) {
	output[["mutcat_SBS_plot"]] = renderPlot({
		mat = DB[["mutcat_SBS"]]

		mat = mat[, colSums(mat) > 100]
		mat = mat/rep(colSums(mat), each = nrow(mat))

		ht = Heatmap(mat, show_column_names = FALSE)
		draw(ht)
	})
}

COMPONENT$mutcat_SBS = list(ui = ui, server = server)



####### mutcat_ID

ui = div(
	plotOutput("mutcat_ID_plot", height = 800)
)

server = function(input, output, session) {
	output[["mutcat_ID_plot"]] = renderPlot({
		mat = DB[["mutcat_ID"]]

		mat = mat[, colSums(mat) > 100]
		mat = mat/rep(colSums(mat), each = nrow(mat))

		ht = Heatmap(mat, show_column_names = FALSE)
		draw(ht)
	})
}

COMPONENT$mutcat_ID = list(ui = ui, server = server)



#### gene

ui = div(
	h3("gene"),
	fluidRow(
		column(
			textInput("gene_search", "Gene", value = "TP53"),
			textAreaInput("gene_sample", "Patients or samples (optional)"),
			actionButton("gene_search_submit", "Submit"),
			width = 3
		),
		column(
			htmlOutput("gene_info"),
			width = 6
		)
	),
	igvShinyOutput("gene_igv"),
	htmlOutput("gene_igv_control")
)

server = function(input, output, session) {
	observeEvent(input$gene_search_submit, {

		showNotification("Loading IGV JavaScript browser", duration = 4, type = "message")

		updateTabItems(session, "sidebar_menu", "gene")

		gene_search = input$gene_search
		updateTextInput(session, "gene_search", value = gene_search)

		output[["gene_igv"]] = renderIgvShiny({

			options = parseAndValidateGenomeSpec(genomeName="hg19",  initialLocus=gene_search,
	                                      stockGenome=TRUE, dataMode="stock",
	                                      fasta=NA, fastaIndex=NA, genomeAnnotation=NA)
			igvShiny(options)
			
		})

		output[["gene_igv_control"]] = renderUI({
			experiments = c("snv", "snv_germline", "indel", "indel_germline", "sv", "sv_germline", "fusion")
			names(experiments) = experiments
			div(
				selectInput("gene_igv_track", label = "New track", choices = experiments),
				actionButton("gene_add_track_submit", "Add new track")
			)
		})

		output[["gene_info"]] = renderUI({
			info = list("symbol" = gene_search)
			gene_id = info[["EntrezID"]] = org.Hs.egSYMBOL2EG[["TP53"]]
			info[["EnsemblID"]] = org.Hs.egENSEMBL[[gene_id]]
			info[["Gene name"]] = org.Hs.egGENENAME[[gene_id]]
			
			tb = data.frame(field = names(info), value = sapply(info, paste, collapse = ", "))
			rownames(tb) = NULL
			colnames(tb) = NULL
			HTML(knitr::kable(tb, format = "html"))
		})

	})

	observeEvent(input$gene_add_track_submit, {
		e = input$gene_igv_track

		samples = NULL
		if(!grepl("^\\s*$", input$gene_sample)) {
			s = strsplit(input$gene_sample, "\\s+")[[1]]
			samples = parse_samples(s)
		}

		if(is.null(samples)) {
			gr = unlist(DB[[e]]@assays)
		} else {
			gr = unlist(DB[[e]]@assays[samples])
		}
		gr = gr[gr$Gene == input$gene_search]
		names(gr) = NULL
		tbl = as.data.frame(gr)[, 1:3]
		colnames(tbl) = c("chr", "start", "end")
		tbl[, 1] = as.character(tbl[, 1])
		loadBedTrack(session, id = "gene_igv", trackName = e, tbl = tbl)		
	})
}

COMPONENT$gene = list(ui = ui, server = server)


#### patient

ui = div(
	h3("Patient/sample"),
	textInput("patient_search", "Patient or sample", value = "59G42F"),
	actionButton("patient_search_submit", "Submit"),
	fluidRow(
		column(width = 3,
			tableOutput("patient_search_table")
		),
		column(width = 8,
			htmlOutput("patient_rainfall")
		)
	)
)

server = function(input, output, session) {
	observeEvent(input$patient_search_submit, {

		l = META$PatientID %in% input$patient_search
		if(!any(l)) {
			l = rownames(META) %in% input$patient_search
		}

		ind = which(l)

		output$patient_search_table = renderTable({

			df = as.data.frame(META[ind, ])
			df = cbind(SampleID = rownames(META)[ind], df)

			df = as.matrix(df)
			t(df)
		}, rownames = TRUE, colnames = FALSE)

		s = rownames(META)[ind][1]
		ee = c("snv", "snv_germline", "indel", "indel_germline", "sv", "sv_germline", "cnv", "cnv_germline", "fusion")
		gl = lapply(ee, function(x) {
			DB[[x]]@assays[[s]]
		})
		names(gl) = ee

		l = !sapply(gl, is.null)
		gl = gl[l]

		output$patient_rainfall = renderUI({
			div(
				selectInput("patient_rainfall_experiment", "Select an experiment", structure(names(gl), names = names(gl))),
				textOutput("patient_rainfall_n"),
				plotOutput("patient_rainfall_plot"),
				tags$script(HTML(qq("Shiny.setInputValue('patient_sample', '@{s}')")))
			)
		})
	})

	output$patient_rainfall_n = renderText({
		experiment_id = input$patient_rainfall_experiment
		s = input$patient_sample
		gr = DB[[experiment_id]]@assays[[s]]

		qq("@{length(gr)} @{experiment_id}.")
	})

	output$patient_rainfall_plot = renderPlot({
		experiment_id = input$patient_rainfall_experiment
		s = input$patient_sample
		showNotification(qq("make rainfall plot for @{experiment_id}, sample @{s}."), duration = 4, type = "message")

		gr = DB[[experiment_id]]@assays[[s]]
		make_rainfall(gr)
	})
}

COMPONENT$patient = list(ui = ui, server = server)

#### cohort

experiments = c("indel", "indel_germline", "snv", "snv_germline", "fusion", "rna")
names(experiments) = experiments
ui = div(
	h3("OncoPrint"),
	textAreaInput("cohort_genes", label = "Genes of interest", value = "TP53\nAPC\nRB1\n"),
	textAreaInput("cohort_patients", label = "Patients of interest", value = "BD26XP\nPP8L7S\n8NFHC7\nRX17A3\nFY5ZQW\nBA95UQ\nFKPP\n"),
	radioButtons("cohort_experiments_category", label = "Experiments", choices = c("Mutations" = "mutation", "CNVs" = "cnv"), selected = "mutation", inline = TRUE),
	htmlOutput("cohort_experiments_select"),
	actionButton("cohort_submit", "Submit"),
	originalHeatmapOutput("cohort_oncoprint", response = c("click", "brush-output"), width = 400),
    HeatmapInfoOutput("cohort_oncoprint", output_ui = htmlOutput("cohort_oncoprint_output"), width = "90%")
)

server = function(input, output, session) {
	
	output$cohort_experiments_select = renderUI({
		if(input$cohort_experiments_category == "mutation") {
			selectInput("cohort_experiments", label = "", choices = experiments, selected = setdiff(experiments, "rna"), multiple = TRUE)
		} else {
			cnv_all = c("DUP", "AMP", "DEL", "HDEL", "LOH")
			cnv_all = structure(cnv_all, names = cnv_all)
			selectInput("cohort_experiments", label = "", choices = cnv_all, selected = cnv_all, multiple = TRUE)
		}
	})

	observeEvent(input$cohort_submit, {

		output$cohort_oncoprint_output = renderUI({
			HTML("")
		})

		experiments = input$cohort_experiments

		genes = strsplit(input$cohort_genes, "\\s")[[1]]
		patients = strsplit(input$cohort_patients, "\\s")[[1]]

		samples = rownames(META)[META$PatientID %in% patients]

		ht = make_oncoprint(genes, samples, experiments)

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
				output$cohort_oncoprint_output = renderUI({
					if(value != "") {
						experiments = strsplit(value, ";")[[1]]
						grl = lapply(experiments, function(e) {
							gr = DB[[e]]@assays[[sample]]
							gr[gr$Gene == gene]
						})
						names(grl) = experiments
						grl = lapply(grl, as.data.frame)

						tbl = list()
						for(i in seq_along(experiments)) {
							tb = grl[[i]]
							tb = t(as.matrix(tb))
							tbl[[i]] = box(
								title = experiments[i],
								width = 6,
								HTML(knitr::kable(tb, format = "html"))
							)
						}

						do.call(fluidRow, tbl)
					} else {
						p("No alteration.")
					}
				})	
			}
		}

		makeInteractiveComplexHeatmap(input, output, session, ht, "cohort_oncoprint",
			click_action = click_action)
	})
}

COMPONENT$cohort = list(ui = ui, server = server)



##### signature

ui = div(
	p("signature")
)

server = function(input, output, session) {

}

COMPONENT$signature = list(ui = ui, server = server)



