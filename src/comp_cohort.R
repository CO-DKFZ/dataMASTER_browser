
COMPONENT$cohort = list()

COMPONENT$cohort$ui = div(
	h3("OncoPrint"),
	textAreaInput("cohort_genes", label = "Genes of interest", value = "TP53\nAPC\nRB1\n"),
	textAreaInput("cohort_patients", label = "Patients of interest", value = "BD26XP\nPP8L7S\n8NFHC7\nRX17A3\nFY5ZQW\nBA95UQ\nFKPP\n"),
	radioButtons("cohort_experiments_category", label = "Experiments", choices = c("Mutations" = "mutation", "CNVs" = "cnv"), selected = "mutation", inline = TRUE),
	htmlOutput("cohort_experiments_select"),
	actionButton("cohort_submit", "Submit"),
	originalHeatmapOutput("cohort_oncoprint", response = c("click", "brush-output"), width = 400),
    HeatmapInfoOutput("cohort_oncoprint", output_ui = htmlOutput("cohort_oncoprint_output"), width = "90%")
)

COMPONENT$cohort$server = function(input, output, session) {

	experiments = c("indel", "indel_germline", "snv", "snv_germline", "fusion", "rna")
	names(experiments) = experiments

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
