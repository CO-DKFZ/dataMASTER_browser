
COMPONENT$patient = list()

COMPONENT$patient$ui = div(
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

COMPONENT$patient$server = function(input, output, session) {
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
