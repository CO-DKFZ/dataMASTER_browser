
COMPONENT$patient = list()

COMPONENT$patient$ui = div(
	h3("Patients"),
	div(
		textInput("patient_search", "Patient/Sample", placeholder = "WES.59G42F.metastasis"),
		actionButton("patient_search_submit", "Search"),
		tags$script(HTML(qq('
	$("#patient_search").typeahead({
	    source: @{jsonlite::toJSON(rownames(META))},
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
	),
	br(),
	tabsetPanel(
		type = "tabs",
		tabPanel(qq("Meta table"),
			uiOutput("patient_search_table_ui"),
			style = "padding:16px 0px"
		),
		tabPanel(qq("analysis"),
			uiOutput("patient_analysis_ui"),
			style = "padding:16px 0px"
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

		if(length(ind) == 0) {
			output$patient_search_table_ui = renderUI({
				p("No sample is found.")
			})
			return(NULL)
		}

		output$patient_search_table_ui = renderUI({
			tableOutput("patient_search_table")
		})

		output$patient_analysis_ui = renderUI({
			tableOutput("patient_analysis")
		})

		output$patient_search_table = renderTable({

			df = as.data.frame(META[ind, ])
			df = cbind(SampleID = rownames(META)[ind], df)

			df = as.matrix(df)
			t(df)
		}, rownames = TRUE, colnames = FALSE,
		   sanitize.rownames.function = function(x) paste0("<b>", x, "</b>"))

		samples = rownames(META)[ind]
		ee = c("snv", "snv_germline", "indel", "indel_germline", "sv", "sv_germline", "cnv", "cnv_germline", "fusion")
		
		m = matrix("0", nrow = length(samples), ncol = length(ee))
		rownames(m) = samples
		colnames(m) = ee
		for(s in samples) {
			for(e in ee) {
				gr = DB[[e]]@assays[[s]]
				if(!is.null(gr)) {
					k = length(gr)
				} else {
					k = 0
				}

				if(k == 0) {
					m[s, e] = paste0("no ", e)
				} else {
					m[s, e] = qq("<a href=\"#\" onclick=\"link_to_other_tab($(this), '@{e}', {@{e}_config_sample_list:'@{s}'}, '@{e}_config_submit');false\">@{k} @{e}</a>")
				}
			}
		}

		output[["patient_analysis"]] = renderTable({
			m
		}, rownames = TRUE, colnames = FALSE, sanitize.text.function = function(x) x)

	})
}
