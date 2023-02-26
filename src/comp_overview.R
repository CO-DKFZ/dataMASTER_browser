
### overview
CD = colData(DB)
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

COMPONENT$overview = list()

COMPONENT$overview$ui = div(
	h3("Overview"),
	div(
		h4("The dataMASTER object", style = "padding-top:20px;"),
		div(
			verbatimTextOutput("master_obj"),
			style = "width:600px"
		),
		h4("Overview of meta columns", style = "padding-top:20px;"),
		do.call(fluidRow, box_list)
	)
)

COMPONENT$overview$server = local({

	CD = CD
	function(input, output, session) {

	output$master_obj = renderPrint({
		output = capture.output(print(DB))
		i = grep("^Functionality", output)
		cat(output[seq(1, i-1)], sep = "\n")
	})
	
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
})
