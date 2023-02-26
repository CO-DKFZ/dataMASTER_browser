

experiment_id = "sequencingQC"
experiment_name = "sequencingQC"

COMPONENT[[experiment_id]] = list()

mat = DB[[experiment_id]]

COMPONENT[[experiment_id]]$ui = div(
	id = qq("@{experiment_id}_workspace"),
	h3(qq("Experiment: @{experiment_name}")),
	p("Not implemented yet...")
)

COMPONENT[[experiment_id]]$server = function(input, output, session) {}
