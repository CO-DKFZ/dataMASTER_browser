
function link_to_other_tab(self, target_tab, param, submit) {
	param_name = Object.keys(param)
	for(var i = 0; i < param_name.length; i ++) {
		Shiny.setInputValue(param_name[i], param[[ param_name[i] ]]);
	}
	Shiny.setInputValue(submit, Math.random());
	
	false;
}

