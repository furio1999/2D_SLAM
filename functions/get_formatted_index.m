% return the progressive index of the searched real id inside allIDs
function index_order = get_formatted_index(allIDs,searched_real_id)
	index_order = 0;
	for i=1:length(allIDs)
		if allIDs(1,i) == searched_real_id
			index_order = i;
		endif
	endfor

endfunction