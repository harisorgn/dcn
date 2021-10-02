
struct comp_t
	x ::Float64
	y ::Float64
	z ::Float64
	d ::Float64
	ctype ::Symbol
end

struct reduced_comp_t
	l ::Float64
	d ::Float64
	Ra ::Float64
	n_comp ::Float64
	mrg_n_comp ::Float64
	area ::Float64
	full_area_sum ::Float64
	id ::String
	connect_to ::String
	ctype ::Symbol
end

function resistance_to_soma(comp::String, connect_dict::Dict{String, String}, comp_dict::Dict{String, comp_t})

	R = 0.0 ;
	current_comp = comp

	while current_comp != "soma"

		l_current_comp = sqrt((comp_dict[current_comp].x - comp_dict[connect_dict[current_comp]].x)^2 +
					  (comp_dict[current_comp].y - comp_dict[connect_dict[current_comp]].y)^2 +
					  (comp_dict[current_comp].z - comp_dict[connect_dict[current_comp]].z)^2 ) ;

		R += 4.0 * RA * l_current_comp / (pi * comp_dict[current_comp].d^2) ;

		current_comp = connect_dict[current_comp] ;
	end

	return R 
end

function increment_strahler_chk(v::Array{Int64})

	m = maximum(v) ;
	m_idx = find(x -> x == m, v) ;

	if length(m_idx) > 1
		return true
	else
		return false
	end

end

function find_value_dict(v :: T, d :: Dict{String, T}) where T 

	idx = find(x -> x == v, values(d)) ;
	return collect(keys(d))[idx] 

end

function rec_strahler(connect_dict::Dict{String, String}, n_connect_dict::Dict{String, Int64}, 
					  strahler_dict::Dict{String, Int64}, comp::String)

	if n_connect_dict[comp] == 0 

		return 1 

	else

		connect_idx = find(x -> x == comp, values(connect_dict)) ;
		connect_v = collect(keys(connect_dict))[connect_idx] ;
		strahler_connect_v = Array{Int64}(length(connect_v)) ;

		for i = 1 : length(connect_v)
			strahler_connect_v[i] = rec_strahler(connect_dict, n_connect_dict, strahler_dict, connect_v[i]) ;
			strahler_dict[connect_v[i]] = strahler_connect_v[i] ;
		end

		if length(connect_v) == 1

			return strahler_dict[connect_v[1]]

		elseif increment_strahler_chk(strahler_connect_v)

			return maximum(strahler_connect_v) + 1 ;

		else

			return maximum(strahler_connect_v) ;

		end
	end

end

function rec_reduced_comp(connect_dict::Dict{String, String}, comp_dict::Dict{String, comp_t},
						sect_connect_dict::Dict{String, String}, other_type_cluster_root_v::Array{String,1},
						strahler_dict ::Dict{String, Int64}, strahler_thrs ::Int64, 
						reduced_comp_v ::Array{reduced_comp_t,1}, comp::String, ctype::Symbol)
		
	l_seq = 0.0 ;
	d_seq = 0.0 ;
	Ra_seq = 0.0 ;
	n_comp_seq = 0.0 ;
	full_area_sum = 0.0 ;

	if ctype != :all 
		if strahler_dict[comp] <= strahler_thrs && comp_dict[comp].ctype == ctype
			
			current_comp = comp ;
			while current_comp != connect_dict[sect_connect_dict[comp]]
				
				l_current_comp = sqrt((comp_dict[current_comp].x - comp_dict[connect_dict[current_comp]].x)^2 +
						  (comp_dict[current_comp].y - comp_dict[connect_dict[current_comp]].y)^2 +
						  (comp_dict[current_comp].z - comp_dict[connect_dict[current_comp]].z)^2 ) ;

				l_seq += l_current_comp ;
				Ra_seq += 4.0 * RA * l_current_comp / (pi * comp_dict[current_comp].d^2) ;
				n_comp_seq += 1.0 ;
				full_area_sum += pi * l_current_comp * comp_dict[current_comp].d ;

				current_comp = connect_dict[current_comp] ;
			end

			d_seq = sqrt(4.0 * RA * l_seq / (pi * Ra_seq)) ;

			comp_id = comp ;
			
		elseif comp != "soma" && comp_dict[comp].ctype == ctype
			
			current_comp = comp ;

			while current_comp != connect_dict[sect_connect_dict[comp]]

				l_current_comp = sqrt((comp_dict[current_comp].x - comp_dict[connect_dict[current_comp]].x)^2 +
						  (comp_dict[current_comp].y - comp_dict[connect_dict[current_comp]].y)^2 +
						  (comp_dict[current_comp].z - comp_dict[connect_dict[current_comp]].z)^2 ) ;

				Ra_current_comp = 4.0 * RA * l_current_comp / (pi * comp_dict[current_comp].d^2) ;

				push!(reduced_comp_v, reduced_comp_t(l_current_comp, comp_dict[current_comp].d, Ra_current_comp,
												1.0, 1.0, pi * l_current_comp * comp_dict[current_comp].d,
												pi * l_current_comp * comp_dict[current_comp].d,
												current_comp, connect_dict[current_comp], ctype)) ;

				current_comp = connect_dict[current_comp] ;
			end
			
		end
	else
		if strahler_dict[comp] <= strahler_thrs 
			
			current_comp = comp ;
			while current_comp != connect_dict[sect_connect_dict[comp]]
				
				l_current_comp = sqrt((comp_dict[current_comp].x - comp_dict[connect_dict[current_comp]].x)^2 +
						  (comp_dict[current_comp].y - comp_dict[connect_dict[current_comp]].y)^2 +
						  (comp_dict[current_comp].z - comp_dict[connect_dict[current_comp]].z)^2 ) ;

				l_seq += l_current_comp ;
				Ra_seq += 4.0 * RA * l_current_comp / (pi * comp_dict[current_comp].d^2) ;
				n_comp_seq += 1.0 ;
				full_area_sum += pi * l_current_comp * comp_dict[current_comp].d ;

				current_comp = connect_dict[current_comp] ;
			end

			d_seq = sqrt(4.0 * RA * l_seq / (pi * Ra_seq)) ;

			comp_id = comp ;
		end
	end

	connect_comp_v = find_value_dict(comp, connect_dict) ;

	n_connect_comp_below_thrs = 0 ;

	l_par = 0.0 ;
	Ra_par = 0.0 ;
	la_par = 0.0 ;
	a_par = 0.0 ;
	d_par = 0.0 ;
	Ra_par_sum = 0.0 ;
	Ra_par_prod = 1.0 ;
	n_comp_par = 0.0 ;

	for c in connect_comp_v
		if ctype != :all
			if comp_dict[c].ctype == ctype && strahler_dict[c] <= strahler_thrs
				
				n_connect_comp_below_thrs += 1 ;

				sect_start_connect_comp = find_value_dict(c, sect_connect_dict) ;
				
				reduced_connect_comp = rec_reduced_comp(connect_dict, comp_dict, sect_connect_dict,
														other_type_cluster_root_v, strahler_dict, strahler_thrs,
														reduced_comp_v, sect_start_connect_comp[1], ctype) ;

				a_par += reduced_connect_comp.area ;
				la_par += reduced_connect_comp.l * reduced_connect_comp.area ;
				d_par += reduced_connect_comp.d ^ 2.0 ;
				Ra_par_prod *= reduced_connect_comp.Ra ;
				Ra_par_sum += reduced_connect_comp.Ra ;
				
				n_comp_par += reduced_connect_comp.n_comp ;
				full_area_sum += reduced_connect_comp.full_area_sum ;
				comp_id = c ;

			elseif comp_dict[c].ctype != ctype
				if ~any(x -> x == comp, other_type_cluster_root_v) && ctype == :pdend
					push!(other_type_cluster_root_v, comp) ;
				end
			end
		else

			if strahler_dict[c] <= strahler_thrs
			
				n_connect_comp_below_thrs += 1 ;

				sect_start_connect_comp = find_value_dict(c, sect_connect_dict) ;
				
				reduced_connect_comp = rec_reduced_comp(connect_dict, comp_dict, sect_connect_dict,
														other_type_cluster_root_v, strahler_dict, strahler_thrs,
														reduced_comp_v, sect_start_connect_comp[1], ctype) ;

				a_par += reduced_connect_comp.area ;
				la_par += reduced_connect_comp.l * reduced_connect_comp.area ;
				d_par += reduced_connect_comp.d ^ 2.0 ;
				Ra_par_prod *= reduced_connect_comp.Ra ;
				Ra_par_sum += reduced_connect_comp.Ra ;
				
				n_comp_par += reduced_connect_comp.n_comp ;
				full_area_sum += reduced_connect_comp.full_area_sum ;
				comp_id = c ;
			end
		end
	end

	if n_connect_comp_below_thrs > 1
		d_par = sqrt(d_par) ;
		l_par = la_par / a_par ;
		Ra_par = Ra_par_prod / Ra_par_sum ;
	elseif n_connect_comp_below_thrs == 1
		d_par = sqrt(d_par) ;
		l_par = la_par / a_par ;
		Ra_par = Ra_par_prod ;
	end

	if ctype != :all 
		if strahler_dict[comp] <= strahler_thrs && comp_dict[comp].ctype == ctype 
			l_eq = l_seq + l_par ;
			Ra_eq = Ra_seq + Ra_par ;
			d_eq = sqrt(4.0 * RA * l_eq / (pi * Ra_eq)) ;
			n_comp = n_comp_seq + n_comp_par ;
			mrg_n_comp = 1.0 ;

			connect_to = connect_dict[sect_connect_dict[comp]] ;

			return reduced_comp_t(l_eq, d_eq, Ra_eq, n_comp, mrg_n_comp, pi*d_eq*l_eq, full_area_sum, comp_id, connect_to, ctype)

		elseif (strahler_dict[comp] > strahler_thrs || comp_dict[comp].ctype != ctype) && 
				n_connect_comp_below_thrs > 0
			Ra_eq = Ra_par_sum ;
			d_eq = d_par ; 
			n_comp = n_comp_seq + n_comp_par ;
			mrg_n_comp = n_connect_comp_below_thrs ;
			if ctype == :pdend
				l_eq = l_par ;
				reduced_area = pi * l_eq * d_eq ; 
			elseif ctype == :ddend
				l_eq = la_par ;
				reduced_area = a_par ;
			end

			connect_to = comp ;

			return reduced_comp_t(l_eq, d_eq, Ra_eq, n_comp, mrg_n_comp, a_par, full_area_sum, comp_id, connect_to, ctype)
		end
	else
		if strahler_dict[comp] <= strahler_thrs 
			l_eq = l_seq + l_par ;
			Ra_eq = Ra_seq + Ra_par ;
			d_eq = sqrt(4.0 * RA * l_eq / (pi * Ra_eq)) ;
			n_comp = n_comp_seq + n_comp_par ;
			mrg_n_comp = 1.0 ;

			connect_to = connect_dict[sect_connect_dict[comp]] ;

			return reduced_comp_t(l_eq, d_eq, Ra_eq, n_comp, mrg_n_comp, pi*d_eq*l_eq, full_area_sum, comp_id, connect_to, ctype)

		elseif strahler_dict[comp] > strahler_thrs  && n_connect_comp_below_thrs > 0	
			Ra_eq = Ra_par_sum / n_connect_comp_below_thrs ;
			d_eq = d_par ; 
			n_comp = n_comp_seq + n_comp_par ;
			mrg_n_comp = n_connect_comp_below_thrs ;
			l_eq = l_par ;

			connect_to = comp ;

			return reduced_comp_t(l_eq, d_eq, Ra_eq, n_comp, mrg_n_comp, a_par, full_area_sum, comp_id, connect_to, ctype)
		end

	end
end

function strahler(strahler_thrs ::Int64)


	reduced_comp_v = Array{reduced_comp_t, 1}() ;

	file = open("./data/pDend.txt", "r")
	pdend = readdlm(file) ;
	close(file)

	file = open("./data/dDend.txt", "r")
	ddend = readdlm(file) ;
	close(file)

	file = open("./data/axIS.txt", "r")
	axis = readdlm(file) ;
	close(file)

	file = open("./data/axIN.txt", "r")
	axin = readdlm(file) ;
	close(file)

	strahler_dict = Dict{String, Int64}() ;
	sizehint!(strahler_dict, size(pdend,1) + size(ddend,1) + 1) ;

	comp_dict = Dict{String, comp_t}() ;
	sizehint!(comp_dict, size(pdend,1) + size(ddend,1) + 1) ;
	comp_dict["soma"] = comp_t(0.0, 0.0, 0.0, 21.597e-6, :s) ;

	connect_dict = Dict{String, String}() ;
	sizehint!(connect_dict, size(pdend,1) + size(ddend,1)) ;

	n_connect_dict = Dict{String, Int64}() ;
	sizehint!(n_connect_dict, size(pdend,1) + size(ddend,1) + 1) ;

	push!(reduced_comp_v, reduced_comp_t(0.0, comp_dict["soma"].d, 8.0 * RA / (pi * comp_dict["soma"].d),
										1.0, 1.0, pi * comp_dict["soma"].d^2.0, pi * comp_dict["soma"].d^2.0,
										"soma", "", :s)) ;
	
	l_axis = 0.0 ;
	Ra_axis = 0.0 ;
	n_comp = 0.0 ;
	full_area_sum = 0.0 ;
	for i = 2 : size(axis,1)

		current_l_axis = sqrt((axis[i,1] - axis[i - 1,1])^2.0 + 
						 (axis[i,2] - axis[i - 1,2])^2.0 +
						 (axis[i,3] - axis[i - 1,3])^2.0) ;

		l_axis += current_l_axis * 1e-6 ;
		Ra_axis += 4.0 * RA * current_l_axis * 1e-6 / (pi * (axis[i,4] * 1e-6)^2.0) ;
		n_comp += 1.0 ;
		full_area_sum += pi * current_l_axis * 1e-6 * axis[i,4] * 1e-6 ; 
	end

	d_axis = sqrt(4.0 * RA * l_axis / (pi * Ra_axis)) ;

	push!(reduced_comp_v, reduced_comp_t(l_axis, d_axis, Ra_axis, n_comp, 1.0, pi*l_axis*d_axis, 
										full_area_sum, "axis", "axhill", :axis)) ;

	l_axin = 0.0 ;
	Ra_axin = 0.0 ;
	n_comp = 0.0 ;
	full_area_sum = 0.0 ;
	for i = 2 : size(axin,1)

		current_l_axin = sqrt((axin[i,1] - axin[i - 1,1])^2.0 + 
						 (axin[i,2] - axin[i - 1,2])^2.0 +
						 (axin[i,3] - axin[i - 1,3])^2.0) ;

		l_axin += current_l_axin * 1e-6 ;
		Ra_axin += 4.0 * RA * current_l_axin * 1e-6 / (pi * (axin[i,4] * 1e-6)^2.0) ;
		n_comp += 1.0 ;
		full_area_sum += pi * current_l_axin * 1e-6 * axin[i,4] * 1e-6 ; 
	end

	d_axin = sqrt(4.0 * RA * l_axin / (pi * Ra_axin)) ;

	push!(reduced_comp_v, reduced_comp_t(l_axin, d_axin, Ra_axin, n_comp, 1.0, pi*l_axin*d_axin, 
										full_area_sum, "axin", "axis", :axin)) ;

	push!(reduced_comp_v, reduced_comp_t(5.0e-6, 4.75e-6, 4.0 * RA * 5.0e-6 / (pi * (4.75e-6)^2.0), 1.0, 1.0,
										pi * 5.0e-6 * 4.75e-6, pi * 5.0e-6 * 4.75e-6,
										"axhill", "soma", :axhill)) ;
	
	for i = 1 : size(pdend, 1)

		comp_dict[pdend[i,1]] = comp_t(pdend[i,3]*1e-6, pdend[i,4]*1e-6, pdend[i,5]*1e-6, pdend[i,6]*1e-6, 
										:pdend) ;

		connect_dict[pdend[i,1]] = pdend[i,2] ;

		get!(n_connect_dict, pdend[i,1], 0);
		n_connect_dict[pdend[i,2]] = get(n_connect_dict, pdend[i,2], 0) + 1 ;

	end

	for i = 1 : size(ddend, 1)

		comp_dict[ddend[i,1]] = comp_t(ddend[i,3]*1e-6, ddend[i,4]*1e-6, ddend[i,5]*1e-6, ddend[i,6]*1e-6, 
										:ddend) ;

		connect_dict[ddend[i,1]] = ddend[i,2] ;

		get!(n_connect_dict, ddend[i,1], 0);
		n_connect_dict[ddend[i,2]] = get(n_connect_dict, ddend[i,2], 0) + 1 ;

	end

	strahler_dict["soma"] = rec_strahler(connect_dict, n_connect_dict, strahler_dict, "soma") ;

	sect_connect_dict = Dict{String, String}() ;
	cluster_root_v = Array{String,1}() ;

	open_v = find_value_dict(0, n_connect_dict) ;

	closed_v = Array{String,1}() ;
	push!(closed_v, "soma") ;
	
	while ~isempty(open_v)

		current_comp = open_v[1] ;
		ctype = comp_dict[open_v[1]].ctype ;
		
		while n_connect_dict[connect_dict[current_comp]] == 1 && 
			comp_dict[connect_dict[current_comp]].ctype == ctype

				current_comp = connect_dict[current_comp] ;
		end

		sect_connect_dict[open_v[1]] = current_comp ;

		if strahler_dict[connect_dict[current_comp]] > strahler_thrs &&
			~any(x -> x == connect_dict[current_comp], cluster_root_v) 

			push!(cluster_root_v, connect_dict[current_comp])
		end

		if ~any(x -> x == connect_dict[current_comp], closed_v) &&
		   ~any(x -> x == connect_dict[current_comp], open_v)			   
				push!(open_v, connect_dict[current_comp]) ;
		end

		push!(closed_v, open_v[1]) ;
		shift!(open_v) ;

	end
	
	for comp in cluster_root_v
		
		other_type_cluster_root_v = Array{String,1}() ; 
		joint_reduced_comp = rec_reduced_comp(connect_dict, comp_dict, sect_connect_dict,
												other_type_cluster_root_v, strahler_dict, strahler_thrs,
												reduced_comp_v, comp, :all) ;
		
		pdend_reduced_comp = rec_reduced_comp(connect_dict, comp_dict, sect_connect_dict,
												other_type_cluster_root_v, strahler_dict, strahler_thrs,
												reduced_comp_v, comp, :pdend) ;
		
		la_sum = 0.0 ;
		area_sum = 0.0 ;
		d_eq = 0.0 ;
		Ra_sum = 0.0 ;
		n_comp = 0.0 ;
		mrg_n_comp = 0.0 ;
		comp_id = " ";
		full_area_sum = 0.0 ;
		for other_type_comp in other_type_cluster_root_v

			chk_cluster_root_v = Array{String,1}() ;

			ddend_part_reduced_comp = rec_reduced_comp(connect_dict, comp_dict, sect_connect_dict,
														chk_cluster_root_v, strahler_dict, strahler_thrs,
														reduced_comp_v, other_type_comp, :ddend) ;
			if ddend_part_reduced_comp != nothing

				area_sum += ddend_part_reduced_comp.area ;
				la_sum += ddend_part_reduced_comp.l ;
				d_eq += ddend_part_reduced_comp.d ^ 2.0 ;
				n_comp += ddend_part_reduced_comp.n_comp ;
				
				comp_id = ddend_part_reduced_comp.id ;
				full_area_sum += ddend_part_reduced_comp.full_area_sum ;

				Ra_sum += ddend_part_reduced_comp.Ra ;
				mrg_n_comp += ddend_part_reduced_comp.mrg_n_comp ;
				
			end
		end
		l_eq = la_sum / area_sum ;
		d_eq = sqrt(d_eq) ;
		Ra_eq = Ra_sum ;
		reduced_area = area_sum ;

		if ~isempty(other_type_cluster_root_v)
			if pdend_reduced_comp != nothing
				ddend_reduced_comp = reduced_comp_t(l_eq, d_eq, Ra_eq, n_comp, mrg_n_comp, reduced_area, 
														full_area_sum, comp_id, reduced_comp_v[end].id, :ddend) ;
			else
				idx = find(x -> x.id == comp, reduced_comp_v) ;
				ddend_reduced_comp = reduced_comp_t(l_eq, d_eq, Ra_eq, n_comp, mrg_n_comp, reduced_area, 
														full_area_sum, comp_id, reduced_comp_v[idx[1]].id, :ddend) ;
			end
		end
		
		
		ddend_ratio = ddend_reduced_comp.l / joint_reduced_comp.l ;
		ddend_l = ddend_ratio * joint_reduced_comp.l ;
		ddend_Ra = ddend_ratio * joint_reduced_comp.Ra ;

		pdend_l = joint_reduced_comp.l - ddend_l ;
		pdend_Ra = joint_reduced_comp.Ra - ddend_Ra ;

		push!(reduced_comp_v, reduced_comp_t(pdend_reduced_comp.l, joint_reduced_comp.d,
											pdend_reduced_comp.Ra/pdend_reduced_comp.mrg_n_comp, 
											pdend_reduced_comp.n_comp,
											1.0, pdend_reduced_comp.area , 
											0.5*pdend_reduced_comp.full_area_sum, pdend_reduced_comp.id, 
											pdend_reduced_comp.connect_to, :pdend)) ;

		push!(reduced_comp_v, reduced_comp_t(joint_reduced_comp.l - pdend_reduced_comp.l, joint_reduced_comp.d,
											joint_reduced_comp.Ra - (pdend_reduced_comp.Ra/pdend_reduced_comp.mrg_n_comp), 
											ddend_reduced_comp.n_comp,
											1.0, ddend_reduced_comp.area ,
											0.5*ddend_reduced_comp.full_area_sum, ddend_reduced_comp.id, 
											ddend_reduced_comp.connect_to, :ddend)) ;
	end



	nc = 0.0 ;
	for i = 1 : length(reduced_comp_v)
		nc += reduced_comp_v[i].n_comp ;
	end 
	nc -= reduced_comp_v[1].n_comp + reduced_comp_v[2].n_comp + reduced_comp_v[3].n_comp + reduced_comp_v[4].n_comp;
	println("dendritic compartments included in reduced model: ", nc)
	println("actual number of dendritic compartments: ", size(pdend,1) + size(ddend,1))
	println("reduced model compartments: ",length(reduced_comp_v))
	
	return reduced_comp_v
end

function full_model()


	reduced_comp_v = Array{reduced_comp_t, 1}() ;

	file = open("./data/pDend.txt", "r")
	pdend = readdlm(file) ;
	close(file)

	file = open("./data/dDend.txt", "r")
	ddend = readdlm(file) ;
	close(file)

	file = open("./data/axIS.txt", "r")
	axis = readdlm(file) ;
	close(file)

	file = open("./data/axIN.txt", "r")
	axin = readdlm(file) ;
	close(file)
	
	strahler_dict = Dict{String, Int64}() ;
	sizehint!(strahler_dict, size(pdend,1) + size(ddend,1) + 1) ;

	comp_dict = Dict{String, comp_t}() ;
	sizehint!(comp_dict, size(pdend,1) + size(ddend,1) + 1) ;
	comp_dict["soma"] = comp_t(0.0, 0.0, 0.0, 21.597e-6, :s) ;

	connect_dict = Dict{String, String}() ;
	sizehint!(connect_dict, size(pdend,1) + size(ddend,1)) ;

	n_connect_dict = Dict{String, Int64}() ;
	sizehint!(n_connect_dict, size(pdend,1) + size(ddend,1) + 1) ;

	push!(reduced_comp_v, reduced_comp_t(0.0, comp_dict["soma"].d, 8.0 * RA / (pi * comp_dict["soma"].d),
										1.0, 1.0, pi * comp_dict["soma"].d^2.0, pi * comp_dict["soma"].d^2.0,
										"soma", "", :s)) ;

	l_axis = 0.0 ;
	Ra_axis = 0.0 ;
	n_comp = 0.0 ;
	full_area_sum = 0.0 ;
	for i = 2 : size(axis,1)

		current_l_axis = sqrt((axis[i,1] - axis[i - 1,1])^2.0 + 
						 (axis[i,2] - axis[i - 1,2])^2.0 +
						 (axis[i,3] - axis[i - 1,3])^2.0) ;

		l_axis += current_l_axis * 1e-6 ;
		Ra_axis += 4.0 * RA * current_l_axis * 1e-6 / (pi * (axis[i,4] * 1e-6)^2.0) ;
		n_comp += 1.0 ;
		full_area_sum += pi * current_l_axis * 1e-6 * axis[i,4] * 1e-6 ; 
	end

	d_axis = sqrt(4.0 * RA * l_axis / (pi * Ra_axis)) ;

	push!(reduced_comp_v, reduced_comp_t(l_axis, d_axis, Ra_axis, n_comp, 1.0, pi*l_axis*d_axis, 
										full_area_sum, "axis", "axhill", :axis)) ;

	l_axin = 0.0 ;
	Ra_axin = 0.0 ;
	n_comp = 0.0 ;
	full_area_sum = 0.0 ;
	for i = 2 : size(axin,1)

		current_l_axin = sqrt((axin[i,1] - axin[i - 1,1])^2.0 + 
						 (axin[i,2] - axin[i - 1,2])^2.0 +
						 (axin[i,3] - axin[i - 1,3])^2.0) ;

		l_axin += current_l_axin * 1e-6 ;
		Ra_axin += 4.0 * RA * current_l_axin * 1e-6 / (pi * (axin[i,4] * 1e-6)^2.0) ;
		n_comp += 1.0 ;
		full_area_sum += pi * current_l_axin * 1e-6 * axin[i,4] * 1e-6 ; 
	end

	d_axin = sqrt(4.0 * RA * l_axin / (pi * Ra_axin)) ;

	push!(reduced_comp_v, reduced_comp_t(l_axin, d_axin, Ra_axin, n_comp, 1.0, pi*l_axin*d_axin, 
										full_area_sum, "axin", "axis", :axin)) ;

	push!(reduced_comp_v, reduced_comp_t(5.0e-6, 4.75e-6, 4.0 * RA * 5.0e-6 / (pi * (4.75e-6)^2.0), 1.0, 1.0,
										pi * 5.0e-6 * 4.75e-6, pi * 5.0e-6 * 4.75e-6,
										"axhill", "soma", :axhill)) ;

	for i = 1 : size(pdend, 1)

		comp_dict[pdend[i,1]] = comp_t(pdend[i,3]*1e-6, pdend[i,4]*1e-6, pdend[i,5]*1e-6, pdend[i,6]*1e-6, 
										:pdend) ;

		connect_dict[pdend[i,1]] = pdend[i,2] ;

		get!(n_connect_dict, pdend[i,1], 0);
		n_connect_dict[pdend[i,2]] = get(n_connect_dict, pdend[i,2], 0) + 1 ;


	end

	for i = 1 : size(ddend, 1)

		comp_dict[ddend[i,1]] = comp_t(ddend[i,3]*1e-6, ddend[i,4]*1e-6, ddend[i,5]*1e-6, ddend[i,6]*1e-6, 
										:ddend) ;

		connect_dict[ddend[i,1]] = ddend[i,2] ;

		get!(n_connect_dict, ddend[i,1], 0);
		n_connect_dict[ddend[i,2]] = get(n_connect_dict, ddend[i,2], 0) + 1 ;

	end


	for comp in keys(comp_dict)
		if comp != "soma" 
		l_current_comp = sqrt((comp_dict[comp].x - comp_dict[connect_dict[comp]].x)^2 +
					  (comp_dict[comp].y - comp_dict[connect_dict[comp]].y)^2 +
					  (comp_dict[comp].z - comp_dict[connect_dict[comp]].z)^2 ) ;

		push!(reduced_comp_v, reduced_comp_t(l_current_comp, comp_dict[comp].d, 4.0 * RA * l_current_comp / (pi * comp_dict[comp].d^2.0),
											1.0, 1.0, pi * l_current_comp * comp_dict[comp].d, 
											pi * l_current_comp * comp_dict[comp].d, comp, connect_dict[comp], comp_dict[comp].ctype)) ;
		end
	end
	
	return reduced_comp_v
end
