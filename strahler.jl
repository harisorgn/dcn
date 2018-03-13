
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
	area_sum ::Float64
	id ::String
	connect_to ::String
	ctype ::Symbol
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

	below_thrs_comp_chk = false ;
	
	l_sect = 0.0 ;
	Ra_sect = 0.0 ;
	n_comp_sect = 0.0 ;
	area_sum = 0.0 ;

	if strahler_dict[comp] <= strahler_thrs && comp_dict[comp].ctype == ctype
		
		below_thrs_comp_chk = true ;	

		current_comp = comp ;
		while current_comp != connect_dict[sect_connect_dict[comp]]

			l_current_comp = sqrt((comp_dict[current_comp].x - comp_dict[connect_dict[current_comp]].x)^2 +
					  (comp_dict[current_comp].y - comp_dict[connect_dict[current_comp]].y)^2 +
					  (comp_dict[current_comp].z - comp_dict[connect_dict[current_comp]].z)^2 ) ;

			l_sect += l_current_comp ;
			Ra_sect += 4.0 * RA * l_current_comp / (pi * comp_dict[current_comp].d^2) ;
			n_comp_sect += 1.0 ;
			area_sum += pi * l_sect * comp_dict[current_comp].d ;

			current_comp = connect_dict[current_comp] ;
		end

		d_sect = sqrt(RA * l_sect / (pi * Ra_sect)) ;

		comp_id = comp ;
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
		if comp_dict[c].ctype == ctype && strahler_dict[c] <= strahler_thrs

			n_connect_comp_below_thrs += 1 ;

			sect_start_connect_comp = find_value_dict(c, sect_connect_dict) ;

			reduced_connect_comp = rec_reduced_comp(connect_dict, comp_dict, sect_connect_dict,
													other_type_cluster_root_v, strahler_dict, strahler_thrs,
													reduced_comp_v, sect_start_connect_comp[1], ctype) ;

			a_par += pi * reduced_connect_comp.l * reduced_connect_comp.d ;
			la_par += reduced_connect_comp.l * pi * reduced_connect_comp.l * reduced_connect_comp.d ;
			d_par += reduced_connect_comp.d ^ 2.0 ;
			Ra_par_prod *= reduced_connect_comp.Ra ;
			Ra_par_sum += reduced_connect_comp.Ra ;
			
			n_comp_par += reduced_connect_comp.n_comp ;
			area_sum += reduced_connect_comp.area_sum ;
			comp_id = c ;
		elseif comp_dict[c].ctype == ctype && strahler_dict[c] > strahler_thrs
			current_comp = find_value_dict(c, sect_connect_dict)[1] ;

			while current_comp != connect_dict[c]

				l_current_comp = sqrt((comp_dict[current_comp].x - comp_dict[connect_dict[current_comp]].x)^2 +
						  (comp_dict[current_comp].y - comp_dict[connect_dict[current_comp]].y)^2 +
						  (comp_dict[current_comp].z - comp_dict[connect_dict[current_comp]].z)^2 ) ;

				Ra_current_comp = 4.0 * RA * l_current_comp / (pi * comp_dict[current_comp].d^2) ;

				push!(reduced_comp_v, reduced_comp_t(l_current_comp, comp_dict[current_comp].d, Ra_current_comp,
												1.0, pi * l_current_comp * comp_dict[current_comp].d,
												current_comp, connect_dict[current_comp], ctype)) ;

				current_comp = connect_dict[current_comp] ;
			end
		else
			if ~any(x -> x == comp, other_type_cluster_root_v) && ctype == :pdend
				push!(other_type_cluster_root_v, comp) ;
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

	if below_thrs_comp_chk || n_connect_comp_below_thrs > 0 
		l_eq = l_sect + l_par ;
		Ra_eq = Ra_sect + Ra_par ;
		d_eq = sqrt(RA * l_eq / (pi * Ra_eq)) ;
		n_comp = n_comp_sect + n_comp_par ;

		if strahler_dict[comp] > strahler_thrs
			connect_to = comp ;
		else
			connect_to = connect_dict[sect_connect_dict[comp]] ;
		end
		
		return reduced_comp_t(l_eq, d_eq, Ra_eq, n_comp, area_sum, comp_id, connect_to, ctype)
	end

	#=
	for c in connect_comp_v
		if comp_dict[c].ctype == ctype #&& strahler_dict[c] <= strahler_thrs
			cluster_end_chk = false ;
		else
			if ~any(x -> x == comp, other_type_cluster_root_v) && ctype == :pdend
				push!(other_type_cluster_root_v, comp) ;
			end
		end
	end

	if cluster_end_chk && strahler_dict[comp] <= strahler_thrs

		l_eq = 0.0 ;
		d_eq = 0.0 ;
		Ra_eq = 0.0 ;
		n_comp = 0.0 ;

		current_comp = comp ;
		while current_comp != connect_dict[sect_connect_dict[comp]]

			l_current_comp = sqrt((comp_dict[current_comp].x - comp_dict[connect_dict[current_comp]].x)^2 +
					  (comp_dict[current_comp].y - comp_dict[connect_dict[current_comp]].y)^2 +
					  (comp_dict[current_comp].z - comp_dict[connect_dict[current_comp]].z)^2 ) ;

			l_eq += l_current_comp ;
			Ra_eq += 4.0 * RA * l_current_comp / (pi * comp_dict[current_comp].d^2)
			n_comp += 1.0 ;

			current_comp = connect_dict[current_comp] ;
		end

		d_eq = sqrt(RA * l_eq / (pi * Ra_eq)) ;
		Ra_eq = Ra_eq / n_comp ;

		return reduced_comp_t(l_eq, d_eq, Ra_eq, n_comp, comp, current_comp, ctype) ;
	else
	=#

	#=
		la_sum = 0.0 ;
		a_sum = 0.0 ;
		d_eq = 0.0 ;
		Ra_prod = 1.0 ;
		Ra_sum = 0.0 ;
		n_comp = 0.0 ;
		Y_merge_chk = false ;
		for c in connect_comp_v
			if comp_dict[c].ctype == ctype && strahler_dict[c] <= strahler_thrs
				sect_start_connect_comp = find_value_dict(c, sect_connect_dict) ;

				reduced_connect_comp = rec_reduced_comp(connect_dict, n_connect_dict, sect_connect_dict,
														other_type_cluster_root_v, strahler_dict, strahler_thrs,
														reduced_comp_v, sect_start_connect_comp[1], ctype) ;
				#if reduced_connect_comp != nothing
					a_sum += pi * reduced_connect_comp.l * reduced_connect_comp.d ;
					la_sum += reduced_connect_comp.l * pi * reduced_connect_comp.l * reduced_connect_comp.d ;
					d_eq += reduced_connect_comp.d ^ 2.0 ;
					Ra_prod *= reduced_connect_comp.Ra ;
					Ra_sum += reduced_connect_comp.Ra ;
					n_comp += reduced_connect_comp.n_comp ;

					Y_merge_chk = true ;
				#end

			elseif comp_dict[c].ctype == ctype && strahler_dict[c] > strahler_thrs
				current_comp = find_value_dict(c, sect_connect_dict) ;

				while current_comp != connect_dict[c]

					l_current_comp = sqrt((comp_dict[current_comp].x - comp_dict[connect_dict[current_comp]].x)^2 +
							  (comp_dict[current_comp].y - comp_dict[connect_dict[current_comp]].y)^2 +
							  (comp_dict[current_comp].z - comp_dict[connect_dict[current_comp]].z)^2 ) ;

					Ra_current_comp = 4.0 * RA * l_current_comp / (pi * comp_dict[current_comp].d^2) ;

					push!(reduced_comp_v, reduced_comp_t(l_current_comp, comp_dict[current_comp].d, Ra_current_comp,
													1.0, current_comp, connect_dict[current_comp], ctype)) ;

					current_comp = connect_dict[current_comp] ;
				end
			else
				if ~any(x -> x == comp, other_type_cluster_root_v) && ctype == :pdend
					push!(other_type_cluster_root_v, comp) ;
				end
			end
		end

		if strahler_dict[comp] <= strahler_thrs && reduced_connect_comp.ctype == comp_dict[comp].ctype
			if Y_merge_chk
				l_eq = la_sum / a_sum ;
				Ra_eq = Ra_prod / Ra_sum ;
				d_eq = sqrt(d_eq) ;
			else
				l_eq = 0.0 ;
				Ra_eq = 0.0 ;
				d_eq = 0.0 ;
			end

			current_comp = comp ;
		#if reduced_connect_comp.ctype == comp_dict[comp].ctype
			l_eq_sect = 0.0 ;
			Ra_eq_sect = 0.0 ;

			while current_comp != connect_dict[sect_connect_dict[comp]]

				l_current_comp = sqrt((comp_dict[current_comp].x - comp_dict[connect_dict[current_comp]].x)^2 +
						  (comp_dict[current_comp].y - comp_dict[connect_dict[current_comp]].y)^2 +
						  (comp_dict[current_comp].z - comp_dict[connect_dict[current_comp]].z)^2 ) ;

				l_eq_sect += l_current_comp ;
				Ra_eq_sect += 4.0 * RA * l_current_comp / (pi * comp_dict[current_comp].d^2)
				n_comp += 1.0 ;

				current_comp = connect_dict[current_comp] ;
			end

			d_eq_sect = sqrt(RA * l_eq_sect / (pi * Ra_eq_sect)) ;
			l_eq += l_eq_sect ;
			d_eq = sqrt(d_eq^2.0 + d_eq_sect^2.0) ;
			Ra_eq += Ra_eq_sect ;
		#end

			return reduced_comp_t(l_eq, d_eq, Ra_eq, n_comp, comp, current_comp, ctype) 
		end
	end
	=#


end

function strahler(strahler_thrs ::Int64)


	reduced_comp_v = Array{reduced_comp_t, 1}() ;

	file = open("pdend2.txt", "r")
	pdend = readdlm(file) ;
	close(file)

	file = open("ddend2.txt", "r")
	ddend = readdlm(file) ;
	close(file)

	file = open("axis.txt", "r")
	axis = readdlm(file) ;
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
										1.0, pi * comp_dict["soma"].d^2.0, "soma", "", :s)) ;

	l_axis = 0.0 ;
	Ra_axis = 0.0 ;
	n_comp = 0.0 ;
	area_sum = 0.0 ;
	for i = 1 : size(axis,1)-1

		current_l_axis = sqrt((axis[end - i + 1,1] - axis[end - i,1])^2 + 
						 (axis[end - i + 1,2] - axis[end - i,2])^2 +
						 (axis[end - i + 1,3] - axis[end - i,3])^2) ;

		l_axis += current_l_axis * 1e-6 ;
		Ra_axis += 4.0 * RA * current_l_axis * 1e-6 / (pi * (axis[i,4] * 1e-6)^2.0) ;
		n_comp += 1.0 ;
		area_sum += pi * current_l_axis * 1e-6 * axis[i,4] * 1e-6 ; 
	end

	d_axis = sqrt(RA * l_axis / (pi * Ra_axis)) ;

	push!(reduced_comp_v, reduced_comp_t(l_axis, d_axis, Ra_axis, n_comp, area_sum, "axis", "axhill", :axis)) ;

	push!(reduced_comp_v, reduced_comp_t(5.0e-6, 4.75e-6, 4.0 * RA * 5.0e-6 / (pi * (4.75e-6)^2.0), 1.0, 
										pi * 5.0e-6 * 4.75e-6, "axhill", "soma", :axhill)) ;

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
		next_comp = current_comp ;

		while n_connect_dict[connect_dict[next_comp]] == 1 && 
			comp_dict[connect_dict[next_comp]].ctype == comp_dict[current_comp].ctype

				next_comp = connect_dict[next_comp] ;
		end

		sect_connect_dict[current_comp] = next_comp ;

		if (strahler_dict[connect_dict[next_comp]] > strahler_thrs) &&
			~any(x -> x == connect_dict[next_comp], cluster_root_v)

			push!(cluster_root_v, connect_dict[next_comp])
		end

		if ~any(x -> x == connect_dict[next_comp], closed_v) &&
		   ~any(x -> x == connect_dict[next_comp], open_v)			   
				push!(open_v, connect_dict[next_comp]) ;
		end

		push!(closed_v, current_comp) ;
		shift!(open_v) ;

	end

	for comp in cluster_root_v 

		other_type_cluster_root_v = Array{String,1}() ; 
		reduced_connect_comp = rec_reduced_comp(connect_dict, comp_dict, sect_connect_dict,
												other_type_cluster_root_v, strahler_dict, strahler_thrs,
												reduced_comp_v, comp, :pdend) ;
		if reduced_connect_comp != nothing
			push!(reduced_comp_v, reduced_connect_comp) ;
		end

		la_sum = 0.0 ;
		a_sum = 0.0 ;
		d_eq = 0.0 ;
		Ra_prod = 1.0 ;
		Ra_sum = 0.0 ;
		n_comp = 0.0 ;
		reduced_comp_chk = false ;
		comp_id = " ";
		for other_type_comp in other_type_cluster_root_v

			chk_cluster_root_v = Array{String,1}() ;

			reduced_connect_comp = rec_reduced_comp(connect_dict, comp_dict, sect_connect_dict,
														chk_cluster_root_v, strahler_dict, strahler_thrs,
														reduced_comp_v, other_type_comp, :ddend) ;
			if reduced_connect_comp != nothing
				a_sum += pi * reduced_connect_comp.l * reduced_connect_comp.d ;
				la_sum += reduced_connect_comp.l * pi * reduced_connect_comp.l * reduced_connect_comp.d ;
				d_eq += reduced_connect_comp.d ^ 2.0 ;
				Ra_prod *= reduced_connect_comp.Ra ;
				Ra_sum += reduced_connect_comp.Ra ;
				n_comp += reduced_connect_comp.n_comp ;

				reduced_comp_chk = true ;
				comp_id = reduced_connect_comp.id ;
			end
		end
		l_eq = la_sum / a_sum ;
		Ra_eq = Ra_prod / Ra_sum ;
		d_eq = sqrt(d_eq) ;

		if reduced_comp_chk
			other_type_reduced_comp = reduced_comp_t(l_eq, d_eq, Ra_eq, n_comp, a_sum, comp_id,
													reduced_comp_v[end].id, reduced_connect_comp.ctype) ;

			push!(reduced_comp_v, other_type_reduced_comp) ;
		end
	end

	nc = 0.0 ;
	for i = 1 : length(reduced_comp_v)
		nc += reduced_comp_v[i].n_comp ;
	end 
	nc -= reduced_comp_v[1].n_comp + reduced_comp_v[2].n_comp + reduced_comp_v[3].n_comp ;

	println(nc)
	println(size(pdend,1) + size(ddend,1))
	println(length(reduced_comp_v))

	return reduced_comp_v

end