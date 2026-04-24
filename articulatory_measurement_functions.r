
measure_blade_root <- function(us_data, speakers=NULL){

	if (is.null(speakers)) speakers=names(us_data)

	for (sp in speakers){
		print(sp)
		plane='sag'
		tonguedata = us_data[[sp]][[plane]]$tongue_traces
		if (nrow(tonguedata)){
			tonguedata$root_advancement = with(tonguedata, (X1+X2+X3)/3)
			tonguedata$blade_angle = with(tonguedata, atan((Y10-Y8)/(X10/X8)))
			us_data[[sp]][[plane]]$tongue_traces = tonguedata
		}
	}
	us_data
}

measure_concavity <- function(us_data, speakers=NULL, indices=1:11){

	if (is.null(speakers)) speakers=names(us_data)

	require(pracma)
	for (sp in speakers){
		print(sp)
		us_data[[sp]]$sag$tongue_traces$absolute_concavity = NA
		us_data[[sp]]$sag$tongue_traces$relative_concavity = NA
		for (r in 1:nrow(us_data[[sp]]$sag$tongue_traces)){

			one_token = us_data[[sp]]$sag$tongue_traces[r,]
			original_polygon = data.frame(X=c(us_data[[sp]]$sag$origin[1], as.numeric(one_token[,paste0('X',indices)])), 
				                          Y=c(us_data[[sp]]$sag$origin[2], as.numeric(one_token[,paste0('Y',indices)])))
			convex_indices <- chull(original_polygon$X, original_polygon$Y)
			tongue_area = round(abs(polyarea(original_polygon$X,original_polygon$Y)),6)
			convex_area = round(abs(polyarea(original_polygon$X[convex_indices],original_polygon$Y[convex_indices])),6)
			# print (c(r,tongue_area,convex_area))
			us_data[[sp]]$sag$tongue_traces$absolute_concavity[r] = convex_area - tongue_area
			us_data[[sp]]$sag$tongue_traces$relative_concavity[r] = sqrt((convex_area-tongue_area)/tongue_area)
			# print (c(r,tongue_area,convex_area,absolute_concavity,relative_concavity))
		}
	}
	us_data
}





measure_velocity <- function(us_data, speakers=NULL, plane='sag', indices=1:11){

	require(matrixStats)

	if (is.null(speakers)) speakers=names(us_data)

	for (sp in speakers){
		print(sp)

		if (plane%in%c('sag','cor')){
			trace_data = us_data[[sp]][[plane]]$tongue_traces
			xcols = paste0('X',indices)
			ycols = paste0('Y',indices)
		}else if (plane=='video'){
			trace_data = us_data[[sp]][[plane]]$lip_traces
			xcols = names(trace_data)[grepl('_x$',names(trace_data))]
			ycols = names(trace_data)[grepl('_y$',names(trace_data))]
		}

		# identify non-sequential rows to exclude
		all_frame_steps = c(NA,diff(trace_data$Time_of_sample_in_recording))
		median_frame_step = median(all_frame_steps, na.rm=TRUE)
		is_a_step = c(all_frame_steps > 0 & all_frame_steps < median_frame_step*3)
		is_a_step[1] = FALSE

		# calculate differences and exclude non-sequential rows
		trace_dX = rbind(NA, colDiffs(as.matrix(trace_data[,xcols])))
		trace_dY = rbind(NA, colDiffs(as.matrix(trace_data[,ycols])))
		trace_dX[!is_a_step] = NA
		trace_dY[!is_a_step] = NA

		# calculate displacement and mean velocity (probably in mm/s)
		trace_displacement = sqrt(trace_dX^2 + trace_dY^2)
		trace_velocity = rowMeans(trace_displacement) / all_frame_steps

		if (plane%in%c('sag','cor')){
			us_data[[sp]][[plane]]$tongue_traces$tongue_velocity = trace_velocity
		}else if (plane=='video'){
			us_data[[sp]][[plane]]$lip_traces$lip_velocity = trace_velocity
		}
	}
	us_data
}
