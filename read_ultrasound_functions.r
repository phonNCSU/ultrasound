
# updated 2023-12-03 to handle the headers created on t-phon-1

colMax <- function (colData, na.rm=FALSE) {
    apply(colData, MARGIN=c(2), max, na.rm=na.rm)
}

colMin <- function (colData, na.rm=FALSE) {
    apply(colData, MARGIN=c(2), min, na.rm=na.rm)
}

colQuantile <- function (colData, prob=0.95, na.rm=TRUE) {
    apply(colData, MARGIN=c(2), quantile, probs=prob, na.rm=na.rm)
}

# make.cartesian <- function(tr, origin=c(0,0), flip=TRUE){    
#   X <- apply(tr, 1, function(x,y) origin[1]-x[2]*cos(x[1]))
#   Y <- apply(tr, 1, function(x,y) x[2]*sin(x[1])-origin[2])
#   # this was added 12/13/2021 to prevent positive unflipped values from coming back with negative Y values
#   if (flip){
#     xy <- cbind(X, Y)
#   }else{
#     xy <- cbind(X, -Y)
#   }
#   return(xy)
# }


aaa2ssanova <- function(data_wide, thetas){

	# take data in wide format and put it in long format for SSANOVA 

	data_long <- c()

	for (row in 1:nrow(data_wide)){
	    one_row <- data_wide[row,]

	    R_cols <- grepl('R[0-9]+',names(one_row))
	    C_cols <- grepl('C[0-9]+',names(one_row))

	    other_colnames <- names(one_row)[!(R_cols|C_cols)]

	    one_token <- data.frame(T=thetas, R=as.numeric(one_row[,R_cols]))

	    if (length(C_cols)){
	    	one_token$Confidence <- as.numeric(one_row[,C_cols])
	    }
	    one_token$token <- paste(one_row$Annotation_Label, one_row$Sequence_number, sep='_')
	    for (colname in other_colnames){
	    	one_token[[colname]] <- one_row[[colname]]
	    }
	    data_long <- rbind(data_long, one_token)
    }

	data_long_xy <- make.cartesian(data_long[,c('T','R')])
	data_long$X <- data_long_xy[,2]
	data_long$Y <- -data_long_xy[,1]

	data_long
}

parse_aaa_export <- function(export_filename){

	# read the AAA export file and format in a way that we can work with. 
	# returns a data frame still in wide format

	# read the txt file and separate the header from the data. 
	# The header doesn't have as many columns as there are columns of data
	# because many headers are for 42 or 84 columns of data

	if (file.exists(export_filename)){
		print (paste('reading',export_filename))
		aaa_raw_data <- read.table(export_filename, sep='\t', quote="", fill=TRUE, header=FALSE)
		aaa_raw_header <- paste(aaa_raw_data[1,!is.na(aaa_raw_data[1,])])
		aaa_data <- aaa_raw_data[2:nrow(aaa_raw_data),]

		found_columns <- ncol(aaa_data)
		print(aaa_raw_header)
		# TEMPORARY HACK TO REPLACE THE HEADER WITH THE LAPTOP-STYLE HEADER
		aaa_raw_header = gsub('Client family name', 'Client Surname', aaa_raw_header)
		aaa_raw_header = gsub('Date Time of recording', 'Date and time of recording', aaa_raw_header)
		aaa_raw_header = gsub('Annotation Title', 'Annotation Label', aaa_raw_header)
    if ('Sequence ' %in% aaa_raw_header){
      aaa_raw_header[aaa_raw_header=='Sequence '] = 'Sequence number'
    }
		
		# aaa_raw_header = gsub('Sequence ', 'Sequence number', aaa_raw_header)
		# aaa_raw_header = gsub('DLC_Tongue has no radius(%)', 'R% values of spline \"DLC_Tongue\"', aaa_raw_header)
    # aaa_raw_header = gsub('X0 DLC_Tongue	Y0 DLC_Tongue	X1 DLC_Tongue	Y1 DLC_Tongue	X2 DLC_Tongue	Y2 DLC_Tongue	X3 DLC_Tongue	Y3 DLC_Tongue	X4 DLC_Tongue	Y4 DLC_Tongue	X5 DLC_Tongue	Y5 DLC_Tongue	X6 DLC_Tongue	Y6 DLC_Tongue	X7 DLC_Tongue', 'X,Y values of spline \"DLC_Tongue\"', aaa_raw_header)
		# aaa_raw_header = c(aaa_raw_header, 'R% values of spline \"DLC_Tongue\"', 'X,Y values of spline \"DLC_Tongue\"', 'Confidence values of spline \"DLC_Tongue\"')
		# aaa_raw_header = c(aaa_raw_header, 'X,Y values of spline \"DLC_Tongue\"', 'Confidence values of spline \"DLC_Tongue\"')
		aaa_raw_header = c(aaa_raw_header, 'X,Y values of spline \"DLC_Tongue\"')
		print('RAW HEADER')
		print(aaa_raw_header)
		
		# if there are columns in the header that aren't listed here, they will be assumed to have only one column of data
		# addional column names should be added here as needed
		expected_columns_for_label <- list()	
		expected_columns_for_label[["Client Surname"]] = 1
		expected_columns_for_label[["Time of sample in recording"]] = 1
		expected_columns_for_label[["Date and time of recording"]] = 1
		expected_columns_for_label[["Annotation Label"]] = 1
		expected_columns_for_label[["Sequence number"]] = 1
		# expected_columns_for_label[["TD"]] = 1
		# expected_columns_for_label[["TR"]] = 1
		# expected_columns_for_label[["TT"]] = 1
		# expected_columns_for_label[["LT"]] = 1
		# expected_columns_for_label[["MT"]] = 1
		# expected_columns_for_label[["RT"]] = 1
		# expected_columns_for_label[["Prompt"]] = 1
		# expected_columns_for_label[["US Frame number of ultrasonic data"]] = 1
		# expected_columns_for_label[["Image Frame of ultrasonic data"]] = 1
		expected_columns_for_label[["X,Y values of spline \"DLC_Tongue\""]] = 22
		# expected_columns_for_label[["Angle values of spline \"Tongue\""]] = 42
		# expected_columns_for_label[["R mm values of spline \"Tongue\""]] = 42
		expected_columns_for_label[["R% values of spline \"DLC_Tongue\""]] = 1
		expected_columns_for_label[["Confidence values of spline \"DLC_Tongue\""]] = 11

		# make a table of the columns found in the header and what columns they appear to refer to in the data
		aaa_columns <- data.frame(label=aaa_raw_header, columns=NA, first=NA, last=NA)

		# build a header to match the data
		# note that R mm and R % will both be labeled as just "R"
		expected_columns = 0
		aaa_header = c()
		for (header_column in aaa_raw_header){
			if (header_column %in% names(expected_columns_for_label)){
				aaa_columns[aaa_columns$label==header_column,'columns'] <- expected_columns_for_label[[header_column]]
				aaa_columns[aaa_columns$label==header_column,'first'] <- expected_columns + 1
				aaa_columns[aaa_columns$label==header_column,'last'] <- expected_columns + expected_columns_for_label[[header_column]]
				expected_columns = expected_columns + expected_columns_for_label[[header_column]]

				if (header_column %in% c("R mm values of spline \"DLC_Tongue\"", "R% values of spline \"DLC_Tongue\"")){
					aaa_header <- c(aaa_header, paste0('R',0))
					print ("WHAT!!!!")
				}else if (header_column %in% c("X,Y values of spline \"DLC_Tongue\"", "R% values of spline \"DLC_Tongue\"")){
					aaa_header <- c(aaa_header, paste0(c('X','Y'),sort(rep(seq(0:10),2))))
				}else if (header_column %in% c("X,Y values of spline \"DLC_Tongue\"", "Confidence values of spline \"DLC_Tongue\"")){			
					aaa_header <- c(aaa_header, paste0('C',0:10))
				}else if (header_column %in% c("Angle values of spline \"DLC_Tongue\"")){			
					aaa_header <- c(aaa_header, paste0('A',0:10))
				}else{
					aaa_header <- c(aaa_header, header_column)
				}
			}else{
				warning(paste('found unexpected column',header_column,'assuming it only has one column (might be wrong!)'))
				aaa_columns[aaa_columns$label==header_column,'columns'] <- 1
				aaa_columns[aaa_columns$label==header_column,'first'] <- expected_columns + 1
				aaa_columns[aaa_columns$label==header_column,'last'] <- expected_columns + 1
				expected_columns = expected_columns + 1
				aaa_header <- c(aaa_header, header_column)
			}
		}

		aaa_header <- gsub(' ', '_', aaa_header)

		print (aaa_columns)
		print (paste('found rows:', nrow(aaa_data)))
		print (paste('found columns:', found_columns))
		print (paste('expected columns:', expected_columns))
    	print(aaa_header)
		colnames(aaa_data) <- aaa_header

		# make sure the multi-column data is numeric
		for (dataprefix in c('R','C','X','Y','A')){
			for (i in 0:41){
				if (paste0(dataprefix,i) %in% aaa_header){
					aaa_data[[paste0(dataprefix,i)]] <- as.numeric(paste(aaa_data[[paste0(dataprefix,i)]]))
				}
			}
		}

		# HACK
		if ('R0' %in% names(aaa_data)){
		  aaa_data$R0 = NULL
		}
		
		# format the date and time columns
		aaa_data$Time_of_sample_in_recording <- as.numeric(paste(aaa_data$Time_of_sample_in_recording))
		aaa_data$date_time <- strptime(paste(aaa_data$Date_and_time_of_recording), format='%m/%d/%Y %H:%M:%S %p')

		# sort by recording time and then by time in recording
		aaa_data <- aaa_data[order(aaa_data$Time_of_sample_in_recording),]
		aaa_data <- aaa_data[order(aaa_data$date_time),]

	}else{
		print (paste('did not find',export_filename))
		print (' (this is only a problem if you expected to have this type of data)')
		aaa_data <- NULL
	}

	aaa_data
}


parse_aaa_export_old <- function(export_filename){

	# read the AAA export file and format in a way that we can work with. 
	# returns a data frame still in wide format

	aaa_raw_data <- read.table(export_filename, sep='\t', quote="", fill=TRUE, header=FALSE)
	aaa_raw_header <- aaa_raw_data[1,]
	aaa_data <- aaa_raw_data[2:nrow(aaa_raw_data),]
	# first_radius <- which(aaa_raw_header=='R mm values of spline "Tongue"')
	first_radius <- which(aaa_raw_header=='R% values of spline "Tongue"')

	first_confidence <- which(aaa_raw_header=='Confidence values of spline "Tongue"')
	# first_na_col <- min(which(is.na(aaa_raw_header=='R mm values of spline Tongue')))
	first_na_col <- min(which(is.na(aaa_raw_header=='R% values of spline Tongue')))

	aaa_header <- c()
	for (i in 1:(first_radius-1)){
		aaa_header <- c(aaa_header, paste(aaa_raw_header[1,i]))
	}
	aaa_header <- c(aaa_header, paste0('R',0:41))
	aaa_header <- c(aaa_header, paste0(c('X','Y'),rep(seq(0:41),2)))
	aaa_header <- c(aaa_header, paste0('C',0:41))
	aaa_header <- gsub(' ', '_', aaa_header)
	colnames(aaa_data) <- aaa_header

	for (i in 0:41){
		aaa_data[[paste0('R',i)]] <- as.numeric(paste(aaa_data[[paste0('R',i)]]))
	}

	# aaa_data <- separate(aaa_data, 'Annotation_Label', c('phone','word'), sep='_', remove=FALSE)
	# aaa_data$phone <- factor(aaa_data$phone)
	# aaa_data$word <- factor(aaa_data$word)
	aaa_data$Time_of_sample_in_recording <- as.numeric(paste(aaa_data$Time_of_sample_in_recording))

	aaa_data$date_time <- strptime(paste(aaa_data$Date_and_time_of_recording), format='%m/%d/%Y %H:%M:%S %p')
	aaa_data <- aaa_data[order(aaa_data$Time_of_sample_in_recording),]
	aaa_data <- aaa_data[order(aaa_data$date_time),]

	aaa_data
}

parse_aaa_pal_occ <- function(pal_occ_filepath){

	if (file.exists(pal_occ_filepath)){
		print (paste('reading',pal_occ_filepath))
		pal_occ <- read.csv(pal_occ_filepath, sep='\t', quote="", fill=TRUE, header=FALSE)
		all_rows <- 1:nrow(pal_occ)
		pal_label_row <- grep('pal',pal_occ[,1], ignore.case = TRUE)
		occ_label_row <- grep('occ',pal_occ[,1], ignore.case = TRUE)
		
		if (!sum(pal_label_row)){
			pal_label_row <- grep('swallow',pal_occ[,1], ignore.case = TRUE)
		}
		if (!sum(occ_label_row)){
			occ_label_row <- grep('bite',pal_occ[,1], ignore.case = TRUE)
		}
		
		if (length(occ_label_row)==0){
			palate_trace <- pal_occ[all_rows>pal_label_row & pal_occ[all_rows,1]!='*',]
			occlusal_trace <- data.frame(X=NA,Y=NA)
		}else if (length(pal_label_row)==0){
			occlusal_trace <- pal_occ[all_rows>occ_label_row & pal_occ[all_rows,1]!='*',]
			palate_trace <- data.frame(X=NA,Y=NA)
		}else if (pal_label_row < occ_label_row){
			palate_trace <- pal_occ[all_rows>pal_label_row & all_rows<occ_label_row & pal_occ[all_rows,1]!='*',]
			occlusal_trace <- pal_occ[all_rows>occ_label_row & pal_occ[all_rows,1]!='*',]
		}else{
			occlusal_trace <- pal_occ[all_rows>occ_label_row & all_rows<pal_label_row & pal_occ[all_rows,1]!='*',]
			palate_trace <- pal_occ[all_rows>pal_label_row & pal_occ[all_rows,1]!='*',]	
		}

		palate_trace = na.omit(data.frame(X=palate_trace[,1], Y=palate_trace[,2]))
		occlusal_trace = na.omit(data.frame(X=occlusal_trace[,1], Y=occlusal_trace[,2]))

		occlusal_trace$X <- as.numeric(occlusal_trace$X)
		occlusal_trace$Y <- as.numeric(occlusal_trace$Y)
		palate_trace$X <- as.numeric(palate_trace$X)
		palate_trace$Y <- as.numeric(palate_trace$Y)

		if (length(occ_label_row)>0){
			occlusal_lm <- lm(Y~X,occlusal_trace)
			occlusal_angle <- 180*atan(coef(occlusal_lm)[2])/pi
		}else{
			occlusal_lm = NA
			occlusal_angle = -20
		}
		list(palate_trace=palate_trace, occlusal_trace=occlusal_trace, occlusal_lm=occlusal_lm, occlusal_angle=occlusal_angle)
	}else{
		print (paste('did not find',pal_occ_filepath))
		print (' (this is only a problem if you expected to have this type of data)')
		list(palate_trace=NULL, occlusal_trace=NULL, occlusal_lm=NULL, occlusal_angle=NULL)
	}
}


add_manual_pal_occ <- function(aaa_data, palate_filepath=NULL, occlusal_filepath=NULL, palate_scale=1, flip=TRUE, xmax=320, ymax=240){

	if (!is.null(palate_filepath)){
		print(paste('preparing to read from',palate_filepath))
		all_palate_traces = read.csv(palate_filepath)

		all_palate_traces$x = all_palate_traces$x*palate_scale
		all_palate_traces$y = all_palate_traces$y*palate_scale
		if (flip){
			all_palate_traces$y = ymax - all_palate_traces$y*palate_scale
		}

	}

	if (!is.null(occlusal_filepath)){
		print(paste('preparing to read from',occlusal_filepath))
			all_occlusal_angles = read.csv(occlusal_filepath)
	}

	rotate_center = c(xmax/2,ymax/2)

	for (sp in names(aaa_data)){
		aaa_data[[sp]]$sag$occlusal_angle = all_occlusal_angles[all_occlusal_angles$speaker==sp,'occlusal_angle']
		palate_trace = all_palate_traces[all_palate_traces$speaker==sp,c('x','y')]
		names(palate_trace) = c('X','Y')
		aaa_data[[sp]]$sag$palate_trace = palate_trace
		print(paste0(sp,': rotating counterclockwise ',aaa_data[[sp]]$sag$occlusal_angle,' degrees around (',rotate_center[1],',',rotate_center[2],')'))
		aaa_data[[sp]]$sag$palate_trace_rotated = rotateXY(aaa_data[[sp]]$sag$palate_trace, aaa_data[[sp]]$sag$occlusal_angle, center=rotate_center)
	}

	aaa_data
}


update_pal_occ <- function(aaa_data, sp, plane, pal_occ_filepath, center=c(0,0)){
	pal_occ = parse_aaa_pal_occ(pal_occ_filepath)


	if (!is.null(pal_occ$occlusal_trace)){
		aaa_data[[sp]][[plane]][['occlusal_trace']] <- pal_occ$occlusal_trace
		aaa_data[[sp]][[plane]][['occlusal_lm']] <- pal_occ$occlusal_lm
		aaa_data[[sp]][[plane]][['occlusal_angle']] <- pal_occ$occlusal_angle
	}
	if (!is.null(pal_occ$palate_trace)){
		aaa_data[[sp]][[plane]][['palate_trace']] <- pal_occ$palate_trace
	}

	if (is.null(aaa_data[[sp]][[plane]]$palate_trace)){
		aaa_data[[sp]][[plane]][['palate_trace_rotated']] = NULL
	}else{
		aaa_data[[sp]][[plane]][['palate_trace_rotated']] = rotateXY(aaa_data[[sp]][[plane]]$palate_trace, aaa_data[[sp]][[plane]]$occlusal_angle, center=center)
	}

	if (is.null(aaa_data[[sp]][[plane]]$occlusal_trace)){
		aaa_data[[sp]][[plane]][['occlusal_trace_rotated']] = NULL
	}else{
		aaa_data[[sp]][[plane]][['occlusal_trace_rotated']] = rotateXY(aaa_data[[sp]][[plane]]$occlusal_trace, aaa_data[[sp]][[plane]]$occlusal_angle, center=center)
	}
	aaa_data
}

read_cl_export <- function(){

	# Make a table of all the exported txt and wav files and their textgrids
	clip_info = data.frame()

	subdirs = list.files('/home/jeff/covart/cl/', pattern='^export_cl.*')
	for (subdir in subdirs){
		print (subdir)
		path_to_files = paste0('/home/jeff/covart/cl/', subdir, '/')
		txt_filelist = list.files(path_to_files, pattern='.*\\.txt')
		# wav_filelist = list.files(paste0('/home/jeff/covart/cl/',subdir), pattern='.*\\.wav') 
		# tg_filelist = list.files(paste0('/home/jeff/covart/cl/',subdir), pattern='.*\\.TextGrid') 
		for (filename in txt_filelist){

			txt_file = paste0(path_to_files, filename)
			txt_content = readChar(txt_file, file.info(txt_file)$size)
			txt_split = unlist(strsplit(txt_content,'\\r\\n'))
			prompt_text = txt_split[1]
			date_time = txt_split[2]
			client_name = unlist(strsplit(txt_split[3],','))[1]

			# look for the textgrid
			expected_tg_filename = gsub('.txt', '.TextGrid', filename)
			expected_tg_path = paste0(path_to_files, expected_tg_filename)
			if (file.exists(expected_tg_path)){
				textgrid_filename = expected_tg_filename
			}else{
				textgrid_filename = NA
			}

			# look for wav files
			expected_hq_wav_filename = gsub('.txt', '_Track0_hq.wav', filename)
			expected_mq_wav_filename = gsub('.txt', '_Track0_mq.wav', filename)
			expected_lq_wav_filename = gsub('.txt', '_Track0.wav', filename)

			expected_hq_wav_path = paste0(path_to_files, expected_hq_wav_filename)
			expected_mq_wav_path = paste0(path_to_files, expected_mq_wav_filename)
			expected_lq_wav_path = paste0(path_to_files, expected_lq_wav_filename)

			if (file.exists(expected_hq_wav_path)){
				wav_filename = expected_hq_wav_filename
			}else if (file.exists(expected_mq_wav_path)){
				wav_filename = expected_mq_wav_filename
			}else if (file.exists(expected_lq_wav_path)){
				wav_filename = expected_lq_wav_filename
			}else{
				wav_filename = NA
			}
			clip_info <- rbind(clip_info, data.frame(client=client_name, 
				                                     date_time=date_time, 
				                                     export_path=path_to_files,
				                                     wav_filename=wav_filename, 
				                                     textgrid_filename=textgrid_filename, 
				                                     prompt_text=prompt_text))
		}
	}
	clip_info
}

#########################################################################
# CL-specific functions I put here when simplifying starting 2-9-22
#########################################################################

read_all_aaa_data <- function(speakers, aaa_data=list(), planes=c('sag'), offset=NULL, center=c(0,0)){

	for (sp in speakers){
		print(sp)
		for (plane in planes){
			print(plane)
			aaa_data[[sp]][[plane]] = list()

			# read and parse the exported data
			export_filepath <- paste0(readwritepath, sp, '_', plane, '_tongues.txt')
			pal_occ_filepath <- paste0(readwritepath,sp,'_', plane, '_pal_occ.txt')

			aaa_data[[sp]][[plane]][['tongue_traces']] = parse_aaa_export(export_filepath)

			pal_occ = parse_aaa_pal_occ(pal_occ_filepath)

			aaa_data[[sp]][[plane]][['palate_trace']] <- pal_occ$palate_trace
			aaa_data[[sp]][[plane]][['occlusal_trace']] <- pal_occ$occlusal_trace
			aaa_data[[sp]][[plane]][['occlusal_lm']] <- pal_occ$occlusal_lm
			aaa_data[[sp]][[plane]][['occlusal_angle']] <- pal_occ$occlusal_angle
			
			if (!is.null(offset)){
				aaa_data[[sp]][[plane]][['tongue_traces']][,'Time_of_sample_in_recording'] = aaa_data[[sp]][[plane]][['tongue_traces']][,'Time_of_sample_in_recording'] + offset[[sp]]
			}
			valid_data_points <- aaa_data[[sp]][[plane]]$tongue_traces[,grepl('C[0-9]+',names(aaa_data[[sp]][[plane]]$tongue_traces))] > 0

			C_names <- which(grepl('C[0-9]+',names(aaa_data[[sp]][[plane]]$tongue_traces)))

			if (is.null(aaa_data[[sp]][[plane]]$palate_trace)){
				aaa_data[[sp]][[plane]][['palate_trace_rotated']] = NULL
			}else{
				aaa_data[[sp]][[plane]][['palate_trace_rotated']] = rotateXY(aaa_data[[sp]][[plane]]$palate_trace, aaa_data[[sp]][[plane]]$occlusal_angle, center=center)
			}

			if (is.null(aaa_data[[sp]][[plane]]$occlusal_trace)){
				aaa_data[[sp]][[plane]][['occlusal_trace_rotated']] = NULL
			}else{
				aaa_data[[sp]][[plane]][['occlusal_trace_rotated']] = rotateXY(aaa_data[[sp]][[plane]]$occlusal_trace, aaa_data[[sp]][[plane]]$occlusal_angle, center=center)
			}

			if (length(C_names)>0){
				aaa_data[[sp]][[plane]]$tongue_traces[,grepl('X[0-9]+',names(aaa_data[[sp]][[plane]]$tongue_traces))][!valid_data_points] <- NA
				aaa_data[[sp]][[plane]]$tongue_traces[,grepl('Y[0-9]+',names(aaa_data[[sp]][[plane]]$tongue_traces))][!valid_data_points] <- NA
				# aaa_data[[sp]][[plane]]$tongue_traces[,grepl('R[0-9]+',names(aaa_data[[sp]][[plane]]$tongue_traces))][!valid_data_points] <- NA
			}else{

			}
		
			aaa_data[[sp]][[plane]]$tongue_traces$token <- NA
			for (unique_label in unique(aaa_data[[sp]][[plane]]$tongue_traces$Annotation_Label)){
				print(unique_label)
				date_times <- subset(aaa_data[[sp]][[plane]]$tongue_traces, Annotation_Label==unique_label)$Date_and_time_of_recording
				aaa_data[[sp]][[plane]]$tongue_traces[aaa_data[[sp]][[plane]]$tongue_traces$Annotation_Label==unique_label,'token'] <- paste0(unique_label, '_', plane, as.numeric(factor(date_times)))
			}

			write.csv(aaa_data[[sp]][[plane]][['tongue_traces']], paste0(readwritepath, sp, '_', plane, '_tongues_parsed.csv'))
		}	
	}
	aaa_data
}


add_analysis_angles <- function(filepath, aaa_data, speakers=NULL, angles_only=FALSE){

	analysis_value_radii = read.csv(filepath, na.strings='NA')

	if (is.null(speakers)) speakers=names(aaa_data)

	universal_multiplier = 0.6953337

	for (sp in speakers){
		print(sp)
		for (plane in names(aaa_data[[sp]])){
			if (plane=='sag'){
			  if (angles_only){
			    aaa_data[[sp]][[plane]]$analysis_value_radii = analysis_value_radii[analysis_value_radii$speaker==sp,grepl('angle',names(analysis_value_radii))]
			  }else{
  				aaa_data[[sp]][[plane]]$analysis_value_radii = analysis_value_radii[analysis_value_radii$speaker==sp,c('TT','TD','TR','sag_multiplier','TTangle','TDangle','TRangle')]
  				if ('TT'%in%names(aaa_data[[sp]][[plane]]$tongue_traces)){
  					aaa_data[[sp]][[plane]]$tongue_traces$TT2 = aaa_data[[sp]][[plane]]$tongue_traces[,paste0('R',aaa_data[[sp]][[plane]]$analysis_value_radii$TT)]/universal_multiplier
  					aaa_data[[sp]][[plane]]$tongue_traces$TD2 = aaa_data[[sp]][[plane]]$tongue_traces[,paste0('R',aaa_data[[sp]][[plane]]$analysis_value_radii$TD)]/universal_multiplier
  					aaa_data[[sp]][[plane]]$tongue_traces$TR2 = aaa_data[[sp]][[plane]]$tongue_traces[,paste0('R',aaa_data[[sp]][[plane]]$analysis_value_radii$TR)]/universal_multiplier
  				}else{
  					aaa_data[[sp]][[plane]]$tongue_traces$TT = aaa_data[[sp]][[plane]]$tongue_traces[,paste0('R',aaa_data[[sp]][[plane]]$analysis_value_radii$TT)]/universal_multiplier
  					aaa_data[[sp]][[plane]]$tongue_traces$TD = aaa_data[[sp]][[plane]]$tongue_traces[,paste0('R',aaa_data[[sp]][[plane]]$analysis_value_radii$TD)]/universal_multiplier
  					aaa_data[[sp]][[plane]]$tongue_traces$TR = aaa_data[[sp]][[plane]]$tongue_traces[,paste0('R',aaa_data[[sp]][[plane]]$analysis_value_radii$TR)]/universal_multiplier
  				}
			  }
			}else{
			  if (angles_only){
			    aaa_data[[sp]][[plane]]$analysis_value_radii = analysis_value_radii[analysis_value_radii$speaker==sp,grepl('angle',names(analysis_value_radii))]
			  }else{
  				aaa_data[[sp]][[plane]]$analysis_value_radii = analysis_value_radii[analysis_value_radii$speaker==sp,c('LT','MT','RT','cor_multiplier')]
  				if ('LT'%in%names(aaa_data[[sp]][[plane]]$tongue_traces)){
  					aaa_data[[sp]][[plane]]$tongue_traces$LT2 = aaa_data[[sp]][[plane]]$tongue_traces[,paste0('R',aaa_data[[sp]][[plane]]$analysis_value_radii$LT)]/universal_multiplier
  					aaa_data[[sp]][[plane]]$tongue_traces$MT2 = aaa_data[[sp]][[plane]]$tongue_traces[,paste0('R',aaa_data[[sp]][[plane]]$analysis_value_radii$MT)]/universal_multiplier
  					aaa_data[[sp]][[plane]]$tongue_traces$RT2 = aaa_data[[sp]][[plane]]$tongue_traces[,paste0('R',aaa_data[[sp]][[plane]]$analysis_value_radii$RT)]/universal_multiplier
  				}
			  }
			}
		}
	}
	aaa_data
}


add_textgrid_segmentation <- function(aaa_data, speakers=NULL, merge_vl=FALSE, planes=c('sag'), 
	tiers=c('phone','word'), match_words_to_phrase=TRUE, middle_is_quarter=c(), simple_tiers=c('YPR'), one_folder=NULL){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		print(sp)
		if (is.null(one_folder)){
			print(paste0('looking for ','export_',sp,' within the directory ',readwritepath))
			export_path <- paste0(readwritepath,'/export_',sp)
			textgrid_filelist <- list.files(export_path, '*TextGrid')
		}else{
			print(paste0('looking for textgrid in ',one_folder,' matching ',paste0(sp,'.*TextGrid')))
			export_path <- one_folder
			textgrid_filelist <- list.files(export_path, paste0(sp,'.*TextGrid'))
		}
		

		for (plane in planes){

			
			dfname = ifelse(plane=='video', 'lip_traces', 'tongue_traces')
			# if('Annotation_Label'%in%names(aaa_data[[sp]][[plane]][[dfname]])){
			if (nrow(aaa_data[[sp]][[plane]][[dfname]])){
				aaa_data[[sp]][[plane]][[dfname]]$word <- NA
				aaa_data[[sp]][[plane]][[dfname]]$phone <- NA
				aaa_data[[sp]][[plane]][[dfname]]$middle_frame <- FALSE
				aaa_data[[sp]][[plane]][[dfname]]$token_id <- NA
				aaa_data[[sp]][[plane]][[dfname]]$exclude <- FALSE
				aaa_data[[sp]][[plane]][[dfname]]$word_start <- NA
				aaa_data[[sp]][[plane]][[dfname]]$word_end <- NA
				aaa_data[[sp]][[plane]][[dfname]]$phone_start <- NA
				aaa_data[[sp]][[plane]][[dfname]]$phone_end <- NA
				aaa_data[[sp]][[plane]][[dfname]]$phone_time <- NA
			}
		}
		if(!length(textgrid_filelist)){
			warning("No textgrid files found. Please check your folder name.")
		}
		for (fn in textgrid_filelist){
			print(paste('reading', fn))
			tg_path <- paste0(export_path,'/',fn)
			# print(tg_path)
			txt_path <- paste0(export_path,'/',gsub('TextGrid','txt',fn))
			# print('x')
			# textgrid  <- tg2df(read.TextGrid(tg_path),tiers)
			textgrid = tgfile2list(tg_path)
			# print('xx')
			names(textgrid) = gsub('words','word',names(textgrid))
			names(textgrid) = gsub('phones','phone',names(textgrid))
			
			textgrid$phone$text <- trimws(textgrid$phone$text)
			textgrid$word$text <- trimws(textgrid$word$text)

			if (file.exists(txt_path)){
				txt <- readChar(txt_path, file.info(txt_path)$size)
				txt_info <- unlist(strsplit(txt,'\r\n'))
				txt_prompts <- txt_info[1]
				txt_date_time <- txt_info[2]
				txt_subject <- txt_info[3]
			}else{
				txt_date_time = NA
			}

			if (merge_vl){
				word_final = which(textgrid$phone$xmax %in% textgrid$word$xmax)
				vowel = which(grepl('[012]', textgrid$phone$text))
				preliquid = which(textgrid$phone$text %in% c('L','R')) - 1

				preliquid_vowels = setdiff(intersect(vowel,preliquid),word_final)

				postvocalic_liquids = preliquid_vowels+1

				if (length(postvocalic_liquids)){
					textgrid$phone[preliquid_vowels,'text'] = paste0(textgrid$phone[preliquid_vowels,'text'], textgrid$phone[postvocalic_liquids,'text'])
					textgrid$phone[preliquid_vowels,'xmax'] = textgrid$phone[postvocalic_liquids,'xmax']
					textgrid$phone = textgrid$phone[-postvocalic_liquids,]
				}
			}

			# print (txt_info)
			for (plane in planes){
				dfname = ifelse(plane=='video', 'lip_traces', 'tongue_traces')
				if (nrow(aaa_data[[sp]][[plane]][[dfname]])){

					if(is.na(txt_date_time)){
						rownumbers_for_textgrid = 1:nrow(aaa_data[[sp]][[plane]][[dfname]])
					}else{
						rownumbers_for_textgrid = which(aaa_data[[sp]][[plane]][[dfname]]$Date_and_time_of_recording==txt_date_time)
					}

					rows_for_textgrid = aaa_data[[sp]][[plane]][[dfname]][rownumbers_for_textgrid,]

					for (tier in tiers){

						for (r in 1:nrow(textgrid[[tier]])){

							if (!r%%1000){
								print(paste0(plane,': processing ',tier,' tier interval ',r,' of ',nrow(textgrid[[tier]])))
							}

							interval_label <- textgrid[[tier]][r,'text']
							interval_xmin <- textgrid[[tier]][r,'xmin']
							interval_xmax <- textgrid[[tier]][r,'xmax']

							is_a_match = which(with(rows_for_textgrid, 
												Time_of_sample_in_recording>=interval_xmin
												 & Time_of_sample_in_recording<interval_xmax))

							matching_rows = rownumbers_for_textgrid[is_a_match]

							if (length(matching_rows)){
								aaa_data[[sp]][[plane]][[dfname]][matching_rows,tier] = interval_label
								if (!tier %in% simple_tiers){
									aaa_data[[sp]][[plane]][[dfname]][matching_rows,paste(tier,'start',sep='_')] <- interval_xmin
									aaa_data[[sp]][[plane]][[dfname]][matching_rows,paste(tier,'end',sep='_')] <- interval_xmax
									aaa_data[[sp]][[plane]][[dfname]][matching_rows,paste(tier,'id',sep='_')] <- paste(gsub('\\.TextGrid', '', fn), 1, interval_label, round(interval_xmin,3),sep='_')
								}
							}

							if (tier=='word'){
								# make sure this word is really in the target phrase. sometimes we get an adjacent word overlapping 
								# the target phrase interval by < 1 ms for some reason and we want to exclude those frames instead
								# of labeling them as something else

								# this is so that words like "a" don't happen to match a substring of a word in the phrase
								current_phrase = unlist(strsplit(tolower(unique(aaa_data[[sp]][[plane]][[dfname]][matching_rows,'Annotation_Label'])),'_'))
								is_in_phrase = tolower(gsub("'","",interval_label)) %in% current_phrase | interval_label == 'sp'

								# if (sum(na.omit(matches_word))){
								if (length(matching_rows)){
									if (match_words_to_phrase){
										if (!is_in_phrase){
											print (paste('excluding', sum(matching_rows), 'frame(s) of ', tier, 
												interval_label, 'from phrase', current_phrase))
											aaa_data[[sp]][[plane]][[dfname]][matching_rows & !is_in_phrase,'exclude'] = TRUE
										}
									}
								}
							

							}else if (tier=='phone'){

								if (length(matching_rows)){

									time_from_mid <- abs(aaa_data[[sp]][[plane]][[dfname]][matching_rows,'Time_of_sample_in_recording'] - mean(c(interval_xmin,interval_xmax)))
									time_from_quarter <- abs(aaa_data[[sp]][[plane]][[dfname]][matching_rows,'Time_of_sample_in_recording'] - (interval_xmin + 0.25*(interval_xmax-interval_xmin)))
									# if (sum(is_this_phone)){
									if (length(matching_rows)){
										# print (aaa_data[[sp]][[plane]]$tongue_traces[is_this_phone,'Time_of_sample_in_recording'])
										# print(time_from_mid)
										# print(c(phone_xmin,phone_xmax))
										# print(c(mean(c(phone_xmin,phone_xmax)), phone_xmin + 0.5*(phone_xmax-phone_xmin)))
										# print('***')
										if (interval_label%in%middle_is_quarter){
											# mid_frame <- which(is_this_phone)[which(time_from_quarter==min(time_from_quarter))[1]]
											mid_frame <- matching_rows[which(time_from_quarter==min(time_from_quarter))[1]]
										}else{
											# mid_frame <- which(is_this_phone)[which(time_from_mid==min(time_from_mid))[1]]
											mid_frame <- matching_rows[which(time_from_mid==min(time_from_mid))[1]]
										}
										aaa_data[[sp]][[plane]][[dfname]][mid_frame,'middle_frame'] <- TRUE
									}

									if (r>3){
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'left2'] = textgrid$phone[r-3,'text']
									}else{
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'left2'] = NA								
									}
									if (r>2){
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'left1'] = textgrid$phone[r-2,'text']
									}else{
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'left1'] = NA								
									}
									if (r>1){
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'left'] = textgrid$phone[r-1,'text']
									}else{
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'left'] = NA								
									}
									
									# aaa_data[[sp]][[plane]][[dfname]][matching_rows,'phone'] <- phone_label

									if (r<nrow(textgrid$phone)){
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'right'] = textgrid$phone[r+1,'text']
									}else{
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'right'] = NA								
									}
									if (r<(nrow(textgrid$phone)-1)){
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'right1'] = textgrid$phone[r+2,'text']
									}else{
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'right1'] = NA								
									}
									if (r<(nrow(textgrid$phone)-2)){
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'right2'] = textgrid$phone[r+3,'text']
									}else{
										aaa_data[[sp]][[plane]][[dfname]][matching_rows,'right2'] = NA								
									}

									# aaa_data[[sp]][[plane]][[dfname]][matching_rows,'phone_start'] <- phone_xmin
									# aaa_data[[sp]][[plane]][[dfname]][matching_rows,'phone_end'] <- phone_xmax
									aaa_data[[sp]][[plane]][[dfname]][matching_rows,'phone_start'] <- interval_xmin
									aaa_data[[sp]][[plane]][[dfname]][matching_rows,'phone_end'] <- interval_xmax
									word_label <- aaa_data[[sp]][[plane]][[dfname]][matching_rows,'word'][1]
									# print(paste(word_label,phone_label))
									aaa_data[[sp]][[plane]][[dfname]][matching_rows,'token_id'] <- paste(gsub('\\.TextGrid', '', fn), 1, word_label, interval_label, round(interval_xmin,3),sep='_')
								}
							}
							# print(names(aaa_data[[sp]][[plane]]$tongue_traces))
							# print('a')
							aaa_data[[sp]][[plane]][[dfname]]$phone_time = with(aaa_data[[sp]][[plane]][[dfname]], (Time_of_sample_in_recording-phone_start)/(phone_end-phone_start))
							# print('b')
						}
					}
				}
			}
		}
	}
	aaa_data
}

traces_wide_to_long <- function(tongue_data_wide, occlusal_angle=NULL, center=c(0,0), factors_to_retain=c('left','phone','right','word'), fast=FALSE, polar=FALSE){
	tongue_data_long <- c()
	columns_to_retain <- c('Annotation_Label','token','token_id','Time_of_sample_in_recording', factors_to_retain)
	
	if (fast){
		if (polar){
			tongue_data_long = data.frame(X=as.numeric(unlist(tongue_data_wide[,grepl('^T[0-9]',names(tongue_data_wide))])), 
					                      Y=as.numeric(unlist(tongue_data_wide[,grepl('^R[0-9]',names(tongue_data_wide))])))
		}else{
			tongue_data_long = data.frame(X=as.numeric(unlist(tongue_data_wide[,grepl('^X[0-9]',names(tongue_data_wide))])), 
					                      Y=as.numeric(unlist(tongue_data_wide[,grepl('^Y[0-9]',names(tongue_data_wide))])))
		}
	}else{
		for(i in 1:nrow(tongue_data_wide)){
			if (i/100==round(i/100)){
				print(paste(i,'of',nrow(tongue_data_wide)))
			}
			token_data_frame <- tongue_data_wide[i,intersect(columns_to_retain, colnames(tongue_data_wide))]
			rownames(token_data_frame) <- NULL
			if (polar){
				token_data_frame_XY <- data.frame(X=as.numeric(paste(tongue_data_wide[i,grepl('^T[0-9]',names(tongue_data_wide))])), 
					                              Y=as.numeric(paste(tongue_data_wide[i,grepl('^R[0-9]',names(tongue_data_wide))])))
			}else{
				token_data_frame_XY <- data.frame(X=as.numeric(paste(tongue_data_wide[i,grepl('^X[0-9]',names(tongue_data_wide))])), 
					                              Y=as.numeric(paste(tongue_data_wide[i,grepl('^Y[0-9]',names(tongue_data_wide))])))
			}
			
			tongue_data_long <- rbind(tongue_data_long, cbind(token_data_frame, token_data_frame_XY))
		}

		for (ftr in intersect(factors_to_retain, colnames(tongue_data_long))){
			tongue_data_long[,ftr] = factor(tongue_data_long[,ftr])
		}
	}
	tongue_data_long <- tongue_data_long[!is.na(tongue_data_long$X),]
	if (!is.null(occlusal_angle)){
		tongue_data_long <- rotateXY(tongue_data_long, occlusal_angle, center=center)
	}
	tongue_data_long
}

# put_middle_frames_in_long_format <- function(aaa_data, speakers=NULL, planes=c('sag')){

# 	if (is.null(speakers)) speakers=names(aaa_data)

# 	for (sp in speakers){
# 		for (plane in planes){
# 			print(paste(sp,plane))
# 			if ('Annotation_Label' %in% names(aaa_data[[sp]][[plane]]$tongue_traces)){
# 				# transform into long format
# 				tongue_data_wide <- subset(aaa_data[[sp]][[plane]]$tongue_traces, middle_frame==TRUE)
# 				occlusal_angle = aaa_data[[sp]][[plane]]$occlusal_angle
# 				tongue_data_long = traces_wide_to_long(tongue_data_wide, occlusal_angle)
# 				aaa_data[[sp]][[plane]]$tongue_traces_long_rotated <- tongue_data_long
# 			}else{
# 				print('no data (only a problem if you expected data of this type)')
# 			}
# 		}
# 	}
# 	aaa_data
# }




choose_polar_origin <- function(aaa_data, method='xmid_ymid', speakers=NULL, planes=c('sag'), flip=FALSE){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		for (plane in planes){
			# if ('Annotation_Label' %in% names(aaa_data[[sp]][[plane]]$tongue_traces)){
			#use the same polar origin for all comparisons and the same axis ranges for all the plots
			# tongue_data_wide <- subset(aaa_data[[sp]][[plane]]$tongue_traces, middle_frame==TRUE)
			tongue_data_wide <- subset(aaa_data[[sp]][[plane]]$tongue_traces)
			# occlusal_angle = aaa_data[[sp]][[plane]]$occlusal_angle
			tongue_data_long = traces_wide_to_long(tongue_data_wide, fast=TRUE)
			polar_origin <- select.origin(tongue_data_long$X, tongue_data_long$Y, 
										  tongue_data_long$token_id, method=method, flip=flip)
			aaa_data[[sp]][[plane]]$origin <- polar_origin
			# }
		}
	}
	aaa_data
}

xy2occlusal <- function(aaa_data, speakers=NULL, replaceXY=TRUE, planes=c('sag'), flip=FALSE, center=c(0,0)){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		print(sp)
		for (plane in planes){
			point_colnames = names(aaa_data[[sp]][[plane]]$tongue_traces)[grepl('_[xy]$',names(aaa_data[[sp]][[plane]]$tongue_traces))|grepl('^[XY][0-9]',names(aaa_data[[sp]][[plane]]$tongue_traces))]
			tongue_data_wide = aaa_data[[sp]][[plane]]$tongue_traces[,point_colnames]

			for (i in seq(1,length(point_colnames),2)){
				one_radius_raw = data.frame(X=tongue_data_wide[,i], Y=tongue_data_wide[,i+1])
				one_radius_rotated <- rotateXY(one_radius_raw, aaa_data[[sp]][[plane]]$occlusal_angle, center=center)
				if (i==1){
					aaa_data[[sp]][[plane]]$tongue_traces_rotated = one_radius_rotated
					names(aaa_data[[sp]][[plane]]$tongue_traces_rotated) = point_colnames[1:2]
				}else{
					aaa_data[[sp]][[plane]]$tongue_traces_rotated[,point_colnames[i]] = one_radius_rotated$X
					aaa_data[[sp]][[plane]]$tongue_traces_rotated[,point_colnames[i+1]] = one_radius_rotated$Y
				}			
			}	
			if (replaceXY){
				aaa_data[[sp]][[plane]]$tongue_traces[,point_colnames] = aaa_data[[sp]][[plane]]$tongue_traces_rotated
				aaa_data[[sp]][[plane]]$tongue_traces_rotated = NULL
			}		
		}
	}
	aaa_data
}

xy2polar <- function(aaa_data, speakers=NULL, replaceXY=FALSE, addTR=TRUE, planes=c('sag'), flip=FALSE){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		print(sp)
		for (plane in planes){
			point_colnames = names(aaa_data[[sp]][[plane]]$tongue_traces)[grepl('^[XY][0-9]',names(aaa_data[[sp]][[plane]]$tongue_traces))]
			tongue_data_wide = aaa_data[[sp]][[plane]]$tongue_traces[,point_colnames]
			polar_origin = aaa_data[[sp]][[plane]]$origin
			
			for (i in seq(1,length(point_colnames),2)){
				one_radius_raw = data.frame(X=tongue_data_wide[,i], Y=tongue_data_wide[,i+1])
				one_radius_polar <- make.polar(one_radius_raw, origin=polar_origin, flip=flip)

				if (i==1){
					aaa_data[[sp]][[plane]]$tongue_traces_polar = one_radius_polar
					names(aaa_data[[sp]][[plane]]$tongue_traces_polar) = point_colnames[1:2]
				}else{
					aaa_data[[sp]][[plane]]$tongue_traces_polar[,point_colnames[i]] = one_radius_polar$X
					aaa_data[[sp]][[plane]]$tongue_traces_polar[,point_colnames[i+1]] = one_radius_polar$Y
				}			
			}
			if (replaceXY){
				aaa_data[[sp]][[plane]]$tongue_traces[,point_colnames] = aaa_data[[sp]][[plane]]$tongue_traces_polar
				aaa_data[[sp]][[plane]]$tongue_traces_polar = NULL
			}	
			if (addTR){
				names(aaa_data[[sp]][[plane]]$tongue_traces_polar) = gsub('X','T',names(aaa_data[[sp]][[plane]]$tongue_traces_polar))
				names(aaa_data[[sp]][[plane]]$tongue_traces_polar) = gsub('Y','R',names(aaa_data[[sp]][[plane]]$tongue_traces_polar))

				aaa_data[[sp]][[plane]]$tongue_traces = aaa_data[[sp]][[plane]]$tongue_traces[,!grepl('[TR][0-9]',names(aaa_data[[sp]][[plane]]$tongue_traces))]

				aaa_data[[sp]][[plane]]$tongue_traces = cbind(aaa_data[[sp]][[plane]]$tongue_traces, aaa_data[[sp]][[plane]]$tongue_traces_polar)
				aaa_data[[sp]][[plane]]$tongue_traces_polar = NULL
			}
		}
	}		
	aaa_data
}

# xy2polar <- function(aaa_data, speakers=NULL, replaceXY=FALSE, planes=c('sag'), flip=FALSE){

# 	if (is.null(speakers)) speakers=names(aaa_data)

# 	for (sp in speakers){
# 		print(sp)
# 		for (plane in planes){
# 			# if ('Annotation_Label' %in% names(aaa_data[[sp]][[plane]]$tongue_traces)){
# 			tongue_data_wide <- aaa_data[[sp]][[plane]]$tongue_traces[,grepl('^[XY][0-9]',names(aaa_data[[sp]][[plane]]$tongue_traces))]
# 			xycols = ncol(tongue_data_wide)/2

# 			polar_origin = aaa_data[[sp]][[plane]]$origin
# 			if (flip) polar_origin[2] = -polar_origin[2]

# 			  for (i in 1:xycols){
# 				one_radius_raw = tongue_data_wide[,c(-1,0)+2*i]
# 				names(one_radius_raw) = c('X','Y') 
# 				# one_radius_rotated <- rotateXY(one_radius_raw, aaa_data[[sp]][[plane]]$occlusal_angle, center=c(0,0))

# 				one_radius_polar <- make.polar(one_radius_raw, origin=polar_origin, flip=flip)

# 				if (i==1){
# 					names(one_radius_polar) = paste0(c('T','R'),i)
# 					aaa_data[[sp]][[plane]]$tongue_traces_polar <- one_radius_polar
# 					# aaa_data[[sp]][[plane]]$tongue_traces_rotated <- one_radius_rotated
# 				}else{
# 					aaa_data[[sp]][[plane]]$tongue_traces_polar[,paste0(c('T','R'),i)] = one_radius_polar
# 					# aaa_data[[sp]][[plane]]$tongue_traces_rotated[,paste0(c('X','Y'),i)] = one_radius_rotated
# 				}			
# 			}
# 			aaa_data[[sp]][[plane]]$tongue_traces_polar[,paste0('T',1:xycols)] = -aaa_data[[sp]][[plane]]$tongue_traces_polar[,paste0('T',1:xycols)]
# 			# }
# 		}
# 	}		
# 	if (replaceXY){
# 	  aaa_data[[sp]][[plane]]$tongue_traces[,colnames(aaa_data[[sp]][[plane]]$tongue_traces_polar)] = aaa_data[[sp]][[plane]]$tongue_traces_polar
# 	  aaa_data[[sp]][[plane]]$tongue_traces_polar = NULL
# 	}
# 	aaa_data
# }


measure_traces_at_angles <- function(aaa_data, speakers=NULL, plotting=FALSE){

	if (is.null(speakers)) speakers=names(aaa_data)


	for (sp in speakers){
		print(sp)
		plane='sag'
		Ts = aaa_data[[sp]][[plane]]$tongue_traces[,grepl('^T[0-9]',names(aaa_data[[sp]][[plane]]$tongue_traces))]
		Rs = aaa_data[[sp]][[plane]]$tongue_traces[,grepl('^R[0-9]',names(aaa_data[[sp]][[plane]]$tongue_traces))]

		selected_angles = aaa_data[[sp]][[plane]]$analysis_value_radii[,c('TTangle','TDangle','TRangle')]*pi/180
		for (angle_name in names(selected_angles)){
			aaa_data[[sp]][[plane]]$tongue_traces[,angle_name] = NA
		}
		
		if (plotting){
			cairo_pdf(paste0('measure_traces_at_angles_',sp,'.pdf'), onefile=TRUE)
		}
		for (i in 1:nrow(aaa_data[[sp]][[plane]]$tongue_traces)){
			if (i/1000==round(i/1000)){
				print(paste(i,'of',nrow(aaa_data[[sp]][[plane]]$tongue_traces)))
			}
			one_trace_long = data.frame(T=as.numeric(Ts[i,]), R=as.numeric(Rs[i,]))
			selected_radii = approx(x=one_trace_long$T, y=one_trace_long$R, xout=selected_angles, rule=2)
			# print(one_trace_long)
			# print(selected_radii)
			if (plotting){
				plot(one_trace_long, type='o')		
				points(selected_radii, col='red')
			}
			aaa_data[[sp]][[plane]]$tongue_traces[i,names(selected_angles)] = selected_radii$y
		}
		if (plotting){
			dev.off()
		}
	}
	aaa_data
}

liparea <- function(lipdatarow){
	require(pracma)
	polyarea(lipdatarow[1:8], lipdatarow[9:16])
}

measure_lips <- function(aaa_data, speakers=NULL, plotting=FALSE){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		print(sp)
		plane='video'
		lipdata = aaa_data[[sp]][[plane]]$lip_traces
		lipdata$lips_horiz = with(lipdata, rightLip_x-leftLip_x)
		lipdata$lips_vert = with(lipdata, topmidinner_y-bottommidinner_y)
		lipdata$lips_purse = with(lipdata, lips_vert / lips_horiz)



		lipdata$lips_area = apply(lipdata[,c('leftLip_x', 'bottomleftinner_x', 'bottommidinner_x', 'bottomrightinner_x', 'rightLip_x', 'toprightinner_x', 'topmidinner_x', 'topleftinner_x',
			                                  'leftLip_y', 'bottomleftinner_y', 'bottommidinner_y', 'bottomrightinner_y', 'rightLip_y', 'toprightinner_y', 'topmidinner_y', 'topleftinner_y')], 
								  1, liparea)
		
		aaa_data[[sp]][[plane]]$lip_traces = lipdata
	}
	aaa_data
}



#### PLOTTING ####

plot_pal_occ <- function(aaa_data, planes=c('sag'), xlim=NULL, ylim=NULL, center=c(0,0)){

	pdf('palate_and_occlusal_traces_rotated.pdf', h=5,w=5*length(planes), onefile=T)

	par(mfrow=c(1,length(planes)))

	for (sp in speakers){
		for (plane in planes){
			if (is.null(aaa_data[[sp]][[plane]]$palate_trace)){
				plot(0,0, type='n', xlim=xlim, ylim=ylim, 
					main=paste(sp, plane, 'occlusal angle is unknown'))
			}else{
				all_XY = rbind(aaa_data[[sp]][[plane]]$occlusal_trace, aaa_data[[sp]][[plane]]$palate_trace, center)
				rotated_and_unrotated = rbind(all_XY, rotateXY(all_XY, aaa_data[[sp]][[plane]]$occlusal_angle, center=center))
				
				if (is.null(xlim) | is.null(ylim)){
					xlim = range(rotated_and_unrotated$X, na.rm=T)
					ylim = range(rotated_and_unrotated$Y, na.rm=T)
					if (diff(ylim) > diff(xlim)){
						xlim = xlim + c(-0.5,0.5)*(diff(ylim)-diff(xlim))
					}else{
						ylim = ylim + c(-0.5,0.5)*(diff(xlim)-diff(ylim))
					}
				}

				plot(aaa_data[[sp]][[plane]]$palate_trace, xlim=xlim, ylim=ylim, type='o',
					main=paste(sp, plane, '\npalate/occlusal traces\nocclusal angle is', round(aaa_data[[sp]][[plane]]$occlusal_angle,1), 'degrees'))
				if (sum(!is.na(aaa_data[[sp]][[plane]]$occlusal_trace))){
					points(aaa_data[[sp]][[plane]]$occlusal_trace)
					abline(coef(aaa_data[[sp]][[plane]]$occlusal_lm))
					points(aaa_data[[sp]][[plane]]$occlusal_trace_rotated, col='red')
					abline(coef(lm(Y~X,aaa_data[[sp]][[plane]]$occlusal_trace_rotated)), col='red')
				}
				points(aaa_data[[sp]][[plane]]$palate_trace_rotated, type='o', col='red')
				legend('bottomleft',pch=1,c('raw','rotated'),col=c('black','red'))
			}
		}
	}
	dev.off()
}



plot_traces <- function(aaa_data, speakers=NULL, planes=c('sag'), token=NULL, phone_column='phone', show_middle=FALSE, center=c(0,0)){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		for (plane in planes){
			if (length(aaa_data[[sp]][[plane]]$tongue_traces$token)>1){

				all_XY = data.frame(X=as.numeric(unlist(aaa_data[[sp]][[plane]]$tongue_traces[grepl('X',names(aaa_data[[sp]][[plane]]$tongue_traces))])),
					                Y=as.numeric(unlist(aaa_data[[sp]][[plane]]$tongue_traces[grepl('Y',names(aaa_data[[sp]][[plane]]$tongue_traces))])))
				all_XY = rbind(all_XY, aaa_data[[sp]][[plane]]$palate_trace)
				all_XY = all_XY[!is.na(all_XY$X),]
				rotated_convex_hull = rotateXY(all_XY[chull(all_XY),], aaa_data[[sp]][[plane]]$occlusal_angle, center=center)
				xlim = range(rotated_convex_hull$X)
				ylim = range(rotated_convex_hull$Y)
				if (diff(ylim) > diff(xlim)){
					xlim = xlim + c(-0.5,0.5)*(diff(ylim)-diff(xlim))
				}else{
					ylim = ylim + c(-0.5,0.5)*(diff(xlim)-diff(ylim))
				}
				aaa_data[[sp]][[plane]]$xlim <- xlim
				aaa_data[[sp]][[plane]]$ylim <- ylim

				all_tokens <- unique(aaa_data[[sp]][[plane]]$tongue_traces$token)
				tokens_sorted <- all_tokens[order(all_tokens)]

				if (is.null(token)){
					cairo_pdf(paste0('all_token_traces_',sp,'_',plane,'.pdf'), onefile=T)
				}else{
					tokens_sorted = token
				}

				for (token_label in tokens_sorted){
					print(token_label)

					plot(0, 0, type='n', xlab='X', ylab='Y', xlim=xlim, ylim=ylim, main=paste(sp, token_label))

					polygon(c(xlim+c(-5,5),rev(xlim+c(-5,5))),
						    c(rep(min(ylim)+0.02*diff(ylim),2),rep(max(ylim)+5,2)),
						    border=NA, col="#F3F3F3")

					mean_tongue_XY_rotated = rotateXY(aaa_data[[sp]][[plane]]$mean_tongue, aaa_data[[sp]][[plane]]$occlusal_angle, center=center)
					min_tongue_XY_rotated = rotateXY(aaa_data[[sp]][[plane]]$min_tongue, aaa_data[[sp]][[plane]]$occlusal_angle, center=center)
					max_tongue_XY_rotated = rotateXY(aaa_data[[sp]][[plane]]$max_tongue, aaa_data[[sp]][[plane]]$occlusal_angle, center=center)

					polygon(na.omit(c(min_tongue_XY_rotated$X,rev(max_tongue_XY_rotated$X))),
						    na.omit(c(min_tongue_XY_rotated$Y,rev(max_tongue_XY_rotated$Y))),
						    col='white',border=NA)
					box()

					# points(min_tongue_XY_rotated, type='l', lty=3, lwd=3, col='darkgray')
					points(mean_tongue_XY_rotated, type='l', lty=1, lwd=3, col='darkgray')
					# points(max_tongue_XY_rotated, type='l', lty=3, lwd=3, col='darkgray')

					if (plane=='sag'){
						polar_origin = aaa_data[[sp]][[plane]]$origin
						# print('ORIGIN!')
						# print(polar_origin)
						draw_angles = seq(-10,170,10)
						# points(polar_origin[1], polar_origin[2], cex=5)
						segments(rep(polar_origin[1],length(draw_angles)), rep(polar_origin[2],length(draw_angles)), 
					             polar_origin[1]-100*cos(draw_angles*pi/180), polar_origin[2]+100*sin(draw_angles*pi/180),
							     col='lightgray')
						if ('analysis_value_radii'%in%names(aaa_data[[sp]][[plane]])){
							selected_angles = as.numeric(aaa_data[[sp]][[plane]]$analysis_value_radii[,c('TTangle','TDangle','TRangle')])

							segments(rep(polar_origin[1],length(selected_angles)), rep(polar_origin[2],length(draw_angles)), 
						             polar_origin[1]-100*cos(selected_angles*pi/180), polar_origin[2]+100*sin(selected_angles*pi/180),
								     col='gray', lwd=2)
						# }else{
						# 	print('skipping analysis_value_radii because none were found')
						}

						points(polar_origin[1]-10*cos(draw_angles*pi/180), polar_origin[2]+10*sin(draw_angles*pi/180), pch=19,col='white',cex=2)
						text(polar_origin[1]-10*cos(draw_angles*pi/180), polar_origin[2]+10*sin(draw_angles*pi/180), labels=draw_angles, cex=0.5, col='gray')
					}

					first_frame = min(which(aaa_data[[sp]][[plane]]$tongue_traces$token==token_label))
					last_frame = max(which(aaa_data[[sp]][[plane]]$tongue_traces$token==token_label))
					n_frames = 1+last_frame-first_frame

					#draw the intensity
					tg = substr(aaa_data[[sp]][[plane]]$tongue_traces[first_frame,'token_id'],1,7)
					start_time = aaa_data[[sp]][[plane]]$tongue_traces[first_frame,'Time_of_sample_in_recording']
					end_time = aaa_data[[sp]][[plane]]$tongue_traces[last_frame,'Time_of_sample_in_recording']
					if (exists('all_intensities')){
						token_intensities = subset(all_intensities, speaker==sp & textgrid==tg & time>=start_time & time<=end_time)[,c('time','intensity')]
						token_intensities$time_scaled_for_x_axis = (token_intensities$time - min(token_intensities$time))/diff(range(token_intensities$time))*diff(xlim) + min(xlim)
						token_intensities$intensity_scaled_for_y_axis = (token_intensities$intensity - min(token_intensities$intensity))/diff(range(token_intensities$intensity))*0.12*diff(ylim) + min(ylim)+0.02*diff(ylim)
						polygon(c(token_intensities$time_scaled_for_x_axis, rev(token_intensities$time_scaled_for_x_axis)),
							    c(token_intensities$intensity_scaled_for_y_axis, rep(min(ylim)+0.02*diff(ylim), nrow(token_intensities))),
							    border=NA, col='gray')
					# }else{
					# 	print('skipping intensities because none were found')
					}
					points(aaa_data[[sp]][[plane]]$palate_trace_rotated, type='l', col='gray', lwd=3)

					intersecting_plane = c()
					for (i in first_frame:last_frame){
						one_tongue_trace_info = aaa_data[[sp]][[plane]]$tongue_traces[i,]

						raw_tongue_trace = data.frame(X=as.numeric(one_tongue_trace_info[grepl('X',names(one_tongue_trace_info))]),
							                          Y=as.numeric(one_tongue_trace_info[grepl('Y',names(one_tongue_trace_info))]),
							                          C=as.numeric(one_tongue_trace_info[grepl('^C[0-9]',names(one_tongue_trace_info))]))

						rotated_tongue_trace = rotateXY(raw_tongue_trace, aaa_data[[sp]][[plane]]$occlusal_angle, center=center)

						# where does the other plane intersect this one?
						if ('intersects_plane'%in%names(aaa_data[[sp]][[plane]])){
							if (!is.na(aaa_data[[sp]][[plane]]$intersects_plane)){
								intersecting_plane = rbind(intersecting_plane, rotated_tongue_trace[aaa_data[[sp]][[plane]]$intersects_plane,])
							}
						}
						if (show_middle){
							# print(one_tongue_trace_info$middle_frame)
							points(rotated_tongue_trace, type='l', lwd=ifelse(one_tongue_trace_info$middle_frame,3,0.5),
								col=rainbow(round(1.1*(n_frames)), v=0.85)[i-first_frame+1])
						}else{
							points(rotated_tongue_trace, type='l', col=rainbow(round(1.1*(n_frames)), v=0.85)[i-first_frame+1])
						}

						if (plane=='sag'){
							if (!is.na(aaa_data[[sp]][[plane]]$analysis_value_radii[,c('TT','TD','TR')][1,1])){
								points(rotated_tongue_trace[as.numeric(aaa_data[[sp]][[plane]]$analysis_value_radii[,c('TT','TD','TR')]),], pch=19, cex=0.5,
									col=rainbow(round(1.1*(n_frames)), v=0.85)[i-first_frame+1])
							}
						}else{
							if (!is.na(aaa_data[[sp]][[plane]]$analysis_value_radii[,c('LT','MT','RT')][1,1])){
								points(rotated_tongue_trace[as.numeric(aaa_data[[sp]][[plane]]$analysis_value_radii[,c('LT','MT','RT')]),], pch=19, cex=0.5,
									col=rainbow(round(1.1*(n_frames)), v=0.85)[i-first_frame+1])
							}
						}
					}
					# print(intersecting_plane)
					if (length(intersecting_plane)){
						intersecting_plane_extremes = rbind(intersecting_plane[intersecting_plane[,2]==min(intersecting_plane[,2]),],
															intersecting_plane[intersecting_plane[,2]==max(intersecting_plane[,2]),])
						points(intersecting_plane_extremes, type='l', lty=2)
					}
					text(seq(xlim[1], xlim[2], length.out=n_frames), 
						 rep(ylim[1],n_frames), 
						 aaa_data[[sp]][[plane]]$tongue_traces[first_frame:last_frame,phone_column], 
						 cex=0.45, col=rainbow(round(1.1*(n_frames)), v=0.85)[1:n_frames])
				}

				if (is.null(token)){
					dev.off()
				}
			}
		}
	}
}



plot_some_signals <- function(data=data, sp='', plane='sag', unique_label=NULL, signals_to_plot=c('TT','TD','TR'), ylim=c(0,100),
	time_warp_functions=NULL, all_phone_boundaries=NULL, all_phone_label_text=NULL, token.col='token', chunk.col='Annotation_Label'){

	full_name_of = list(sag='sagittal', cor='coronal')

	plot(0,0, type='n', ylim=ylim, xlim=c(0,1), 
		xlab='relative time', ylab='tongue distance (mm)', main=paste(sp, full_name_of[[plane]], unique_label))
	if (chunk.col %in% names(data[[sp]][[plane]]$tongue_traces)){
		matching_tokens <- unique(data[[sp]][[plane]]$tongue_traces[data[[sp]][[plane]]$tongue_traces[,chunk.col]==unique_label,][,token.col])
		for (j in 1:length(matching_tokens)){
			mt=matching_tokens[j]
			subdata <- data[[sp]][[plane]]$tongue_traces[data[[sp]][[plane]]$tongue_traces[,token.col]==mt,]
			if(nrow(subdata)){
				if (!is.null(time_warp_functions)){
					warp_time = time_warp_functions[[mt]]
					relative_time = warp_time(subdata$Time_of_sample_in_recording)
				}else{
					relative_time <- (subdata$Time_of_sample_in_recording - min(subdata$Time_of_sample_in_recording)) / diff(range(subdata$Time_of_sample_in_recording))
				}

				data[[sp]][[plane]]$tongue_traces[data[[sp]][[plane]]$tongue_traces[,token.col]==mt,'relative_time'] = relative_time

				is_new_phone = which((c(subdata$phone[1],subdata$phone) != c(subdata$phone,NA))[1:nrow(subdata)])
				phone_boundaries_x=warp_time(all_phone_boundaries[mt,])
				phone_boundaries_x = phone_boundaries_x[2:(length(phone_boundaries_x)-1)]
				phone_boundaries_ys = c()
				for (i in 1:3){
					signal <- signals_to_plot[i]
					# print (subdata[,signal])
					if (sum(!is.na(subdata[,signal]))>1){
						phone_boundaries_y = approx(subdata$Time_of_sample_in_recording, subdata[,signal], all_phone_boundaries[mt,], rule=2)$y
						phone_boundaries_y = phone_boundaries_y[2:(length(phone_boundaries_y)-1)]

						phone_boundaries_ys = rbind(phone_boundaries_ys, phone_boundaries_y)
						points(relative_time, as.numeric(subdata[,signal]), type='l', col=rainbow(3)[i])
						segments(phone_boundaries_x, phone_boundaries_y-diff(ylim)/100, 
							     phone_boundaries_x, phone_boundaries_y+diff(ylim)/100, col=rainbow(3)[i])
					}
				}
				if (j==1){
					phone_label_time = all_phone_boundaries[mt,][1:length(all_phone_boundaries[mt,])-1] + diff(all_phone_boundaries[mt,])/2
					text(warp_time(phone_label_time), rep(min(ylim), length(all_phone_label_text[[mt]])), all_phone_label_text[[mt]])
				}

				max_boundary_ys = apply(phone_boundaries_ys,2,max)
				min_boundary_ys = apply(phone_boundaries_ys,2,min)
				if (length(phone_boundaries_x) & length(phone_boundaries_y) & length(is_new_phone)){
					segments(phone_boundaries_x, rep(min(ylim), length(is_new_phone)), 
						     phone_boundaries_x, min_boundary_ys, col='gray', lwd=0.5)
					segments(phone_boundaries_x, min_boundary_ys, 
						     phone_boundaries_x, max_boundary_ys, col='gray', lwd=0.5)
				}
			}
		}
		legend('topleft', signals_to_plot, lty=1, col=rainbow(3))		
	}
	data[[sp]][[plane]]
}


plot_trajectories <- function(aaa_data, speakers=NULL, planes=c('sag'), coordinate_system='polar', normalization='average_phones', target=NULL, plot_formants=NULL, token.col='token', chunk.col='Annotation_Label'){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		print (sp)
		all_annotation_labels <- unique(c(aaa_data[[sp]][['sag']]$tongue_traces[,chunk.col]))         
		if (chunk.col %in% names(aaa_data[[sp]]$cor$tongue_traces)){
		  all_annotation_labels = unique(c(all_annotation_labels, aaa_data[[sp]][['cor']]$tongue_traces[,chunk.col]))
		}            
		annotation_labels_sorted <- all_annotation_labels[order(all_annotation_labels)]

		aaa_data[[sp]][['sag']]$tongue_traces$relative_time = NA
		aaa_data[[sp]][['cor']]$tongue_traces$relative_time = NA

		if (is.null(target)){
			cairo_pdf(paste0('all_token_trajectories_',sp,'.pdf'), onefile=T, h=10, w=5)
		}else{
			annotation_labels_sorted = target
		}
		# par(mfrow=c(2,1))

		if (is.null(plot_formants)){
			layout(mat = matrix(c(1,2,3), nrow = 3, ncol = 1),
			       heights = c(3, 3, 2),
			       widths = c(1))
		}else{
			layout(mat = matrix(c(1,2,3,4), nrow = 4, ncol = 1),
			       heights = c(3, 3, 2, 3),
			       widths = c(1))
		}

		for (unique_label in annotation_labels_sorted){
			print(unique_label)
			matching_rows = aaa_data[[sp]][['sag']]$tongue_traces[aaa_data[[sp]][['sag']]$tongue_traces[,chunk.col]==unique_label,][,c('token','token_id','Time_of_sample_in_recording','phone','word','word_id')]
			if ('Annotation_Label'%in%names(aaa_data[[sp]][['cor']]$tongue_traces)){
				matching_rows = rbind(matching_rows, aaa_data[[sp]][['cor']]$tongue_traces[aaa_data[[sp]][['cor']]$tongue_traces[,chunk.col]==unique_label,][,c('token','token_id','Time_of_sample_in_recording','phone','word','word_id')])
			}

			#investigate phone durations
			all_phone_boundaries = c()
			all_phone_durations = c()
			all_token_times = list()
			all_phone_label_text = list()
			time_warp_functions = list()
			matching_tokens = unique(matching_rows[,token.col])

			for (i in 1:length(matching_tokens)){	
				matching_token=matching_tokens[i]
				subdata = matching_rows[matching_rows[,token.col]==matching_token,]
				tg = substr(subdata[1,'token_id'],1,7)
				start_time = subdata[1,'Time_of_sample_in_recording']
				end_time = subdata[nrow(subdata),'Time_of_sample_in_recording']
				all_token_times[[matching_token]] = subdata$Time_of_sample_in_recording
				# print (subdata$phone)
				is_new_phone = which((c(subdata$phone[1],subdata$phone) != c(subdata$phone,NA))[1:nrow(subdata)])

				half_frame = mean(diff(subdata$Time_of_sample_in_recording))/2
				phone_boundaries = c(subdata$Time_of_sample_in_recording[1] - half_frame,
									 subdata$Time_of_sample_in_recording[is_new_phone] - half_frame, 
									 subdata$Time_of_sample_in_recording[nrow(subdata)] + half_frame)

				# find word boundaries
				is_new_word = which((c(subdata$word[1],subdata$word) != c(subdata$word,NA))[1:nrow(subdata)])

				# find consecutive words that are not sp
				add_sp_at = c()
				for (j in is_new_word){
					if (!'sp'%in%subdata$word[j-c(1,0)]){
						# need to add a zero-length sp
						# print (paste('add sp at frame', j))
						add_sp_at = c(add_sp_at, j)
					}
				}
				phone_label_text_at = sort(c(1, is_new_phone, add_sp_at))
				phone_label_text = subdata$phone[phone_label_text_at]
				phone_label_text[diff(phone_label_text_at)==0] = 'sp'
				all_phone_label_text[[matching_token]] = phone_label_text

				phone_boundaries = sort(c(phone_boundaries, subdata$Time_of_sample_in_recording[add_sp_at] - half_frame))

				# REPEAT THE LAST TIME POINT FOR ANY MISSING PHONES
				# THIS ONLY WORKS IF THE FIRST TOKEN IS NOT MISSING PHONES!
				if (!is.null(all_phone_boundaries)){
					if (length(phone_boundaries)<ncol(all_phone_boundaries)){
						phone_boundaries = c(phone_boundaries, rep(phone_boundaries[length(phone_boundaries)], ncol(all_phone_boundaries)-length(phone_boundaries)))
					}
				}
				all_phone_boundaries = rbind(all_phone_boundaries, phone_boundaries)
				phone_durations = diff(phone_boundaries)
				all_phone_durations = rbind(all_phone_durations, phone_durations)
			}
			rownames(all_phone_boundaries) = matching_tokens

			# print(all_phone_label_text)
			# print(all_phone_boundaries)

			mean_phone_durations = colMeans(all_phone_durations)
			mean_phone_times = cumsum(mean_phone_durations)
			scaled_mean_phone_times = c(0,mean_phone_times/mean_phone_times[length(mean_phone_times)])

			# make the time warping functions
			for (i in 1:length(matching_tokens)){
			
				matching_token = matching_tokens[i]	
				print(matching_token)
				frames_in_token = length(all_token_times[[matching_token]])
				phone_boundaries = all_phone_boundaries[matching_token,]

				# print (frames_in_token)
				# print (phone_boundaries)
				# print (phone_boundaries[c(1,length(phone_boundaries))])
				# print (0:1 + c(-1,1)*(((frames_in_token+1)/frames_in_token)-1))
				# warp_time_equal_phrase_old = approxfun(all_token_times[[i]][c(1,frames_in_token)], 0:1)
				warp_time_equal_phrase = approxfun(phone_boundaries[c(1,length(phone_boundaries))], 0:1 + c(-1,1)*(((frames_in_token+1)/frames_in_token)-1))
				
				# warp_time_equal_phones = approxfun(phone_boundaries, seq(0,1,length.out=length(phone_boundaries)))
				warp_time_average_phones = approxfun(phone_boundaries, scaled_mean_phone_times)
				if (normalization=='average_phones'){
					time_warp_functions[[matching_token]] = warp_time_average_phones
				}else{
					time_warp_functions[[matching_token]] = warp_time_equal_phrase
				}
			}
			# print('!')
			plane='sag'
			if (coordinate_system=='aaa'){
				if ('TT2' %in% names(aaa_data[[sp]][[plane]]$tongue_traces)){
					signals_to_plot <- c('TT2','TD2','TR2')
				}else{
					signals_to_plot <- c('TT','TD','TR')
				}
			}else{
				signals_to_plot <- c('TTangle','TDangle','TRangle')
			}
			# print(signals_to_plot)
			ylim = range(as.numeric(unlist(aaa_data[[sp]][[plane]]$tongue_traces[,signals_to_plot])), na.rm=T)
			ylim = ylim + c(0,0.3*diff(ylim))
			
			print(unique_label)
			sp_data = plot_some_signals(aaa_data, sp, 'sag', unique_label, signals_to_plot, ylim, time_warp_functions, all_phone_boundaries, all_phone_label_text, token.col=token.col)
			# add relative time
			aaa_data[[sp]][['sag']] = sp_data
			######

			if ('cor' %in% planes){
				plane='cor'
				signals_to_plot <- c('LT2','MT2','RT2')
				if (chunk.col %in% names(aaa_data[[sp]][[plane]]$tongue_traces)){
					ylim = range(as.numeric(unlist(aaa_data[[sp]][[plane]]$tongue_traces[,signals_to_plot])), na.rm=T)
					ylim = ylim + c(0,0.3*diff(ylim))
				}else{
					ylim=c(NA,NA)
				}
				sp_data = plot_some_signals(aaa_data, sp, 'cor', unique_label, signals_to_plot, ylim, time_warp_functions, all_phone_boundaries, all_phone_label_text)
				# add relative time
				aaa_data[[sp]][['cor']] = sp_data
			
			}else{
				plot(0,0,type='n',axes=F,xlab='',ylab='')
			}

			######

			if (exists('all_intensities')){
				plot(0,0,type='n', xlim=c(0,1), ylim=c(6,0), axes=F, xlab='relative time', ylab='', main='intensity by token')
				axis(1)
				box()
				for (i in 1:length(matching_tokens)){

					mt = matching_tokens[i]
					subdata = matching_rows[matching_rows[,token.col]==mt,]
					tg = substr(subdata[1,'token_id'],1,7)
					start_time = subdata[1,'Time_of_sample_in_recording']
					end_time = subdata[nrow(subdata),'Time_of_sample_in_recording']

					phone_boundaries = all_phone_boundaries[mt,]
					warp_time = time_warp_functions[[matching_tokens[i]]]
					relative_time = warp_time(subdata$Time_of_sample_in_recording)
					phone_boundaries_x = warp_time(phone_boundaries[2:(length(phone_boundaries)-1)])
					phone_label_time = phone_boundaries[1:length(phone_boundaries)-1] + diff(phone_boundaries)/2
					phone_label_text = all_phone_label_text[[mt]]
					token_intensities = subset(all_intensities, speaker==sp & textgrid==tg & time>=start_time & time<=end_time)[,c('time','intensity')]
					token_intensities$time_scaled_for_x_axis = warp_time(token_intensities$time)
					token_intensities$intensity_scaled_for_y_axis = (token_intensities$intensity - min(token_intensities$intensity))/diff(range(token_intensities$intensity))*0.9
				
					phone_starts_x = c(0,phone_boundaries_x)
					phone_ends_x = c(phone_boundaries_x,1)

					mean_scaled_durations = diff(scaled_mean_phone_times)
					current_scaled_durations = diff((phone_boundaries - min(phone_boundaries))/diff(range(phone_boundaries)))
					current_scaled_durations / mean_scaled_durations
					gray_levels = round(128+127*(log2(current_scaled_durations / mean_scaled_durations)))
					gray_levels[gray_levels>255] = 255
					gray_levels[gray_levels<10] = 10
					
					for (k in 1:length(phone_starts_x)){
						within_phone = token_intensities$time_scaled_for_x_axis>phone_starts_x[k] & token_intensities$time_scaled_for_x_axis<=phone_ends_x[k]
						sub_x = token_intensities$time_scaled_for_x_axis[within_phone]
						sub_y = token_intensities$intensity_scaled_for_y_axis[within_phone]

						polygon(c(sub_x, rev(sub_x)), i-c(sub_y, rep(0, length(sub_y))),
							    border=NA, col=gray.colors(256)[gray_levels[k]])	
					}

					text(phone_label_time, rep(i-0.25, length(phone_label_text)), labels=phone_label_text, cex=0.5)

					segments(phone_boundaries_x, rep(i,length(phone_boundaries_x)), 
						     phone_boundaries_x, rep(i-1,length(phone_boundaries_x)))
				}
				text(rep(0.05,length(matching_tokens)),-0.75+1:length(matching_tokens),matching_tokens, cex=0.5)
			}else{
				plot(0,0,type='n',axes=F,xlab='',ylab='')
			}

			if (!is.null(plot_formants)){

				plot(0,0,type='n', xlim=c(0,1), ylim=c(0,5000), xlab='relative time', ylab='formant frequency (Hz)', 
					main='formant frequency by token\n(questionable values: + frequency, x articulation, o bandwidth')

				for (i in 1:length(matching_tokens)){
					matching_token = matching_tokens[i]
					plane=ifelse(grepl('cor[0-9]',matching_token), 'cor','sag')
					# formant_rows = subset(all_formants, token_id%in%unique(subset(aaa_data[[sp]][[plane]]$tongue_traces, token==matching_token)$token_id))
					# formant_data = c()
					# for (row in 1:nrow(formant_rows)){
					# 	formants_long = long_format_formants(formant_rows[row,], c('F','B'), 1:3)
					# 	formant_data = rbind(formant_data, formants_long)
					# 	formant_data = rbind(formant_data, NA)
					# }
					  
					formant_data = subset(aaa_data[[sp]][[plane]]$tongue_traces, token==matching_token)
					warp_time = time_warp_functions[[matching_token]]
					formant_data$tnorm = warp_time(formant_data$Time_of_sample_in_recording)
					is_vowel_formant = formant_data$formant_time>0 & formant_data$formant_time<1
					if (plot_formants=='all'){					
						points(formant_data$tnorm, formant_data$F1, type='l', lty=2, col=rgb(1,0,0))
						points(formant_data$tnorm, formant_data$F2, type='l', lty=2, col=rgb(0,1,0))
						points(formant_data$tnorm, formant_data$F3, type='l', lty=2, col=rgb(0,0,1))
						vowel_formant_data_to_plot = formant_data
					}else{
						vowel_formant_data_to_plot = formant_data[is_vowel_formant,]
					}
					points(formant_data[is_vowel_formant,]$tnorm, formant_data[is_vowel_formant,]$F1, type='l', col=rgb(1,0,0))
					points(formant_data[is_vowel_formant,]$tnorm, formant_data[is_vowel_formant,]$F2, type='l', col=rgb(0,1,0))
					points(formant_data[is_vowel_formant,]$tnorm, formant_data[is_vowel_formant,]$F3, type='l', col=rgb(0,0,1))

					points(vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F1_freq_res_sd)>2,]$tnorm, 
						   vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F1_freq_res_sd)>2,]$F1, pch=3)
					points(vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F2_freq_res_sd)>2,]$tnorm, 
						   vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F2_freq_res_sd)>2,]$F2, pch=3)
					points(vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F3_freq_res_sd)>2,]$tnorm, 
						   vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F3_freq_res_sd)>2,]$F3, pch=3)

					if(plane=='sag'){
						points(vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F1_res_sd)>2,]$tnorm, 
							   vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F1_res_sd)>2,]$F1, pch=4)
						points(vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F2_res_sd)>2,]$tnorm, 
							   vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F2_res_sd)>2,]$F2, pch=4)
						points(vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F3_res_sd)>2,]$tnorm, 
							   vowel_formant_data_to_plot[abs(vowel_formant_data_to_plot$F3_res_sd)>2,]$F3, pch=4)
					}
					points(vowel_formant_data_to_plot[vowel_formant_data_to_plot$logB1_sd>2,]$tnorm, 
						   vowel_formant_data_to_plot[vowel_formant_data_to_plot$logB1_sd>2,]$F1, pch=1)
					points(vowel_formant_data_to_plot[vowel_formant_data_to_plot$logB2_sd>2,]$tnorm, 
						   vowel_formant_data_to_plot[vowel_formant_data_to_plot$logB2_sd>2,]$F2, pch=1)
					points(vowel_formant_data_to_plot[vowel_formant_data_to_plot$logB3_sd>2,]$tnorm, 
						   vowel_formant_data_to_plot[vowel_formant_data_to_plot$logB3_sd>2,]$F3, pch=1)
				}

				# legend('topright',legend=c('surprising frequency','surprising given articulation','large bandwidth'),pch=c(3,4,1))
			}


		}
		if (is.null(target)){
			dev.off()
		}
	}
	aaa_data
}

# "F1_res_sd"                  
# [197] "F2_res_sd"                   "F3_res_sd"                  
# [199] "logB1_sd"                    "logB2_sd"                   
# [201] "logB3_sd"                    "badness1"                   
# [203] "badness1_sd"                 "badness2"                   
# [205] "badness2_sd"                 "formant_phone"              
# [207] "formant_time"                "F1_freq_res"                
# [209] "F2_freq_res"                 "F3_freq_res"                
# [211] "F1_freq_res_sd"              "F2_freq_res_sd"             
# [213] "F3_freq_res_sd"      

set_basic_ssanova_targets <- function(){

	targets <- list(after_IY1=list(list(Annotation_Label=c("peep"), phone=c('IY1')),
					               list(Annotation_Label=c("peel"), phone=c('L')),
					               list(Annotation_Label=c("teeth"), phone=c('TH')),
					               list(Annotation_Label=c("beak"), phone=c('K'))))
	targets$after_IH1=list(list(Annotation_Label=c("pip"), phone=c('IH1')),
				                   list(Annotation_Label=c("pill"), phone=c('L')),
								   list(Annotation_Label=c("myth"), phone=c('TH')))		                                  
	targets$after_EH1=list(list(Annotation_Label=c("beb"), phone=c('EH1')),
				                   list(Annotation_Label=c("bell"), phone=c('L')),
				                   list(Annotation_Label=c("beth"), phone=c('TH')))		                                  
	targets$after_EY1=list(list(Annotation_Label=c("babe"), phone=c('EY1')),
				                   list(Annotation_Label=c("bale"), phone=c('L')),
				                   list(Annotation_Label=c("faith"), phone=c('TH')))		                                  
	targets$after_UW1=list(list(Annotation_Label=c("poop"), phone=c('UW1')),
				                   list(Annotation_Label=c("pool"), phone=c('L')),
				                   list(Annotation_Label=c("booth"), phone=c('TH')))		                                  
	targets$after_AO1=list(list(Annotation_Label=c("saw_magic"), phone=c('AO1')),
				                   list(Annotation_Label=c("fall_again"), phone=c('L')),
				                   list(Annotation_Label=c("hawk"), phone=c('K')))    	                                  
	targets$after_OW1=list(list(Annotation_Label=c("bo_again"), phone=c('OW1')),
				                   list(Annotation_Label=c("pole"), phone=c('L')),
				                   list(Annotation_Label=c("both"), phone=c('TH')))                                  
	targets$in_CALL=list(list(Annotation_Label=c("call_someone"), phone=c('L')),
				                 list(Annotation_Label=c("call_someone"), phone=c('AO1')),
				                 list(Annotation_Label=c("call_someone"), phone=c('K')))
	targets$in_HAL_KEPT=list(list(Annotation_Label=c("hal_kept"), phone=c('L')),
				                     list(Annotation_Label=c("hal_kept"), phone=c('AE1')),
				                     list(Annotation_Label=c("hal_kept"), phone=c('EH1')),
				                     list(Annotation_Label=c("hal_kept"), phone=c('K')))			                                  
	targets$in_KIM=list(list(Annotation_Label=c("saw_kim"), phone=c('IH1')),
				                list(Annotation_Label=c("saw_kim"), phone=c('M')),
				                list(Annotation_Label=c("saw_kim"), phone=c('K')))
	targets$before_EH1=list(list(Annotation_Label=c("hefty"), phone=c('EH1')),
				                    list(Annotation_Label=c("left"), phone=c('L')),
				                    list(Annotation_Label=c("theft"), phone=c('TH')),
				                    list(Annotation_Label=c("wept"), phone=c('W')))		                                  
	targets$AO1_EH1=list(list(Annotation_Label=c("saw_magic"), phone=c('AO1')),
				                list(Annotation_Label=c("see_magic"), phone=c('IY1')),
				                list(Annotation_Label=c("beb"), phone=c('EH1')))		                                  
	targets$in_BO_AGAIN=list(list(Annotation_Label=c("bo_again"), phone=c('EH1')),
				                list(Annotation_Label=c("bo_again"), phone=c('G')),
				                list(Annotation_Label=c("bo_again"), phone=c('AH0')))			                                  
	targets$AA1_IY1=list(list(Annotation_Label=c("mom"), phone=c('AA1')),
				                list(Annotation_Label=c("peep"), phone=c('IY1')),
				                list(Annotation_Label=c("beb"), phone=c('EH1')))			                                  
	targets$UW1_IY1=list(list(Annotation_Label=c("poop"), phone=c('UW1')),
				                list(Annotation_Label=c("peep"), phone=c('IY1')),
				                list(Annotation_Label=c("beb"), phone=c('EH1')))			                                  
	targets$all_Ls=list(list(Annotation_Label=c("peel"), phone=c('L')),
				                list(Annotation_Label=c("pill"), phone=c('L')),
				                list(Annotation_Label=c("bale"), phone=c('L')),
				                list(Annotation_Label=c("bell"), phone=c('L')),
				                list(Annotation_Label=c("pool"), phone=c('L')),
				                list(Annotation_Label=c("pole"), phone=c('L')))			                                  
	targets$all_Vs=list(list(Annotation_Label=c("peep"), phone=c('IY1')),
				                list(Annotation_Label=c("pip"), phone=c('IH1')),
				                list(Annotation_Label=c("babe"), phone=c('EY1')),
				                list(Annotation_Label=c("beb"), phone=c('EH1')),
				                list(Annotation_Label=c("poop"), phone=c('UW1')),
				                list(Annotation_Label=c("bo again"), phone=c('OW1')))			                                  
	targets$all_THs=list(list(Annotation_Label=c("teeth"), phone=c('TH')),
				                list(Annotation_Label=c("myth"), phone=c('TH')),
				                list(Annotation_Label=c("faith"), phone=c('TH')),
				                list(Annotation_Label=c("beth"), phone=c('TH')),
				                list(Annotation_Label=c("booth"), phone=c('TH')),
				                list(Annotation_Label=c("both"), phone=c('TH')))
	targets
}

set_basic_comparison_levels <- function(){

	diff_levels_of  <- list(after_IY1=list(c('IY1','L'), c('IY1','TH')),
		                    after_IH1=list(c('IH1','L'), c('IH1','TH')),
							after_EH1=list(c('EH1','L'), c('EH1','TH')),
							after_EY1=list(c('EY1','L'), c('EY1','TH')),
							after_UW1=list(c('UW1','L'), c('UW1','TH')),
							after_OW1=list(c('OW1','L'), c('OW1','TH')),
							after_AO1=list(c('AO1','L'), c('AO1','K')),
							in_CALL=list(c('AO1','L'), c('AO1','K')),
							in_KIM=list(c('IH1','M'), c('IH1','K')),
							in_HAL_KEPT=list(c('AE1','L'), c('EH1','K')),
							before_EH1=list(c('EH1','L'), c('EH1','W')),
							AO1_EH1=list(c('EH1','AO1'), c('EH1','IY1')),
							in_BO_AGAIN=list(c('EH1','G'), c('AH0','G')),
							AA1_IY1=list(c('EH1','AA1'), c('EH1','IY1')),
							UW1_IY1=list(c('EH1','UW1'), c('EH1','IY1')))
	diff_levels_of	                                  
}


do_cl_ssanova_comparisons <- function(aaa_data, speakers=NULL, aaa_ssanovas=list(), alpha=1.4, SDs=2.58){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		print(sp)
		aaa_ssanovas[[sp]] = list()
		plane = 'sag'
		aaa_ssanovas[[sp]][[plane]] <- list()

		tongue_data_long <- aaa_data[[sp]][[plane]]$tongue_traces_long_rotated
		polar_origin <- aaa_data[[sp]][[plane]]$origin
		xlim=aaa_data[[sp]][[plane]]$xlim
		ylim=aaa_data[[sp]][[plane]]$ylim

		# aaa_ssanovas[[sp]][['sag']][['all_phones']] <- polar.ssanova(tongue_data_long, data.cat='phone', 
		#   token.label='token_id', main=paste(sp, 'phones', collapse=' '), crop=TRUE, xlim=xlim, ylim=ylim, CI.fill=TRUE, 
		#   origin=polar_origin, flip=FALSE, legend.pos='bottomleft', SDs=2.58, drop_levels=TRUE)

		# if('Annotation_Label'%in%names(aaa_data[[sp]][[plane]]$tongue_traces)){
		# 	aaa_ssanovas[[sp]][['cor']] <- list()
		# 	aaa_ssanovas[[sp]][['cor']][['all_phones']] <- polar.ssanova(aaa_data[[sp]][['cor']]$tongue_traces_long_rotated, data.cat='phone', 
		# 	  token.label='token_id', main=paste(sp, 'phones', collapse=' '), crop=TRUE, xlim=xlim, ylim=ylim, CI.fill=TRUE, 
		# 	  origin=polar_origin, flip=FALSE, legend.pos='bottomleft', SDs=2.58, drop_levels=TRUE)
		# }

		for (one_target_set_name in setdiff(names(targets), names(aaa_ssanovas[[sp]][[plane]]))){
			one_target_set = targets[[one_target_set_name]]
			print(c(sp,plane,one_target_set_name))
			tongue_data_targets <- c()
			for (item in one_target_set){
				# subset(tongue_data_long, Annotation_Label==item$Annotation_Label & phone==item$phone)

				more_targets <- tongue_data_long
				for (item_column in names(item)){
					more_targets <- more_targets[more_targets[,item_column]%in%item[[item_column]],]
				}
			  	tongue_data_targets <- rbind(tongue_data_targets, more_targets)
			}

			# dc=ifelse(one_target_set_name%in%c('all_Ls','all_THs'), 'Annotation_Label', 'phone')
			dc = names(one_target_set[[1]])[1]
			print(dc)
			tongue_data_targets[,dc] = factor(tongue_data_targets[,dc])
			print(unique(tongue_data_targets[,dc]))
			# print (head(tongue_data_targets))
			aaa_ssanovas[[sp]][[plane]][[one_target_set_name]] <- polar.ssanova(tongue_data_targets, data.cat=dc, 
			  token.label='token_id', main=paste(sp, one_target_set_name, collapse=' '), crop=TRUE, xlim=xlim, ylim=ylim, CI.fill=TRUE, 
			  origin=polar_origin, alpha=alpha, flip=FALSE, legend.pos='bottomleft', SDs=SDs, drop_levels=TRUE)
		}
	}
	aaa_ssanovas
}


get_target_mean_signals <- function(aaa_ssanovas, targets, speakers=NULL){

	require(plyr)

	if (is.null(speakers)) speakers = names(aaa_ssanovas)

	targets_means_by_plane = list(sag=c(),cor=c())

	# get trajectory values at middle frames
	for (sp in speakers){
		# aaa_ssanovas[[sp]] = list()
		for (plane in c('sag','cor')){
			print(c(sp,plane))
			# aaa_ssanovas[[sp]][[plane]] <- list()


			if ('middle_frame' %in% names(aaa_data[[sp]][[plane]]$tongue_traces)){
				middle_frames = subset(aaa_data[[sp]][[plane]]$tongue_traces, middle_frame==TRUE)

				for (one_target_set_name in names(targets)){
					one_target_set = targets[[one_target_set_name]]
					tongue_data_targets <- c()
					for (item in one_target_set){
						# subset(tongue_data_long, Annotation_Label==item$Annotation_Label & phone==item$phone)
						more_targets <- middle_frames
						for (item_column in names(item)){
							more_targets <- more_targets[more_targets[,item_column]%in%item[[item_column]],]
						}
					  	tongue_data_targets <- rbind(tongue_data_targets, more_targets)
					}
					if (plane=='sag'){
						targets_for_ddply=data.frame(phone=tongue_data_targets$phone, word=tongue_data_targets$word,
							TR=as.numeric(tongue_data_targets$TR),
							TD=as.numeric(tongue_data_targets$TD),
							TT=as.numeric(tongue_data_targets$TT))
						targets_means = ddply(targets_for_ddply, .(phone, word), summarize, TR=mean(TR), TD=mean(TD), TT=mean(TT))
					}else{
						targets_for_ddply=data.frame(phone=tongue_data_targets$phone, word=tongue_data_targets$word,
							LT=as.numeric(tongue_data_targets$LT),
							MT=as.numeric(tongue_data_targets$MT),
							RT=as.numeric(tongue_data_targets$RT))
						targets_means = ddply(targets_for_ddply, .(phone, word), summarize, LT=mean(LT), MT=mean(MT), RT=mean(RT))
					}

					if (nrow(targets_means)){
						targets_means=data.frame(speaker = sp, model_name=one_target_set_name, targets_means)
					}else{
						targets_means=c()
					}
					targets_means_by_plane[[plane]] = rbind(targets_means_by_plane[[plane]], targets_means)

				}
			}
		}
	}
	if (is.null(targets_means_by_plane[['cor']])){
		target_mean_signals = targets_means_by_plane[['sag']]
	}else{
		target_mean_signals = merge(targets_means_by_plane[['sag']], targets_means_by_plane[['cor']], all=T)
	}
	target_mean_signals
}

make_ssanovas_plots_and_calculate_displacements <- function(aaa_ssanova, targets){
	make_ssanovas_plots(aaa_ssanova, targets, calculate_displacements=TRUE)
}

make_ssanovas_plots <- function(aaa_ssanova, targets, calculate_displacements=FALSE, legend.pos='bottomright'){

	if (calculate_displacements){
		pdf('ssanova_figures.pdf', h=6, w=20, onefile=TRUE)
		par(mfrow=c(1,4))
		displacements <- c()
	}else{
		pdf('ssanova_figures.pdf', h=6, w=10, onefile=TRUE)
		par(mfrow=c(1,2))		
	}

	for (one_target_set_name in names(targets)){
		print (one_target_set_name)
		if (is.null(speakers)) speakers=names(aaa_ssanovas)

		for (sp in speakers){
			plane = 'sag'
			main = paste(sp, one_target_set_name, collapse=' ')
			# dc=ifelse(one_target_set_name%in%c('all_Ls','all_THs'), 'Annotation_Label', 'phone')
			dc = names(targets[[one_target_set_name]][[1]])[1]
			replot.tongue.ss(aaa_ssanovas[[sp]][[plane]][[one_target_set_name]], CI.fill=TRUE, legend.pos=legend.pos)
			show.traces(aaa_ssanovas[[sp]][[plane]][[one_target_set_name]]$data, data.cat=dc, 
		  	token.label='token_id', main=main, xlim=aaa_data[[sp]][[plane]]$xlim, ylim=aaa_data[[sp]][[plane]]$ylim, origin=polar_origin, flip=FALSE, legend.pos=legend.pos, drop_levels=TRUE)

			if (calculate_displacements){
				if (one_target_set_name %in% names(diff_levels_of)){
					for (diff_levels in diff_levels_of[[one_target_set_name]]){
						# diff_levels = c('EH1', 'L')
						if (length(intersect(diff_levels, unique(aaa_ssanovas[[sp]]$sag[[one_target_set_name]]$pol.cart$phone)))==length(diff_levels)){
							diff_return = difference_plot(aaa_ssanovas[[sp]]$sag[[one_target_set_name]], diff_levels[1], diff_levels[2], main=main, flip=FALSE)
							new_displacement = data.frame(speaker=sp, plane=plane, model_name=one_target_set_name,
														  phone1=diff_levels[1], phone2=diff_levels[2],
														  total_area=diff_return$total_area, 
														  displacement_angle=diff_return$displacement_angle, 
														  displacement_distance=diff_return$displacement_distance,
												          contour_displacement_angle=diff_return$contour_displacement_angle,
												          contour_displacement_magnitude=diff_return$contour_displacement_magnitude,
												          polygon_displacement_angle=diff_return$polygon_displacement_angle,
												          polygon_displacement_magnitude=diff_return$polygon_displacement_magnitude)

							displacements = rbind(displacements, new_displacement)

						}else{
							plot(0,0,type='n',main=sp)
						}
					}
				}else{
					plot(0,0,type='n',main=sp)
				}	
			}
		}
	}

	dev.off()
	if (calculate_displacements){
		displacements
	}
}


plot_displacement_measures <- function(displacements, targets, plane='sag'){
	pdf('displacement_measures.pdf', h=8, w=12, onefile=TRUE)

	par(mfrow=c(2,3))

	for (one_target_set_name in names(targets)){
		for (diff_levels in diff_levels_of[[one_target_set_name]]){

			subdata <- subset(displacements, model_name==one_target_set_name & phone2==diff_levels[2])

			if (nrow(subdata)){
				plot(subdata$total_area, subdata$displacement_angle, type='n', xlab='magnitude', ylab='angle', 
					main = paste(one_target_set_name, plane, diff_levels[2], 'vs.', diff_levels[1], collapse=' '),
					ylim=c(-pi,pi), xlim=range(displacements$total_area))
				text(subdata$total_area, subdata$displacement_angle, subdata$sp, col='red')

				plot(subdata$contour_displacement_magnitude, subdata$contour_displacement_angle, type='n', xlab='magnitude (contour)', ylab='angle (contour)', 
					main = paste(one_target_set_name, plane, diff_levels[2], 'vs.', diff_levels[1], collapse=' '),
					ylim=c(-pi,pi), xlim=range(displacements$contour_displacement_magnitude))
				text(subdata$contour_displacement_magnitude, subdata$contour_displacement_angle, subdata$sp, col='black')
			
				plot(subdata$polygon_displacement_magnitude, subdata$polygon_displacement_angle, type='n', xlab='magnitude (polygon)', ylab='angle (polygon)', 
					main = paste(one_target_set_name, plane, diff_levels[2], 'vs.', diff_levels[1], collapse=' '),
					ylim=c(-pi,pi), xlim=range(displacements$polygon_displacement_magnitude))
				text(subdata$polygon_displacement_magnitude, subdata$polygon_displacement_angle, subdata$sp, col='blue')
			}
		}
	}
	dev.off()
}


plot_analysis_values <- function(target_mean_signals, targets){

	pdf('aaa_analysis_values.pdf', h=8, w=12, onefile=TRUE)

	par(mfrow=c(2,3))

	for (one_target_set_name in names(targets)){
		for (diff_levels in diff_levels_of[[one_target_set_name]]){

			subdata1 <- subset(target_mean_signals, model_name==one_target_set_name & phone==diff_levels[1])
			subdata2 <- subset(target_mean_signals, model_name==one_target_set_name & phone==diff_levels[2])

			names(subdata1)[5:10] = paste0(names(subdata1)[5:10], '1')
			names(subdata2)[5:10] = paste0(names(subdata2)[5:10], '2')

			subdata_both = merge(subdata1[,c(1,5:10)], subdata2[,c(1,5:10)], all=T)
			subdata_both$TR_diff = subdata_both$TR2 - subdata_both$TR1
			subdata_both$TD_diff = subdata_both$TD2 - subdata_both$TD1
			subdata_both$TT_diff = subdata_both$TT2 - subdata_both$TT1
			subdata_both$LT_diff = subdata_both$LT2 - subdata_both$LT1
			subdata_both$MT_diff = subdata_both$MT2 - subdata_both$MT1
			subdata_both$RT_diff = subdata_both$RT2 - subdata_both$RT1

			subdata_both$latT1 = (subdata_both$LT1 + subdata_both$RT1)/2
			subdata_both$latT2 = (subdata_both$LT2 + subdata_both$RT2)/2
			subdata_both$latT_diff = subdata_both$latT2 - subdata_both$latT1

			# xylim=range(subdata_both[,grep('diff',names(subdata_both))], na.rm=T)
			xylim=c(-20,20)

			if (nrow(subdata_both)){
				plot(subdata_both$TR_diff, subdata_both$TD_diff, type='n', xlab='TR difference', ylab='TD difference', 
					main = paste(one_target_set_name, '', diff_levels[2], 'vs.', diff_levels[1], collapse=' '),
					xlim=xylim, ylim=xylim)
				abline(0,1, col='gray') ; abline(v=0, col='gray') ; abline(h=0, col='gray')
				text(subdata_both$TR_diff, subdata_both$TD_diff, subdata_both$sp)

				plot(subdata_both$TT_diff, subdata_both$TD_diff, type='n', xlab='TT difference', ylab='TD difference', 
					main = paste(one_target_set_name, '', diff_levels[2], 'vs.', diff_levels[1], collapse=' '),
					xlim=xylim, ylim=xylim)
				abline(0,1, col='gray') ; abline(v=0, col='gray') ; abline(h=0, col='gray')
				text(subdata_both$TT_diff, subdata_both$TD_diff, subdata_both$sp)

				plot(subdata_both$MT_diff, subdata_both$latT_diff, type='n', xlab='MT difference', ylab='mean(LT,RT) difference', 
					main = paste(one_target_set_name, '', diff_levels[2], 'vs.', diff_levels[1], collapse=' '),
					xlim=xylim, ylim=xylim)
				abline(0,1, col='gray') ; abline(v=0, col='gray') ; abline(h=0, col='gray')
				text(subdata_both$MT_diff, subdata_both$latT_diff, subdata_both$sp)

			}
		}
	}
	dev.off()
}


add_interpolated_formants <- function(aaa_data, all_formants, speakers=NULL){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		print(sp)
		speaker_formants = subset(all_formants, speaker==sp)
		for (param in c('formant_phone','F1','F2','F3','B1','B2','B3')){
			for (plane in c('sag','cor')){
				if('Annotation_Label'%in%names(aaa_data[[sp]][[plane]]$tongue_traces)){
					aaa_data[[sp]][[plane]]$tongue_traces[,param] = NA
					# print (paste('added column',param))
				}
			}
		}

		for (tg in unique(speaker_formants$textgrid)){
			print(tg)
			textgrid_formants = subset(speaker_formants, textgrid==tg)

			formant_data = c()
			formant_phones = c()
			for (row in 1:nrow(textgrid_formants)){
				formants_long = long_format_formants(textgrid_formants[row,], c('F','B'), 1:3)
				formant_data = rbind(formant_data, formants_long)
				# formant_data = rbind(formant_data, data.frame(formants_long, formant_phone=textgrid_formants[row,'phone']))
				formant_phones = rbind(formant_phones, data.frame(formant_phone=textgrid_formants[row,'phone'], 
					start=formants_long[1,'t'], end=formants_long[nrow(formants_long),'t'],
					phonestart=textgrid_formants[row,'phonestart'], phoneend=textgrid_formants[row,'phoneend']))
				if (row < nrow(textgrid_formants)) formant_data = rbind(formant_data, NA)
			}
			formant_data[which(is.na(formant_data$t)),'t'] = (formant_data[which(is.na(formant_data$t))-1,'t'] + formant_data[which(is.na(formant_data$t))+1,'t'])/2
			# print('*')
			# print (formant_phones)
			for (plane in c('sag','cor')){
				if('Annotation_Label'%in%names(aaa_data[[sp]][[plane]]$tongue_traces)){
					rows_needing_formants = which(grepl(tg, aaa_data[[sp]][[plane]]$tongue_traces$token_id))
					# print (length(rows_needing_formants))
					if (length(rows_needing_formants)){
						formant_times = aaa_data[[sp]][[plane]]$tongue_traces[rows_needing_formants,'Time_of_sample_in_recording']
						# aaa_data[[sp]][[plane]]$tongue_traces[rows_needing_formants,'formant_phone'] = formant_data$formant_phone
						for (param in c('F1','F2','F3','B1','B2','B3')){
							aaa_data[[sp]][[plane]]$tongue_traces[rows_needing_formants,param] = approx(formant_data$t,formant_data[,param],formant_times, na.rm=FALSE)$y
							# print (paste('added data for',param))
						}
						for (row in 1:nrow(formant_phones)){
							rows_this_phone = intersect(rows_needing_formants, which(aaa_data[[sp]][[plane]]$tongue_traces$Time_of_sample_in_recording >= formant_phones[row,'start'] & aaa_data[[sp]][[plane]]$tongue_traces$Time_of_sample_in_recording <= formant_phones[row,'end']))
							aaa_data[[sp]][[plane]]$tongue_traces[rows_this_phone,'formant_phone'] = formant_phones[row,'formant_phone']
							aaa_data[[sp]][[plane]]$tongue_traces[rows_this_phone,'formant_time'] = (aaa_data[[sp]][[plane]]$tongue_traces[rows_this_phone,'Time_of_sample_in_recording']-formant_phones[row,'phonestart'])/(formant_phones[row,'phoneend']-formant_phones[row,'phonestart'])
						
						}
					}
				}
			}
		}
	}
	aaa_data
}




measure_formant_badness <- function(aaa_data, speakers=NULL){

	if (is.null(speakers)) speakers=names(aaa_data)

	require(lme4)

	for (sp in speakers){
		print(sp)
		# get residuals for sagittal tokens by tongue measurements
		is_vowel_sag = grepl('[012]',aaa_data[[sp]]$sag$tongue_traces$phone)
		F1_sag_lmer = lmer(F1~TRangle+TDangle+TTangle+(1|phone), aaa_data[[sp]]$sag$tongue_traces)
		F2_sag_lmer = lmer(F2~TRangle+TDangle+TTangle+(1|phone), aaa_data[[sp]]$sag$tongue_traces)
		F3_sag_lmer = lmer(F3~TRangle+TDangle+TTangle+(1|phone), aaa_data[[sp]]$sag$tongue_traces)

		sag_residuals = data.frame(aaa_data[[sp]]$sag$tongue_traces[,c('phone','F1','F2','F3')], F1_res=NA, F2_res=NA, F3_res=NA)
		sag_residuals[!is.na(sag_residuals$F1),]$F1_res = residuals(F1_sag_lmer)
		sag_residuals[!is.na(sag_residuals$F2),]$F2_res = residuals(F2_sag_lmer)
		sag_residuals[!is.na(sag_residuals$F3),]$F3_res = residuals(F3_sag_lmer)

		aaa_data[[sp]]$sag$tongue_traces$F1_res_sd = (sag_residuals$F1_res-mean(sag_residuals[is_vowel_sag,'F1_res'], na.rm=T))/sd(sag_residuals[is_vowel_sag,'F1_res'], na.rm=T)
		aaa_data[[sp]]$sag$tongue_traces$F2_res_sd = (sag_residuals$F2_res-mean(sag_residuals[is_vowel_sag,'F2_res'], na.rm=T))/sd(sag_residuals[is_vowel_sag,'F2_res'], na.rm=T)
		aaa_data[[sp]]$sag$tongue_traces$F3_res_sd = (sag_residuals$F3_res-mean(sag_residuals[is_vowel_sag,'F3_res'], na.rm=T))/sd(sag_residuals[is_vowel_sag,'F3_res'], na.rm=T)

		if ('Annotation_Label'%in%names(aaa_data[[sp]]$cor$tongue_traces)){

			# put all bandwidth data (from sag and cor tokens) together and calculate sd all at once
			bandwidth_data = rbind(aaa_data[[sp]]$sag$tongue_traces[,c('phone','B1','B2','B3')], 
				                   aaa_data[[sp]]$cor$tongue_traces[,c('phone','B1','B2','B3')])
			is_vowel_bandwidth = grepl('[012]',bandwidth_data$phone)
			aaa_data[[sp]]$sag$tongue_traces$logB1_sd = (log(aaa_data[[sp]]$sag$tongue_traces$B1)-mean(log(bandwidth_data$B1[is_vowel_bandwidth]), na.rm=T))/sd(log(bandwidth_data$B1[is_vowel_bandwidth]), na.rm=T)
			aaa_data[[sp]]$sag$tongue_traces$logB2_sd = (log(aaa_data[[sp]]$sag$tongue_traces$B2)-mean(log(bandwidth_data$B2[is_vowel_bandwidth]), na.rm=T))/sd(log(bandwidth_data$B2[is_vowel_bandwidth]), na.rm=T)
			aaa_data[[sp]]$sag$tongue_traces$logB3_sd = (log(aaa_data[[sp]]$sag$tongue_traces$B3)-mean(log(bandwidth_data$B3[is_vowel_bandwidth]), na.rm=T))/sd(log(bandwidth_data$B3[is_vowel_bandwidth]), na.rm=T)
			aaa_data[[sp]]$cor$tongue_traces$logB1_sd = (log(aaa_data[[sp]]$cor$tongue_traces$B1)-mean(log(bandwidth_data$B1[is_vowel_bandwidth]), na.rm=T))/sd(log(bandwidth_data$B1[is_vowel_bandwidth]), na.rm=T)
			aaa_data[[sp]]$cor$tongue_traces$logB2_sd = (log(aaa_data[[sp]]$cor$tongue_traces$B2)-mean(log(bandwidth_data$B2[is_vowel_bandwidth]), na.rm=T))/sd(log(bandwidth_data$B2[is_vowel_bandwidth]), na.rm=T)
			aaa_data[[sp]]$cor$tongue_traces$logB3_sd = (log(aaa_data[[sp]]$cor$tongue_traces$B3)-mean(log(bandwidth_data$B3[is_vowel_bandwidth]), na.rm=T))/sd(log(bandwidth_data$B3[is_vowel_bandwidth]), na.rm=T)

			formant_frequency_data = rbind(aaa_data[[sp]]$sag$tongue_traces[,c('formant_phone','formant_time','F1','F2','F3')], 
				                           aaa_data[[sp]]$cor$tongue_traces[,c('formant_phone','formant_time','F1','F2','F3')])
			is_vowel_formant = formant_frequency_data$formant_time>0 & formant_frequency_data$formant_time<1
			phone_levels = unique(formant_frequency_data[is_vowel_formant,'formant_phone'])
			F1_lm = lm(F1~formant_phone*formant_time, formant_frequency_data, subset=is_vowel_formant)
			F2_lm = lm(F2~formant_phone*formant_time, formant_frequency_data, subset=is_vowel_formant)
			F3_lm = lm(F3~formant_phone*formant_time, formant_frequency_data, subset=is_vowel_formant)

			formant_frequency_data$F1_freq_res = NA
			formant_frequency_data$F2_freq_res = NA
			formant_frequency_data$F3_freq_res = NA
			formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F1_freq_res'] = formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F1'] - predict(F1_lm, formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,])
			formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F2_freq_res'] = formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F2'] - predict(F2_lm, formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,])
			formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F3_freq_res'] = formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F3'] - predict(F3_lm, formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,])

			aaa_data[[sp]]$sag$tongue_traces$F1_freq_res = NA
			aaa_data[[sp]]$sag$tongue_traces$F2_freq_res = NA
			aaa_data[[sp]]$sag$tongue_traces$F3_freq_res = NA
			aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F1_freq_res']=aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F1']-predict(F1_lm, aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,])
			aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F2_freq_res']=aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F2']-predict(F2_lm, aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,])
			aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F3_freq_res']=aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F3']-predict(F3_lm, aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,])

			aaa_data[[sp]]$cor$tongue_traces$F1_freq_res = NA
			aaa_data[[sp]]$cor$tongue_traces$F2_freq_res = NA
			aaa_data[[sp]]$cor$tongue_traces$F3_freq_res = NA
			aaa_data[[sp]]$cor$tongue_traces[aaa_data[[sp]]$cor$tongue_traces$formant_phone%in%phone_levels,'F1_freq_res']=aaa_data[[sp]]$cor$tongue_traces[aaa_data[[sp]]$cor$tongue_traces$formant_phone%in%phone_levels,'F1']-predict(F1_lm, aaa_data[[sp]]$cor$tongue_traces[aaa_data[[sp]]$cor$tongue_traces$formant_phone%in%phone_levels,])
			aaa_data[[sp]]$cor$tongue_traces[aaa_data[[sp]]$cor$tongue_traces$formant_phone%in%phone_levels,'F2_freq_res']=aaa_data[[sp]]$cor$tongue_traces[aaa_data[[sp]]$cor$tongue_traces$formant_phone%in%phone_levels,'F2']-predict(F2_lm, aaa_data[[sp]]$cor$tongue_traces[aaa_data[[sp]]$cor$tongue_traces$formant_phone%in%phone_levels,])
			aaa_data[[sp]]$cor$tongue_traces[aaa_data[[sp]]$cor$tongue_traces$formant_phone%in%phone_levels,'F3_freq_res']=aaa_data[[sp]]$cor$tongue_traces[aaa_data[[sp]]$cor$tongue_traces$formant_phone%in%phone_levels,'F3']-predict(F3_lm, aaa_data[[sp]]$cor$tongue_traces[aaa_data[[sp]]$cor$tongue_traces$formant_phone%in%phone_levels,])

			aaa_data[[sp]]$sag$tongue_traces$F1_freq_res_sd = (aaa_data[[sp]]$sag$tongue_traces$F1_freq_res-mean(formant_frequency_data[is_vowel_formant,'F1_freq_res'], na.rm=T))/sd(formant_frequency_data[is_vowel_formant,'F1_freq_res'], na.rm=T)
			aaa_data[[sp]]$sag$tongue_traces$F2_freq_res_sd = (aaa_data[[sp]]$sag$tongue_traces$F2_freq_res-mean(formant_frequency_data[is_vowel_formant,'F2_freq_res'], na.rm=T))/sd(formant_frequency_data[is_vowel_formant,'F2_freq_res'], na.rm=T)
			aaa_data[[sp]]$sag$tongue_traces$F3_freq_res_sd = (aaa_data[[sp]]$sag$tongue_traces$F3_freq_res-mean(formant_frequency_data[is_vowel_formant,'F3_freq_res'], na.rm=T))/sd(formant_frequency_data[is_vowel_formant,'F3_freq_res'], na.rm=T)

			aaa_data[[sp]]$cor$tongue_traces$F1_freq_res_sd = (aaa_data[[sp]]$cor$tongue_traces$F1_freq_res-mean(formant_frequency_data[is_vowel_formant,'F1_freq_res'], na.rm=T))/sd(formant_frequency_data[is_vowel_formant,'F1_freq_res'], na.rm=T)
			aaa_data[[sp]]$cor$tongue_traces$F2_freq_res_sd = (aaa_data[[sp]]$cor$tongue_traces$F2_freq_res-mean(formant_frequency_data[is_vowel_formant,'F2_freq_res'], na.rm=T))/sd(formant_frequency_data[is_vowel_formant,'F2_freq_res'], na.rm=T)
			aaa_data[[sp]]$cor$tongue_traces$F3_freq_res_sd = (aaa_data[[sp]]$cor$tongue_traces$F3_freq_res-mean(formant_frequency_data[is_vowel_formant,'F3_freq_res'], na.rm=T))/sd(formant_frequency_data[is_vowel_formant,'F3_freq_res'], na.rm=T)

		}else{
			# just calculate sd from sag bandwidths because cor doesn't exist
			bandwidth_data = aaa_data[[sp]]$sag$tongue_traces[,c('B1','B2','B3')]
			is_vowel_bandwidth = grepl('[012]',bandwidth_data$phone)
			aaa_data[[sp]]$sag$tongue_traces$logB1_sd = (log(aaa_data[[sp]]$sag$tongue_traces$B1)-mean(log(bandwidth_data$B1[is_vowel_bandwidth]), na.rm=T))/sd(log(bandwidth_data$B1[is_vowel_bandwidth]), na.rm=T)
			aaa_data[[sp]]$sag$tongue_traces$logB2_sd = (log(aaa_data[[sp]]$sag$tongue_traces$B2)-mean(log(bandwidth_data$B2[is_vowel_bandwidth]), na.rm=T))/sd(log(bandwidth_data$B2[is_vowel_bandwidth]), na.rm=T)
			aaa_data[[sp]]$sag$tongue_traces$logB3_sd = (log(aaa_data[[sp]]$sag$tongue_traces$B3)-mean(log(bandwidth_data$B3[is_vowel_bandwidth]), na.rm=T))/sd(log(bandwidth_data$B3[is_vowel_bandwidth]), na.rm=T)

			formant_frequency_data = aaa_data[[sp]]$sag$tongue_traces[,c('formant_phone','formant_time','F1','F2','F3')]
			is_vowel_formant = formant_frequency_data$formant_time>0 & formant_frequency_data$formant_time<1
			phone_levels = unique(formant_frequency_data[is_vowel_formant,'formant_phone'])
			F1_lm = lm(F1~formant_phone*formant_time, formant_frequency_data, subset=is_vowel_formant)
			F2_lm = lm(F2~formant_phone*formant_time, formant_frequency_data, subset=is_vowel_formant)
			F3_lm = lm(F3~formant_phone*formant_time, formant_frequency_data, subset=is_vowel_formant)

			formant_frequency_data$F1_freq_res = NA
			formant_frequency_data$F2_freq_res = NA
			formant_frequency_data$F3_freq_res = NA
			formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F1_freq_res'] = formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F1'] - predict(F1_lm, formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,])
			formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F2_freq_res'] = formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F2'] - predict(F2_lm, formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,])
			formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F3_freq_res'] = formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,'F3'] - predict(F3_lm, formant_frequency_data[formant_frequency_data$formant_phone%in%phone_levels,])

			aaa_data[[sp]]$sag$tongue_traces$F1_freq_res = NA
			aaa_data[[sp]]$sag$tongue_traces$F2_freq_res = NA
			aaa_data[[sp]]$sag$tongue_traces$F3_freq_res = NA
			aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F1_freq_res']=aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F1']-predict(F1_lm, aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,])
			aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F2_freq_res']=aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F2']-predict(F2_lm, aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,])
			aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F3_freq_res']=aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,'F3']-predict(F3_lm, aaa_data[[sp]]$sag$tongue_traces[aaa_data[[sp]]$sag$tongue_traces$formant_phone%in%phone_levels,])

			aaa_data[[sp]]$sag$tongue_traces$F1_freq_res_sd = (aaa_data[[sp]]$sag$tongue_traces$F1_freq_res-mean(formant_frequency_data[is_vowel_formant,'F1_freq_res'], na.rm=T))/sd(formant_frequency_data[is_vowel_formant,'F1_freq_res'], na.rm=T)
			aaa_data[[sp]]$sag$tongue_traces$F2_freq_res_sd = (aaa_data[[sp]]$sag$tongue_traces$F2_freq_res-mean(formant_frequency_data[is_vowel_formant,'F2_freq_res'], na.rm=T))/sd(formant_frequency_data[is_vowel_formant,'F2_freq_res'], na.rm=T)
			aaa_data[[sp]]$sag$tongue_traces$F3_freq_res_sd = (aaa_data[[sp]]$sag$tongue_traces$F3_freq_res-mean(formant_frequency_data[is_vowel_formant,'F3_freq_res'], na.rm=T))/sd(formant_frequency_data[is_vowel_formant,'F3_freq_res'], na.rm=T)

			sag_residuals = data.frame(aaa_data[[sp]]$sag$tongue_traces[,c('phone','F1','F2','F3')], F1_res=NA, F2_res=NA, F3_res=NA)

			sag_residuals[!is.na(sag_residuals$F1),]$F1_res = residuals(F1_sag_lmer)
			sag_residuals[!is.na(sag_residuals$F2),]$F2_res = residuals(F2_sag_lmer)
			sag_residuals[!is.na(sag_residuals$F3),]$F3_res = residuals(F3_sag_lmer)

			aaa_data[[sp]]$sag$tongue_traces$F1_res_sd = (sag_residuals$F1_res-mean(sag_residuals[is_vowel_sag,'F1_res'], na.rm=T))/sd(sag_residuals[is_vowel_sag,'F1_res'], na.rm=T)
			aaa_data[[sp]]$sag$tongue_traces$F2_res_sd = (sag_residuals$F2_res-mean(sag_residuals[is_vowel_sag,'F2_res'], na.rm=T))/sd(sag_residuals[is_vowel_sag,'F2_res'], na.rm=T)
			aaa_data[[sp]]$sag$tongue_traces$F3_res_sd = (sag_residuals$F3_res-mean(sag_residuals[is_vowel_sag,'F3_res'], na.rm=T))/sd(sag_residuals[is_vowel_sag,'F3_res'], na.rm=T)

		}


# aaa_data[[sp]]$sag$tongue_traces$badness1 = with(aaa_data[[sp]]$sag$tongue_traces, sqrt(logB1_sd^2+logB2_sd^2+logB3_sd^2+F1_res_sd^2+F2_res_sd^2+F3_res_sd^2))
# aaa_data[[sp]]$sag$tongue_traces$badness1_sd = (aaa_data[[sp]]$sag$tongue_traces$badness1-mean(aaa_data[[sp]]$sag$tongue_traces[is_vowel_sag,'badness1'], na.rm=T))/sd(aaa_data[[sp]]$sag$tongue_traces[is_vowel_sag,'badness1'], na.rm=T)

# aaa_data[[sp]]$sag$tongue_traces$badness2 = with(aaa_data[[sp]]$sag$tongue_traces, sqrt(logB1_sd^2+logB2_sd^2+logB3_sd^2))
# aaa_data[[sp]]$sag$tongue_traces$badness2_sd = (aaa_data[[sp]]$sag$tongue_traces$badness2-mean(aaa_data[[sp]]$sag$tongue_traces[is_vowel_sag,'badness2'], na.rm=T))/sd(aaa_data[[sp]]$sag$tongue_traces[is_vowel_sag,'badness2'], na.rm=T)

# if ('Annotation_Label'%in%names(aaa_data[[sp]]$cor$tongue_traces)){
# 	is_vowel_cor = grepl('[012]',aaa_data[[sp]]$cor$tongue_traces$phone)
# 	aaa_data[[sp]]$cor$tongue_traces$badness2 = with(aaa_data[[sp]]$cor$tongue_traces, sqrt(logB1_sd^2+logB2_sd^2+logB3_sd^2))
# 	aaa_data[[sp]]$cor$tongue_traces$badness2_sd = (aaa_data[[sp]]$cor$tongue_traces$badness2-mean(aaa_data[[sp]]$cor$tongue_traces[is_vowel_cor,'badness2'], na.rm=T))/sd(aaa_data[[sp]]$cor$tongue_traces[is_vowel_cor,'badness2'], na.rm=T)
# }

	}
	aaa_data
}


find_mean_tongue <- function(aaa_data, speakers=NULL){
	require(dplyr)
	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		print(sp)
		for (plane in c('sag','cor')){
			if('Annotation_Label'%in%names(aaa_data[[sp]][[plane]]$tongue_traces)){

				phone_aaa = aaa_data[[sp]][[plane]]$tongue_traces
				phone_aaa_means <- phone_aaa %>% group_by(phone) %>% summarize(across(matches("^[XYR][0-9]"), ~ mean(.x, na.rm=T)))
				grand_aaa_means = colMeans(phone_aaa_means[,-1], na.rm=T)
				grand_aaa_means_XYR = data.frame(X=grand_aaa_means[grepl("^X[0-9]",names(grand_aaa_means))],
					                             Y=grand_aaa_means[grepl("^Y[0-9]",names(grand_aaa_means))],
					                             R=grand_aaa_means[grepl("^R[0-9]",names(grand_aaa_means))])
				# grand_aaa_means_XY_rotated = rotateXY(grand_aaa_means_XY, aaa_data[[sp]][[plane]]$occlusal_angle, center=c(0,0))
				# print(names(phone_aaa))
				aaa_data[[sp]][[plane]]$mean_tongue = grand_aaa_means_XYR
				aaa_data[[sp]][[plane]]$max_tongue = data.frame(X=c(colQuantile(phone_aaa[,paste0("X",1:21)],0.999),
					                                                colQuantile(phone_aaa[,paste0("X",22:42)],0.001)),
					                                            Y=colQuantile(phone_aaa[,paste0("Y",1:42)],0.999),
					                                            R=colQuantile(phone_aaa[,paste0("R",0:41)],0.999))
				aaa_data[[sp]][[plane]]$min_tongue = data.frame(X=c(colQuantile(phone_aaa[,paste0("X",1:21)],0.001),
					                                                colQuantile(phone_aaa[,paste0("X",22:42)],0.999)),
					                                            Y=colQuantile(phone_aaa[,paste0("Y",1:42)],0.001),
					                                            R=colQuantile(phone_aaa[,paste0("R",0:41)],0.001))

			}
		}
	}
	aaa_data
}

# this was to do something similar with polar coordinates but it doesn't work very well:
# all_x = c()
# all_y = c()
# for (ph in unique(aaa_data[[sp]][[plane]]$tongue_traces$phone)){
# 	print(ph)
# 	subdata = aaa_data[[sp]][[plane]]$tongue_traces_polar[aaa_data[[sp]][[plane]]$tongue_traces$phone==ph,]
# 	X = subdata[,grepl("^T[0-9]", names(subdata))]
# 	Y = subdata[,grepl("^R[0-9]", names(subdata))]
# 	plot(unlist(X),unlist(Y))
# 	smsp=smooth.spline(na.omit(unlist(X)),na.omit(unlist(Y)),spar=1)
# 	points(smsp,type='l')
# 	all_x = c(all_x, smsp$x)
# 	all_y = c(all_y, smsp$y)
# }
# smsp_all=smooth.spline(all_x,all_y,spar=0.9)
# smsp_all_cart = make.cartesian(cbind(-smsp_all$x, smsp_all$y), origin=aaa_data[[sp]][[plane]]$origin, flip=FALSE)
# # smsp_all_cart[,2] = -smsp_all_cart[,2]
# points(smsp_all_cart[,1], smsp_all_cart[,2], lty=1, lwd=5, col='darkgray')


find_displacement = function(scaled_radii){
	scaled_radii_rle = rle(as.numeric(scaled_radii/abs(scaled_radii)))
	# print(scaled_radii_rle)
	if (1 %in% scaled_radii_rle$values){

		sorted_spanlengths = c(rev(sort(scaled_radii_rle$lengths[scaled_radii_rle$values==1])),NA,NA,NA)
		maxlength = sorted_spanlengths[1]
		maxlength2 = sorted_spanlengths[2]
		maxlength3 = sorted_spanlengths[3]

		if (is.na(maxlength)){
			which_maxlength=NA
			which_maxlength2=NA
		}else if (is.na(maxlength2)){
			which_maxlength = which(scaled_radii_rle$lengths==maxlength & scaled_radii_rle$values==1)
			which_maxlength2=NA
		}else if (maxlength==maxlength2){
			which_maxlength = which(scaled_radii_rle$lengths==maxlength & scaled_radii_rle$values==1)[2]
			which_maxlength2 = which(scaled_radii_rle$lengths==maxlength & scaled_radii_rle$values==1)[1]
		}else{
			which_maxlength = which(scaled_radii_rle$lengths==maxlength & scaled_radii_rle$values==1)
			which_maxlength2 = which(scaled_radii_rle$lengths==maxlength2 & scaled_radii_rle$values==1)[1]
		}
		# maxlength = max(scaled_radii_rle$lengths[scaled_radii_rle$values==1],na.rm=T)[1]
		# which_maxlength = which(scaled_radii_rle$lengths==maxlength & scaled_radii_rle$values==1)[1]

		maxspan_end = sum(scaled_radii_rle$lengths[1:(which_maxlength)])
		maxspan_start = maxspan_end-(scaled_radii_rle$lengths[which_maxlength]-1)
		span_vals = scaled_radii[maxspan_start:maxspan_end]
		displacement_cog = sum((span_vals * maxspan_start:maxspan_end)/sum(span_vals))
		displacement_mean = mean(as.numeric(span_vals))
		displacement_max = max(as.numeric(span_vals))
		displacement_mass = sum(span_vals)

		if (is.na(maxlength2)){
			c(displacement_cog, displacement_mean, displacement_max, displacement_mass,
			  NA, 0, 0, 0,
			  maxlength,maxlength2,maxlength3)
		}else{
			maxspan2_end = sum(scaled_radii_rle$lengths[1:(which_maxlength2)])
			maxspan2_start = maxspan2_end-(scaled_radii_rle$lengths[which_maxlength2]-1)
			span2_vals = scaled_radii[maxspan2_start:maxspan2_end]
			displacement_cog2 = sum((span2_vals * maxspan2_start:maxspan2_end)/sum(span2_vals))
			displacement_mean2 = mean(as.numeric(span2_vals))
			displacement_max2 = max(as.numeric(span2_vals))
			displacement_mass2 = sum(span2_vals)

			c(displacement_cog, displacement_mean, displacement_max, displacement_mass,
			  displacement_cog2, displacement_mean2, displacement_max2, displacement_mass2,
			  maxlength,maxlength2,maxlength3)
		}
	}else{
		c(NA,0,0,0,NA,0,0,0,NA,NA,NA)
	}
}


find_all_displacements = function(aaa_data, speakers=NULL, planes=c('sag','cor')){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		for (plane in planes){
			if('Annotation_Label'%in%names(aaa_data[[sp]][[plane]]$tongue_traces)){
				print(paste(sp,plane))

				mean_tongue_R = aaa_data[[sp]][[plane]]$mean_tongue$R
				max_tongue_R = aaa_data[[sp]][[plane]]$max_tongue$R

				all_scaled_radii = data.frame(mapply(`/`,data.frame(mapply(`-`,aaa_data[[sp]][[plane]]$tongue_traces[,paste0('R',0:41)],mean_tongue_R)),max_tongue_R-mean_tongue_R))
				# print('x')
				displacement_info = data.frame(t(as.matrix(apply(all_scaled_radii, 1, find_displacement))))
				# print('y')
				names(displacement_info)=c('displacement_cog', 'displacement_mean', 'displacement_max', 'displacement_mass',
					'displacement_cog2', 'displacement_mean2', 'displacement_max2', 'displacement_mass2',
					'span1','span2','span3')

				displacement_info$norm_displacement_cog = (displacement_info$displacement_cog-mean(displacement_info$displacement_cog,na.rm=T))/sd(displacement_info$displacement_cog,na.rm=T)
				displacement_info$norm_displacement_max = (displacement_info$displacement_max-mean(displacement_info$displacement_max,na.rm=T))/sd(displacement_info$displacement_max,na.rm=T)
				displacement_info$delta = c(NA,sqrt(diff(displacement_info$norm_displacement_cog)^2+diff(displacement_info$norm_displacement_max)^2))

				displacement_info$norm_displacement_cog2 = (displacement_info$displacement_cog2-mean(displacement_info$displacement_cog2,na.rm=T))/sd(displacement_info$displacement_cog2,na.rm=T)
				displacement_info$norm_displacement_max2 = (displacement_info$displacement_max2-mean(displacement_info$displacement_max2,na.rm=T))/sd(displacement_info$displacement_max2,na.rm=T)
				displacement_info$delta2 = c(NA,sqrt(diff(displacement_info$norm_displacement_cog2)^2+diff(displacement_info$norm_displacement_max2)^2))

				for (colname in colnames(displacement_info)){
					aaa_data[[sp]][[plane]]$tongue_traces[,colname] = NULL
				}
				aaa_data[[sp]][[plane]]$tongue_traces = cbind(aaa_data[[sp]][[plane]]$tongue_traces,displacement_info)
			}
		}
	}
	aaa_data
}

phone_locations_plot_nongeneric <- function(aaa_data, speakers=NULL){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		plane='sag'

		phone_locations_static = subset(aaa_data[[sp]][[plane]]$phone_locations, phone%in%c('K','G','L','W','NG','R'))

		phone_locations_dynamic = subset(aaa_data[[sp]][[plane]]$phone_locations, grepl('1',phone)&!grepl('1R',phone))

		# cairo_pdf('test_displacement_again.pdf', onefile=T)
		plot(0,0, type='n', main=paste(sp,'dorsal sounds'), xlab='displacement_cog (anterior to posterior)', ylab='displacement_max (0=neutral, 1=99.9th percentile)',
			xlim=range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_cog, na.rm=T), 
			ylim=range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_max, na.rm=T))

		# phone_locations_l = subset(phone_locations_dynamic, grepl('L',phone))
		# lchull = chull(phone_locations_l[,c('x3','y3')])

		# polygon(phone_locations_l[lchull,'x3'],phone_locations_l[lchull,'y3'],col='gray',border=NA)

		arrows(phone_locations_dynamic$x1,phone_locations_dynamic$y1,
			   phone_locations_dynamic$x3,phone_locations_dynamic$y3, length=0.1)
		points(phone_locations_dynamic$x1,phone_locations_dynamic$y1,pch=21,cex=1.75,col='black',
			bg=ifelse(phone_locations_dynamic$phone%in%c('AA1','AO1'),'pink','white'))
		text(phone_locations_dynamic$x1,phone_locations_dynamic$y1,labels=phone_locations_dynamic$phone, cex=0.25)

		points(phone_locations_static$x,phone_locations_static$y,pch=21,cex=1.75,col='black',
			bg=ifelse(phone_locations_static$phone%in%c('K','G','NG','W'),'cyan','white'))
		text(phone_locations_static$x,phone_locations_static$y,labels=phone_locations_static$phone, cex=0.25)

	}
}

phone_locations_plot <- function(aaa_data, speakers=NULL, phone_table, main_label='dorsal sounds'){

	require(scales)

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		plane='sag'

		phone_table$fillcolor = alpha(phone_table$color,0.3)
		phone_locations = merge(aaa_data[[sp]][[plane]]$phone_locations, phone_table)

		# cairo_pdf('test_displacement_again.pdf', onefile=T)
		plot(0,0, type='n', main=paste(sp,main_label), xlab='displacement_cog (posterior to anterior)', ylab='displacement_max (0=neutral, 1=99.9th percentile)',
			xlim=rev(range(aaa_data[[sp]][[plane]]$phone_locations[,grepl('^x',names(aaa_data[[sp]][[plane]]$phone_locations))], na.rm=T)), 
			ylim=range(aaa_data[[sp]][[plane]]$phone_locations[,grepl('^y',names(aaa_data[[sp]][[plane]]$phone_locations))], na.rm=T))

		phone_locations_l = phone_locations[phone_locations$dynamic & phone_locations$for_chull2,]
		lchull = chull(phone_locations_l[,c('x3','y3')])
		polygon(phone_locations_l[lchull,'x3'],phone_locations_l[lchull,'y3'],col='lightgray',border=NA)

		phone_locations_l = phone_locations[phone_locations$dynamic & phone_locations$for_chull1,]
		lchull = chull(phone_locations_l[,c('x3','y3')])
		polygon(phone_locations_l[lchull,'x3'],phone_locations_l[lchull,'y3'],col='darkgray',border=NA)

		if (sum(phone_locations$dynamic)){
			arrows(phone_locations[phone_locations$dynamic,'x1'],phone_locations[phone_locations$dynamic,'y1'],
				   phone_locations[phone_locations$dynamic,'x3'],phone_locations[phone_locations$dynamic,'y3'], 
				   col=phone_locations[phone_locations$dynamic,'color'], length=0.1)
			points(phone_locations[phone_locations$dynamic,'x1'],phone_locations[phone_locations$dynamic,'y1'],
				   col=NA, 
				   bg='white', pch=phone_locations[phone_locations$dynamic,'pch'], cex=1.75)
			points(phone_locations[phone_locations$dynamic,'x1'],phone_locations[phone_locations$dynamic,'y1'],
				   col=phone_locations[phone_locations$dynamic,'color'], 
				   bg=phone_locations[phone_locations$dynamic,'fillcolor'], pch=phone_locations[phone_locations$dynamic,'pch'],cex=1.75)
			text(phone_locations[phone_locations$dynamic,'x1'],phone_locations[phone_locations$dynamic,'y1'],
				 labels=phone_locations[phone_locations$dynamic,'phone'], 
				 col=phone_locations[phone_locations$dynamic,'color'], cex=0.25)
		}
		points(phone_locations[!phone_locations$dynamic,'x'],phone_locations[!phone_locations$dynamic,'y'],
			   col=phone_locations[!phone_locations$dynamic,'color'], 
			   bg=phone_locations[!phone_locations$dynamic,'fillcolor'],
			   pch=phone_locations[!phone_locations$dynamic,'pch'],cex=1.75)
		text(phone_locations[!phone_locations$dynamic,'x'],phone_locations[!phone_locations$dynamic,'y'],
			labels=phone_locations[!phone_locations$dynamic,'phone'], 
			col=phone_locations[!phone_locations$dynamic,'color'], cex=0.25)

	}
}

find_displacement_density <- function(aaa_data, speakers=NULL, kde_n=50){
	# require(ggplot2)
	require(MASS)

	name_of=list(P='Posterior', A='Anterior')

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		plane='sag'
		cairo_pdf(paste0('displacement_density_',sp,'.pdf'), onefile=T, h=6, w=6)

		phone_locations = c()
		# xlim=range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_cog, na.rm=T)
		xlim=range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_cog, na.rm=T)
		ylim=range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_max, na.rm=T)
		# print(xlim)
		# print(ylim)
		for (ph in sort(unique(aaa_data[[sp]][[plane]]$tongue_traces$phone))){
		# for (ph in 'AH1'){
			print(paste(sp,ph))
			new_phone_locations=data.frame(phone=ph,x=NA, y=NA, x1=NA, y1=NA, x2=NA, y2=NA, x3=NA, y3=NA,
				                                    xP=NA,yP=NA,xP1=NA,yP1=NA,xP2=NA,yP2=NA,xP3=NA,yP3=NA,
				                                    xA=NA,yA=NA,xA1=NA,yA1=NA,xA2=NA,yA2=NA,xA3=NA,yA3=NA)
			subdata = subset(aaa_data[[sp]][[plane]]$tongue_traces, phone==ph)

			subdata$displacement_cogA = subdata$displacement_cog
			subdata$displacement_maxA = subdata$displacement_max
			subdata[which(subdata$largest_is_posterior),'displacement_cogA'] = subdata[which(subdata$largest_is_posterior),'displacement_cog2']
			subdata[which(subdata$largest_is_posterior),'displacement_maxA'] = subdata[which(subdata$largest_is_posterior),'displacement_max2']
			subdata$displacement_cogP = subdata$displacement_cog
			subdata$displacement_maxP = subdata$displacement_max
			subdata[which(!subdata$largest_is_posterior),'displacement_cogP'] = subdata[which(!subdata$largest_is_posterior),'displacement_cog2']
			subdata[which(!subdata$largest_is_posterior),'displacement_maxP'] = subdata[which(!subdata$largest_is_posterior),'displacement_max2']

			# print(subdata[,c('displacement_cog','displacement_cogA','displacement_cogP')])
			# print(subdata[,c('displacement_max','displacement_maxA','displacement_maxP')])
			third = ceiling(subdata$phone_time*3)
			
			# print(subdata$phone_time)
			# xlim = range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_cog,na.rm=T)
			# ylim = range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_max,na.rm=T)

			kde_results = list()
			for (PA_specifier in c('','P','A')){
				# print (PA_specifier)
				what_cog = paste0('displacement_cog',PA_specifier)
				what_max = paste0('displacement_max',PA_specifier)
				kdedata=subdata[!is.na(subdata[,what_cog])&!is.na(subdata[,what_max]),c(what_cog,what_max)]
				if (nrow(kdedata)>1){
					# print (nrow(kdedata))
					kde_results[[paste0('kde_',PA_specifier)]]=kde2d(x=kdedata[,what_cog], y=kdedata[,what_max], n=kde_n, lims=c(xlim,ylim))
					which_max_z = which(kde_results[[paste0('kde_',PA_specifier)]]$z==max(kde_results[[paste0('kde_',PA_specifier)]]$z))[1]
					which_max_x = (which_max_z-1)%%kde_n+1
					which_max_y = floor((which_max_z-1)/kde_n)+1
					new_phone_locations[,paste0('x',PA_specifier)]=kde_results[[paste0('kde_',PA_specifier)]]$x[which_max_x]
					new_phone_locations[,paste0('y',PA_specifier)]=kde_results[[paste0('kde_',PA_specifier)]]$y[which_max_y]
				}

				kdedata1=subdata[third==1&!is.na(subdata[,what_cog])&!is.na(subdata[,what_max]),c(what_cog,what_max)]
				# print(nrow(kdedata1))
				if (nrow(kdedata1)>1){
					kde_results[[paste0('kde_',PA_specifier,'1')]]=kde2d(x=kdedata1[,what_cog], y=kdedata1$displacement_max, n=kde_n, lims=c(xlim,ylim))
					which_max_z1 = which(kde_results[[paste0('kde_',PA_specifier,'1')]]$z==max(kde_results[[paste0('kde_',PA_specifier,'1')]]$z))[1]
					which_max_x1 = (which_max_z1-1)%%kde_n+1
					which_max_y1 = floor((which_max_z1-1)/kde_n)+1
					new_phone_locations[,paste0('x',PA_specifier,'1')]=kde_results[[paste0('kde_',PA_specifier,'1')]]$x[which_max_x1]
					new_phone_locations[,paste0('y',PA_specifier,'1')]=kde_results[[paste0('kde_',PA_specifier,'1')]]$y[which_max_y1]
				}

				kdedata2=subdata[third==2&!is.na(subdata[,what_cog])&!is.na(subdata[,what_max]),c(what_cog,what_max)]
				if (nrow(kdedata2)>1){
					kde_results[[paste0('kde_',PA_specifier,'2')]]=kde2d(x=kdedata2[,what_cog], y=kdedata2$displacement_max, n=kde_n, lims=c(xlim,ylim))
					which_max_z2 = which(kde_results[[paste0('kde_',PA_specifier,'2')]]$z==max(kde_results[[paste0('kde_',PA_specifier,'2')]]$z))[1]
					which_max_x2 = (which_max_z2-1)%%kde_n+1
					which_max_y2 = floor((which_max_z2-1)/kde_n)+1
					new_phone_locations[,paste0('x',PA_specifier,'2')]=kde_results[[paste0('kde_',PA_specifier,'2')]]$x[which_max_x2]
					new_phone_locations[,paste0('y',PA_specifier,'2')]=kde_results[[paste0('kde_',PA_specifier,'2')]]$y[which_max_y2]
				}

				kdedata3=subdata[third==3&!is.na(subdata[,what_cog])&!is.na(subdata[,what_max]),c(what_cog,what_max)]
				if (nrow(kdedata3)>1){
					kde_results[[paste0('kde_',PA_specifier,'3')]]=kde2d(x=kdedata3[,what_cog], y=kdedata3$displacement_max, n=kde_n, lims=c(xlim,ylim))
					which_max_z3 = which(kde_results[[paste0('kde_',PA_specifier,'3')]]$z==max(kde_results[[paste0('kde_',PA_specifier,'3')]]$z))[1]
					which_max_x3 = (which_max_z3-1)%%kde_n+1
					which_max_y3 = floor((which_max_z3-1)/kde_n)+1
					new_phone_locations[,paste0('x',PA_specifier,'3')]=kde_results[[paste0('kde_',PA_specifier,'3')]]$x[which_max_x3]
					new_phone_locations[,paste0('y',PA_specifier,'3')]=kde_results[[paste0('kde_',PA_specifier,'3')]]$y[which_max_y3]
				}
			
			}

			for (PA_specifier in c('','P','A')){
				main_base = ph  
				if (PA_specifier %in% names(name_of)){
					main_base = paste(main_base, name_of[[PA_specifier]])
				}
				plot(0,0, type='n', main=main_base, xlim=xlim, ylim=ylim, 
					xlab='displacement_cog (anterior to posterior)', ylab='displacement_max (0=neutral, 1=99.9th percentile)')
				for (tk in unique(subdata$token)){
					subsubdata = subset(subdata, token==tk)
					points(subsubdata[,paste0('displacement_cog',PA_specifier)][1], 
						   subsubdata[,paste0('displacement_max',PA_specifier)][1], 
						   col=rainbow(1.25*nrow(subsubdata))[1],pch=19)
					if (nrow(subsubdata)>1){
						for (i in 1:(nrow(subsubdata)-1)){
							if (i<(nrow(subsubdata)-1)){
								points(subsubdata[,paste0('displacement_cog',PA_specifier)][i:(i+1)], subsubdata[,paste0('displacement_max',PA_specifier)][i:(i+1)], 
									   type='l', col=rainbow(1.25*nrow(subsubdata))[i],lwd=2)
							}else{
								arrows(subsubdata[,paste0('displacement_cog',PA_specifier)][i], subsubdata[,paste0('displacement_max',PA_specifier)][i], 
									   subsubdata[,paste0('displacement_cog',PA_specifier)][i+1], subsubdata[,paste0('displacement_max',PA_specifier)][i+1], 
									   type='l', col=rainbow(1.25*nrow(subsubdata))[i], length=0.1,lwd=2)
							}
						}
					}
				}

				# print(subdata$phone_time)
				plot(0,0, type='n', main=main_base, xlim=xlim, ylim=ylim, 
					xlab='displacement_cog (anterior to posterior)', ylab='displacement_max (0=neutral, 1=99.9th percentile)')
				this_kde_result = kde_results[[paste0('kde_',PA_specifier)]]
				if (!is.null(this_kde_result$z)){
					image(this_kde_result$x, this_kde_result$y, this_kde_result$z, 
						useRaster=TRUE, add=TRUE, col = hcl.colors(36, "YlOrRd", rev = TRUE))
				}
				points(new_phone_locations[,paste0('x',PA_specifier)],new_phone_locations[,paste0('y',PA_specifier)],col='white', cex=2)
				points(new_phone_locations[,paste0('x',PA_specifier,'1')],new_phone_locations[,paste0('y',PA_specifier,'1')],col='red', cex=2)
				points(new_phone_locations[,paste0('x',PA_specifier,'2')],new_phone_locations[,paste0('y',PA_specifier,'2')],col='green', cex=2)
				points(new_phone_locations[,paste0('x',PA_specifier,'3')],new_phone_locations[,paste0('y',PA_specifier,'3')],col='blue', cex=2)

				segments(new_phone_locations[,paste0('x',PA_specifier,'1')],new_phone_locations[,paste0('y',PA_specifier,'1')],
					     new_phone_locations[,paste0('x',PA_specifier,'2')],new_phone_locations[,paste0('y',PA_specifier,'2')])
				arrows(new_phone_locations[,paste0('x',PA_specifier,'2')],new_phone_locations[,paste0('y',PA_specifier,'2')],
					   new_phone_locations[,paste0('x',PA_specifier,'3')],new_phone_locations[,paste0('y',PA_specifier,'3')], length=0.1)
						
				plot(0,0, type='n', main=paste(main_base,'first third'), 
					# xlim=range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_cog, na.rm=T), 
					# ylim=range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_max, na.rm=T), 
					xlim=xlim, ylim=ylim, 
					xlab='displacement_cog (anterior to posterior)', ylab='displacement_max (0=neutral, 1=99.9th percentile)')
				# print (names(kde_results))
				this_kde_result = kde_results[[paste0('kde_',PA_specifier,'1')]]
				if (!is.null(this_kde_result$z)){
					# print(paste0('kde_1',PA_specifier))
					# print (this_kde_result$x)
					image(this_kde_result$x, this_kde_result$y, this_kde_result$z, 
						useRaster=TRUE, add=TRUE, col = hcl.colors(36, "YlOrRd", rev = TRUE))
				}
				points(new_phone_locations[,paste0('x',PA_specifier,'1')],new_phone_locations[,paste0('y',PA_specifier,'1')],col='red', cex=2)

				plot(0,0, type='n', main=paste(main_base,'middle third'), 
					# xlim=range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_cog, na.rm=T), 
					# ylim=range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_max, na.rm=T), 
					xlim=xlim, ylim=ylim, 
					xlab='displacement_cog (anterior to posterior)', ylab='displacement_max (0=neutral, 1=99.9th percentile)')

				this_kde_result = kde_results[[paste0('kde_',PA_specifier,'2')]]
				if (!is.null(this_kde_result$z)){
					image(this_kde_result$x, this_kde_result$y, this_kde_result$z, 
						useRaster=TRUE, add=TRUE, col = hcl.colors(36, "YlOrRd", rev = TRUE))
				}
				points(new_phone_locations[,paste0('x',PA_specifier,'2')],new_phone_locations[,paste0('y',PA_specifier,'2')],col='green', cex=2)

				plot(0,0, type='n', main=paste(main_base,'last third'), 
					# xlim=range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_cog, na.rm=T), 
					# ylim=range(aaa_data[[sp]][[plane]]$tongue_traces$displacement_max, na.rm=T), 
					xlim=xlim, ylim=ylim, 
					xlab='displacement_cog (anterior to posterior)', ylab='displacement_max (0=neutral, 1=99.9th percentile)')

				this_kde_result = kde_results[[paste0('kde_',PA_specifier,'3')]]
				if (!is.null(this_kde_result$z)){
					image(this_kde_result$x, this_kde_result$y, this_kde_result$z, 
						useRaster=TRUE, add=TRUE, col = hcl.colors(36, "YlOrRd", rev = TRUE))
				}
				points(new_phone_locations[,paste0('x',PA_specifier,'3')],new_phone_locations[,paste0('y',PA_specifier,'3')],col='blue', cex=2)

			}
			phone_locations = rbind(phone_locations, new_phone_locations)
		}
		aaa_data[[sp]][[plane]]$phone_locations = phone_locations

		# phone_locations_plot_nongeneric(aaa_data, sp)

		dev.off()
	}
	aaa_data
}


add_target_phone_time <- function(aaa_data, speakers=NULL, planes=c('sag','cor')){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		for (plane in planes){

			print(paste(sp,plane))

			if('Annotation_Label'%in%names(aaa_data[[sp]][[plane]]$tongue_traces)){


				aaa_data[[sp]][[plane]]$tongue_traces$target_phone = FALSE
				for (tk in unique(aaa_data[[sp]][[plane]]$tongue_traces$token)){
					which_token = aaa_data[[sp]][[plane]]$tongue_traces$token==tk

					# find target phone
					phones = rle(aaa_data[[sp]][[plane]]$tongue_traces[which_token,'phone'])$values
					which_target_phone = which(grepl('1',phones))[1]
					target_phone = phones[which_target_phone]

					which_token_target_phone = aaa_data[[sp]][[plane]]$tongue_traces$token==tk & aaa_data[[sp]][[plane]]$tongue_traces$phone==target_phone
					aaa_data[[sp]][[plane]]$tongue_traces[which_token_target_phone,'target_phone'] = TRUE

					# phones[which_target_phone] = paste0('*',target_phone,'*')
					# print(paste(c(tk, length(target_phone), paste(phones,collapse='-'))))

					# relative_time scaled so that 0-1 is target phone interval
					phone_rel = coef(lm(phone_time~relative_time, aaa_data[[sp]][[plane]]$tongue_traces[which_token_target_phone,]))
					aaa_data[[sp]][[plane]]$tongue_traces[which_token,'target_phone_time'] = aaa_data[[sp]][[plane]]$tongue_traces[which_token,'relative_time'] * phone_rel[2] + phone_rel[1]
					
				}
			}
		}
	}
	aaa_data
}





make_trajectory_comparisons <- function(aaa_data, trajectory_comparisons, trajectory_ssanovas=list(), sp, plane='sag'){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		print (sp)

		if (!sp%in%names(trajectory_ssanovas)){
			trajectory_ssanovas[[sp]] = list()
		}

		for (tc in names(trajectory_comparisons)){
			tokens_data = c()
			labels = trajectory_comparisons[[tc]]$labels
			ref_phone = trajectory_comparisons[[tc]]$ref_phone
			ref_range = trajectory_comparisons[[tc]]$ref_range
			if (is.null(trajectory_comparisons[[tc]]$phonex)){
				phonex = round(subset(aaa_data[[sp]][[plane]]$phone_locations, phone==ref_phone)[1,'x']) + (-ref_range:ref_range)
			}else{
				phonex = trajectory_comparisons[[tc]]$phonex
			}
			for (label in labels){
				matching_tokens = unique(subset(aaa_data[[sp]][[plane]]$tongue_traces, Annotation_Label==label)$token)
				for (tid in matching_tokens){
					print(tid)
					one_token = subset(aaa_data[[sp]][[plane]]$tongue_traces, token==tid)
					scaled_radii = data.frame(mapply(`/`,data.frame(mapply(`-`,one_token[,paste0('R',0:41)],mean_tongue_R)),max_tongue_R-mean_tongue_R))
					max_loc = phonex[1]-1 + apply(scaled_radii[,phonex], 1, function(x) which(x==max(x, na.rm=T))[1])
					max_val = apply(scaled_radii[,phonex], 1, function(x) max(x, na.rm=T))
					tokens_data = rbind(tokens_data, data.frame(label=label,token=tid,max_loc=max_loc,max_val=max_val,target_phone_time=one_token$target_phone_time, scaled_radii))
				}
			}
			 
			ylim=c(-1,1)

			tokens_data$X = tokens_data$target_phone_time
			tokens_data$label = factor(tokens_data$label)
			for (l in rev(labels)){
				if (l %in% levels(tokens_data$label)){
					tokens_data$label = relevel(tokens_data$label, l)
				}
			}
			# 
			
			# main = paste(unique(tokens_data$label), collapse=' ')
			main = sp

			trajectory_ssanovas[[sp]][[tc]] = list()

			trajectory_ssanovas[[sp]][[tc]][['tokens_data']] = tokens_data
			
			tokens_data$Y = (tokens_data$max_loc-gx)/(ax-gx)
			trajectory_ssanovas[[sp]][[tc]][['max_loc_ssanova']] = cart.ssanova(tokens_data, data.cat='label', 
			                         CI.fill=TRUE, flip=FALSE, plot.labels=c(main,'relative time','constriction location'), 
			                         xlim=c(0,1), ylim=ylim, legend.pos='bottomright', crop=TRUE, plotting=FALSE)

			tokens_data$Y = tokens_data$max_val
			trajectory_ssanovas[[sp]][[tc]][['max_val_ssanova']] = cart.ssanova(tokens_data, data.cat='label', 
			                         CI.fill=TRUE, flip=FALSE, plot.labels=c(main,'relative time','constriction degree (open to closed)'), 
			                         xlim=c(0,1), ylim=ylim, legend.pos='bottomright', crop=TRUE, overplot=T)
			
		}
	}
	trajectory_ssanovas
}



plot_trajectory_comparisons <- function(aaa_data, trajectory_comparisons, trajectory_ssanovas, speaker, plane='sag'){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		print (sp)

		# plane='sag'
		# sp='cl14'
		gx = subset(aaa_data[[sp]][[plane]]$phone_locations, phone=='G')[1,'x']
		ax = subset(aaa_data[[sp]][[plane]]$phone_locations, phone=='AA1')[1,'x']
		normalized_angles = (seq(0,41)-gx)/(ax-gx)
		mean_tongue_R = aaa_data[[sp]][[plane]]$mean_tongue$R
		max_tongue_R = aaa_data[[sp]][[plane]]$max_tongue$R

		par(mfrow=c(1,2))

		for (tc in names(trajectory_comparisons)){
			ylim=c(-1,1)
			plot(0,0,type='n',xlim=c(-1,2), ylim=ylim, main=sp, xlab='relative time', ylab='constriction degree (open to closed)')
			abline(v=c(0,1), col='gray')

			tokens_data = trajectory_ssanovas[[sp]][[tc]][['tokens_data']]

			max_val_ssanova = trajectory_ssanovas[[sp]][[tc]][['max_val_ssanova']]
			replot.tongue.ss(max_val_ssanova, overplot=T)

			# max_val_ssanova = cart.ssanova(tokens_data, data.cat='label', 
			#                          CI.fill=TRUE, flip=FALSE, plot.labels=c(main,'relative time','constriction degree (open to closed)'), 
			#                          xlim=c(0,1), ylim=ylim, legend.pos='bottomright', crop=TRUE, overplot=T)

			max_loc_ssanova = trajectory_ssanovas[[sp]][[tc]][['max_loc_ssanova']]
			# max_loc_ssanova = cart.ssanova(tokens_data, data.cat='label', 
			#                          CI.fill=TRUE, flip=FALSE, plot.labels=c(main,'relative time','constriction location'), 
			#                          xlim=c(0,1), ylim=ylim, legend.pos='bottomright', crop=TRUE, plotting=FALSE)

			plot(0,0,type='n',xlim=c(1.5,-1),ylim=c(-1,1),main=main,
				xlab='constriction location (posterior to anterior)', ylab='constriction degree (open to closed)')
			segments(c(0,1), rep(-0.925,2), c(0,1), rep(1,2), col='gray', lwd=2)
			abline(v=(range(phonex)-gx)/(ax-gx), col='gray')
			text(c(0,1), rep(-1,2), labels=c('ɡ','ɑ'),col='gray')
			for (i in 1:length(unique(tokens_data$label))){
				annotation_label = unique(tokens_data$label)[i]
				subdata = subset(tokens_data, label==annotation_label)

				for (tid in unique(subdata$token)){
					token_data = subset(subdata, token==tid)
					for (tr in 1:nrow(token_data)){
						token_row = token_data[tr,]
						points((seq(1,42)-gx)/(ax-gx), token_row[,paste0('R',0:41)], type='l', col=rainbow(length(unique(tokens_data$label)),v=0.7,a=0.05)[i])
					}
				}

				# max_val_spline = smooth.spline(subdata$relative_time, subdata$max_val)
				# max_loc_spline = smooth.spline(subdata$relative_time, (subdata$max_loc-gx)/(ax-gx))

				max_val_spline = data.frame(x=subset(max_val_ssanova$ss.cart, label==annotation_label)$X,
					                          y=subset(max_val_ssanova$ss.cart, label==annotation_label)$ss.Fit)
				max_loc_spline = data.frame(x=subset(max_loc_ssanova$ss.cart, label==annotation_label)$X,
					                          y=subset(max_loc_ssanova$ss.cart, label==annotation_label)$ss.Fit)

				max_val_spline_target = max_val_spline[max_val_spline$x>=0&max_val_spline$x<=1,]
				max_loc_spline_target = max_loc_spline[max_loc_spline$x>=0&max_loc_spline$x<=1,]

				points(max_loc_spline_target$y,max_val_spline_target$y, type='l',
					col=rainbow(length(unique(tokens_data$label)),v=0.7)[i], lwd=3)
				arrows(max_loc_spline_target$y[length(max_loc_spline_target$y)-2], max_val_spline_target$y[length(max_val_spline_target$y)-2],
					     max_loc_spline_target$y[length(max_loc_spline_target$y)], max_val_spline_target$y[length(max_val_spline_target$y)],
					col=rainbow(length(unique(tokens_data$label)),v=0.7)[i], lwd=2, length=0.1)

			}
			legend('bottomright',legend=unique(tokens_data$label),lty=1,lwd=2,col=rainbow(length(unique(tokens_data$label)),v=0.7))
		}
	}
}



tgfile2list <- function(textgrid_filepath, tiers=NULL){
  
  tglines = readLines(textgrid_filepath)
  tglist <- list()
  
  for (i in 1:length(tglines)){
    line = tglines[i]
    # DOES IT INTRODUCE A NEW TIER?
    if (grepl('name = ',line)){
      tiername = strsplit(line,'"')[[1]][2]
      tglist[[tiername]] = c()
    # IS IT A TIER INTERVAL?
    }else if (grepl('intervals \\[[0-9]*\\]',line)){
      tglist[[tiername]] = rbind(tglist[[tiername]], 
                               data.frame(xmin = as.numeric(strsplit(tglines[i+1],'=')[[1]][2]),
                                          xmax = as.numeric(strsplit(tglines[i+2],'=')[[1]][2]),
                                          text = strsplit(tglines[i+3],'"')[[1]][2]))
    }
  }

  for (tiername in names(tglist)){
    newtier = tglist[[tiername]]
    if (nrow(newtier)==1){
      tglist[[paste0('last',tiername)]] <- data.frame(xmin=NA, xmax=NA, text=NA)
      tglist[[paste0('next',tiername)]] <- data.frame(xmin=NA, xmax=NA, text=NA)
      tglist[[paste0('lastlast',tiername)]] <- data.frame(xmin=NA, xmax=NA, text=NA)
      tglist[[paste0('nextnext',tiername)]] <- data.frame(xmin=NA, xmax=NA, text=NA)
    }else{
      tglist[[paste0('last',tiername)]] <- data.frame(xmin=as.numeric(newtier[,1]), xmax=as.numeric(newtier[,2]), text=c(NA,newtier[1:(nrow(newtier)-1),3]))
      tglist[[paste0('next',tiername)]] <- data.frame(xmin=as.numeric(newtier[,1]), xmax=as.numeric(newtier[,2]), text=c(newtier[2:nrow(newtier),3], NA))
      if (nrow(newtier)==2){
        tglist[[paste0('lastlast',tiername)]] <- data.frame(xmin=c(NA,NA), xmax=c(NA,NA), text=c(NA,NA))
        tglist[[paste0('nextnext',tiername)]] <- data.frame(xmin=c(NA,NA), xmax=c(NA,NA), text=c(NA,NA))
      }else{
        tglist[[paste0('lastlast',tiername)]] <- data.frame(xmin=as.numeric(newtier[,1]), xmax=as.numeric(newtier[,2]), text=c(NA,NA,newtier[1:(nrow(newtier)-2),3]))
        tglist[[paste0('nextnext',tiername)]] <- data.frame(xmin=as.numeric(newtier[,1]), xmax=as.numeric(newtier[,2]), text=c(newtier[3:nrow(newtier),3], NA, NA))
      }
    }  
  }
  tglist
}


compare_trajectories <- function(aaa_data, sp, plane='sag', signal='TRangle', words=c(), token.col='word_id', main=''){
  
  name_of = list(TRangle='tongue root retraction', TDangle='tongue dorsum height', TTangle='tongue tip height')
  how_many_words = length(words)
  xlab=c(0,1)
  ylab = ifelse(signal%in%names(name_of), name_of[[signal]], signal)
  
  if (plane=='video'){
  	plotdata = aaa_data[[sp]][[plane]]$lip_traces[aaa_data[[sp]][[plane]]$lip_traces$word%in%words,]
  }else{
	plotdata = aaa_data[[sp]][[plane]]$tongue_traces[aaa_data[[sp]][[plane]]$tongue_traces$word%in%words,]
  }
  # print(plotdata)
  word_palette = rainbow(length(words), v=0.6, a=0.8)
  
  # ylim=range(plotdata[,signal]) + c(-5, 0)
  ylim=range(plotdata[,signal]) + diff(range(plotdata[,signal]))*c(-0.1,0)
  xlim=c(0,1)
  
  plot(0, 0, type='n', xlim=xlim, ylim=ylim, xlab='relative time', ylab=ylab, main=main)
  
  for (i in 1:how_many_words){
    w = words[i]
    subdata = plotdata[plotdata[,'word']==w,]
    for (tk in unique(subdata[,token.col])){
      tokendata = subdata[subdata[,token.col]==tk,]
      
      tokendata$relative_time = (tokendata$Time_of_sample_in_recording - min(tokendata$Time_of_sample_in_recording)) / diff(range(tokendata$Time_of_sample_in_recording))
      
      points(tokendata$relative_time, tokendata[,signal], type='l', col=word_palette[i])
    }
  }
  legend('bottomleft', legend=words, col=word_palette[1:how_many_words], lwd=1, ncol=how_many_words, bg='white')
}



add_radial_grid <- function(aaa_data, sp, plane='sag', from=-10, to=170, length=100, flip=FALSE){
  polar_origin = aaa_data[[sp]][[plane]]$origin

  if (flip) polar_origin[2] = -polar_origin[2]

  # print('ORIGIN!')
  # print(polar_origin)
  draw_angles = seq(from,to,10)
  # points(polar_origin[1], polar_origin[2], cex=5)
  segments(rep(polar_origin[1],length(draw_angles)), rep(polar_origin[2],length(draw_angles)), 
           polar_origin[1]-length*cos(draw_angles*pi/180), polar_origin[2]+length*sin(draw_angles*pi/180),
           col='lightgray')
  if ('analysis_value_radii'%in%names(aaa_data[[sp]][[plane]])){
    selected_angles = as.numeric(aaa_data[[sp]][[plane]]$analysis_value_radii[,c('TTangle','TDangle','TRangle')])
    
    segments(rep(polar_origin[1],length(selected_angles)), rep(polar_origin[2],length(draw_angles)), 
             polar_origin[1]-length*cos(selected_angles*pi/180), polar_origin[2]+length*sin(selected_angles*pi/180),
             col='gray', lwd=2)
    # }else{
    # 	print('skipping analysis_value_radii because none were found')
  }
  
  points(polar_origin[1]-10*cos(draw_angles*pi/180), polar_origin[2]+10*sin(draw_angles*pi/180), pch=19,col='white',cex=2)
  text(polar_origin[1]-10*cos(draw_angles*pi/180), polar_origin[2]+10*sin(draw_angles*pi/180), labels=draw_angles, cex=0.5, col='gray')

}



# MAKE A FUNCTION TO CONVERT ARPABET SYMBOLS INTO IPA SYMBOLS FOR PLOTTING. FOR NOW JUST PASTE THIS IN.
# LATER YOU CAN ADD MORE ENTRIES AND REPASTE, IF YOU USE OTHER SOUNDS
arpabet2ipa = function(arpabet){
  ipa_of = list(K='k', B='b', P='p', T='t', S='s', SH='ʃ', CH='tʃ', HH='h', R='ɹ', 
                N='n', T='t', D='d', Z='z', G='ɡ', AA1='ɑ', IH1='ɪ', IY1='i', AE1='æ', 
                UW1='u', ER0='ɚ', ER1='ɜ˞', EH1='ɛ', EY1='e', AH0='ə')
  ipa = c()
  for (a in arpabet){
  	# print(a)
    if (a %in% names(ipa_of)){
      ipa = c(ipa, paste0('/',ipa_of[[a]],'/'))
    }else{
      ipa = c(ipa, a)
      if (is.na(a)){
      	print('missing value for phone')
      }else if (a!=''){
	    print (paste('ipa_of has no entry for',a))
	  }
    }
  }
  ipa
}



read_all_dlc_data <- function(speakers, readwritepath=NULL, ultrasound_model_name=NULL,	lips_model_name=NULL){

	require(tidyr)

	aaa_data = list()
	for (speaker in speakers){
		# print(speaker)
		aaa_data[[speaker]] = list(sag=list(),cor=list(),video=list())
	}

	if (!is.null(ultrasound_model_name)){
		read_csv = paste0(readwritepath, '/', ultrasound_model_name, '.csv')
		print(paste('Reading ultrasound data from',read_csv))
		all_ultrasound_data = read.csv(read_csv, header=FALSE)
		ultrasound_data_names = all_ultrasound_data[2:3,-1]
		ultrasound_filenames = all_ultrasound_data[-c(1:3),1]
		ultrasound_speaker_times_x = separate(all_ultrasound_data[-c(1:3),], col=1, into=c('speaker','ms','x'), sep='_')
		ultrasound_speaker = ultrasound_speaker_times_x[,1]
		ultrasound_times = as.numeric(ultrasound_speaker_times_x[,2])/1000
		ultrasound_datapoints = all_ultrasound_data[-c(1:3),-1]
	}
	###########

	if (!is.null(lips_model_name)){
		read_csv = paste0(readwritepath, '/', lips_model_name, '.csv')
		print(paste('Reading lips data from',read_csv))
		all_video_data = read.csv(read_csv, header=FALSE)
		video_data_names = all_video_data[2:3,-1]
		video_filenames = all_video_data[-c(1:3),1]
		video_speaker_times_x = separate(all_video_data[-c(1:3),], col=1, into=c('speaker','ms','x'), sep='_')
		video_speaker = video_speaker_times_x[,1]
		video_times = as.numeric(video_speaker_times_x[,2])/1000
		video_datapoints = all_video_data[-c(1:3),-1]
	}

	if (!is.null(ultrasound_model_name)){
		for (speaker in speakers){
			print(paste('Processing ultrasound data for',speaker))
			# aaa_data[[speaker]] = list(sag=list(),cor=list(),video=list())

			us_colnames = c()
			tongue_traces = data.frame(image=ultrasound_filenames[ultrasound_speaker==speaker], 
				                       Time_of_sample_in_recording=ultrasound_times[ultrasound_speaker==speaker])
			for (i in 1:ncol(ultrasound_datapoints)){
				data_name = paste(ultrasound_data_names[1,i],ultrasound_data_names[2,i],sep='_')
				us_colnames = c(us_colnames, data_name)
				tongue_traces[,data_name] = as.numeric(paste(ultrasound_datapoints[ultrasound_speaker==speaker,i]))
			}
			# number the tongue points (not doing this for video data)
			tongue_columns = which(grepl('vallecula',names(tongue_traces))|grepl('tongue',names(tongue_traces)))
			tongue_points = length(tongue_columns)/3
			tongue_point_numbers = paste0(rep(c('X','Y','C'),tongue_points), sort(rep(c(1:tongue_points),3)))
			names(tongue_traces)[tongue_columns] = tongue_point_numbers
			
			aaa_data[[speaker]]$sag$tongue_traces = tongue_traces
			aaa_data[[speaker]]$sag$us_colnames = us_colnames
		}
	}

	if (!is.null(lips_model_name)){
		for (speaker in speakers){
			print(paste('Processing lip data for',speaker))
			video_colnames = c()
			lip_traces = data.frame(image=video_filenames[video_speaker==speaker], 
				                    Time_of_sample_in_recording=video_times[video_speaker==speaker])
			for (i in 1:ncol(video_datapoints)){
				data_name = paste(video_data_names[1,i],video_data_names[2,i],sep='_')
				video_colnames = c(video_colnames, data_name)
				lip_traces[,data_name] = as.numeric(paste(video_datapoints[video_speaker==speaker,i]))
			}

			aaa_data[[speaker]]$video$lip_traces = lip_traces
			aaa_data[[speaker]]$video$video_colnames = video_colnames
		}
	}

	aaa_data
}

plot_sample_frames <- function(aaa_data, speakers=NULL, label='', xlim=NULL, ylim=NULL, polar=FALSE, palate=NULL, img=NULL){

	if (is.null(speakers)) speakers=names(aaa_data)

	pdf(paste0('sample_frames_',label,'.pdf'), height=5, width=6, onefile=TRUE)
	for (speaker in speakers){
		print(speaker)
		if ('tongue_traces' %in% names(aaa_data[[speaker]]$sag)){
			if (!length(intersect(img,aaa_data[[speaker]]$sag$tongue_traces$image))){
				sample_trace = aaa_data[[speaker]]$sag$tongue_traces[round(nrow(aaa_data[[speaker]]$sag$tongue_traces)/2),]
			}else{
				sample_trace = aaa_data[[speaker]]$sag$tongue_traces[aaa_data[[speaker]]$sag$tongue_traces$image%in%img,][1,]
			}
			if (polar){
				Xs = as.numeric(sample_trace[,grepl('T[0-9]', names(sample_trace))])
				Ys = as.numeric(sample_trace[,grepl('R[0-9]', names(sample_trace))])
			}else{
				Xs = as.numeric(sample_trace[,grepl('X[0-9]', names(sample_trace))|grepl('_x', names(sample_trace))])
				Ys = as.numeric(sample_trace[,grepl('Y[0-9]', names(sample_trace))|grepl('_y', names(sample_trace))])
			}
			plot(Xs,Ys,col=rainbow(length(Xs)), pch=19, main=paste(speaker, 'sagittal', label), xlim=xlim, ylim=ylim)
			if ('origin' %in% names(aaa_data[[speaker]]$sag)){
				points(aaa_data[[speaker]]$sag$origin[1], aaa_data[[speaker]]$sag$origin[2])
			}
			if (!is.null(palate)){
				points(aaa_data[[speaker]]$sag[[palate]], type='l')
			}
		}

		if ('tongue_traces' %in% names(aaa_data[[speaker]]$cor)){
			if (!length(intersect(img,aaa_data[[speaker]]$cor$tongue_traces$image))){
				sample_trace = aaa_data[[speaker]]$cor$tongue_traces[round(nrow(aaa_data[[speaker]]$cor$tongue_traces)/2),]
			}else{
				sample_trace = aaa_data[[speaker]]$cor$tongue_traces[aaa_data[[speaker]]$cor$tongue_traces$image%in%img,][1,]
			}
			if (polar){
				Xs = as.numeric(sample_trace[,grepl('T[0-9]', names(sample_trace))])
				Ys = as.numeric(sample_trace[,grepl('R[0-9]', names(sample_trace))])
			}else{
				Xs = as.numeric(sample_trace[,grepl('X[0-9]', names(sample_trace))])
				Ys = as.numeric(sample_trace[,grepl('Y[0-9]', names(sample_trace))])
			}
			plot(Xs,Ys,col=rainbow(length(Xs)), pch=19, main=paste(speaker, 'coronal', label), xlim=xlim, ylim=ylim)
			if ('origin' %in% names(aaa_data[[speaker]]$cor)){
				points(aaa_data[[speaker]]$cor$origin[1], aaa_data[[speaker]]$cor$origin[2])
			}
			if (!is.null(palate)){
				points(aaa_data[[speaker]]$cor[[palate]], type='l')
			}
		}
		
		if ('lip_traces' %in% names(aaa_data[[speaker]]$video)){
			if (!length(intersect(img,aaa_data[[speaker]]$video$lip_traces$image))){
				sample_trace = aaa_data[[speaker]]$video$lip_traces[round(nrow(aaa_data[[speaker]]$video$lip_traces)/2),]
			}else{
				sample_trace = aaa_data[[speaker]]$video$lip_traces[aaa_data[[speaker]]$video$lip_traces$image%in%img,][1,]
			}
			Xs = as.numeric(sample_trace[,grepl('_x', names(sample_trace))])
			Ys = as.numeric(sample_trace[,grepl('_y', names(sample_trace))])
			if(nrow(aaa_data[[speaker]]$video$lip_traces)){
				plot(Xs,Ys,col=rainbow(length(Xs)), pch=19, main=paste(speaker, 'lips', label), xlim=xlim, ylim=ylim)
				if ('origin' %in% names(aaa_data[[speaker]]$video)){
					points(aaa_data[[speaker]]$video$origin[1], aaa_data[[speaker]]$video$origin[2])
				}
			}else{
				plot(0,0,type='n', main=paste(speaker, 'lips', label), xlim=xlim, ylim=ylim)				
			}
		}
	}
	dev.off()
}



flip_once <- function(aaa_data, speakers=NULL, planes=NULL, ymax=0){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (speaker in speakers){
		if ('sag' %in% planes & 'tongue_traces' %in% names(aaa_data[[speaker]]$sag)){
			ycols = grepl('Y', names(aaa_data[[speaker]]$sag$tongue_traces))|grepl('_y', names(aaa_data[[speaker]]$sag$tongue_traces))
			aaa_data[[speaker]]$sag$tongue_traces[,ycols] = ymax - aaa_data[[speaker]]$sag$tongue_traces[,ycols]
		}

		if ('cor' %in% planes & 'tongue_traces' %in% names(aaa_data[[speaker]]$cor)){
			ycols = grepl('Y', names(aaa_data[[speaker]]$cor$tongue_traces))|grepl('_y', names(aaa_data[[speaker]]$sag$tongue_traces))
			aaa_data[[speaker]]$cor$tongue_traces[,ycols] = ymax - aaa_data[[speaker]]$cor$tongue_traces[,ycols]
		}
		
		if ('video' %in% planes & 'lip_traces' %in% names(aaa_data[[speaker]]$video)){
			ycols = grepl('_y', names(aaa_data[[speaker]]$video$lip_traces))
			aaa_data[[speaker]]$video$lip_traces[,ycols] = ymax - aaa_data[[speaker]]$video$lip_traces[,ycols]
		}
	}
	aaa_data
}


find_tongue_point_density <- function(speakers=NULL, us_data, xn, yn, scale=1, phones=NULL){

	require(MASS)

	tongue_point_density = list()

	for (sp in speakers){
		traces = us_data[[sp]]$sag$tongue_traces
		if (!is.null(phones)){
			traces = subset(traces, phone %in% phones)
		}
		print(paste(sp, nrow(traces)))
		x=c()
		y=c()
		for (i in 1:11){
		    x = c(x, traces[,paste0('X',i)])
		    y = c(y, traces[,paste0('Y',i)])
		}
		tongue_point_density[[sp]] = kde2d(x*scale, y*scale, n = c(xn*scale, yn*scale))
	}
	tongue_point_density
}
