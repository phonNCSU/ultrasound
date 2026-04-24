###########################################################################################################################################
# HERE IS WHAT YOU NEED IN YOUR WORKING DIRECTORY TO GET STARTED (WHERE speakername IS A LABEL YOU CHOSE FOR YOUR RECORDING):
# YOU MAY NEED TO RENAME SOME OF YOUR DATA FILES
#
# export_SESSIONNAME               a folder containing all your AAA-exported txt files and your corresponding MFA-created textgrids
#                                  the textgrid should have word and phone tiers (i.e. they are NOT the minimal textgrids reimported into AAA)
# SESSIONNAME_sag_pal_occ.txt      your palate and occlusal plane traces, exported from AAA
# SESSIONNAME_sag_tongues.txt      your tongue traces, exported from AAA
#
# THESE ARE THE R SCRIPTS INVOLVED (SEE scriptpath BELOW)
#
# read_aaa_export.r                this file, containing commands for reading your export
# plot_aaa_export_functions.r      commands for working with this data after importing it
# read_aaa_export_functions.r      R functions for working with data exported from AAA
# tongue_ssanova.r                 R functions for making SSANOVA comparisons in polar coordinates
#
############################################################################################################################################

###########################################################################################################################################
# PREPARE R FOR YOUR ANALYSIS
###########################################################################################################################################

# IF YOU DON't HAVE gss OR plyr LIBRARIES ALREADY, INSTALL THEM WITH THESE COMMANDS (MINUS THE INITIAL #)
# install.packages('gss')
# install.packages('plyr')

# CHANGE TO YOUR WORKING DIRECTORY
# setwd('C:\\Users\\jrlang2\\Downloads\\phon-main\\Ultrasound Exports_Jacox Lab\\OSP_152\\NEW Pre-Op_ OSP_152_5_17_2022_Pre_op_ultrasound')
# WE'LL READ AND WRITE FROM THE WORKING DIRECTORY (NOT A SUBDIRECTORY)
readwritepath = '/home/jimielke/covart/cr/dlc/data_october_2024'
# readwritepath = '/phon/covart/cr/dlc/'
# scriptpath = '/home/jimielke/scripts/phon/AAA_DLC_R/'
# scriptpath = '/phon/scripts/phon/AAA_DLC_R/'
scriptpath = '~/scripts/phonNCSU/ultrasound/'
# LOAD NECESSARY LIBRARIES AND FUNCTIONS

source(paste0(scriptpath,'read_ultrasound_functions.r'))
source(paste0(scriptpath,'tongue_ssanova.r'))

###########################################################################################################################################
# IMPORT AND ORGANIZE YOUR DATA
###########################################################################################################################################

# CHANGE TO YOUR speakername
# speakers = c('cr01')
speakers = c('cr01','cr02','cr03','cr05','cr06','cr07','cr08','cr09','cr10','cr11',
	         'cr13','cr14','cr15','cr18','cr19','cr23','cr24','cr28','cr29','cr30',
	         'cr32','cr33','cr26','cr04','cr21')

# radii_filename = c('cr_analysis_value_radii.csv')

image_xlim=c(0,320)
image_ylim=c(0,240)

# LOOP THROUGH ALL YOUR (ONE) SPEAKER(S) AND LOAD THEIR DATA
# THE RESULT IS A COMPLEX LIST STRUCTURE, WITH AN ENTRY FOR EACH SPEAKER, AND EACH PLANE (sagittal, coronal, video, if applicable).

# us_data = read_all_us_data(speakers)

# for (sp in names(us_data)){
# 	us_data[[sp]]$sag$tongue_traces$Annotation_Label = us_data[[sp]]$sag$tongue_traces$word
# 	us_data[[sp]]$sag$tongue_traces$token = us_data[[sp]]$sag$tongue_traces$token_id
# }

###########################################################################################################################################
# STEP 1: read the raw coordinates from the file
# us_data = read_all_dlc_data(speakers, readwritepath=readwritepath, ultrasound_model_name = 'all_small_us_pngsDLC_mobnet_100_SpeechProductionFeb12shuffle2_800000',
	# lips_model_name = 'all_video_pngsDLC_mobnet_100_Tal_LipsJan28shuffle2_800000')
us_data = read_all_dlc_data(speakers, readwritepath=readwritepath, ultrasound_model_name = 'all_small_us_pngsDLC_mobnet_100_SpeechProductionFeb12shuffle1_1030000',
	lips_model_name = 'all_video_pngsDLC_mobnet_100_Tal_LipsJan28shuffle1_1030000')
# save(us_data, file='cr_dlc_2025_01_29.RData')
# load('cr_dlc_2025_01_29.RData')
plot_sample_frames(us_data, label='read_all_dlc_data', xlim=image_xlim, ylim=image_ylim)
###########################################################################################################################################

###########################################################################################################################################
# READ TEXTGRIDS AND ADD TO THE ARTICULATORY DATA
us_data = add_textgrid_segmentation(speakers=speakers, us_data, planes=c('sag','video'), tiers=c('phone','word','phrase','quality'), simple_tiers=c('quality'), match_words_to_phrase=FALSE, 
	one_folder='~/covart/cr/dlc/data_october_2024/cr_textgrids')

save(us_data, file='cr_dlc_2025_03_10_step_one_with_textgrids.RData')

for (sp in speakers){
	print(sp)
	quality_phrase = us_data[[sp]]$sag$tongue_traces$quality=='' & us_data[[sp]]$sag$tongue_traces$phrase!=''
	us_data[[sp]]$sag$tongue_traces = us_data[[sp]]$sag$tongue_traces[quality_phrase,]
	quality_phrase = us_data[[sp]]$video$lip_traces$quality=='' & us_data[[sp]]$video$lip_traces$phrase!=''
	us_data[[sp]]$video$lip_traces = us_data[[sp]]$video$lip_traces[quality_phrase,]
	us_data[[sp]]$sag$tongue_traces$quality = NULL
	us_data[[sp]]$video$lip_traces$quality = NULL
}

save(us_data, file='cr_dlc_2025_03_10_step_one_with_textgrids_quality_phrases.RData')
# load('cr_dlc_2025_03_10_step_one_with_textgrids_quality_phrases.RData')
###########################################################################################################################################



###########################################################################################################################################
# STEP xxx: find tongue density to help palate tracing
# stops = c('D','G','K','N','NG','T')
# fricatives = c('DH','S','SH','TH','Z','ZH')
# high_vowels = c('ER1','IY1','R','UW1','W','Y')
# tongue_point_density = find_tongue_point_density(setdiff(names(us_data),'cr15'), us_data, xn=xlim[2], yn=ylim[2], phones=stops)
# save(tongue_point_density, file='CR_all_tongue_point_density_stops_only.RData')
# tongue_point_density = find_tongue_point_density(names(us_data), us_data, xn=xlim[2], yn=ylim[2], phones=c(stops,fricatives))
# save(tongue_point_density, file='CR_all_tongue_point_density_stops_and_fricatives_only.RData')
# tongue_point_density = find_tongue_point_density(names(us_data), us_data, xn=xlim[2], yn=ylim[2], phones=c(stops,fricatives,high_vowels))
# save(tongue_point_density, file='CR_all_tongue_point_density_stops_fricatives_and_high_vowels_only.RData')
###########################################################################################################################################


###########################################################################################################################################
# NEW STEP: flip tongue and lip traces
us_data = flip_once(us_data, planes=c('sag','video'), ymax=max(image_ylim))
plot_sample_frames(us_data, label='flip_once', xlim=image_xlim, ylim=image_ylim)
###########################################################################################################################################


###########################################################################################################################################
# READ PALATE TRACES AND OCCLUSAL ANGLES (CREATED MANUALLY SINCE DLC DOESN'T DO THESE WELL)
# NEW STEP: scale and flip palate traces
# STEP 2: rotate only the palate and occlusal plane
# us_data = add_manual_pal_occ(us_data, palate_filepath=paste0(readwritepath,'/','shiny_palate_traces_2025Feb20_14h27s49.csv'), 
# 	occlusal_filepath='cr_occlusal_angles.csv', palate_scale=1)
source(paste0(scriptpath,'read_ultrasound_functions.r'))
us_data = add_manual_pal_occ(us_data, palate_filepath=paste0(readwritepath,'/','shiny_palate_traces_2025Feb20_14h27s49.csv'), 
	occlusal_filepath=paste0(readwritepath,'/','shiny_occlusal_traces_2025Mar13_16h43s43.csv'), measure_occlusal=TRUE, palate_scale=1)

	
# plot_pal_occ(us_data, planes='sag', xlim=image_xlim, ylim=image_ylim, center=c(mean(image_xlim),mean(image_ylim)))
plot_sample_frames(us_data, label='read_palate', xlim=image_xlim, ylim=image_ylim, palate='palate_trace')
###########################################################################################################################################


###########################################################################################################################################
# SCALE TO CM, ALL X AND Y COLUMNS AND PALATE
us_data = scale_to_cm(us_data, depths_filepath='/home/jimielke/covart/cr/dlc/data_october_2024/cr_depth.csv')

xlim_ultrasound = c(-5, 10)
ylim_ultrasound = c(-3.1, 6.8)
xlim_video = c(-5.9, 6.4)
ylim_video = c(-4.6, 3.6)


source(paste0(scriptpath,'read_ultrasound_functions.r'))
plot_sample_frames(us_data, label='scaled_to_cm', xlim_ultrasound=xlim_ultrasound, ylim_ultrasound=ylim_ultrasound, xlim_video=xlim_video, ylim_video=ylim_video, palate='palate_trace')
###########################################################################################################################################
       


###########################################################################################################################################
# STEP 3: make polar (storing rotated XY traces in the tongue_traces data frame as TR instead of XY)
us_data = choose_polar_origin(us_data, method='shortTendon')
us_data = xy2polar(us_data)
plot_sample_frames(us_data, label='xy2polar', xlim_ultrasound=xlim_ultrasound, ylim_ultrasound=ylim_ultrasound, xlim_video=xlim_video, ylim_video=ylim_video, polar=TRUE)
###########################################################################################################################################

###########################################################################################################################################
# STEP 2a/3a: read in the selected angles and measure traces
# us_data = add_analysis_angles('cr_analysis_value_radii.csv', us_data, angles_only=TRUE)
# us_data = measure_traces_at_angles(us_data)
###########################################################################################################################################

###########################################################################################################################################
# STEP xxx: make lip signals
us_data = measure_lips(us_data)
###########################################################################################################################################



save(us_data, file='cr_dlc_2025_04_02.RData')

###########################################################################################################################################
# DONE IMPORTING! NOW CONTINUE WITH analyze_ultrasound_covart_cr.r
###########################################################################################################################################

# widths = c()
# for (sp in names(us_data)){
# 	widths = c(widths, as.numeric(quantile(us_data[[sp]]$video$lip_traces$lips_horiz,0.5,na.rm=TRUE)))
# }
# mean(widths, na.rm=TRUE)
# mean 117.517 pixels

# mean 50.1 mm
# 50.1 / 117.517 mmperpixels 0.4263213

# 95% 135.3286
# 3D Stereophotogrammetry Quantitative Lip Analysis Sawyer et al 2009






















###########################################################################################################################################
# STEP 2.9: rotate traces to occlusal plane (overwriting raw traces), and choose a polar origin
# us_data = xy2occlusal(us_data, center=c(mean(image_xlim),mean(image_ylim)), rotate_palate=TRUE)

# plot_sample_frames(us_data, label='rotate_to_occlusal', xlim=image_xlim, ylim=image_ylim, palate='palate_trace')	
   
# surprisingly low: cr01, cr07, cr08, cr28, cr30, cr04, cr21
# surprisingly high: cr03, cr10, cr14, cr23, cr24, cr29, cr26

# us_data = choose_polar_origin(us_data, method='xmean_y025')
# us_data = choose_polar_origin(us_data, method='shortTendon')
# plot_sample_frames(us_data, label='choose_polar_origin', xlim_ultrasound=xlim_ultrasound, ylim_ultrasound=ylim_ultrasound, xlim_video=xlim_video, ylim_video=ylim_video, palate='palate_trace_rotated')
###########################################################################################################################################



# > mean(us_data$cr03$sag$tongue_traces$X6) - mean(us_data$cr03$sag$tongue_traces$hyoid_x)
# [1] 81.75715


# SCALE TO CM HERE, ALL X AND Y COLUMNS AND PALATE
# depths = read.csv('/home/jimielke/covart/cr/dlc/data_october_2024/cr_depth.csv')
# for (sp in speakers){
# 	print(sp)

# 	tongue_xcols = names(us_data[[sp]]$sag$tongue_traces)[grepl('X[0-9]', names(us_data[[sp]]$sag$tongue_traces))|grepl('_x', names(us_data[[sp]]$sag$tongue_traces))]
# 	tongue_ycols = names(us_data[[sp]]$sag$tongue_traces)[grepl('Y[0-9]', names(us_data[[sp]]$sag$tongue_traces))|grepl('_y', names(us_data[[sp]]$sag$tongue_traces))]

# 	lips_xcols   = names(us_data[[sp]]$video$lip_traces)[grepl('X[0-9]', names(us_data[[sp]]$video$lip_traces))|grepl('_x', names(us_data[[sp]]$video$lip_traces))]
# 	lips_ycols   = names(us_data[[sp]]$video$lip_traces)[grepl('Y[0-9]', names(us_data[[sp]]$video$lip_traces))|grepl('_y', names(us_data[[sp]]$video$lip_traces))]

# 	usmmpermixel = depths[depths$speaker==sp,'usmmpermixel']
# 	videommpermixel = depths[depths$speaker==sp,'videommpermixel']

# 	# ULTRASOUND ORIGIN: MEAN HYOID POSITION
# 	us_data[[sp]]$sag$tongue_traces[,tongue_xcols] = (us_data[[sp]]$sag$tongue_traces[,tongue_xcols] - mean(us_data[[sp]]$sag$tongue_traces[,'hyoid_x'])) * usmmpermixel/10
# 	us_data[[sp]]$sag$tongue_traces[,tongue_ycols] = (us_data[[sp]]$sag$tongue_traces[,tongue_ycols] - mean(us_data[[sp]]$sag$tongue_traces[,'hyoid_y'])) * usmmpermixel/10

# 	us_data[[sp]]$sag$palate_trace$X = (us_data[[sp]]$sag$palate_trace$X - mean(us_data[[sp]]$sag$tongue_traces[,'hyoid_x'])) * usmmpermixel/10
# 	us_data[[sp]]$sag$palate_trace$Y = (us_data[[sp]]$sag$palate_trace$Y - mean(us_data[[sp]]$sag$tongue_traces[,'hyoid_y'])) * usmmpermixel/10

# 	us_data[[sp]]$sag$palate_trace_rotated$X = (us_data[[sp]]$sag$palate_trace_rotated$X - mean(us_data[[sp]]$sag$tongue_traces[,'hyoid_x'])) * usmmpermixel/10
# 	us_data[[sp]]$sag$palate_trace_rotated$Y = (us_data[[sp]]$sag$palate_trace_rotated$Y - mean(us_data[[sp]]$sag$tongue_traces[,'hyoid_y'])) * usmmpermixel/10

# 	# VIDEO ORIGIN: MEAN LIP DOT POSITION
# 	us_data[[sp]]$video$lip_traces[,lips_xcols] = (us_data[[sp]]$video$lip_traces[,lips_xcols] - mean(unlist(us_data[[sp]]$video$lip_traces[,lips_xcols]))) * videommpermixel/10
# 	us_data[[sp]]$video$lip_traces[,lips_ycols] = (us_data[[sp]]$video$lip_traces[,lips_ycols] - mean(unlist(us_data[[sp]]$video$lip_traces[,lips_ycols]))) * videommpermixel/10

# }

# xlim_ultrasound = c(-4,9.5)
# ylim_ultrasound = c(-2,7)
# xlim_video = c(-4.5,4.5)
# ylim_video = c(-3,3)

# # source(paste0(scriptpath,'read_ultrasound_functions.r'))
# plot_sample_frames(us_data, label='scaled_to_cm', xlim_ultrasound=xlim_ultrasound, ylim_ultrasound=ylim_ultrasound, xlim_video=xlim_video, ylim_video=ylim_video,show_grid=FALSE, palate='palate_trace_rotated')
       
# for (sp in speakers){
# 	print(sp)

# 	tongue_xcols = names(us_data[[sp]]$sag$tongue_traces)[grepl('X[0-9]', names(us_data[[sp]]$sag$tongue_traces))|grepl('_x', names(us_data[[sp]]$sag$tongue_traces))]
# 	tongue_ycols = names(us_data[[sp]]$sag$tongue_traces)[grepl('Y[0-9]', names(us_data[[sp]]$sag$tongue_traces))|grepl('_y', names(us_data[[sp]]$sag$tongue_traces))]

# 	lips_xcols   = names(us_data[[sp]]$video$lip_traces)[grepl('X[0-9]', names(us_data[[sp]]$video$lip_traces))|grepl('_x', names(us_data[[sp]]$video$lip_traces))]
# 	lips_ycols   = names(us_data[[sp]]$video$lip_traces)[grepl('Y[0-9]', names(us_data[[sp]]$video$lip_traces))|grepl('_y', names(us_data[[sp]]$video$lip_traces))]

# 	# print(range(unlist(us_data[[sp]]$sag$tongue_traces[,tongue_xcols])))
# 	# print(range(unlist(us_data[[sp]]$sag$tongue_traces[,tongue_ycols])))
# 	# print(range(unlist(us_data[[sp]]$video$lip_traces[,lips_xcols])))
# 	# print(range(unlist(us_data[[sp]]$video$lip_traces[,lips_ycols])))
# }