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

# CHANGE TO YOUR WORKING DIRECTORY
# setwd('C:\\Users\\jrlang2\\Downloads\\phon-main\\Ultrasound Exports_Jacox Lab\\OSP_152\\NEW Pre-Op_ OSP_152_5_17_2022_Pre_op_ultrasound')
# WE'LL READ AND WRITE FROM THE WORKING DIRECTORY (NOT A SUBDIRECTORY)
readwritepath = '/home/jimielke/covart/cr/dlc/'
# readwritepath = '/phon/covart/cr/dlc/'
# scriptpath = '/home/jimielke/scripts/phon/AAA_DLC_R/'
# scriptpath = '/phon/scripts/phon/AAA_DLC_R/'
scriptpath = '~/scripts/phonNCSU/ultrasound/'
# LOAD NECESSARY LIBRARIES AND FUNCTIONS

# IF YOU DON't HAVE gss OR plyr LIBRARIES ALREADY, INSTALL THEM WITH THESE COMMANDS (MINUS THE INITIAL #)
# install.packages('gss')
# install.packages('plyr')

source(paste0(scriptpath,'read_ultrasound_functions.r'))
source(paste0(scriptpath,'tongue_ssanova.r'))

###########################################################################################################################################
# IMPORT AND ORGANIZE YOUR DATA
###########################################################################################################################################

# CHANGE TO YOUR speakername
speakers = c('cr01')
# speakers = c('cr01','cr02','cr03','cr05','cr06','cr07','cr08','cr09','cr10','cr11',
# 	         'cr13','cr14','cr15','cr18','cr19','cr23','cr24','cr28','cr29','cr30',
# 	         'cr32','cr33','cr26','cr04','cr21')
radii_filename = c('cr_analysis_value_radii.csv')

xlim=c(0,320)
ylim=c(0,240)

# LOOP THROUGH ALL YOUR (ONE) SPEAKER(S) AND LOAD THEIR DATA
# THE RESULT IS A COMPLEX LIST STRUCTURE, WITH AN ENTRY FOR EACH SPEAKER, AND EACH PLANE (sagittal, coronal, video, if applicable).

# us_data = read_all_us_data(speakers)

# for (sp in names(us_data)){
# 	us_data[[sp]]$sag$tongue_traces$Annotation_Label = us_data[[sp]]$sag$tongue_traces$word
# 	us_data[[sp]]$sag$tongue_traces$token = us_data[[sp]]$sag$tongue_traces$token_id
# }

###########################################################################################################################################
# STEP 1: read the raw coordinates from the file
us_data = read_all_dlc_data(speakers, readwritepath=readwritepath, ultrasound_model_name = 'all_small_us_pngsDLC_mobnet_100_SpeechProductionFeb12shuffle2_800000',
	lips_model_name = 'all_video_pngsDLC_mobnet_100_Tal_LipsJan28shuffle2_800000')
save(us_data, file='cr_dlc_2024_06_xx.RData')
# load('cr_dlc_2024_06_07.RData')
plot_sample_frames(us_data, label='read_all_dlc_data', xlim=xlim, ylim=ylim)
###########################################################################################################################################

###########################################################################################################################################
# NEW STEP: flip tongue and lip traces
us_data = flip_once(us_data, planes=c('sag','video'), ymax=max(ylim))
plot_sample_frames(us_data, label='flip_once', xlim=xlim, ylim=ylim)
###########################################################################################################################################

###########################################################################################################################################
# READ PALATE TRACES AND OCCLUSAL ANGLES (CREATED MANUALLY SINCE DLC DOESN'T DO THESE WELL)
# NEW STEP: scale and flip palate traces
# STEP 2: rotate only the palate and occlusal plane
us_data = add_manual_pal_occ(us_data, palate_filepath='shiny_palate_traces_2024May20_15h56s52_modified.csv', occlusal_filepath='cr_occlusal_angles.csv',
	palate_scale=0.5)
plot_pal_occ(us_data, planes='sag', xlim=xlim, ylim=ylim, center=c(mean(xlim),mean(ylim)))
plot_sample_frames(us_data, label='read_palate', xlim=xlim, ylim=ylim, palate='palate_trace')
###########################################################################################################################################

###########################################################################################################################################
# STEP 2.9: rotate traces to occlusal plane (overwriting raw traces), and choose a polar origin
us_data = xy2occlusal(us_data, center=c(mean(xlim),mean(ylim)))
us_data = choose_polar_origin(us_data, method='xmean_y025')
plot_sample_frames(us_data, label='choose_polar_origin', xlim=xlim, ylim=ylim, palate='palate_trace_rotated')
###########################################################################################################################################

###########################################################################################################################################
# STEP 3: make polar (storing rotated XY traces in the tongue_traces data frame as TR instead of XY)
us_data = xy2polar(us_data)
plot_sample_frames(us_data, label='xy2polar', polar=TRUE)
###########################################################################################################################################

###########################################################################################################################################
# STEP 2a/3a: read in the selected angles and measure traces
us_data = add_analysis_angles('cr_analysis_value_radii.csv', us_data, angles_only=TRUE)
us_data = measure_traces_at_angles(us_data)
###########################################################################################################################################

###########################################################################################################################################
# READ TEXTGRIDS AND ADD TO THE ARTICULATORY DATA
us_data = add_textgrid_segmentation(us_data, planes=c('sag','video'), tiers=c('phone','word','YPR'), match_words_to_phrase=FALSE)
save(us_data, file='cr_dlc_2024_06_12.RData')
###########################################################################################################################################

###########################################################################################################################################
# DONE IMPORTING! NOW CONTINUE WITH analyze_ultrasound_covart_cr.r
###########################################################################################################################################
