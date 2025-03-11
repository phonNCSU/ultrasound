###########################################################################################################################################
# HERE IS WHAT YOU NEED IN YOUR WORKING DIRECTORY TO GET STARTED (WHERE speakername IS A LABEL YOU CHOSE FOR YOUR RECORDING):
# YOU MAY NEED TO RENAME SOME OF YOUR DATA FILES
#





# export_cr01					a folder containing the textgrid for cr01
# all_small_us_pngsDLC_mobnet_100_SpeechProductionFeb12shuffle1_1030000.csv
# all_video_pngsDLC_mobnet_100_Tal_LipsJan28shuffle1_1030000.csv
# K25_analysis_value_radii.csv
# K25_occlusal_angles.csv
# shiny_palate_traces_fake.csv
#
# THESE ARE THE R SCRIPTS INVOLVED (SEE scriptpath BELOW)
#
# read_dlc_export.r                this file, containing commands for reading your export
# read_ultrasound_functions.r      R functions for working with data exported from DLC
# tongue_ssanova.r                 R functions for making SSANOVA comparisons in polar coordinates
#
############################################################################################################################################

###########################################################################################################################################
# PREPARE R FOR YOUR ANALYSIS
###########################################################################################################################################

# CHANGE TO YOUR WORKING DIRECTORY
# setwd('C:\\Users\\jrlang2\\Downloads\\phon-main\\Ultrasound Exports_Jacox Lab\\OSP_152\\NEW Pre-Op_ OSP_152_5_17_2022_Pre_op_ultrasound')
# WE'LL READ AND WRITE FROM THE WORKING DIRECTORY (NOT A SUBDIRECTORY)
readwritepath = '/home/jimielke/Kalasha/analysis_2025/dlc'
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
# speakers = c('Kal1')
speakers = c('Dame1', 'Dame2', 'Dame4', 'Dame5', 'Dame6', 'Kal1', 'Kal4', 'Kal5', 'Kal8', 
	'KalBi1', 'KalBi2', 'KalBi3', 'KalBi4', 'KalBi5', 'KalBi9', 'KalBu1', 'KalBu2', 'KalBu3', 
	'KalBu4', 'KalBu5', 'KalBu6', 'KalBu8', 'KalBu12', 'KalBu13', 'KalRu1', 'KalRu2', 
	'KalRu3', 'KalRu4', 'KalRu5', 'KalRu6', 'KalRu7', 'KalRu9', 'KalRu10', 'Kami1', 'Kami2', 
	'Kami3', 'Kati1', 'Kati2', 'Kati3', 'Kati4', 'Kati7', 'Khow1', 'Khow2', 'Khow3', 'Khow4', 
	'Palu1', 'Palu2', 'Palu3', 'Palu4', 'Shin1', 'Shin2', 'Shin3', 'Shin4', 'Shin5')

radii_filename = c('K25_analysis_value_radii.csv')

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
# us_data = read_all_dlc_data(speakers, readwritepath=readwritepath, ultrasound_model_name = 'all_small_us_pngsDLC_mobnet_100_SpeechProductionFeb12shuffle2_800000',
	# lips_model_name = 'all_video_pngsDLC_mobnet_100_Tal_LipsJan28shuffle2_800000')
us_data = read_all_dlc_data(speakers, readwritepath=readwritepath, ultrasound_model_name = 'all_small_us_pngsDLC_mobnet_100_SpeechProductionFeb12shuffle1_1030000',
	lips_model_name = 'all_video_pngsDLC_mobnet_100_Tal_LipsJan28shuffle1_1030000')
save(us_data, file='K25_dlc_2025_01_20.RData')
# load('K25_dlc_2025_01_20.RData')
plot_sample_frames(us_data, label='read_all_dlc_data', xlim=xlim, ylim=ylim)
###########################################################################################################################################


tongue_point_density = find_tongue_point_density(speakers, us_data, xn=xlim[2], yn=ylim[2])

save(tongue_point_density, file='K25_tongue_point_density.RData')
# plot(0, 0, type='n', xlim=c(0,320), ylim=c(0,240))
# contour(z, lwd = 2, add = TRUE, col = hcl.colors(10, "Spectral"), nlevels=50)
# dev.off()

# filled.contour(tongue_point_density$Kal1)

###########################################################################################################################################
# NEW STEP: flip tongue and lip traces
us_data = flip_once(us_data, planes=c('sag','video'), ymax=max(ylim))
plot_sample_frames(us_data, label='flip_once', xlim=xlim, ylim=ylim)
###########################################################################################################################################



###########################################################################################################################################
# READ PALATE TRACES AND OCCLUSAL ANGLES (CREATED MANUALLY SINCE DLC DOESN'T DO THESE WELL)
# NEW STEP: scale and flip palate traces
# STEP 2: rotate only the palate and occlusal plane
us_data = add_manual_pal_occ(us_data, palate_filepath='shiny_palate_traces_fake.csv', occlusal_filepath='K25_occlusal_angles.csv',
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
us_data = add_analysis_angles('K25_analysis_value_radii.csv', us_data, angles_only=TRUE)
us_data = measure_traces_at_angles(us_data)
###########################################################################################################################################

###########################################################################################################################################
# READ TEXTGRIDS AND ADD TO THE ARTICULATORY DATA
source(paste0(scriptpath,'read_ultrasound_functions.r'))
us_data = add_textgrid_segmentation(us_data, planes=c('sag','video'), tiers=c('phone','word'), match_words_to_phrase=FALSE, one_folder='textgrids')




save(us_data, file='K25_dlc_2025_01_20b.RData')
###########################################################################################################################################

###########################################################################################################################################
# DONE IMPORTING! NOW CONTINUE WITH analyze_ultrasound_covart_cr.r
###########################################################################################################################################



#  LOOK FOR PALATE SHINY APP, OCCLUSAL PLANE?
#  x find palate app
#    add occlusal angle


# process all of the Kalasha data and put in shiny app
# x rotate in app


plot(0,0,type='n',xlim=xlim,ylim=ylim)
for (i in 1:11){
        points(us_data$Kal5$sag$tongue_traces[,paste0(c('X','Y'),i)],cex=0.1,col=rainbow(11,a=0.2,v=0.5)[i])
}

plot(0,0,type='n',xlim=c(pi/2,7/4*pi),ylim=c(0,150))
for (i in 1:11){
        points(us_data$Kal1$sag$tongue_traces[,paste0(c('T','R'),i)],cex=0.1,col=rainbow(11,a=0.2,v=0.5)[i])
}

png('tongue_migraine.png')
# install.packages("MASS")
library(MASS)

# Data
x=c()
y=c()
for (i in 1:11){
    x = c(x, us_data$Kal1$sag$tongue_traces[,paste0('X',i)])
    y = c(y, us_data$Kal1$sag$tongue_traces[,paste0('Y',i)])
}

z <- kde2d(x, y, n = 50)
z <- kde2d(x, y, n = c(320,240))
plot(0, 0, type='n', xlim=c(0,320), ylim=c(0,240))
contour(z, lwd = 2, add = TRUE, col = hcl.colors(10, "Spectral"), nlevels=50)


dev.off()



mkdir textgrids
cp /phon/Kalasha/manual_correction/Qandeel_corrected/Kal1.TextGrid textgrids
cp /phon/Kalasha/Kalasha_Aero_Ultrasound/Kal4/Kal4_Ultrasound/Kal4.TextGrid textgrids
cp /phon/Kalasha/Kalasha_Aero_Ultrasound/Kal5/Kal5_Ultrasound/Kal5.TextGrid textgrids
cp /phon/Kalasha/Kalasha_Aero_Ultrasound/Kal8/Kal8_Ultrasound/Kal8.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu1.TextGrid textgrids
cp /phon/Kalasha/Align_fieldwork_2018/Qandeel_corrected/KalBu2.TextGrid textgrids
cp /phon/Kalasha/Align_fieldwork_2018/KalBu3.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu4.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu5.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu6.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu8.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu10.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu12.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu13.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu14.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu17.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu18.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu19.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu20.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu21.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu22.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBu23.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalRu1.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalRu2.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalRu3.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalRu4.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalRu5.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalRu6.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalRu7.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalRu9.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalRu10.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBi1.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBi2.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBi3.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBi4.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBi5.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Aligned_Kalasha/KalBi9.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Nuristani/Ultrasound/Kami1/Kami1.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Nuristani/Ultrasound/Kami2/Kami2.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Nuristani/Ultrasound/Kami3/Kami3.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Nuristani/Ultrasound/Kati1/Kati1.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Nuristani/Ultrasound/Kati2/Kati2.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Nuristani/Ultrasound/Kati3/Kati3.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Nuristani/Ultrasound/Kati4/Kati4.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Nuristani/Ultrasound/Kati7/Kati7.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Dameli/Dame1_ultrasound/Dame1.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Dameli/Dame2_ultrasound/Dame2.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Dameli/Dame4_ultrasound/Dame4.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Dameli/Dame5_ultrasound/Dame5.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Dameli/Dame6_ultrasound/Dame6.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Palula/Ultrasound/Palu1_ultrasound/Palu1.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Palula/Ultrasound/Palu2_ultrasound/Palu2.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Palula/Ultrasound/Palu3_ultrasound/Palu3.TextGrid textgrids
cp /phon/Kalasha/Fieldwork2_Pakistan_2018/Palula/Ultrasound/Palu4_ultrasound/Palu4.TextGrid textgrids
cp /phon/Kalasha/manual_correction/Qandeel_corrected/Khow1.TextGrid textgrids
cp /phon/Kalasha/manual_correction/Qandeel_corrected/Khow2.TextGrid textgrids
cp /phon/Kalasha/manual_correction/Qandeel_corrected/Khow3.TextGrid textgrids
cp /phon/Kalasha/manual_correction/Qandeel_corrected/Khow4.TextGrid textgrids
cp /phon/Kalasha/manual_correction/Qandeel_corrected/Shin1.TextGrid textgrids
cp /phon/Kalasha/manual_correction/Qandeel_corrected/Shin2.TextGrid textgrids
cp /phon/Kalasha/manual_correction/Qandeel_corrected/Shin3.TextGrid textgrids
cp /phon/Kalasha/manual_correction/Qandeel_corrected/Shin4.TextGrid textgrids
cp /phon/Kalasha/manual_correction/Qandeel_corrected/Shin5.TextGrid textgrids
