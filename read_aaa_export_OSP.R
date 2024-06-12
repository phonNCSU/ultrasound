###########################################################################################################################################
# HERE IS WHAT YOU NEED IN YOUR WORKING DIRECTORY TO GET STARTED (WHERE speakername IS A LABEL YOU CHOSE FOR YOUR RECORDING):
# YOU MAY NEED TO RENAME SOME OF YOUR DATA FILES
#
# export_speakername               a folder containing all your AAA-exported txt files and your corresponding MFA-created textgrids
# speakername_sag_pal_occ.txt      your palate and occlusal plane traces, exported from AAA
# speakername_sag_tongues.txt      your tongue traces, exported from AAA
# read_aaa_export_for_lab6.r       this file
# read_aaa_export_functions.r      R functions for working with data exported from AAA
# tongue_ssanova.r                 R functions for making SSANOVA comparisons in polar coordinates
###########################################################################################################################################

###########################################################################################################################################
# PREPARE R FOR YOUR ANALYSIS
###########################################################################################################################################

# CHANGE TO YOUR WORKING DIRECTORY
# setwd('C:/Users/griff/jacox_speech_new/')
# WE'LL READ AND WRITE FROM THE WORKING DIRECTORY (NOT A SUBDIRECTORY)
readwritepath = '/home/jimielke/orthognathic/OSP_ultrasound/jacox_speech_new/data/'
scriptpath = '/home/jimielke/scripts/phonNCSU/ultrasound/'

# IF YOU DON't HAVE gss OR plyr LIBRARIES ALREADY, INSTALL THEM WITH THESE COMMANDS (MINUS THE INITIAL #)
# install.packages('gss')
# install.packages('plyr')

source(paste0(scriptpath,'read_ultrasound_functions.r'))
source(paste0(scriptpath,'tongue_ssanova.r'))

###########################################################################################################################################
# IMPORT AND ORGANIZE YOUR DATA
###########################################################################################################################################

# CHANGE TO YOUR speakernames
speakers = c('OSP_161_10_11_2022_2_Mo_Post_Op', 
             'OSP_163_8_23_2022_pre_op', 
             'OSP_163_12_20_2022_2_mo_Take_2',
             'OSP_163_8_22_2023_1_Yr_Post_Op', 
             'OSP_164_8_30_2022_Pre-Op', 
             'OSP_167_09_30_2022_Pre_Op_Take2',
             'OSP_167_3_24_2023_3_mo', 
             'OSP_167_11_16_2023_1_Yr_Post_Op', 
             'OSP_170_1_23_2024_1_Yr_Post_Op')#, 
             # 'OSP_170_10_11_2022_Pre_Op')

xlim=c(0,160)
ylim=c(0,120)

###########################################################################################################################################
# STEP 1: read the tongue, palate, and occlusal traces from files
# LOOP THROUGH ALL YOUR (ONE) SPEAKER(S) AND LOAD THEIR EXPORTED DATA
# THE RESULT IS KIND OF A COMPLEX LIST STRUCTURE, WITH AN ENTRY FOR EACH SPEAKER, AND EACH PLANE (sagittal and coronal, if applicable).
us_data = read_all_aaa_data(speakers, center=c(mean(xlim),mean(ylim)))
plot_sample_frames(us_data, label='read_all_aaa_data', xlim=xlim, ylim=ylim)
###########################################################################################################################################

###########################################################################################################################################
# PLOT THE OCCLUSAL ANGLE AND SHOW THE ROTATION OF THE PALATE TRACE
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
us_data = add_analysis_angles('OSP_analysis_value_radii.csv', us_data, angles_only=TRUE)
us_data = measure_traces_at_angles(us_data)
###########################################################################################################################################

###########################################################################################################################################
# READ TEXTGRIDS AND ADD TO THE ARTICULATORY DATA
us_data = add_textgrid_segmentation(us_data, planes=c('sag'), tiers=c('phone','word'), match_words_to_phrase=FALSE, middle_is_quarter=c('T','K','CH'))
# ADD THE IPA COLUMN
for (speaker in speakers){
  us_data[[speaker]]$sag$tongue_traces$ipa = factor(arpabet2ipa(us_data[[speaker]]$sag$tongue_traces$phone))
}
save(us_data, file='OSP_2024_06_12.RData')
###########################################################################################################################################
