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
# readwritepath = '/home/jimielke/orthognathic/OSP_ultrasound/data/'
readwritepath = '/home/jimielke/orthognathic/OSP_ultrasound/data/'
scriptpath = '/home/jimielke/scripts/phon/AAA_DLC_R/'

# LOAD NECESSARY LIBRARIES AND FUNCTIONS

# IF YOU DON't HAVE gss OR plyr LIBRARIES ALREADY, INSTALL THEM WITH THESE COMMANDS (MINUS THE INITIAL #)
# install.packages('gss')
# install.packages('plyr')

source(paste0(scriptpath,'read_aaa_export_functions.r'))
source(paste0(scriptpath,'tongue_ssanova.r'))

###########################################################################################################################################
# IMPORT AND ORGANIZE YOUR DATA
###########################################################################################################################################

# CHANGE TO YOUR speakername
# speakers = c('OSP_152_5_17_2022_Pre_op_ultrasound', 'OSP_140_11_9_2021_ultrasound_pre-op')
speakers = c('OSP_140_3_8_2022_3_mo')
# speakers = c('OSP_161_10_11_2022_2_Mo_Post_Op')
# speakers = c('OSP_164_8_30_2022_Pre-Op')

radii_filename = c('OSP_analysis_value_radii.csv')

# LOOP THROUGH ALL YOUR (ONE) SPEAKER(S) AND LOAD THEIR EXPORTED DATA
# THE RESULT IS KIND OF A COMPLEX LIST STRUCTURE, WITH AN ENTRY FOR EACH SPEAKER, AND EACH PLANE (sagittal and coronal, if applicable).
aaa_data = read_all_aaa_data(speakers)

# PLOT THE OCCLUSAL ANGLE AND SHOW THE ROTATION OF THE PALATE TRACE
plot_pal_occ(aaa_data, planes='sag')

# ADD YOUR TEXTGRID SEGMENTATION TO THE DATA FRAME CONTAINING THE TONGUE TRACES
# TEXTGRIDS ARE EXPECTED TO BE IN A FOLDER NAMED export_speakername (WITH YOUR ACTUAL SPEAKER NAME) WITHIN THE WORKING DIRECTORY
aaa_data = add_textgrid_segmentation(aaa_data, merge_vl=FALSE, match_words_to_phrase=FALSE, middle_is_quarter=c('P','T','K','CH'))

# PREPARE THE MIDDLE FRAME OF EACH SEGMENT FOR SSANOVA COMPARISON
# aaa_data = put_middle_frames_in_long_format(aaa_data)

# CHOOSE A GOOD ORIGIN FOR OUR POLAR COORDINATE SYSTEM
aaa_data = choose_polar_origin(aaa_data, 'xmid_ymax', flip=FALSE)

# MEASURE TRACES AT SELECTED ANGLES
aaa_data = add_analysis_angles(radii_filename, aaa_data, angles_only=TRUE)
aaa_data = xy2polar(aaa_data)
aaa_data = measure_traces_at_angles(aaa_data)

###########################################################################################################################################
# DONE IMPORTING! NOW CONTINUE WITH plot_aaa_export.r
###########################################################################################################################################


