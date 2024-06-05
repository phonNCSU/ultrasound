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
# readwritepath = '/home/jimielke/covart/cr/dlc/'
readwritepath = '/phon/covart/cr/dlc/'
# scriptpath = '/home/jimielke/scripts/phon/AAA_DLC_R/'
scriptpath = '/phon/scripts/phon/AAA_DLC_R/'

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
# speakers = c('cr01')
speakers = c('cr01','cr02','cr03','cr05','cr06','cr07','cr08','cr09','cr10','cr11',
	         'cr13','cr14','cr15','cr18','cr19','cr23','cr24','cr28','cr29','cr30',
	         'cr32','cr33','cr26','cr04','cr21')
radii_filename = c('cr_analysis_value_radii.csv')

# LOOP THROUGH ALL YOUR (ONE) SPEAKER(S) AND LOAD THEIR EXPORTED DATA
# THE RESULT IS KIND OF A COMPLEX LIST STRUCTURE, WITH AN ENTRY FOR EACH SPEAKER, AND EACH PLANE (sagittal and coronal, if applicable).

# aaa_data = read_all_aaa_data(speakers)


library(tidyr)

ultrasound_model_name = 'all_small_us_pngsDLC_mobnet_100_SpeechProductionFeb12shuffle2_800000'
lips_model_name = 'all_video_pngsDLC_mobnet_100_Tal_LipsJan28shuffle2_800000'
all_ultrasound_data = read.csv(paste0(readwritepath, '/', ultrasound_model_name, '.csv'), header=FALSE)
all_video_data = read.csv(paste0(readwritepath, '/', lips_model_name, '.csv'), header=FALSE)

aaa_data = list()


ultrasound_data_names = all_ultrasound_data[2:3,-1]
ultrasound_filenames = all_ultrasound_data[-c(1:3),1]
ultrasound_speaker_times_x = separate(all_ultrasound_data[-c(1:3),], col=1, into=c('speaker','ms','x'), sep='_')
ultrasound_speaker = ultrasound_speaker_times_x[,1]
ultrasound_times = as.numeric(ultrasound_speaker_times_x[,2])/1000
ultrasound_datapoints = all_ultrasound_data[-c(1:3),-1]


###########

video_data_names = all_video_data[2:3,-1]
video_filenames = all_video_data[-c(1:3),1]
video_speaker_times_x = separate(all_video_data[-c(1:3),], col=1, into=c('speaker','ms','x'), sep='_')
video_speaker = video_speaker_times_x[,1]
video_times = as.numeric(video_speaker_times_x[,2])/1000
video_datapoints = all_video_data[-c(1:3),-1]



for (speaker in speakers){
	print(speaker)
	aaa_data[[speaker]] = list(sag=list(),cor=list(),video=list())

	us_colnames = c()
	tongue_traces = data.frame(image=ultrasound_filenames[ultrasound_speaker==speaker], Time_of_sample_in_recording=ultrasound_times[ultrasound_speaker==speaker])
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

	video_colnames = c()
	lip_traces = data.frame(image=video_filenames[video_speaker==speaker], Time_of_sample_in_recording=video_times[video_speaker==speaker])
	for (i in 1:ncol(video_datapoints)){
		data_name = paste(video_data_names[1,i],video_data_names[2,i],sep='_')
		video_colnames = c(video_colnames, data_name)
		lip_traces[,data_name] = as.numeric(paste(video_datapoints[video_speaker==speaker,i]))
	}

	aaa_data[[speaker]]$sag$tongue_traces = tongue_traces
	aaa_data[[speaker]]$sag$us_colnames = us_colnames

	aaa_data[[speaker]]$video$lip_traces = lip_traces
	aaa_data[[speaker]]$video$video_colnames = video_colnames
}

# save(aaa_data, file='cr_dlc.RData')

aaa_data = add_textgrid_segmentation(aaa_data, planes=c('sag','video'), tiers=c('phone','word','YPR'), match_words_to_phrase=FALSE)

save(aaa_data, file='cr_dlc_2024_05_20.RData')


######################

readwritepath = '~/covart/cr/dlc/'
scriptpath = '~/scripts/phon/AAA_DLC_R/'

source(paste0(scriptpath,'read_aaa_export_functions.r'))
source(paste0(scriptpath,'tongue_ssanova.r'))

load('cr_dlc_2024_05_20.RData')

speakers=names(aaa_data)

all_occlusal_angles = read.csv('cr_occlusal_angles.csv')
all_palate_traces = read.csv('shiny_palate_traces_2024May20_15h56s52.csv')

for (sp in names(aaa_data)){
	print(sp)
	aaa_data[[sp]]$sag$occlusal_angle = all_occlusal_angles[all_occlusal_angles$speaker==sp,'occlusal_angle']
	palate_trace = all_palate_traces[all_palate_traces$speaker==sp,c('x','y')] * 0.5
	names(palate_trace) = c('X','Y')
	aaa_data[[sp]]$sag$palate_trace = palate_trace
	aaa_data[[sp]]$sag$palate_trace_rotated = rotateXY(aaa_data[[sp]]$sag$palate_trace, aaa_data[[sp]]$sag$occlusal_angle, center=c(0,0))
	aaa_data[[sp]]$sag$tongue_traces$Annotation_Label = aaa_data[[sp]]$sag$tongue_traces$word
	aaa_data[[sp]]$sag$tongue_traces$token = aaa_data[[sp]]$sag$tongue_traces$token_id
}


aaa_data = choose_polar_origin(aaa_data, 'xmean_ymin')

aaa_data = add_analysis_angles('cr_analysis_value_radii.csv', aaa_data, angles_only=TRUE)
aaa_data = xy2polar(aaa_data)
aaa_data = measure_traces_at_angles(aaa_data)

save(aaa_data, file='cr_dlc_2024_05_22.RData')

speaker = names(aaa_data)[1]

data = traces_wide_to_long(subset(aaa_data[[speaker]]$sag$tongue_traces, 
        middle_frame==TRUE & phone%in%c('T','K') & grepl('1',right)), 
        aaa_data[[speaker]]$sag$occlusal_angle, factors_to_retain=c('left','phone','right','word'))#,'ipa'))


cairo_pdf('trajectories_sample_words.pdf', height = 8, width = 5, onefile = T)
par(mfrow=c(2,1))
for (speaker in speakers){
  compare_trajectories(aaa_data, speaker, word=c('sock','cop','see','key'), signal='TTangle', main=paste(speaker,'\n','tongue tip'))
  compare_trajectories(aaa_data, speaker, word=c('sock','cop','see','key'), signal='TDangle', main=paste(speaker,'\n','tongue dorsum'))
  # compare_trajectories(aaa_data, speaker, word=c('sock','cop','see','key'), signal='TRangle', main=paste(speaker,'tongue root'))
}
dev.off()

  compare_trajectories(aaa_data, speaker, word=c('STREAM','SCREAM'), signal='TDangle')#, main=paste(speaker,'\n','tongue tip'))



X "tongue_traces"
X "palate_trace"
no "occlusal_trace"        
no "occlusal_lm"
X "occlusal_angle"
X "palate_trace_rotated"  
no "occlusal_trace_rotated"

# don't need:
 # "Client_Surname"              
 # "Annotation_Label"           
 # "Date_and_time_of_recording"  
 # "date_time"                  
 # "token"    
 # "Sequence_number"        

 "Time_of_sample_in_recording"
 X1-11  Y1-11












# PLOT THE OCCLUSAL ANGLE AND SHOW THE ROTATION OF THE PALATE TRACE
plot_pal_occ(aaa_data, planes='sag')

# ADD YOUR TEXTGRID SEGMENTATION TO THE DATA FRAME CONTAINING THE TONGUE TRACES
# TEXTGRIDS ARE EXPECTED TO BE IN A FOLDER NAMED export_speakername (WITH YOUR ACTUAL SPEAKER NAME) WITHIN THE WORKING DIRECTORY
aaa_data = add_textgrid_segmentation(aaa_data, merge_vl=FALSE, tiers=c('phones','words'), match_words_to_phrase=FALSE)

# PREPARE THE MIDDLE FRAME OF EACH SEGMENT FOR SSANOVA COMPARISON
aaa_data = put_middle_frames_in_long_format(aaa_data)

# CHOOSE A GOOD ORIGIN FOR OUR POLAR COORDINATE SYSTEM
aaa_data = choose_polar_origin(aaa_data, 'xmid_ymid')

# MEASURE TRACES AT SELECTED ANGLES
aaa_data = add_analysis_angles(radii_filename, aaa_data, angles_only=TRUE)
aaa_data = xy2polar(aaa_data)
aaa_data = measure_traces_at_angles(aaa_data)

###########################################################################################################################################
# DONE IMPORTING! NOW CONTINUE WITH plot_aaa_export.r
###########################################################################################################################################


