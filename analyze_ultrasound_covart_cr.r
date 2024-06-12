
###########################################################################################################################################
# analyze_ultrasound_covart_cr.r
###########################################################################################################################################

readwritepath = '/home/jimielke/covart/cr/dlc/'
scriptpath = '~/scripts/phonNCSU/ultrasound/'

source(paste0(scriptpath,'read_ultrasound_functions.r'))
source(paste0(scriptpath,'tongue_ssanova.r'))

# do this to create us_data from the original files:
#  source(paste0(scriptpath,'read_dlc_export_covart_cr.r'))

# --or-- 

# do something like this to reload us_data later:
load('cr_dlc_2024_06_12.RData')

###########################################################################################################################################
# sample plotting
###########################################################################################################################################

speaker = names(us_data)[1]

# STEP 3b: compare trajectories only uses analysis values
compare_trajectories(us_data, speaker, word=c('STREAM','SCREAM'), signal='TDangle')#, main=paste(speaker,'\n','tongue tip'))

# STEP 3c: put polar data in long format for ssanova
# consonant_data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('T','K') & word%in%c('STREAM','SCREAM')))
consonant_data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('T','K') & word%in%c('STREAM','SCREAM')), polar=TRUE)

# STEP 3d: convert long to polar (AGAIN) and interpolate for ssanova
show.traces(consonant_data, data.cat='phone', token.label='token_id', origin=us_data[[speaker]]$sag$origin, 
	drop_levels=TRUE, palate=us_data[[speaker]]$sag$palate_trace_rotated, 
	main=speaker, is.polar=TRUE)
add_radial_grid(us_data, speaker, from=0, to=140, length=150)

# xx = polar.ssanova(consonant_data, data.cat='phone', token.label='token_id', origin=us_data[[speaker]]$sag$origin, 
# 	drop_levels=TRUE, palate=us_data[[speaker]]$sag$palate_trace_rotated, crop=TRUE, interpolate=21, main='scream stream')

xx = polar.ssanova(consonant_data, data.cat='phone', token.label='token_id', origin=us_data[[speaker]]$sag$origin, 
	drop_levels=TRUE, palate=us_data[[speaker]]$sag$palate_trace_rotated, crop=TRUE, interpolate=21, main='scream stream', is.polar=TRUE)

