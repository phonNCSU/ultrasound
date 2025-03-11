
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
# load('cr_dlc_2025_02_10.RData')
load('cr_dlc_2025_02_20_quality_phrases.RData')
###########################################################################################################################################
# sample plotting
###########################################################################################################################################

speaker = names(us_data)[1]

us_data$cr01$video$lip_traces$lips_area_cm2 = us_data$cr01$video$lip_traces$lips_area*(0.4263213^2)/100
# STEP 3b: compare trajectories only uses analysis values

cairo_pdf('sample_lip_trajectories.pdf', height=5,width=6)
compare_trajectories(us_data, 'cr01', word=c('KEEP','CREEP'), plane='video', signal='lips_area_cm2')#, main=paste(speaker,'\n','tongue tip'))
dev.off()



# STEP 3c: put polar data in long format for ssanova
# consonant_data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('T','K') & word%in%c('STREAM','SCREAM')))
consonant_data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('AH1','AA1','IY1') & word%in%c('RUN','ROCK','REED')), polar=TRUE)
consonant_data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('S','K','G') & word%in%c('SEAM','KEEP','GOAT')), polar=TRUE)
consonant_data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('R','AA1','K') & word%in%c('ROCK')), polar=TRUE)

# STEP 3d: convert long to polar (AGAIN) and interpolate for ssanova
show.traces(consonant_data, data.cat='phone', token.label='token_id', origin=us_data[[speaker]]$sag$origin, 
	drop_levels=TRUE, palate=us_data[[speaker]]$sag$palate_trace_rotated, 
	main=speaker, is.polar=TRUE)
add_radial_grid(us_data, speaker, from=0, to=140, length=150)


compare_trajectories(us_data, 'cr01', word=c('ROCK'), signal='TDangle')#, main=paste(speaker,'\n','tongue tip'))
compare_trajectories(us_data, 'cr01', word=c('SEAM','KEEP','GOAT'), signal='Y5')#, main=paste(speaker,'\n','tongue tip'))
compare_trajectories(us_data, 'cr01', plane='video', word=c('SEAM','KEEP','GOAT'), signal='rightLip_x')#, main=paste(speaker,'\n','tongue tip'))

# xx = polar.ssanova(consonant_data, data.cat='phone', token.label='token_id', origin=us_data[[speaker]]$sag$origin, 
# 	drop_levels=TRUE, palate=us_data[[speaker]]$sag$palate_trace_rotated, crop=TRUE, interpolate=21, main='scream stream')

xx = polar.ssanova(consonant_data, data.cat='phone', token.label='token_id', origin=us_data[[speaker]]$sag$origin, 
	drop_levels=TRUE, palate=us_data[[speaker]]$sag$palate_trace_rotated, crop=TRUE, interpolate=21, main='scream stream', is.polar=TRUE)


# CL to do

# DOESN'T MATTER: DON'T NEED ANGLES ANYWAY
# . ANGLES DON'T WORK RIGHT NOW
# . universal multiplier?
# 7: 394-45   7cm: mmpermixel = 70/350
# 8: 401-43   8cm: mmpermixel = 80/359
# 9: 401-43   9cm: mmpermixel = 90/359

# x rotate occlusal the other way
# x update cr02 palate trace

# x sort out phrase ID and phrase

# x correct some palates
#   use YPR to exclude tokens
#   correct some occlusal plane angles
#   improve some polar origins
#   convert to mm

# IDENTIFY PHRASES AND FIND TARGETS

# x use prompt tier
# x match the prompt tier to the target list and isolate the target words
# x correct some mislabeled phrases
# x find target words for each phrase
# x flag YPR words

#   shiny app with phrase chunking

# x make some lip signals

#   compare creep-keep-reed lips and tongue trajectories

measure_traces_at_angles()

plot_pal_occ()
plot_traces()
plot_some_signals()
plot_trajectories()




set_basic_ssanova_targets()
set_basic_comparison_levels()
do_cl_ssanova_comparisons()



get_target_mean_signals()
make_ssanovas_plots_and_calculate_displacements()

make_ssanovas_plots()

plot_displacement_measures()

plot_analysis_values()

add_interpolated_formants()
measure_formant_badness()



find_mean_tongue()
find_displacement()
find_all_displacements()


phone_locations_plot_nongeneric()

phone_locations_plot()

find_displacement_density()

add_target_phone_time()


make_trajectory_comparisons()

plot_trajectory_comparisons()

compare_trajectories()
add_radial_grid()


plot_sample_frames()



#################################################################


target_phrases = read.csv('~/covart/cr/dlc/data_october_2024/cr_target_phrases.csv')


all_phrases = c()
for (sp in names(us_data)){
	print(sp)
	data = us_data[[sp]]$sag$tongue_traces

	# data$new_phrase = c(TRUE, diff(data$Time_of_sample_in_recording)>0.1)
	data$new_phrase = c(TRUE, data$prompt[2:nrow(data)] != data$prompt[1:(nrow(data)-1)])
	data$new_word = c(TRUE, data$word_id[2:nrow(data)] != data$word_id[1:(nrow(data)-1)])
	# data$phrase = 

	data$phrase = trimws(toupper(data$prompt))
	data[data$phrase%in%c('','SP','1'),'phrase'] = NA
	data[data$phrase%in%c('THEY HAVE A DWEEB','THEY MET A DWEEB AGAIN'),'phrase'] = 'THEY HAVE A DWEEB AGAIN'
	data[data$phrase%in%c('THEY HIT A HUGE REEF AGAIN'),'phrase'] = 'THEY HIT SLUDGE REEF AGAIN'
	data[data$phrase%in%c('I MET HIS WEED AGAIN'),'phrase'] = 'I MET ED WEED AGAIN'
	data[data$phrase%in%c('SHE SAW THE CAR SITTING THERE'),'phrase'] = 'SHE SAW THE CAR SITTING IDLE'
	data[data$phrase%in%c('THEY SAW THE RECORD TODAY'),'phrase'] = 'THEY SAW THE SCROD TODAY'
	data[data$phrase%in%c('THEY SAW A GEEK TODAY'),'phrase'] = 'THEY SAW A SOCK TODAY'




595	cr02	537.387661682187	NA	AGAIN THEY SAW A SPROCKET sp TODAY 
They saw    ;a sprocket      ;today  

660	cr02	744.495900868991	NA	sp SHOW ME THE DYE
Show me     ;the dye         ;again   

2135	cr06	575.020462121914	NA	HE DID sp sp HE DID A SQUAT TODAY sp
HE DID A SQUAT TODAY  (repeated with YPR problems)

4784	cr11	723.637741256	NA	THEY SAW A SOCK TODAY sp
They saw    ;a sock          ;today     

4939	cr11	1266.13413566	NA	sp SHOW ME A THRUSH SOMETIME sp SHOW ME A THRUSH SOMETIME sp
Show me     ;a thrush        ;sometime  (twice, first with YPR issues)

11792	cr21	263.935183624	NA	SAW HIS sp PA SITTING sp AT HOME sp
He saw      ;his pa sitting  ;at home

library(plyr)
phrase_log = read.csv('data_october_2024/textgrid_phrase_log.csv')
phone_summary = ddply(subset(phrase_log, matched_target==1), .(target_phrase, all_phones), summarize, n=length(phrase_id))
quality_summary = ddply(subset(phrase_log, matched_target==1), .(target_phrase, speaker, quality_label), summarize, n=length(phrase_id)) 
write.csv(quality_summary, 'quality_summary.csv', row.names=TRUE)
write.csv(phone_summary, 'phone_summary.csv', row.names=TRUE)
	
check_phones = phone_summary[phone_summary$n<5,'all_phones']
phones_to_check = phrase_log[phrase_log$all_phones%in%check_phones,]
phones_to_check = subset(phones_to_check, !(target_phrase=='A SWIPE' & all_phones=='AH0 S W AY1 P Ph'))
phones_to_check = subset(phones_to_check, !(target_phrase=='A WIPE' & all_phones=='AH0 W AY1 P Ph'))
write.csv(phones_to_check, 'phones_to_check.csv', row.names=TRUE)


and what if they advanced the prompt while speaking?
      
	data$phrase_id = NA
	data$target = NA
	data$target_id = NA


	phrase_start_indices = which(data$new_phrase)


	for (i in 1:length(phrase_start_indices)){
	# 	print(i)
		if (i == length(phrase_start_indices)){
			phrase_rows = phrase_start_indices[i]:nrow(data)
		}else{
			phrase_rows = phrase_start_indices[i]:(phrase_start_indices[i+1]-1)
		}
		phrase_data = data[phrase_rows,]
	# 	word_starts = phrase_data[phrase_data$new_word,]
	# 	the_words = word_starts$word
	# 	narrow_words = the_words[!the_words%in%c('','sp')]
	# 	data[phrase_rows,'phrase_id'] = paste(sp,'1',paste(narrow_words, collapse='-'),round(data$phone_start[phrase_start_indices[i]],3), sep='_')
	# 	data[phrase_rows,'phrase_words'] = paste(narrow_words, collapse=' ')
		# all_phrases = rbind(all_phrases, data.frame(speaker=sp, phrase=data$phrase))#,
			# phrase_words=paste(narrow_words, collapse=' '), 
			# phrase_id=paste(sp,'1',paste(narrow_words, collapse='-'), round(data$phone_start[phrase_start_indices[i]],3), sep='_')))
	# }
		# if (is.na(data[phrase_start_indices[i],'phrase'])){
			# print(data$word)
			word_starts = phrase_data[phrase_data$new_word,]
			the_words = word_starts$word
			phrase_words=paste(the_words, collapse=' ')
		# }
		newdata = data.frame(speaker=sp, start=phrase_data$phone_start[1], phrase=phrase_data$phrase[1], words=phrase_words)
		all_phrases = rbind(all_phrases, newdata)
	}
}

all_phrases$matching = all_phrases$phrase %in% target_phrases$phrase
write.csv(all_phrases, 'cl_all_phrases_02_11_25.csv')
# cr01_1_TO_15.675
