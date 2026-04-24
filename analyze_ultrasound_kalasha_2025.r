
###########################################################################################################################################
# analyze_ultrasound_covart_cr.r
###########################################################################################################################################

# CHANGE THESE PATHS

readwritepath = '/home/jimielke/Kalasha/analysis_2025/dlc'
scriptpath = '~/scripts/phonNCSU/ultrasound/'

source(paste0(scriptpath,'read_ultrasound_functions.r'))
source(paste0(scriptpath,'tongue_ssanova.r'))

load('K25_dlc_2025_04_02_processed.RData')
###########################################################################################################################################
# sample plotting
###########################################################################################################################################

speakers = c('Dame1', 'Dame2', 'Dame4', 'Dame5', 'Dame6', 'Kal1', 'Kal4', 'Kal5', 'Kal8', 
	'KalBi1', 'KalBi2', 'KalBi3', 'KalBi4', 'KalBi5', 'KalBi9', 'KalBu1', 'KalBu2', 'KalBu3', 
	'KalBu4', 'KalBu5', 'KalBu6', 'KalBu8', 'KalBu12', 'KalBu13', 'KalRu1', 'KalRu2', 
	'KalRu3', 'KalRu4', 'KalRu5', 'KalRu6', 'KalRu7', 'KalRu9', 'KalRu10', 'Kami1', 'Kami2', 
	'Kami3', 'Kati1', 'Kati2', 'Kati3', 'Kati4', 'Kati7', 'Khow1', 'Khow2', 'Khow3', 'Khow4', 
	'Palu1', 'Palu2', 'Palu3', 'Palu4', 'Shin1', 'Shin2', 'Shin3', 'Shin4', 'Shin5')

kalasha_speakers = c('Kal1', 'Kal4', 'Kal5', 'Kal8', 
	'KalBi1', 'KalBi2', 'KalBi3', 'KalBi4', 'KalBi5', 'KalBi9', 'KalBu1', 'KalBu2', 'KalBu3', 
	'KalBu4', 'KalBu5', 'KalBu6', 'KalBu8', 'KalBu12', 'KalBu13', 'KalRu1', 'KalRu2', 
	'KalRu3', 'KalRu4', 'KalRu5', 'KalRu6', 'KalRu7', 'KalRu9', 'KalRu10')


cairo_pdf('sample_lip_trajectories.pdf', height=5,width=6,onefile=TRUE)
for (sp in setdiff(kalasha_speakers, 'KalRu5')){
	compare_trajectories(us_data, sp, word=c("ba'a'","paa","tap"), plane='video', signal='lips_area', main=paste(sp,'\n','lips area'))
}
dev.off()


cairo_pdf('sample_lip_hyoid_trajectories.pdf', height=12,width=6,onefile=TRUE)
	par(mfrow=c(3,1))
for (sp in setdiff(kalasha_speakers, 'KalRu5')){
	compare_trajectories(us_data, sp, word=c("ba'a'","paa","tap"), plane='video', signal='lips_area', main=paste(sp,'\n','lips area'))
	compare_trajectories(us_data, sp, word=c("ba'a'","paa","tap"), plane='sag', signal='hyoid_x', main=paste(sp,'\n','hyoid X'))
	compare_trajectories(us_data, sp, word=c("ba'a'","paa","tap"), plane='sag', signal='hyoid_y', main=paste(sp,'\n','hyoid Y'))
}
dev.off()

# 1 vallecula
# 2 root1
# 3 root2
# 4 dorsum1
# 5 dorsum2
# 6 body1
# 7 body2
# 8 blade1
# 9 blade2
# 10 tip1
# 11 tip2

#   NAs in token ids
#   prepare for formant measurement

data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & word %in% c("ba", "pa") & phone %in% c('p', 'b')),
                               polar=TRUE, factors_to_retain=c('left','phone','right','word'))
    Kal_ssanovas[[speaker]][['all_consonants']] = ssanova_wrapper(data, main=paste(speaker,'\n','consonants in all contexts'))



    Kal_ssanovas = list()
  for (speaker in kalasha_speakers){
    Kal_ssanovas[[speaker]] = list()
  }

  cairo_pdf('test.pdf', onefile=TRUE)
  # COLLECT THE DATA AND RUN THE SSANOVAS
  for (speaker in setdiff(kalasha_speakers,c('KalRu2'))){
    print(speaker)
    data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & word %in% c("ba", "pa") & phone %in% c('p', 'b')),
                               polar=TRUE, factors_to_retain=c('left','phone','right','word'))
    # show.traces(data)
    Kal_ssanovas[[speaker]][['all_consonants']] = ssanova_wrapper(data, main=paste(speaker,'\n','consonants in all contexts'))
    
  }

  dev.off()


    cairo_pdf('test2.pdf', onefile=TRUE)
  # COLLECT THE DATA AND RUN THE SSANOVAS
  for (speaker in setdiff(kalasha_speakers,c('KalRu2'))){
    print(speaker)
    data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & word %in% c("ba", "pa") & phone %in% c('p', 'b')),
                               polar=TRUE, factors_to_retain=c('left','phone','right','word'))
    reshow.traces(Kal_ssanovas[[speaker]][['all_consonants']])
    
  }

  dev.off()


    cairo_pdf('3_31_validate_tongues.pdf', height = 5, width = 12, onefile = T)
  par(mfrow=c(1,2))
  for (speaker in setdiff(kalasha_speakers,c('KalRu2'))){
  	print(speaker)
    show.traces(Kal_ssanovas[[speaker]][['all_consonants']]$data, data.cat='phone', token.label='token_id', origin=us_data[[speaker]]$sag$origin, 
                interpolate=21, drop_levels=TRUE, is.polar=TRUE, main=Kal_ssanovas[[speaker]][['all_consonants']]$main)
    replot.tongue.ss(Kal_ssanovas[[speaker]][['all_consonants']])
  }
  dev.off()