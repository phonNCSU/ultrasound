
###########################################################################################################################################
# analyze_ultrasound_OSP.r
###########################################################################################################################################

readwritepath = '/home/jimielke/orthognathic/OSP_ultrasound/jacox_speech_new/data/'
scriptpath = '/home/jimielke/scripts/phonNCSU/ultrasound/'

source(paste0(scriptpath,'read_ultrasound_functions.r'))
source(paste0(scriptpath,'tongue_ssanova.r'))

target_words = c('sock','shock','chop','cop','top','sack','shack','chap','cap','tap','see','she','cheap','key','tea','sue','shoe','chew','coo','too')

# this lets us simplifiy our ssanova commands below by filling in the same values for options
ssanova_wrapper <- function(data, data.cat='ipa', token.label='token_id', main){
  polar.ssanova(data, data.cat=data.cat, token.label=token.label, origin=us_data[[speaker]]$sag$origin, 
      palate=us_data[[speaker]]$sag$palate_trace_rotated, is.polar=TRUE, drop_levels=TRUE, crop=TRUE, main=main)
}

# PREPARE THE SSANOVA LIST
osp_ssanovas = list()
for (speaker in speakers){
  osp_ssanovas[[speaker]] = list()
}

###########################################################################################################################################
# PLOT TRAJECTORIES FOR SELECTED WORDS
cairo_pdf('trajectories_sample_words.pdf', height = 8, width = 5, onefile = T)
par(mfrow=c(2,1))
for (speaker in speakers){
  compare_trajectories(us_data, speaker, word=c('sock','cop','see','key'), signal='TTangle', main=paste(speaker,'\n','tongue tip'))
  compare_trajectories(us_data, speaker, word=c('sock','cop','see','key'), signal='TDangle', main=paste(speaker,'\n','tongue dorsum'))
  # compare_trajectories(us_data, speaker, word=c('sock','cop','see','key'), signal='TRangle', main=paste(speaker,'tongue root'))
}
dev.off()
###########################################################################################################################################

###########################################################################################################################################
# COLLECT THE DATA AND RUN THE SSANOVAS
for (speaker in speakers){

  data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & word %in% target_words & phone%in%c('S','SH','CH','T','K') & grepl('1',right)), 
    polar=TRUE, factors_to_retain=c('left','phone','right','word','ipa'))
  osp_ssanovas[[speaker]][['all_consonants']] = ssanova_wrapper(data, main=paste(speaker,'\n','consonants in all contexts'))

  data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & word %in% target_words & grepl('1',phone)), 
    polar=TRUE, factors_to_retain=c('left','phone','right','word','ipa'))
  osp_ssanovas[[speaker]][['all_vowels']] = ssanova_wrapper(data, main=paste(speaker,'\n','vowels in all contexts'))

  for (vowel in c('AA1','AE1','IY1','UW1')){
    data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & word %in% target_words & phone%in%c('S','SH','CH','T','K') & right==vowel), 
      polar=TRUE, factors_to_retain=c('left','phone','right','word','ipa'))
    osp_ssanovas[[speaker]][[vowel]] = ssanova_wrapper(data, main=paste(speaker,'\n','consonants before',arpabet2ipa(vowel)))
  }

  for (consonant in c('S','SH','CH','T','K')){
    data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & word %in% target_words & phone==consonant & grepl('1',right)), 
      polar=TRUE, factors_to_retain=c('left','phone','right','word','ipa'))
    data$right_ipa = factor(arpabet2ipa(data$right))
    osp_ssanovas[[speaker]][[consonant]] = ssanova_wrapper(data, data.cat='right_ipa', main=paste(speaker,'\n',arpabet2ipa(consonant),'by following vowel'))
  }
}
###########################################################################################################################################

###########################################################################################################################################
# PLOT CONSONANT SSANOVAS
cairo_pdf('consonants_all_with_traces.pdf', height = 5, width = 12, onefile = T)
par(mfrow=c(1,2))
for (speaker in speakers){
  show.traces(osp_ssanovas[[speaker]][['all_consonants']]$data, data.cat='ipa', token.label='token_id', origin=us_data[[speaker]]$sag$origin, 
    palate=us_data[[speaker]]$sag$palate_trace_rotated, interpolate=21, drop_levels=TRUE, is.polar=TRUE, main=osp_ssanovas[[speaker]][['all_consonants']]$main)
  replot.tongue.ss(osp_ssanovas[[speaker]][['all_consonants']], palate=us_data[[speaker]]$sag$palate_trace_rotated)
}
dev.off()
###########################################################################################################################################

###########################################################################################################################################
# PLOT VOWEL SSANOVAS
cairo_pdf('vowels_all_with_traces.pdf', height = 5, width = 18, onefile = T)
par(mfrow=c(1,3))
for (speaker in speakers){
  show.traces(osp_ssanovas[[speaker]][['all_vowels']]$data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=us_data[[speaker]]$sag$origin, 
              interpolate=21, palate=us_data[[speaker]]$sag$palate_trace_rotated, drop_levels=TRUE, is.polar=TRUE, main=osp_ssanovas[[speaker]][['all_vowels']]$main)
  add_radial_grid(us_data, speaker, from=-30, to=180)

  show.traces(osp_ssanovas[[speaker]][['all_vowels']]$data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=us_data[[speaker]]$sag$origin, 
              interpolate=21, palate=us_data[[speaker]]$sag$palate_trace_rotated, drop_levels=TRUE, is.polar=TRUE, main=osp_ssanovas[[speaker]][['all_vowels']]$main)
  add_radial_grid(us_data, speaker, from=-30, to=180)

  replot.tongue.ss(osp_ssanovas[[speaker]][['all_vowels']], palate=us_data[[speaker]]$sag$palate_trace_rotated)
  add_radial_grid(us_data, speaker, from=-30, to=180)
}
dev.off()
###########################################################################################################################################

###########################################################################################################################################
# PLOT CONSONANT-BY-CONTEXT SSANOVAS
cairo_pdf('consonants_by_context_with_traces.pdf', height = 20, width = 12, onefile = T)

par(mfrow=c(4,2))
for (speaker in speakers){
  for (vowel in c('AA1','AE1','IY1','UW1')){
    show.traces(osp_ssanovas[[speaker]][[vowel]]$data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=us_data[[speaker]]$sag$origin, 
              interpolate=21, palate=us_data[[speaker]]$sag$palate_trace_rotated, drop_levels=TRUE, is.polar=TRUE, main=osp_ssanovas[[speaker]][[vowel]]$main, 
              xlim=osp_ssanovas[[speaker]][['all_consonants']]$xlim, ylim=osp_ssanovas[[speaker]][['all_consonants']]$ylim)
    add_radial_grid(us_data, speaker, from=-30, to=180)

  }
  for (vowel in c('AA1','AE1','IY1','UW1')){
    replot.tongue.ss(osp_ssanovas[[speaker]][[vowel]], palate=us_data[[speaker]]$sag$palate_trace_rotated, 
              xlim=osp_ssanovas[[speaker]][['all_consonants']]$xlim, ylim=osp_ssanovas[[speaker]][['all_consonants']]$ylim)
  }
}
dev.off()
###########################################################################################################################################

###########################################################################################################################################
# PLOT EACH CONSONANT WITH VOWEL CONTEXT AS FACTOR
cairo_pdf('consonants_separately_by_context_with_traces.pdf', height = 24, width = 12, onefile = T)
par(mfrow=c(5,2))
for (speaker in speakers){
  for (consonant in c('S','SH','CH','T','K')){
    show.traces(osp_ssanovas[[speaker]][[consonant]]$data, data.cat='right_ipa', token.label='token_id', flip='FALSE', origin=us_data[[speaker]]$sag$origin, 
              interpolate=21, palate=us_data[[speaker]]$sag$palate_trace_rotated, drop_levels=TRUE, is.polar=TRUE, main=osp_ssanovas[[speaker]][[consonant]]$main, 
              xlim=osp_ssanovas[[speaker]][['all_consonants']]$xlim, ylim=osp_ssanovas[[speaker]][['all_consonants']]$ylim)
  }
  for (consonant in c('S','SH','CH','T','K')){
    replot.tongue.ss(osp_ssanovas[[speaker]][[consonant]], palate=us_data[[speaker]]$sag$palate_trace_rotated, 
              xlim=osp_ssanovas[[speaker]][['all_consonants']]$xlim, ylim=osp_ssanovas[[speaker]][['all_consonants']]$ylim)
  }
}
dev.off()
###########################################################################################################################################
