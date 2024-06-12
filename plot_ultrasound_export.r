
###########################################################################################################################################
#
# THESE COMMANDS SHOULD WORK IF YOU HAVE ALREADY GONE THROUGH read_aaa_export.r WITH YOUR DATASET
#
###########################################################################################################################################


###########################################################################################################################################
# MAKE SSANOVA COMPARISONS
###########################################################################################################################################


# START BY DEFINING SOME VARIABLES TO MAKE OUR SSANOVA COMMANDS EASIER TO INTERPRET
speaker = speakers[1]
my_origin = aaa_data[[speaker]]$sag$origin

# ADD ipa COLUMNS TO OUR DATA
aaa_data[[speaker]]$sag$tongue_traces$ipa = factor(arpabet2ipa(aaa_data[[speaker]]$sag$tongue_traces$phone))
# aaa_data[[speaker]]$sag$tongue_traces_long_rotated$ipa = factor(arpabet2ipa(aaa_data[[speaker]]$sag$tongue_traces_long_rotated$phone))

# MAKE SOME SUBSETS OF THE DATA TO ANALYZE (CHANGE THE LAST TWO LINES TO MAKE A DATA SUBSET FOR THE COMPARISON YOU ADDED TO THE PROMPTS)

target_words = c('sock','shock','chop','cop','top','sack','shack','chap','cap','tap','see','she','cheap','key','tea','sue','shoe','chew','coo','too')



# consonant_data = subset(aaa_data[[speaker]]$sag$tongue_traces_long_rotated, phone%in%c('S','SH','CH','T','K'))

consonant_data = traces_wide_to_long(subset(aaa_data[[speaker]]$sag$tongue_traces, 
                                            middle_frame==TRUE & word %in% target_words & phone%in%c('S','SH','CH','T','K') & grepl('1',right)), 
                                      aaa_data[[speaker]]$sag$occlusal_angle, factors_to_retain=c('left','phone','right','word','ipa'))










consonant_data$ipa_word = factor(paste(consonant_data$ipa, 'in', consonant_data$word))
for (ph in rev(c('/k/','/t/','/tʃ/','/ʃ/','/s/'))){
    consonant_data$ipa = relevel(consonant_data$ipa, ph)
}
                
consonant_words = c('/s/ in sock',  '/s/ in sack',  '/s/ in see',    '/s/ in sue', 
                    '/ʃ/ in shock', '/ʃ/ in shack', '/ʃ/ in she',    '/ʃ/ in shoe', 
                    '/tʃ/ in chop', '/tʃ/ in chap', '/tʃ/ in cheap', '/tʃ/ in chew', 
                    '/k/ in cop',   '/k/ in cap',   '/k/ in key',    '/k/ in coo',
                    '/t/ in top',   '/t/ in tap',   '/t/ in tea',    '/t/ in too')
consonant_data = subset(consonant_data, ipa_word %in% consonant_words)

consonant_AA_data = subset(consonant_data, word%in%c('sock','shock','chop','cop','top'))
consonant_AE_data = subset(consonant_data, word%in%c('sack','shack','chap','cap','tap'))
consonant_IY_data = subset(consonant_data, word%in%c('see','she','cheap','key','tea'))
consonant_UW_data = subset(consonant_data, word%in%c('sue','shoe','chew','coo','too'))


                          
# vowel_data = subset(aaa_data[[speaker]]$sag$tongue_traces_long_rotated, phone%in%c('AA1','AE1','IY1','UW1'))

vowel_data = traces_wide_to_long(subset(aaa_data[[speaker]]$sag$tongue_traces, 
                                        middle_frame==TRUE & word %in% target_words & phone%in%c('AA1','AE1','IY1','UW1')), 
                                 aaa_data[[speaker]]$sag$occlusal_angle, factors_to_retain=c('left','phone','right','word','ipa'))




vowel_data$ipa_word = factor(paste(vowel_data$ipa, 'in', vowel_data$word))
for (ph in rev(c('/ɑ/','/æ/','/i/','/u/'))){
    vowel_data$ipa = relevel(vowel_data$ipa, ph)
}

vowel_words = c('/ɑ/ in sock',  '/æ/ in sack',  '/i/ in see',   '/u/ in sue', 
                '/ɑ/ in shock', '/æ/ in shack', '/i/ in she',   '/u/ in shoe', 
                '/ɑ/ in chop',  '/æ/ in chap',  '/i/ in cheap', '/u/ in chew', 
                '/ɑ/ in cop',   '/æ/ in cap',   '/i/ in key',   '/u/ in coo',
                '/ɑ/ in top',   '/æ/ in tap',   '/i/ in tea',   '/u/ in too')
vowel_data = subset(vowel_data, ipa_word %in% vowel_words)


library(akima)
XY = subset(consonant_data, token_id=='File001_1_sock_S_1.56')[,c('X','Y')]


# THESE COMMANDS WILL MAKE A 6-PAGE PDF WITH RAW TRACES AND AN SSANOVA FIGURE FOR EACH OF THREE COMPARIONS
# WE USE cairo_pdf() INSTEAD OF pdf() IN ORDER TO MAKE PDFS WITH IPA SYMBOLS IN THEM
# BUT YOU CAN CHANGE TO pdf() TO AVOID ERRORS (AT THE COST OF LOSING IPA UNTIL YOU FIGURE OUT CAIRO)
# IN THE PAST THERE HAVE BEEN PROBLEMS WITH cairo_pdf ON SOME MACS, BUT THE PROBLEMS ARE SOLVABLE


cairo_pdf('OSP_ssanovas_140.pdf', height=5, width=6, onefile=TRUE)

show.traces(vowel_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, interpolate=21, drop_levels=TRUE, main='vowels')
x=polar.ssanova(vowel_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, drop_levels=TRUE, main='vowels', 
              palate=aaa_data[[speaker]]$sag$palate_trace_rotated, crop=TRUE, interpolate=21)

show.traces(consonant_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, interpolate=21, drop_levels=TRUE, main='consonants')
x=polar.ssanova(consonant_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, drop_levels=TRUE, main='consonants',
              palate=aaa_data[[speaker]]$sag$palate_trace_rotated, crop=TRUE, interpolate=21)

show.traces(consonant_AA_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, interpolate=21, drop_levels=TRUE, main='consonants before /ɑ/')
x=polar.ssanova(consonant_AA_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, drop_levels=TRUE, main='consonants before /ɑ/',
              palate=aaa_data[[speaker]]$sag$palate_trace_rotated, crop=TRUE, interpolate=21)

show.traces(consonant_AE_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, interpolate=21, drop_levels=TRUE, main='consonants before /æ/')
x=polar.ssanova(consonant_AE_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, drop_levels=TRUE, main='consonants before /æ/',
              palate=aaa_data[[speaker]]$sag$palate_trace_rotated, crop=TRUE, interpolate=21)

show.traces(consonant_IY_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, interpolate=21, drop_levels=TRUE, main='consonants before /i/')
x=polar.ssanova(consonant_IY_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, drop_levels=TRUE, main='consonants before /i/',
              palate=aaa_data[[speaker]]$sag$palate_trace_rotated, crop=TRUE, interpolate=21)

show.traces(consonant_UW_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, interpolate=21, drop_levels=TRUE, main='consonants before /u/')
x=polar.ssanova(consonant_UW_data, data.cat='ipa', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, drop_levels=TRUE, main='consonants before /u/',
              palate=aaa_data[[speaker]]$sag$palate_trace_rotated, crop=TRUE, interpolate=21)

dev.off()

# SUPERIMPOSE A RADIAL GRID ON A TONGUE TRACE PLOT TO DECIDE WHERE TO MEASURE
show.traces(consonant_data, data.cat='ipa_word', token.label='token_id', flip='FALSE', origin=aaa_data[[speaker]]$sag$origin, drop_levels=TRUE, main='consonant_data by following vowel')
add_radial_grid(aaa_data, speaker, from=-30, to=180)

# PLOT TRAJECTORIES
compare_trajectories(aaa_data, speaker, word=c('cop','cap','key','coo'), signal='TDangle')
compare_trajectories(aaa_data, speaker, word=c('sock','sack','see','sue'), signal='TTangle')


vowel_words = c('/ɑ/ in ',  '/æ/ in ',  '/i/ in ',   '/u/ in ', 
                '/ɑ/ in shock', '/æ/ in shack', '/i/ in she',   '/u/ in shoe', 
                '/ɑ/ in chop',  '/æ/ in chap',  '/i/ in cheap', '/u/ in chew', 
                '/ɑ/ in cop',   '/æ/ in cap',   '/i/ in key',   '/u/ in coo',
                '/ɑ/ in top',   '/æ/ in tap',   '/i/ in tea',   '/u/ in too')

