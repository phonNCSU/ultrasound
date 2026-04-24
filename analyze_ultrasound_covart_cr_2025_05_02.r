
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
load('cr_dlc_2025_04_02.RData')
###########################################################################################################################################
# sample plotting
###########################################################################################################################################
speakers = sort(names(us_data))
speakers = c("cr06", "cr07", "cr08", "cr11", "cr15", "cr18", "cr23", "cr24", "cr26", "cr28", "cr29", "cr32", "cr30", "cr09", "cr13", "cr03", "cr02", "cr05", "cr01", "cr04", "cr19", "cr10", "cr14", "cr21", "cr33")


source(paste0(scriptpath,'read_ultrasound_functions.r'))
source(paste0(scriptpath,'tongue_ssanova.r'))

# sp='cr06'
# xx = us_data[[sp]]$sag$tongue_traces[1:40,]


measure_tongue <- function(aaa_data, speakers=NULL, plotting=FALSE){

	if (is.null(speakers)) speakers=names(aaa_data)

	for (sp in speakers){
		print(sp)
		plane='sag'
		tonguedata = aaa_data[[sp]][[plane]]$tongue_traces
		if (nrow(tonguedata)){
			tonguedata$root_advancement = with(tonguedata, (X1+X2+X3)/3)
			tonguedata$blade_angle = with(tonguedata, atan((Y10-Y8)/(X10/X8)))
			aaa_data[[sp]][[plane]]$tongue_traces = tonguedata
		}
	}
	aaa_data
}

us_data = measure_tongue(us_data)

require(pracma)
for (sp in speakers){
	print(sp)
	us_data[[sp]]$sag$tongue_traces$absolute_concavity = NA
	us_data[[sp]]$sag$tongue_traces$relative_concavity = NA
	for (r in 1:nrow(us_data[[sp]]$sag$tongue_traces)){

		one_token = us_data[[sp]]$sag$tongue_traces[r,]
		original_polygon = data.frame(X=c(us_data[[sp]]$sag$origin[1], as.numeric(one_token[,paste0('X',1:11)])), 
			                          Y=c(us_data[[sp]]$sag$origin[2], as.numeric(one_token[,paste0('Y',1:11)])))
		convex_indices <- chull(original_polygon$X, original_polygon$Y)
		tongue_area = round(abs(polyarea(original_polygon$X,original_polygon$Y)),6)
		convex_area = round(abs(polyarea(original_polygon$X[convex_indices],original_polygon$Y[convex_indices])),6)
		# print (c(r,tongue_area,convex_area))
		us_data[[sp]]$sag$tongue_traces$absolute_concavity[r] = convex_area - tongue_area
		us_data[[sp]]$sag$tongue_traces$relative_concavity[r] = sqrt((convex_area-tongue_area)/tongue_area)
		# print (c(r,tongue_area,convex_area,absolute_concavity,relative_concavity))
	}
}

save(us_data, file='cr_dlc_2025_04_17.RData')
# load('cr_dlc_2025_04_17.RData')

# plot(0,0,xlim=c(-2,7),ylim=c(-2,7))
# convexpolygon <- chull(original_polygon$X, original_polygon$Y)
# polygon(original_polygon$X,original_polygon$Y, border='blue')
# polygon(original_polygon$X[convexpolygon],original_polygon$Y[convexpolygon], border='red')
# library(pracma)
# abs(polyarea(original_polygon$X,original_polygon$Y))
# abs(polyarea(original_polygon$X[convexpolygon],original_polygon$Y[convexpolygon]))



# tongues_measurements$concavity <- with(tongues_measurements, sqrt((convex_area-tongue_area)/tongue_area))





source(paste0(scriptpath,'read_ultrasound_functions.r'))
sp = names(us_data)[1]
# compare_trajectories(us_data, sp, word=c('KEEP','CREEP'), plane='video', signal='lips_area', main=paste(sp,'\n','lips area'))
compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-KEEP','A-CREEP'), 
	plane='video', signal='lips_area', ref_phone='Kh', ref_edge='phone_end', main=paste(sp,'\n','lips area'))


# us_data$cr01$video$lip_traces$lips_area_cm2 = us_data$cr01$video$lip_traces$lips_area*(0.4263213^2)/100
# STEP 3b: compare trajectories only uses analysis values

cairo_pdf('sample_lip_trajectories_sheep_stream_seam_2025_04_23.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-SHEEP','A-STREAM','A-SEAM'), 
		plane='video', signal='lips_area', ref_phone='AH0', ref_edge='phone_end', main=paste(sp,'\n','lips area'))
}
dev.off()
cairo_pdf('sample_lip_trajectories_reef_eave_2025_04_23.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('THIS-REEF','THIS-EAVE'), 
		plane='video', signal='lips_area', ref_phone='S', ref_edge='phone_end', main=paste(sp,'\n','lips area'))
}
dev.off()

cairo_pdf('sample_lip_trajectories_pa_sitting_car_sitting_2025_04_23.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('PA-SITTING','CAR-SITTING'), 
		plane='video', signal='lips_area', ref_phone='S', ref_edge='phone_start', main=paste(sp,'\n','lips area'))
}
dev.off()

cairo_pdf('sample_lip_trajectories_2025_04_09.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-KEEP','A-CREEP'), 
		plane='video', signal='lips_area', ref_phone='Kh', ref_edge='phone_end', main=paste(sp,'\n','lips area'))
}
dev.off()

source(paste0(scriptpath,'read_ultrasound_functions.r'))
cairo_pdf('sample_root_trajectories_2025_04_23.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-KEEP','A-CREEP'), 
		plane='sag', signal='root_advancement', ref_phone='Kh', ref_edge='phone_end', ylim=c(-2,2), main=paste(sp,'\n','root advancement'))
}
dev.off()


source(paste0(scriptpath,'read_ultrasound_functions.r'))
source(paste0(scriptpath,'tongue_ssanova.r'))

blade_crow_coat_comparison = list()

cairo_pdf('sample_blade_trajectories_2025_04_23.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	xxxx = compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-COAT','A-CROW'), 
		plane='sag', signal='blade_angle', ref_phone='Kh', ref_edge='phone_end', ylim=c(-1.5,0.5), main=paste(sp,'\n','blade_angle'))#,
		# show.ssanova=TRUE)
	# blade_creep_keep_comparison[[sp]] = xxxx
}
dev.off()
cairo_pdf('sample_Y9_trajectories_2025_04_23.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	xxxx = compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-COAT','A-CROW'), 
		plane='sag', signal='Y9', ref_phone='Kh', ref_edge='phone_end', ylim=c(1,6), main=paste(sp,'\n','Y9'))#,
		# show.ssanova=TRUE)
	# blade_creep_keep_comparison[[sp]] = xxxx
}
dev.off()

cairo_pdf('sample_Y11_trajectories_2025_04_23.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	xxxx = compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-COAT','A-CROW'), 
		plane='sag', signal='Y11', ref_phone='Kh', ref_edge='phone_end', ylim=c(0,6), main=paste(sp,'\n','Y11'))#,
		# show.ssanova=TRUE)
	# blade_creep_keep_comparison[[sp]] = xxxx
}
dev.off()

root_creep_keep_comparison = list()

cairo_pdf('sample_root_trajectories_2025_04_23.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	xxxx = compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-CREEP','A-KEEP'), 
		plane='sag', signal='root_advancement', ref_phone='Kh', ref_edge='phone_end', ylim=c(-2,2), main=paste(sp,'\n','root advancement'),
		show.ssanova=TRUE)
	root_creep_keep_comparison[[sp]] = xxxx
}
dev.off()


creep_keep_root = c()
for (sp in speakers){
	phrases = c('A-CREEP', 'A-KEEP')
	difference_direction = -1

	sp_model = root_creep_keep_comparison[[sp]]$ssanova$ss.cart
	sp_ref_time = root_creep_keep_comparison[[sp]][[2]]
	phrase1 = subset(sp_model, phrase==phrases[1])
	phrase2 = subset(sp_model, phrase==phrases[2])

	diff_span = phrase1$X[range(which(phrase1$ss.upper.CI.Y < phrase2$ss.lower.CI.Y))] - sp_ref_time
	diff_size = max(difference_direction*(phrase1$ss.Fit-phrase2$ss.Fit))

	diff_peak = phrase1$X[which(difference_direction*(phrase1$ss.Fit-phrase2$ss.Fit) == diff_size)] - sp_ref_time
	creep_keep_root = rbind(creep_keep_root, data.frame(speaker=sp, diff_start=diff_span[1], diff_peak=diff_peak, diff_end=diff_span[2], diff_size=diff_size))
}

lips_creep_keep_comparison = list()

cairo_pdf('sample_lip_trajectories_2025_04_23.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	xxxx = compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-CREEP','A-KEEP'), 
		plane='video', signal='lips_area', ref_phone='Kh', ref_edge='phone_end', main=paste(sp,'\n','lips area'),
		show.ssanova=TRUE)
	lips_creep_keep_comparison[[sp]] = xxxx

}
dev.off()

creep_keep_lips = c()
for (sp in speakers){
	phrases = c('A-CREEP', 'A-KEEP')
	difference_direction = -1

	sp_model = lips_creep_keep_comparison[[sp]]$ssanova$ss.cart
	sp_ref_time = lips_creep_keep_comparison[[sp]][[2]]
	phrase1 = subset(sp_model, phrase==phrases[1])
	phrase2 = subset(sp_model, phrase==phrases[2])

	diff_span = phrase1$X[range(which(phrase1$ss.upper.CI.Y < phrase2$ss.lower.CI.Y))] - sp_ref_time
	diff_size = max(difference_direction*(phrase1$ss.Fit-phrase2$ss.Fit))

	diff_peak = phrase1$X[which(difference_direction*(phrase1$ss.Fit-phrase2$ss.Fit) == diff_size)] - sp_ref_time
	creep_keep_lips = rbind(creep_keep_lips, data.frame(speaker=sp, diff_start=diff_span[1], diff_peak=diff_peak, diff_end=diff_span[2], diff_size=diff_size))
}




cairo_pdf('sample_blade_trajectories_2025_04_16.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-CROW','A-COAT'), 
		plane='sag', signal='blade_angle', ref_phone='Kh', ref_edge='phone_end', ylim=c(-1.5,0.5), main=paste(sp,'\n','blade angle'))
}
dev.off()

cairo_pdf('sample_blade_trajectories_rock_reed_2025_04_16.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-ROCK','A-REED'), 
		plane='sag', signal='blade_angle', ref_phone='R', ref_edge='phone_end', ylim=c(-1.5,0.5), main=paste(sp,'\n','blade angle'))
}
dev.off()


cairo_pdf('sample_Yx_trajectories_rock_reed_2025_04_23.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	for (x in 1:11){
		Yx = paste0('Y',x)
		compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-ROCK','A-REED'), 
			plane='sag', signal=Yx, ref_phone='R', ref_edge='phone_end',  ylim=c(0,10), main=paste(sp,'\n',Yx))
	}
}
dev.off()


cairo_pdf('sample_Yx_trajectories_crow_coat_2025_04_23.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	for (x in 1:11){
		Yx = paste0('Y',x)
		compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-CROW','A-COAT'), 
			plane='sag', signal=Yx, ref_phone='Kh', ref_edge='phone_end',  ylim=c(0,10), main=paste(sp,'\n',Yx))
	}
}
dev.off()


cairo_pdf('sample_concavity_trajectories_keep_creep_2025_04_17.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-KEEP','A-CREEP'), 
		plane='sag', signal='absolute_concavity', ref_phone='Kh', ref_edge='phone_end', ylim=c(-0.5,2), main=paste(sp,'\n','absolute concavity'))
}
dev.off()

cairo_pdf('sample_absolute_concavity_trajectories_rock_reed_2025_04_17.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-ROCK','A-REED'), 
		plane='sag', signal='absolute_concavity', ref_phone='R', ref_edge='phone_end', ylim=c(-0.5,2), main=paste(sp,'\n','absolute concavity'))
}
dev.off()

cairo_pdf('sample_relative_concavity_trajectories_rock_reed_2025_04_17.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	compare_trajectories(us_data, sp, select_col = 'phrase', select_vals=c('A-ROCK','A-REED'), 
		plane='sag', signal='relative_concavity', ref_phone='R', ref_edge='phone_end', ylim=c(-0.1,0.4), main=paste(sp,'\n','relative concavity'))
}
dev.off()

cairo_pdf('rock_reef_2025_04_09.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	consonant_data = traces_wide_to_long(subset(us_data[[sp]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('R') & word%in%c('ROCK','REEF')), polar=TRUE)

	show.traces(consonant_data, data.cat='word', token.label='token_id', origin=us_data[[sp]]$sag$origin, 
		drop_levels=TRUE, #palate=us_data[[speaker]]$sag$palate_trace_rotated, 
		main=sp, is.polar=TRUE, interpolate=100)
}
dev.off()




cairo_pdf('reed_drop_sprocket_spree_stromboli_2025_04_15.pdf', height=5,width=6,onefile=TRUE)
for (sp in speakers){
	consonant_data = traces_wide_to_long(subset(us_data[[sp]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('R')
	 & word%in%c('REED','DROP','SPROCKET','SPREE','STROMBOLI')), polar=TRUE)

	show.traces(consonant_data, data.cat='word', token.label='token_id', origin=us_data[[sp]]$sag$origin, 
		drop_levels=TRUE, #palate=us_data[[speaker]]$sag$palate_trace_rotated, 
		main=sp, is.polar=TRUE, interpolate=100)
}
dev.off()


data_for_coat_goat_scatterplot = c()
for (sp in speakers){
	middle_frame_data = subset(us_data[[sp]]$sag$tongue_traces, middle_frame==TRUE & word %in% c('COAT','GOAT'))
	sp_summary = data.frame(speaker=sp, 
							K_mean_hx=mean(middle_frame_data$hyoid_x[middle_frame_data$phone=='K']), 
							G_mean_hx=mean(middle_frame_data$hyoid_x[middle_frame_data$phone=='G']), 
							K_mean_hy=mean(middle_frame_data$hyoid_y[middle_frame_data$phone=='K']), 
							G_mean_hy=mean(middle_frame_data$hyoid_y[middle_frame_data$phone=='G']))
	data_for_coat_goat_scatterplot = rbind(data_for_coat_goat_scatterplot, sp_summary)
}
data_for_coat_goat_scatterplot$G_K_hx_diff = data_for_coat_goat_scatterplot$G_mean_hx - data_for_coat_goat_scatterplot$K_mean_hx
data_for_coat_goat_scatterplot$G_K_hy_diff = data_for_coat_goat_scatterplot$G_mean_hy - data_for_coat_goat_scatterplot$K_mean_hy

plot(data_for_coat_goat_scatterplot$G_K_hx_diff, data_for_coat_goat_scatterplot$G_K_hy_diff)

# APRIL 15 CONCERNS
#
# some bad tongue traces, especially cr33 tongue tip
# some bad occlusal angles
# some bad stabilization, especially cr11, cr02, cr05


# compare_trajectories(us_data, sp, select_col = 'word', select_vals=c('KEEP','CREEP'), plane='video', signal='lips_area', ref_phone='Kh', ref_edge='phone_end', main=paste(sp,'\n','lips area'))


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

source(paste0(scriptpath,'read_ultrasound_functions.r'))
cairo_pdf('sample_tongue_trajectories_cr01_2025_04_16.pdf', height=5,width=6,onefile=TRUE)
sp = 'cr01'
for (i in 1:11){
	xsignal=paste0('X',i)
	ysignal=paste0('Y',i)	  
	compare_trajectories(us_data, sp, select_col='phrase', select_vals=c('A-COAT','A-CROW'), plane='sag', 
		signal=xsignal, ref_phone='Kh', ref_edge='phone_end', ylim=c(-2,6), main=paste(sp,'\n',xsignal))
	compare_trajectories(us_data, sp, select_col='phrase', select_vals=c('A-COAT','A-CROW'), plane='sag', 
		signal=ysignal, ref_phone='Kh', ref_edge='phone_end', ylim=c(-2,6), main=paste(sp,'\n',ysignal))
}
dev.off()

source(paste0(scriptpath,'read_ultrasound_functions.r'))
cairo_pdf('sample_tongue_trajectories_cr01_gi_2025_04_16.pdf', height=5,width=6,onefile=TRUE)
sp = 'cr01'
for (i in 1:11){
	xsignal=paste0('X',i)
	ysignal=paste0('Y',i)	  
	compare_trajectories(us_data, sp, select_col='phrase', select_vals=c('A-GEEK','THE-GREEN'), plane='sag', 
		signal=xsignal, ref_phone='Gh', ref_edge='phone_end', ylim=c(-2,6), main=paste(sp,'\n',xsignal))
	compare_trajectories(us_data, sp, select_col='phrase', select_vals=c('A-GEEK','THE-GREEN'), plane='sag', 
		signal=ysignal, ref_phone='Gh', ref_edge='phone_end', ylim=c(-2,6), main=paste(sp,'\n',ysignal))
}
dev.off()

cairo_pdf('sample_tongue_trajectories_2025_03_18.pdf', height=12,width=6,onefile=TRUE)
	par(mfrow=c(3,1))
	for (sp in names(us_data)){
		compare_trajectories(us_data, sp, word=c('KEEP','CREEP'), plane='sag', signal='X2', main=paste(sp,'\n','X2'))
		compare_trajectories(us_data, sp, word=c('KEEP','CREEP'), plane='sag', signal='Y6', main=paste(sp,'\n','Y6'))
		compare_trajectories(us_data, sp, word=c('KEEP','CREEP'), plane='video', signal='lips_area', main=paste(sp,'\n','lips area'))
	}
	for (sp in names(us_data)){
		compare_trajectories(us_data, sp, word=c('COAT','CROW'), plane='sag', signal='X2', main=paste(sp,'\n','X2'))
		compare_trajectories(us_data, sp, word=c('COAT','CROW'), plane='sag', signal='Y6', main=paste(sp,'\n','Y6'))
		compare_trajectories(us_data, sp, word=c('COAT','CROW'), plane='video', signal='lips_area', main=paste(sp,'\n','lips area'))
	}
	for (sp in names(us_data)){
		compare_trajectories(us_data, sp, word=c('GOAT','GROAN'), plane='sag', signal='X2', main=paste(sp,'\n','X2'))
		compare_trajectories(us_data, sp, word=c('GOAT','GROAN'), plane='sag', signal='Y6', main=paste(sp,'\n','Y6'))
		compare_trajectories(us_data, sp, word=c('GOAT','GROAN'), plane='video', signal='lips_area', main=paste(sp,'\n','lips area'))
	}
dev.off()



# STEP 3c: put polar data in long format for ssanova
# consonant_data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('T','K') & word%in%c('STREAM','SCREAM')))
consonant_data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('AH1','AA1','IY1') & word%in%c('RUN','ROCK','REED')), polar=TRUE)
consonant_data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('S','K','G') & word%in%c('SEAM','KEEP','GOAT')), polar=TRUE)
consonant_data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('R','AA1','K') & word%in%c('ROCK')), polar=TRUE)


consonant_data = traces_wide_to_long(subset(us_data[[speaker]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('K','Kh','R','OW1') & word%in%c('CROW')), polar=TRUE)

sp = 'cr01'
consonant_data = traces_wide_to_long(subset(us_data[[sp]]$sag$tongue_traces, middle_frame==TRUE & phone%in%c('G') & word%in%c('GEEK','GREEN','GROAN','GOAT')), polar=TRUE)



# STEP 3d: convert long to polar (AGAIN) and interpolate for ssanova
show.traces(consonant_data, data.cat='word', token.label='token_id', origin=us_data[[sp]]$sag$origin, 
	drop_levels=TRUE, #palate=us_data[[speaker]]$sag$palate_trace_rotated, 
	main=sp, is.polar=TRUE, interpolate=100)
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
# 7: 394-45   7cm: mmperpixel = 2*70/350
# 8: 401-43   8cm: mmperpixel = 2*80/359
# 9: 401-43   9cm: mmperpixel = 2*90/359

# these were measured from 640x480 images, but the analysis uses 320x240 images, so I multiplied by 2 to get mmperpixel


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
