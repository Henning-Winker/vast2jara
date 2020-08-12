###############################################################################################################################
## This scripts illustrates the approach of conducting RedList risk assessments
## with VAST - FishLife - JARA for North Sea Thornback skate    
## North Sea IBTS data sourced from DATRAS
###############################################################################################################################
#rm( list = ls() ) 
######## Load required R packages
library( VAST )
library(vast2jara) # has data sets
# devtools::install_github("Henning-Winker/vast2jara")
library(JARA) 
# devtools::install_github("Henning-Winker/JARA")
library(FishLife) 

######## Load sampling data
data(datras2vast)
#RootDir = "F:/Meta_analysis_project/North_Atlantic_allspecies/"
RootDir = "C:/Users/henni/Dropbox/VAST2JARA/NorthSea"
setwd( RootDir )
# Load Survey data
survey_data = nsibts_data 
# Prepare Grid for to use Region = "user" option
grid_extrapolation = nsibts_grid # Define input Grid
Data_Extrap <- as.data.frame( grid_extrapolation$Data_Extrap )
input_grid <- as.data.frame( cbind( Data_Extrap$Lon, Data_Extrap$Lat, as.vector( grid_extrapolation$Area_km2_x ) ) )
names( input_grid ) <- c( "Lon", "Lat", "Area_km2" )

# define species
all.species =unique(survey_data$Category)
nsp <- length( unique( survey_data$Category ) )
# Run check on nominal CPUE


#*********************************************************************
# Select species i here
#********************************************************************
splist= as.data.frame(all.species)
splist
i = c(1,22,9,25)[1] 
species = survey_data$Category[i]
species

# find VAST package
TmbDir = system.file("executables", package="VAST")

test_fit = TRUE # several species hitting the bounds perhaps issue with 2010?
if(i ==25){
  test_fit = FALSE
}
# Non convergens
# e.g. i = 11; Leucoraja naevus
######## Loop over species 
#### Extract data for the species under consideration 
  sampling_data <- survey_data[survey_data$Category == species,]
	OutputDir = paste0( RootDir, '/Output_', species, '/' )
	dir.create( OutputDir,showWarnings = F )
	
  
  #observations_LL <- as.data.frame( cbind( "Lon" = sampling_data$Lon, "Lat" = sampling_data$Lat ) )

	#### Make settings
	settings = make_settings( n_x = 100, Region = "User", purpose = "index2",
		RhoConfig = c( "Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0 ), 
		FieldConfig = matrix (
		  data = c( Omega1 = "IID", Epsilon1 = "IID", Beta1 = "IID", 
		            Omega2 = "0", Epsilon2 = "0", Beta2 = "IID" ), nrow = 3, ncol = 2 ),
		use_anisotropy = TRUE, fine_scale = FALSE )

	#### Run model
	setwd( OutputDir )
	# Running on  test_fit = FALSE due to convergence issues
	fit = fit_model( settings = settings, Lat_i = sampling_data[,'Lat'], Lon_i = sampling_data[,'Lon'], 
		t_i = sampling_data[,'Year'], c_i = rep( 0, nrow( sampling_data ) ), 
		b_i = sampling_data[,'Catch_KG'], a_i = sampling_data[,'AreaSwept_km2'], 
		v_i = sampling_data[,'Vessel'], input_grid = input_grid,
		CompileDir = TmbDir,test_fit=test_fit)
	
	
	Save = list( "Opt" = fit, "Report" = fit$Report, "ParHat" = fit$ParHat, "DataFrame" = fit$data_frame )
	save( Save, file = paste0( OutputDir,species,".RData" ) )

	#### Plot results
	results = plot_results(fit, check_residuals = FALSE, plot_set = 3 )
	MapDetails_List <- results$map_list
	save( MapDetails_List, file = paste0( OutputDir,species, "MapDetails_List.RData" ) )

  #------------------------------------------------------------
	# Get Generation Length from FishLife
	#------------------------------------------------------------
	# get generation time from FishLife
	
	sp = unlist(strsplit(as.character(species)," "))
	taxa = Search_species(Genus=sp[1],Species = sp[2],add_ancestors=TRUE)$match_taxonomy
	png(file = paste0(OutputDir,sp[1],".",sp[2],".traits.png"), width = 6, height = 7, 
	    res = 200, units = "in")
	predfl =Plot_taxa(taxa,mfrow=c(3,2))
	dev.off()
	
	GL = predfl[[5]]$Mean_pred["G"] 
	
	#---------------------------------------------------
	# JARA section
	#---------------------------------------------------
	
	# prepare input data for JARA
	jd = read.csv("Table_for_SS3.csv")
	Year = min(jd$Year):max(jd$Year)
	# Define index and standard error
	index  = data.frame(Year,Index=NA)
	se  = data.frame(Year,Index=0.1)
	# Deal with missing years
	index[Year%in%jd$Year,2] = jd$Estimate_metric_tons
	se[Year%in%jd$Year,2] = jd$SD_log
	
  # build JARA model
	jarainput = build_jara(I=index,se=se,GL=GL,
	                       assessment=species,scenario = paste0("example",i),
	                       pk.prior = c(0.5,0.1) # limit exponentially increasing projections
	                       ) 
  # check input data
	jrplot_indices(jarainput)
	# save
	jrplot_indices(jarainput,as.png = T)
	
	# run jara
 	jara = fit_jara(jarainput,output.dir = getwd(),save.jara = T)
  
 	# check some plots
 	jrplot_fits(jara)
 	jrplot_logfits(jara)
 	# Save all
 	jara_plots(jara)
  
 	# Make multiplot Figure
 	plname = paste0("jara.",sp[1],".",sp[2])
 	pwidth = 8
 	pheight = 8.5
 	res=300
 	windows(width=pwidth,height=pheight) # Sorry Rich
 	jrpar(mfrow=c(3,2),labs=T,plot.cex=1)
 	jrplot_indices(jarainput,add=T)
 	mtext(text=paste0(letters[1],")"), xpd=NA, side=3, adj=adj, cex=1)
 	jrplot_trjfit(jara,add=T)
 	mtext(text=paste0(letters[2],")"), xpd=NA, side=3, adj=adj, cex=1)
 	jrplot_poptrj(jara,add=T)
 	mtext(text=paste0(letters[3],")"), xpd=NA, side=3, adj=adj, cex=1)
 	jrplot_changes(jara,add=T)
 	mtext(text=paste0(letters[4],")"), xpd=NA, side=3, adj=adj, cex=1)
 	jrplot_state(jara,add=T)
 	mtext(text=paste0(letters[5],")"), xpd=NA, side=3, adj=adj, cex=1)
 	jrplot_iucn(jara,add=T)
 	mtext(text=paste0(letters[6],")"), xpd=NA, side=3, adj=adj, cex=1)
 	mtext(text=paste0(letters[6],")"), xpd=NA, side=3, adj=adj, cex=1)
 	mtext(species, side=3, outer=T,line= .5,cex=1.,c(0.5))
 	
 	dev.print(tiff,paste0(getwd(),"/",plname,"_hires.tiff"), width = pwidth, height = pheight, res = res, units = "in")
 	dev.print(jpeg,paste0(getwd(),"/",plname,".jpg"), width = pwidth, height = pheight, res = res, units = "in")
 	
 
  # Do some extra Hindcasting
  hc = jara_hindcast(jarainput,peels = 0:10,save.jarafile = F)
  
  # make plot
  plname = paste0("jarahindcast",sp[1],".",sp[2])
  pwidth = 8
  pheight = 9
  res=300
  windows(width=pwidth,height=pheight) 
  jrpar(mfrow=c(2,1),labs=T,plot.cex=1,mai = c(0.6, 0.6, 0.1, 0.5))
  jrplot_retrobias(hc,add=T)
  mtext(text=paste0(letters[1],")"), xpd=NA, side=3, adj=adj, cex=1)
  jrplot_retroiucn(hc,add=T)
  mtext(text=paste0(letters[2],")"), xpd=NA, side=3, adj=adj, cex=1)
  mtext(paste0("Retrospective Hindcast ",species), side=3, outer=T,line= .5,cex=1.,c(0.5))
  
  dev.print(tiff,paste0(getwd(),"/",plname,"_hires.tiff"), width = pwidth, height = pheight, res = res, units = "in")
  dev.print(jpeg,paste0(getwd(),"/",plname,".jpg"), width = pwidth, height = pheight, res = res, units = "in")
  
  # The end for now
  
