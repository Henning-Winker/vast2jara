###############################################################################################################################
##
## This script runs a VAST model for all species of the North Atlantic 
##
###############################################################################################################################
#rm( list = ls() ) 

######## Load required R packages
library( VAST )
library(vast2jara)
######## Load sampling data
RootDir = "F:/Meta_analysis_project/North_Atlantic_allspecies/"
RootDir = "C:/Users/henni/Dropbox/DemJARA/North_Atlantic_allspecies"
setwd( RootDir )
load( "sampling_data_allspecies.rdata",verbose=T)
sampling_data_allspecies = nsibts_data 
species =unique(sampling_data_allspecies$Category)
ids = c(1,9:18,23:25,29,31)
nsp <- length( unique( sampling_data_allspecies$Category ) )
TmbDir = system.file("executables", package="VAST")
######## Loop over species 
for( i in 10 : nsp ) {
if(i%in%ids){ 
	#### Extract data for the species under consideration 
	sampling_data <- sampling_data_allspecies[sampling_data_allspecies$Category == sampling_data_allspecies$Category[i],]
	OutputDir = paste0( RootDir, '/Output_', sampling_data_allspecies$Category[i], '/' )
	dir.create( OutputDir )
	
	#observations_LL <- as.data.frame( cbind( "Lon" = sampling_data$Lon, "Lat" = sampling_data$Lat ) )

	#### Make settings
	settings = make_settings( n_x = 100, Region = "NorthSeaIBTS", purpose = "index2",
		RhoConfig = c( "Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0 ), 
		FieldConfig = matrix (
		  data = c( Omega1 = "IID", Epsilon1 = "IID", Beta1 = "IID", 
		            Omega2 = "0", Epsilon2 = "0", Beta2 = "IID" ), nrow = 3, ncol = 2 ),
		use_anisotropy = TRUE, fine_scale = FALSE )

	#### Run model
	setwd( OutputDir )
	fit = fit_vast( settings = settings, Lat_i = sampling_data[,'Lat'], Lon_i = sampling_data[,'Lon'], 
		t_i = sampling_data[,'Year'], c_i = rep( 0, nrow( sampling_data ) ), 
		b_i = sampling_data[,'Catch_KG'], a_i = sampling_data[,'AreaSwept_km2'], 
		v_i = sampling_data[,'Vessel'], CompileDir = TmbDir,test_fit=FALSE)
	Save = list( "Opt" = fit, "Report" = fit$Report, "ParHat" = fit$ParHat, "DataFrame" = fit$data_frame )
	save( Save, file = paste0( OutputDir, "Save.RData" ) )

	#### Plot results
	results = plot( fit, check_residuals = FALSE, plot_set = 3 )
	MapDetails_List <- results$map_list
	save( MapDetails_List, file = paste0( OutputDir, "MapDetails_List.RData" ) )
}
######## End the loop over species
}

