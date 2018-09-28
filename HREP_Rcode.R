#### Culvert Evaluation Model ####
#
# Updated:
# 07/26/18, Allison Truhlar
#
# This script is based on the culvert evaluation model developed by Rebecca Marjerison in 2013, and the accompanying
# Python scripts developed by David Gold in 2015.
#
# This script will:
# 1. Determine the runoff peak discharge of given culvert's watershed using the SCS graphical curve number method.
# 2. Calculate the cross sectional area of each culvert and assign c and Y coefficients based on culvert characteristics
#    (FHWA engineering pub HIF12026, appendix A).
# 3. Determine the maximum capacity of a culvert assuming submerged inlet control (Equation A.3, FHWA engineering pub 
#    HIF12026, appendix A), and the headwater at the culvert inlet is at the height of the road surface.
# 4. Determine the maximum return period storm that the culvert can safely pass both current and future rainfall conditions,
#    by comparing the maximum capacity of the culvert to peak flows for the 1, 2, 5, 10, 25, 50, 100, 200 and 500 yr 24-hr 
#    storm events.
#
# Inputs:
# 1. NAACC data: A CSV file containing culvert characteristics downloaded from the North Atlantic Aquatic Connectivity
#    Collaborative's database (https://www.streamcontinuity.org/cdb2/naacc_search_crossing.cfm).
#
# 2. Precipitation data: A CSV file containing extreme precipitation data (in) for the 1, 2, 5, 10, 25, 50, 100, 200 and 500 
#    yr 24-hr storm events, downloaded from the Northeast Regional Climate Center (http://precip.eas.cornell.edu/)
#
# 3. A 10m DEM for the area containing the target culverts, downloaded from the NYS GIS Clearinghouse (https://gis.ny.gov/)
#
# 4. A flow direction raster for the area containing the target culverts, created by applying the Fill and Flow Direction 
#    tools available in ArcMap 10.3.1 to a version of the 10m DEM that has been 'burned' with National Hydrography Dataset
#    flowlines (https://nhd.usgs.gov/data.html)
#    (e.g., http://www.oberlin.edu/OCTET/HowTo/GIS/Tutorial/SpatialAnalysist-Hydrology/Hydrologic%20Modeling%20Using%20GIS.doc)
#
# 5. A flow accumulation raster for the area containing the target culverts, created by applying the Flow Accumulation tool
#    available in ArcMap 10.3.1 to the flow direction raster.
#
# 6. A flow length raster for the the area containing the target culverts, created by applying the Flow Length tool available
#    in ArcMap 10.3.1 to the flow direction raster.
#
# 7. A slope raster for the area containing the target culverts, created by applying the Fill and Slope tools available in
#    ArcMap 10.3.1 to the original 10m DEM.
#
# 8. A curve number raster for the area containing the target culverts, resampled to the resolution of the DEM. A CN raster  
#    created by Rebecca Marjerison using land use data from the 2006 National Land Cover Database and soil characteristics from 
#    the 2010 State Soil Geographic database, is available for NYS.
#
# Outputs:
# 1. Model output: A CSV file containing the peak discharge for each culvert's watershed for the analyzed return period storms 
#    under (1) current and (2) projected 2050 (15% increase in 24-hr storm event precipitation) rainfall conditions, (3) watershed
#    area flowing to each culvert, (4) average CN for the watershed of each culvert, (5) time of concentration for each culvert
#    watershed, (5) the cross-section area of each culvert inlet, (6) the individual capacity of each culvert, (7) the total
#    capacity of each crossing location (equal to the individual culvert capacity if only one culvert at the crossing), and (8) the
#    maximum storm size the crossing can safely pass under current and future rainfall conditions.
#
# 2. Rejected: A CSV file containing crossings that were rejected based solely on the crossing characteristics, prior to the 
#    watershed delineation step.

#####

library(reticulate)
library(BBmisc)

# Clear workspace
rm(list=ls())

#### Import Python modules ####
arcpy <- import("arcpy")
arcpy.sa <- import("arcpy.sa")

#Retrieve Spatial Analyst licencse
arcpy$CheckOutExtension("Spatial")

#Set overwrite to TRUE
arcpy.env <- py_get_attr(arcpy,"env")
py_set_attr(arcpy.env,"overwriteOutput","TRUE")


#### Set working directory ####
dir <- readline(prompt="Enter path to working directory:")
#Z:\Culverts_Allison\For_HREP\R_code\TEST\
#dir <- "C:\\Users\\Public\\Documents\\2015_Culvert\\REN_test" #PC
#dir <- "/Volumes/WRILab/Culverts_Allison/For_HREP/R_code/TEST/" #Mac
setwd(dir)


#### Load NAACC data ####
ws <- readline(prompt="Enter a three letter watershed or county code, in all caps (e.g., REN for Rensselaer):")
#ws <- "REN"
extracted_data <- read.csv(paste(dir,"Data\\",ws,".csv",sep="")) #PC
#NAACC_data <- read.csv(paste(dir,"Data/",ws,".csv",sep="")) #Mac

# colnames <- c("Survey_Id","GPS_Y_Coordinate","GPS_X_Coordinate","Crossing_Type","Inlet_Type",
#               "Number_Of_Culverts","Road_Fill_Height","Crossing_Structure_Length","Inlet_Height",
#               "Inlet_Structure_Type","Inlet_Substrate_Water_Width","Inlet_Water_Depth","Inlet_Width",
#               "Material","Slope_Percent")
# extracted_data <- NAACC_data[,colnames]

#### Create field data and rejected data files ####

# Create modeling notes column
extracted_data[,"Modeling_notes"] <- "NA"

# Create empty "reject" data frame
rejected_data <- data.frame(matrix(ncol=ncol(extracted_data),nrow=1))
colnames(rejected_data) <- colnames(extracted_data)

# Ensure other character variables are read as characters when transfered to rejected data frame
extracted_data[,"Crossing_Type"] <- as.character(extracted_data[,"Crossing_Type"])
extracted_data[,"Inlet_Type"] <- as.character(extracted_data[,"Inlet_Type"])
extracted_data[,"Inlet_Structure_Type"] <- as.character(extracted_data[,"Inlet_Structure_Type"])
extracted_data[,"Material"] <- as.character(extracted_data[,"Material"])

# Replace "No data" with "NA"
extracted_data[extracted_data=="No data"] <- "NA"

# Assume negative slopes are equal to zero
extracted_data$Slope_Percent[extracted_data$Slope_Percent<0] <- 0

# Start indexing for rejected data
n <- 1

for (i in 1:nrow(extracted_data)){
  
  #Remove negative widths, heights, road fill heights, and lengths from extracted_data and save in rejected_data
  
  if (extracted_data$Inlet_Width[i]<0){
    extracted_data$Modeling_notes[i] <- "Negative inlet width"
    rejected_data[n,] <- extracted_data[i,]
    extracted_data$Inlet_Width[i] <- NA
    n <- n+1
  } else if (extracted_data$Inlet_Height[i]<0){
    extracted_data$Modeling_notes[i] <- "Negative inlet height"
    rejected_data[n,] <- extracted_data[i,]
    extracted_data$Inlet_Height[i] <- NA
    n <- n+1
  } else if (extracted_data$Road_Fill_Height[i]<0){
    extracted_data$Modeling_notes[i] <- "Negative road fill height"
    rejected_data[n,] <- extracted_data[i,]
    extracted_data$Road_Fill_Height[i] <- NA
    n <- n+1
  } else if (extracted_data$Crossing_Structure_Length[i]<0){
    extracted_data$Modeling_notes[i] <- "Negative crossing structure length"
    rejected_data[n,] <- extracted_data[i,]
    extracted_data$Crossing_Structure_Length[i] <- NA
    n <- n+1
  } else if (extracted_data$Crossing_Type[i]=="Bridge" & !(extracted_data$Inlet_Structure_Type[i]=="Box/Bridge with Abutments" || extracted_data$Inlet_Structure_Type[i]=="Open Bottom Arch Bridge/Culvert")){
    # If no negative dimensions, then proceed to checking the crossing type and inlet structure type
    extracted_data$Modeling_notes[i] <- "Bridge crossing type AND inlet structure type is NOT Box/Bridge with Abutments OR Open Bottom Arch Bridge/Culvert"
    rejected_data[n,] <- extracted_data[i,]
    n <- n+1
  } else if (extracted_data$Crossing_Type[i]=="Bridge" & (extracted_data$Inlet_Structure_Type[i]=="Box/Bridge with Abutments" || extracted_data$Inlet_Structure_Type[i]== "Open Bottom Arch Bridge/Culvert") & extracted_data$Inlet_Width[i]>=20){
    # If no negative dimensions, then proceed to checking the crossing type and inlet structure type
    extracted_data$Modeling_notes[i] <- "OK Bridge crossing type, but inlet width greater than or equal to 20 ft"
    rejected_data[n,] <- extracted_data[i,]
    n <- n+1
  }else{
    extracted_data$Modeling_notes[i] <- "valid"
  }
  
}

model_data <- extracted_data[extracted_data$Modeling_notes=="valid",]
colnames(model_data) <- colnames(extracted_data)

# Order data by Survey ID. Two culverts located at the same site will have the same survey ID
model_data <- model_data[order(model_data$Survey_Id),]

# Write separate txt files for each row in model_data, to be read separately into ArcGIS later
for (i in 1:nrow(model_data)){
  write.table(model_data[i,],row.names=FALSE,file=paste(dir,"Data\\SurveyID_",model_data$Survey_Id[i],".txt",sep=""),sep=",")
}

write.csv(rejected_data,row.names=FALSE,file=paste(dir,"Data\\",ws,"_rejected.csv",sep=""))

#### Read and format NRCC extreme precip data (24-hr storm) ####
precip_raw <- read.csv(paste(dir,"Data\\",ws,"_precip.csv",sep=""),header=FALSE)
precip_in <- precip_raw[11:19,11]
precip_cm <- as.numeric(as.character(precip_in))*2.54 #convert in to cm
precip_cm <- c(precip_cm, precip_cm*1.15)
P <- precip_cm

#Add columns to factor_values for q_peaks
df <- data.frame(matrix(ncol = 18, nrow = nrow(model_data)))
colnames <- c("1-yr, current","2-yr, current","5-yr, current","10-yr, current", "25-yr, current",
              "50-yr, current","100-yr, current","200-yr, current","500-yr, current",
              "1-yr, future","2-yr, future","5-yr, future","10-yr, future", "25-yr, future",
              "50-yr, future","100-yr, future","200-yr, future","500-yr, future")
colnames(df) <- colnames
model_data <- cbind(model_data,df)

#### Load GIS files ####
flow_dir <- paste(dir,"GIS_files\\Merge_county_DEMs\\",ws,"\\fdir.tif",sep="")
flow_accum <- paste(dir,"GIS_files\\Merge_county_DEMs\\",ws,"\\facc.tif",sep="")
flow_length <- paste(dir,"GIS_files\\Merge_county_NHD\\",ws,"\\flolen",sep="")
slope <- paste(dir,"GIS_files\\Merge_county_DEMS\\",ws,"\\slope",sep="")
CN_file <- paste(dir, "GIS_files\\Merge_county_CN\\",ws,"\\CNres.tif",sep="")
DEM <- paste(dir,"GIS_files\\Merge_county_DEMs\\",ws,"\\DEM.img",sep="")

#### Start loop ####

#### (1) watershed delineation and characterization ####
for(i in 1:nrow(model_data)){

  # Only complete watershed delineation and peak flow for one culvert at any Survey ID
  # Will later complete capacity calculations for each culvert individually to determine maximum storm size
  
  if(i==1 || model_data$Survey_Id[i]!=model_data$Survey_Id[i-1]){
    
    # Create feature from field data coords 
    # create point layer in GCS_WGS_1984 (field data coord system)
    spRef <- arcpy$SpatialReference("WGS 1984")
    temp_pts <- arcpy$MakeXYEventLayer_management(paste(dir,"Data\\SurveyID_",model_data$Survey_Id[i],".txt",sep=""),"GPS_X_Coordinate", "GPS_Y_Coordinate","temp_pts",spRef)
    temp_pts <- arcpy$FeatureClassToFeatureClass_conversion(temp_pts,paste(dir,"GIS_files\\XY_points\\",sep=""),paste("SurveyID_",model_data$Survey_Id[i],"_XY",sep=""))
    # project to UTM
    spRef_utm <- arcpy$SpatialReference("NAD 1983 UTM ZONE 18N")
    Points_UTM <-  arcpy$Project_management(temp_pts,paste(dir,"GIS_files\\XY_points\\SurveyID_",model_data$Survey_Id[i],"_XY_UTM.shp",sep=""),spRef_utm)
    
    #Snap pour point, using 20m snap distance
    snap_code <- arcpy.sa$SnapPourPoint(Points_UTM, flow_accum, 20,"Survey_Id")
    snap_out <- paste(dir,"GIS_files\\Temp\\Snap\\snap_ws",model_data$Survey_Id[i],".img",sep="")
    snap_code$save(snap_out)
    
    #Delineate watershed
    #ws_raster_file <- paste(dir,"ArcGIS_files\\Test_locations\\Temp\\ws_raster.tif",sep="")
    ws_raster <- arcpy.sa$Watershed(flow_dir,snap_code,"VALUE")
    #ws_raster$save(ws_raster_file)
    
    #Convert watershed raster to polygon and dissolve to ensure only one entry in attribute table
    ws_temp_out <- paste(dir,"GIS_files\\Temp\\WS\\ws",model_data$Survey_Id[i],".shp",sep="")
    ws_temp <- arcpy$RasterToPolygon_conversion(ws_raster,ws_temp_out,"SIMPLIFY","Value")
    ws_out <- paste(dir,"GIS_files\\WS_Poly\\ws",model_data$Survey_Id[i],".shp",sep="")
    ws_poly <- arcpy$Dissolve_management(ws_temp,ws_out,"GRIDCODE")
    
    #Add fields Area_sqkm, Tc_hr, CN to the watershed polygon
    arcpy$AddField_management(ws_poly,"Area_sqkm","FLOAT")
    arcpy$AddField_management(ws_poly,"Tc_hr","FLOAT")
    arcpy$AddField_management(ws_poly,"CN","FLOAT")
    
    #Calculate watershed area
    ws_area <- arcpy$CalculateField_management(ws_poly,"Area_sqkm","!SHAPE.area@SQUAREKILOMETERS!","PYTHON")
    cursor <- arcpy$da$SearchCursor(ws_poly,"Area_sqkm")
    model_data$WS_area[i] <- unlist(iterate(cursor))
    
    #Clip CN raster to watershed extent
    CN_clip_out <- paste(dir,"GIS_files\\CN_clip.tif",sep="")
    desc <- arcpy$Describe(CN_file)
    extent <- py_str(py_get_attr(desc,"extent"))
    extent_split <- strsplit(extent, " ")
    extent_split <- unlist(extent_split,use.names=FALSE)
    extent <- paste(extent_split[1], extent_split[2], extent_split[3], extent_split[4])
    CN_clip <- arcpy$Clip_management(CN_file,extent,CN_clip_out,ws_poly,"#","ClippingGeometry","NO_MAINTAIN_EXTENT")
    
    #Calculate mean CN for watershed
    model_data$mean_CN[i] <- as.numeric(py_str(arcpy$GetRasterProperties_management(CN_clip,"MEAN","")))
    #arcpy$CalculateField_management(ws_poly, "CN", CN,"PYTHON_9.3")
    #cursor <- arcpy$da$SearchCursor(ws_poly,"CN")
    #factor_values$CN[i] <-unlist(iterate(cursor))
    
    #Resample CN raster to size of DEM - this is done ahead now
    # CN_clip_res_out <- paste(dir,"GIS_files\\CN_clip_res.tif",sep="")
    # CellSizeX <- arcpy$GetRasterProperties_management(DEM,"CELLSIZEX")
    # CellSizeY <- arcpy$GetRasterProperties_management(DEM,"CELLSIZEY")
    # CN_clip_res <- arcpy$Resample_management(CN_clip,CN_clip_res_out,paste(CellSizeX," ",CellSizeY,sep=""),"NEAREST")
    
    #Clip flow length to watershed area
    flo_out <- paste(dir,"GIS_files\\Temp\\flo_temp.tif",sep="")
    desc <- arcpy$Describe(ws_poly)
    extent <- py_str(py_get_attr(desc,"extent"))
    extent_split <- strsplit(extent, " ")
    extent_split <- unlist(extent_split,use.names=FALSE)
    extent <- paste(extent_split[1], extent_split[2], extent_split[3], extent_split[4])
    arcpy$Clip_management(flow_length, extent, flo_out, ws_poly,"-3.402823e+038", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    
    #Clip slope to watershed area
    slope_out <- paste(dir,"GIS_files\\Temp\\slope_temp.tif",sep="")
    arcpy$Clip_management(slope, extent, slope_out, ws_poly,"-3.402823e+038", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    
    #Get maximum flow length
    MaxLen <- as.numeric(py_str(arcpy$GetRasterProperties_management(flo_out, "MAXIMUM", "")))
    
    #Get average slope
    mean_slope <- as.numeric(py_str(arcpy$GetRasterProperties_management(slope_out, "MEAN", "")))
    
    #Calculate time of concentration
    model_data$tc[i] <- 0.000325*(MaxLen^0.77)*((mean_slope/100)^-0.385)
    
    #### (2) peak flow calculation ####
    
    Tc <- model_data$tc[i]
    CN <- model_data$mean_CN[i]
    WS_area <- model_data$WS_area[i]
    
    # Calculate q_peak
    # Skip any watershed with CN or Tc = 0, or WS_area < 0.01 km
    if(WS_area < 0.01){
      model_data$Modeling_notes[i] <- "Skipped, WS_area < 0.01 km sq"
      model_data[i,colnames] <- rep(NA,length(colnames))
    } else if (Tc == 0){
      model_data$Modeling_notes[i] <- "Skipped, Tc = 0"
      model_data[i,colnames] <- rep(NA,length(colnames))
    } else if (CN == 0){
      model_data$Modeling_notes[i] <- "Skipped, CN = 0"
      model_data[i,colnames] <- rep(NA,length(colnames))
    } else {
      
      # calculate storage, S  and Ia in cm
      Storage = 0.1 * ((25400 / CN) - 254.0) #cm
      Ia = 0.2*Storage #inital abstraction, amount of precip that never has a chance to become runoff (cm)
      
      # calculate depth of runoff from each storm
      # if P < Ia NO runoff is produced
      # Note that P is a vector of the 9 values, so everything hereafter is too.
      Pe = (P - Ia) #cm
      Pe[Pe < 0] <- 0 # get rid of negative Pe's
      Q = (Pe ^ 2) / (P + (Storage - Ia)) #cm
      
      # calculate q_peak, cubic meters per second
      # q_u is an adjustment because these watersheds are very small. It is a function of tc,
      # and constants Const0, Const1, and Const2 which are in turn functions of Ia/P (rain_ratio) and rainfall type
      # We are using rainfall Type II because that is applicable to most of New York State
      # rain_ratio is a vector with one element per input return period
      rain_ratio = Ia / P
      rain_ratio[rain_ratio<0.1] <- 0.1 # keep rain ratio within limits set by TR55
      rain_ratio[rain_ratio>0.5] <- 0.5 # keep rain ratio within limits set by TR55
      #Calculate constants (NRCS, 1986)
      Const0 = (rain_ratio ^ 2) * -2.2349 + (rain_ratio * 0.4759) + 2.5273
      Const1 = (rain_ratio ^ 2) *  1.5555 - (rain_ratio * 0.7081) - 0.5584
      Const2 = (rain_ratio ^ 2) *  0.6041 + (rain_ratio * 0.0437) - 0.1761
      #Calculate q_u, m^3/s per km^2 per cm, equation from Chin (2013)
      qu = 10 ^ (Const0 + Const1 * log10(Tc) + Const2 * (log10(Tc) ^ 2) - 2.366)
      #Caculate q_peak
      q_peak = Q * qu * WS_area #m^3/s
      #Add q_peak tp factor_values dataframe
      model_data[i,colnames] <- q_peak
      
    }
    
    # If survey ID is the same as the previous culvert, assign same values for WS area, CN, Tc, and peak flows
  } else if (model_data$Survey_Id[i]==model_data$Survey_Id[i-1]){
    model_data$WS_area[i] <- model_data$WS_area[i-1]
    model_data$mean_CN[i] <- model_data$mean_CN[i-1]
    model_data$tc[i] <- model_data$tc[i-1]
    model_data[i,colnames] <- model_data[i-1,colnames]
    model_data$Modeling_notes[i] <- model_data$Modeling_notes[i-1]
    
  }
  
  #### (3) culvert capacity calculations ####

  # Modify data needed for culvert capacity calculation
  length <- model_data$Crossing_Structure_Length[i]/3.2808 # converts culvert length from feet to meters
  culvert_slope <- model_data$Slope_Percent[i]/100 # converts slope from percent to meter/meter
  inlet_width <- model_data$Inlet_Width[i]/3.2808 # converts inlet width from feet to meters
  
  if (model_data$Inlet_Structure_Type[i]!="Round Culvert"){
    inlet_height <- model_data$Inlet_Height[i]/3.2808  # if culvert is not round, need inlet height, and convert from feet to meters
  }
  
  # Calculate area of the culvert inlet
  if (model_data$Inlet_Structure_Type[i] == "Round Culvert"){
    model_data$inlet_area_sqm[i] <- (inlet_width/2)^2*pi #Area in m^2
    depth <- inlet_width # if culvert is round, depth is diameter
  } else if (model_data$Inlet_Structure_Type[i] == "Pipe Arch/Elliptical Culvert"){
    model_data$inlet_area_sqm[i] <- (inlet_width/2)*(inlet_height/2)*pi
    depth <- inlet_height # if culvert is eliptical, depth is inlet height
  } else if (model_data$Inlet_Structure_Type[i] == "Box Culvert" | model_data$Inlet_Structure_Type[i] == "Box/Bridge with Abutments"){
    model_data$inlet_area_sqm[i] <- inlet_width*inlet_height
    depth <- inlet_height # if culvert is a box, depth is inlet height
  } else if (model_data$Inlet_Structure_Type[i] == "Open Bottom Arch Bridge/Culvert"){
    model_data$inlet_area_sqm[i] <- ((inlet_width/2)*(inlet_height/2)*pi)/2
    depth <- inlet_height # if culvert is an arch, depth is inlet height
  }
    
  # Calculate head over the culvert invert (bottom of the culvert) by adding road fill height to depth (the distance to the bottom of the culvert)
  HW = model_data$Road_Fill_Height[i]/3.2808 + depth # head in meters
  
  # assign ks (slope coefficient from FHWA engineering pub HIF12026, appendix A)
  if (model_data$Inlet_Structure_Type[i] == 'Mitered to Slope'){
    ks <- 0.7
  } else{
    ks <- -0.5
  }
  
  #Ku constant, unit conversion for SI
  ku <- 1.811
  
  # assign c and y values (coefficients based on shape and material from FHWA engineering pub HIF12026, appendix A) 
  # Calculate capacity, submerged inlet control, Equation A.3 from FHWA engineering pub HIF12026, appendix A
  # Capacity is for the indivual culvert, Qi, assuming an inlet headwater depth equal to the road fill height, m^3/s
  if (model_data$Inlet_Structure_Type[i] == "Open Bottom Arch Bridge/Culvert"){
    
    if (model_data$Material[i] == "Concrete" | model_data$Material[i] == "Stone"){
      
      if (model_data$Inlet_Type[i]=="Headwall" | model_data$Inlet_Type[i]=="Projecting"){
        #Table A.4 - concrete open-bottom arch with headwall, 2:1 span to rise
        c <- 0.041
        Y <- 0.570
        model_data$Culvert_notes[i] <- "Table A.4 - concrete open-bottom arch with headwall, 2:1 span to rise"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else if (model_data$Inlet_Type[i]=="Mitered to Slope"){
        #Table A.4 - concrete open-bottom arch, mitered to slope, 2:1 span to rise
        c <- 0.040
        Y <- 0.48
        model_data$Culvert_notes[i] <- "Table A.4 - concrete open-bottom arch, mitered to slope, 2:1 span to rise"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else if (model_data$Inlet_Type[i]=="Wingwalls" | model_data$Inlet_Type[i]=="Headwall and Wingwalls"){
        #Table A.4 - concrete open-bottom arch, headwall with wingwalls, 2:1 span to rise
        c <- 0.040
        Y <- 0.62
        model_data$Culvert_notes[i] <- "Table A.4 - concrete open-bottom arch, headwall with wingwalls, 2:1 span to rise"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else { 
        c <- NA
        Y <- NA
        model_data$Culvert_notes[i] <- "No c, Y constants currently available for this culvert material/shape combination"
        model_data$Qi_m3_s[i] <- NA
      }
      
    } else if (model_data$Material[i] == "Plastic" | model_data$Material[i] == "Metal"){
      
      if (model_data$Inlet_Type[i]=="Headwall" | model_data$Inlet_Type[i]=="Wingwalls" | model_data$Inlet_Type[i]=="Headwall and Wingwalls"){
        #Table A.6 - embedded elliptical, square headwall
        c <- 0.0431
        Y <- 0.610
        model_data$Culvert_notes[i] <- "Table A.6 - embedded elliptical, square headwall"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else if (model_data$Inlet_Type[i]=="Mitered to Slope"){
        #Table A.6 - embedded elliptical, mitered end
        c <- 0.0541
        Y <- 0.5
        model_data$Culvert_notes[i] <- "Table A.6 - embedded elliptical, mitered end"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else if (model_data$Inlet_Type[i]=="Projecting"){
        #Table A.6 - embedded elliptical, projecting end, channelized (as opposed to ponded)
        c <- 0.649
        Y <- 0.12
        model_data$Culvert_notes[i] <- "Table A.6 - embedded elliptical, projecting end, channelized (as opposed to ponded)"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else { 
        c <- NA
        Y <- NA
        model_data$Culvert_notes[i] <- "No c, Y constants currently available for this culvert material/shape combination"
        model_data$Qi_m3_s[i] <- NA
      }
    } else { 
      c <- NA
      Y <- NA
      model_data$Culvert_notes[i] <- "No c, Y constants currently available for this culvert material/shape combination"
      model_data$Qi_m3_s[i] <- NA
    }
    
  } else if (model_data$Inlet_Structure_Type[i] == "Box Culvert" | model_data$Inlet_Structure_Type[i] == "Box/Bridge with Abutments"){
    
    if (model_data$Material[i] == "Concrete" | model_data$Material[i] == "Stone"){
      #Table A.1 - rectangular concrete, side tapered, more favorable edges
      c <- 0.0378
      Y <- 0.870
      model_data$Culvert_notes[i] <- "Table A.1 - rectangular concrete, side tapered, more favorable edges"
      model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
      
    } else if (model_data$Material[i] == "Plastic" | model_data$Material[i] == "Metal"){
      
      if (model_data$Inlet_Type[i] == "Headwall"){
        #Table A.1 - Box, corrugated metal, headwall
        c <- 0.0379
        Y <- 0.69
        model_data$Culvert_notes[i] <- "Table A.1 - Box, corrugated metal, headwall"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else { 
        c <- NA
        Y <- NA
        model_data$Culvert_notes[i] <- "No c, Y constants currently available for this culvert material/shape combination"
        model_data$Qi_m3_s[i] <- NA
      } 
    } else { 
      c <- NA
      Y <- NA
      model_data$Culvert_notes[i] <- "No c, Y constants currently available for this culvert material/shape combination"
      model_data$Qi_m3_s[i] <- NA
    }  
    
  } else if (model_data$Inlet_Structure_Type[i] == "Pipe Arch/Elliptical Culvert"){
    
    if (model_data$Material[i] == "Plastic"| model_data$Material[i] == "Metal"){
      
      if (model_data$Inlet_Type[i] == "Projecting"){
        #Table A.1 - elliptical face, tapered inlet, thin edge projecting
        c <- 0.0598
        Y <- 0.75
        model_data$Culvert_notes[i] <- "Table A.1 - elliptical face, tapered inlet, thin edge projecting"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else { #Pipe arch/ellipitical culvert, plastic/metal, all other inlet types
        #Table A.1 - Elliptical face, tapered inlet, square edges
        c <- 0.0478
        Y <- 0.80
        model_data$Culvert_notes[i] <- "Table A.1 - Elliptical face, tapered inlet, square edges"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
      }
      
    } else { #Pipe arch/ellipitical culvert, all other materials besides plastic/metal
      #Table A.1 - Elliptical face, tapered inlet, square edges
      c <- 0.0478
      Y <- 0.80
      model_data$Culvert_notes[i] <- "Table A.1 - Elliptical face, tapered inlet, square edges"
      model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
    }
    
  } else if (model_data$Inlet_Structure_Type[i] == "Round Culvert"){
    
    if (model_data$Material[i] == "Concrete" | model_data$Material[i] == "Stone"){
      
      if (model_data$Inlet_Type[i] == "Projecting"){
        #Table A.1 - circular concrete, groove end projecting
        c <- 0.0317
        Y <- 0.69
        model_data$Culvert_notes[i] <- "Table A.1 - circular concrete, groove end projecting"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else if (model_data$Inlet_Type[i] == "Headwall" | model_data$Inlet_Type[i] == "Headwall and Wingwalls"){
        #Table A.1 - circular concrete, groove end with headwall
        c <- 0.0292
        Y <- 0.74
        model_data$Culvert_notes[i] <- "Table A.1 - circular concrete, groove end with headwall"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else { 
        c <- NA
        Y <- NA
        model_data$Culvert_notes[i] <- "No c, Y constants currently available for this culvert material/shape combination"
        model_data$Qi_m3_s[i] <- NA
      }
      
    } else if (model_data$Material[i] == "Plastic"| model_data$Material[i] == "Metal"){
      if (model_data$Inlet_Type[i] == "Projecting"){
        #Table A.1 - circular, corrugated metal, projecting
        c <- 0.0553
        Y <- 0.54
        model_data$Culvert_notes[i] <- "Table A.1 - circular, corrugated metal, projecting"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else if (model_data$Inlet_Type[i] == "Mitered to Slope"){
        #Table A.1 - circular, corrugated metal, mitered to slope
        c <- 0.0463
        Y <- 0.75
        model_data$Culvert_notes[i] <- "Table A.1 - circular, corrugated metal, mitered to slope"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else if (model_data$Inlet_Type[i] == "Headwall" | model_data$Inlet_Type[i] == "Headwall and Wingwalls"){
        #Table A.1 - circular, corrugated metal, headwall
        c <- 0.0379
        Y <- 0.69
        model_data$Culvert_notes[i] <- "Table A.1 - circular, corrugated metal, headwall"
        model_data$Qi_m3_s[i] <- model_data$inlet_area_sqm[i]*sqrt(depth*(((HW/depth)-Y-ks*culvert_slope)/c))/ku
        
      } else { 
        c <- NA
        Y <- NA
        model_data$Culvert_notes[i] <- "No c, Y constants currently available for this culvert material/shape combination"
        model_data$Qi_m3_s[i] <- NA
      }
    }
    
  } else { 
    c <- NA
    Y <- NA
    model_data$Culvert_notes[i] <- "No c, Y constants currently available for this culvert material/shape combination"
    model_data$Qi_m3_s[i] <- NA
  }

  #Write watershed area and culvert capacity to csv before starting on the next culvert
  write.csv(model_data,row.names=FALSE,file=paste(dir,"Data\\",ws,"_model_output.csv",sep=""))
}

#### Determine maximum storm size that can be safely passed ####

# Sum the capacity of culverts located at the same site to determine Q total, m^3/s
for (j in 1:nrow(model_data)){
  model_data$Qtot_m3_s[j] <- sum(model_data$Qi_m3_s[model_data$Survey_Id==model_data$Survey_Id[j]])
}
  
# Detemine current maximum storm size

current <- c("1-yr, current","2-yr, current","5-yr, current","10-yr, current", "25-yr, current",
             "50-yr, current","100-yr, current","200-yr, current","500-yr, current")
current_nos <- c(1,2,5,10,25,50,100,200,500)
future <- c("1-yr, future","2-yr, future","5-yr, future","10-yr, future", "25-yr, future",
            "50-yr, future","100-yr, future","200-yr, future","500-yr, future")
future_nos <- c(1,2,5,10,25,50,100,200,500)

for (k in 1:nrow(model_data)){
  if(is.na(model_data$Qtot_m3_s[k])==TRUE || is.na(model_data$`1-yr, current`[k])==TRUE){
    model_data$Current_max_storm_size[k] <- NA
    model_data$Future_max_storm_size[k] <- NA
  } else{
    model_data$Current_max_storm_size[k] <- c(current_nos[(which.last(model_data[k,current]<model_data$Qtot_m3_s[k]))],0)[1]
    model_data$Future_max_storm_size[k] <- c(future_nos[(which.last(model_data[k,future]<model_data$Qtot_m3_s[k]))],0)[1]
  }
}

write.csv(model_data,row.names=FALSE,file=paste(dir,"Data\\",ws,"_model_output.csv",sep=""))

