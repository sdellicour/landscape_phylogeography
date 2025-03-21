# 1. RRW simulations conducted on an environmental raster impacting the dispersal velocity of lineages
# 2. Simulations conducted along a MCC tree whose branch durations are impacted by environmental distances
# 3. Conducting all the MDS (multi-dimensional scaling) and cartogram transformations for the simulations
# 4. Preparing the XML files for the different continuous phylogeographic analyses to conduct with BEAST
# 5. Performing all the different landscape phylogeographic analyses (post hoc and prior-informed analyses)

library(cartogram)
library(diagram)
library(lubridate)
library(seraphim)
library(sf)

source("~/Dropbox/Temp_MBP3/MDS/All_new_R_functions/cleanResultsTable.r")
source("~/Dropbox/Temp_MBP3/MDS/All_new_R_functions/generateXMLfiles1.r")
source("~/Dropbox/Temp_MBP3/MDS/All_new_R_functions/generateXMLfiles2.r")

# 1. RRW simulations conducted on an environmental raster impacting the dispersal velocity of lineages

simulationDirectory = "All_RRW_simulations"; scalingValue = 2; birthRate = 0.20; samplingRate = 0.20; fourCells = TRUE

nberOfSimulations = 500; envVariable = raster("Elevation_16_k10.asc")
dir.create(file.path(paste0(simulationDirectory)), showWarnings=F)
dir.create(file.path(paste0(simulationDirectory,"/All_the_simulations")), showWarnings=F)
dir.create(file.path(paste0(simulationDirectory,"/Selected_simulations")), showWarnings=F)
mcc_tab = read.csv("RABV_US1_MCC.csv"); ancestID = which(!mcc_tab[,"node1"]%in%mcc_tab[,"node2"])[1]
ancestPosition = c(mcc_tab[ancestID,"startLon"], mcc_tab[ancestID,"startLat"]); resistance = TRUE
startingYear = 1972; samplingWindow = cbind(1982.2,2004.7)
timeSlice = 1/12; timeIntervale = 1/12; showingPlots = FALSE
extractionOfValuesOnMatrix = FALSE
for (i in 1:nberOfSimulations)
	{
		simulation = NULL; worked = FALSE
		while (worked == FALSE)
			{
				trycatch = tryCatch( # simulatorRRW3 performs envVariable[!is.na(envVariable[])] = envVariable[!is.na(envVariable[])]+1
					{
						simulation = simulatorRRW3(envVariable, resistance, scalingValue, ancestPosition, birthRate, samplingRate, startingYear, 
												   samplingWindow, timeSlice, timeIntervale, showingPlots, extractionOfValuesOnMatrix)
					},	error = function(cond) {
					},	finally = {
					})
				if (!is.null(simulation)) worked = TRUE
			}
		write.csv(simulation[[1]], paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".csv"), quote=F, row.names=F)
		write.tree(simulation[[2]], paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".tre"))			
	}
for (i in 1:nberOfSimulations)
	{
		pdf(paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".pdf"), width=8.0, height=4.0) # dev.new(width=10, height=4.0)
		par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(0.3,0,0,0.3), mar=c(1.5,0.5,0.5,0), lwd=0.4, col="gray30")
		tab = read.csv(paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".csv"), head=T)
		for (j in 1:dim(tab)[1])
			{
				if (tab[j,"node1"]%in%tab[,"node2"])
					{
						index = which(tab[,"node2"]==tab[j,"node1"])
						if (tab[index,"endLon"] != tab[j,"startLon"]) print(c(i,j))
					}
			}
		tree = read.tree(paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".tre"))
		col_start = colorRampPalette(brewer.pal(9,"PuBu"))(101)[1]
		cols3 = colorRampPalette(brewer.pal(9,"PuBu"))(101)[(((nodeHeights(tree)[,2])/(samplingWindow[2]-startingYear))*100)+1]
		plot(tree, show.tip.label=F, edge.width=0.4, edge.col="gray30")
		for (j in 1:dim(tree$edge)[1])
			{
				if (j == 1)
					{
						nodelabels(node=tree$edge[j,1], pch=16, cex=0.7, col=col_start)
						nodelabels(node=tree$edge[j,1], pch=1, cex=0.7, col="gray30", lwd=0.25)
					}
				nodelabels(node=tree$edge[j,2], pch=16, cex=0.7, col=cols3[j])
				nodelabels(node=tree$edge[j,2], pch=1, cex=0.7, col="gray30", lwd=0.25)
			}
		axis(side=1, lwd.tick=0, cex.axis=0.65, mgp=c(0,0.1,0), lwd=0.4, tck=0, col.tick="gray30", col.axis="gray30", col="gray30", at=c(0,335), labels=rep(NA,2))
		axis(side=1, lwd.tick=0.4, cex.axis=0.65, mgp=c(0,0.1,0), lwd=0, tck=-0.016, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(-2,38,5), labels=seq(1970,2010,5))
		plot(extent(envVariable), axes=F, ann=F, col="white")
		plot(envVariable, add=T, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(101)[1:101], colNA="gray90", legend=F, useRaster=T, asp=1, axes=F, box=F)
		ext = extent(envVariable); rect(ext@xmin, ext@ymin, ext@xmax, ext@ymax, lwd=0.4, border="gray30", col=NA)
		startingYear = min(tab[,"startYear"]); col_start = colorRampPalette(brewer.pal(9,"PuBu"))(101)[1]
		cols2 = colorRampPalette(brewer.pal(9,"PuBu"))(101)[(((tab[,"endYear"]-startingYear)/(samplingWindow[2]-startingYear))*100)+1]
		for (j in 1:dim(tab)[1])
			{	
				curvedarrow(cbind(tab[j,"startLon"],tab[j,"startLat"]), cbind(tab[j,"endLon"],tab[j,"endLat"]), arr.length=0,
							arr.width=0, lwd=0.3, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
			}
		for (j in dim(tab)[1]:1)
			{
				if (j == 1)
					{
						points(cbind(tab[j,"startLon"],tab[j,"startLat"]), pch=16, cex=0.7, col=col_start)
						points(cbind(tab[j,"startLon"],tab[j,"startLat"]), pch=1, cex=0.7, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
					}
				points(cbind(tab[j,"endLon"],tab[j,"endLat"]), pch=16, cex=0.7, col=cols2[j])
				points(cbind(tab[j,"endLon"],tab[j,"endLat"]), pch=1, cex=0.7, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
			}
		plot(envVariable, legend.only=T, add=T, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(101)[1:101], legend.width=0.5, legend.shrink=0.3,
			 smallplot=c(0.060,0.965,0.079,0.094), legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
			 axis.args=list(cex.axis=0.65, lwd=0, lwd.tick=0.4, tck=-1.0, col.tick="gray30", col.axis="gray30", line=0, mgp=c(0,0.1,0), at=seq(1,10,1)))
		dev.off()
	}
c = 0
for (i in 1:nberOfSimulations)
	{
		tab = read.csv(paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".csv"), head=T); nberOfTips = length(which(!tab[,"node2"]%in%tab[,"node1"])); # print(nberOfTips)
		mean_env = mean(raster::extract(envVariable, tab[which(!tab[,"node2"]%in%tab[,"node1"]),c("endLon","endLat")]))
		sd_env = sd(raster::extract(envVariable, tab[which(!tab[,"node2"]%in%tab[,"node1"]),c("endLon","endLat")]))
		if ((nberOfTips >= 50)&((sd_env/mean_env) > 0.5)&(c < 100))
			{
				c = c+1; print(c(i,c,nberOfTips,sd_env/mean_env))
				file.copy(paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".csv"), paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",c,".csv"), overwrite=T)
				file.copy(paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".tre"), paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",c,".tre"), overwrite=T)
				file.copy(paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".pdf"), paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",c,".pdf"), overwrite=T)
			}
	}
i = 22; timeSlots = c(1985, 1995, 2005) # to plot successice snapshots of a selected RRW simulation
for (h in 1:length(timeSlots))
	{
		pdf(paste0("RRW_simulation_",i,"_",timeSlots[h],".pdf"), width=8.0, height=4.0)
		par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(0.3,0,0,0.3), mar=c(1.5,0.5,0.5,0), lwd=0.4, col="gray30")
		tab = read.csv(paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",i,".csv"), head=T)
		for (j in 1:dim(tab)[1])
			{
				if ((tab[j,"endYear"] <= timeSlots[h])&(tab[j,"node1"]%in%tab[,"node2"]))
					{
						index = which(tab[,"node2"]==tab[j,"node1"])
						if (tab[index,"endLon"] != tab[j,"startLon"]) print(c(i,j))
					}
			}
		tree = read.tree(paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",i,".tre"))
		col_start = colorRampPalette(brewer.pal(9,"PuBu"))(101)[1]
		cols3 = colorRampPalette(brewer.pal(9,"PuBu"))(101)[(((nodeHeights(tree)[,2])/(samplingWindow[2]-startingYear))*100)+1]
		plot(tree, show.tip.label=F, edge.width=0.4, edge.col="gray30")
		for (j in 1:dim(tree$edge)[1])
			{
				if (j == 1)
					{
						nodelabels(node=tree$edge[j,1], pch=16, cex=0.7, col=col_start)
						nodelabels(node=tree$edge[j,1], pch=1, cex=0.7, col="gray30", lwd=0.25)
					}
				nodelabels(node=tree$edge[j,2], pch=16, cex=0.7, col=cols3[j])
				nodelabels(node=tree$edge[j,2], pch=1, cex=0.7, col="gray30", lwd=0.25)
			}
		axis(side=1, lwd.tick=0, cex.axis=0.65, mgp=c(0,0.1,0), lwd=0.4, tck=0, col.tick="gray30", col.axis="gray30", col="gray30", at=c(0,335), labels=rep(NA,2))
		axis(side=1, lwd.tick=0.4, cex.axis=0.65, mgp=c(0,0.1,0), lwd=0, tck=-0.016, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(-2,38,5), labels=seq(1970,2010,5))
		plot(extent(envVariable), axes=F, ann=F, col="white")
		plot(envVariable, add=T, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(101)[1:101], colNA="gray90", legend=F, useRaster=T, asp=1, axes=F, box=F)
		ext = extent(envVariable); rect(ext@xmin, ext@ymin, ext@xmax, ext@ymax, lwd=0.4, border="gray30", col=NA)
		startingYear = min(tab[,"startYear"]); col_start = colorRampPalette(brewer.pal(9,"PuBu"))(101)[1]
		cols2 = colorRampPalette(brewer.pal(9,"PuBu"))(101)[(((tab[,"endYear"]-startingYear)/(samplingWindow[2]-startingYear))*100)+1]
		for (j in 1:dim(tab)[1])
			{
				if (tab[j,"endYear"] <= timeSlots[h])
					{
						curvedarrow(cbind(tab[j,"startLon"],tab[j,"startLat"]), cbind(tab[j,"endLon"],tab[j,"endLat"]), arr.length=0,
									arr.width=0, lwd=0.3, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
					}
			}
		for (j in 1:dim(tab)[1])
			{
				if (tab[j,"endYear"] <= timeSlots[h])
					{
						if (j == 1)
							{
								points(cbind(tab[j,"startLon"],tab[j,"startLat"]), pch=16, cex=0.7, col=col_start)
								points(cbind(tab[j,"startLon"],tab[j,"startLat"]), pch=1, cex=0.7, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
							}
						points(cbind(tab[j,"endLon"],tab[j,"endLat"]), pch=16, cex=0.7, col=cols2[j])
						points(cbind(tab[j,"endLon"],tab[j,"endLat"]), pch=1, cex=0.7, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
					}
			}
		plot(envVariable, legend.only=T, add=T, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(101)[1:101], legend.width=0.5, legend.shrink=0.3,
			 smallplot=c(0.060,0.965,0.079,0.094), legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
			 axis.args=list(cex.axis=0.65, lwd=0, lwd.tick=0.4, tck=-1.0, col.tick="gray30", col.axis="gray30", line=0, mgp=c(0,0.1,0), at=seq(1,10,1)))
		dev.off()
	}
firstExplorations = TRUE
if (firstExplorations)
	{
		wd = getwd(); setwd(simulationDirectory)
		localTreesDirectory = paste0("Selected_simulations"); envVariables = list(envVariable)
		envVariables[[1]][] = envVariables[[1]][]+1; resistances = list(resistance)
		avgResistances = list(resistance); nberOfRandomisations = 0; randomProcedure = 3
		minimumConvexHull = TRUE; showingPlots = FALSE; nberOfCores = 1; OS = "Unix"
		juliaCSImplementation = FALSE; simulations = TRUE; randomisations = FALSE
		dir.create(paste0("Seraphim_analyses"), showWarnings=F)
		dir.create(paste0("Seraphim_analyses"), showWarnings=F)
		for (i in 1:100)
			{
				nberOfExtractionFiles = 1
				tab = read.csv(paste0("Selected_simulations/TreeSimulations_",i,".csv"), head=T)
				dir.create(paste0("Selected_simulations/TreeSimulations_",i,"_ext"), showWarnings=F)
				write.csv(tab, paste0("Selected_simulations/TreeSimulations_",i,"_ext/TreeSimulations_1.csv"), quote=F, row.names=F)
				localTreesDirectory = paste0("Selected_simulations/TreeSimulations_",i,"_ext")
				pathModel = 2; outputName = paste0("TreeSimulations_",i,"_LC")
				spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
							  randomProcedure, outputName, showingPlots, nberOfCores, OS, juliaCSImplementation, simulations, randomisations, minimumConvexHull)
 				file.copy(paste0("TreeSimulations_",i,"_LC_env_distances.txt"), paste0("Seraphim_analyses/TreeSimulations_",i,"_LC_env_distances.txt"), overwrite=T)
 				file.remove(paste0("TreeSimulations_",i,"_LC_env_distances.txt"))
				pathModel = 3; outputName = paste0("TreeSimulations_",i,"_CS")
				spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
							  randomProcedure, outputName, showingPlots, nberOfCores, OS, juliaCSImplementation, simulations, randomisations, minimumConvexHull)
 				file.copy(paste0("TreeSimulations_",i,"_CS_env_distances.txt"), paste0("Seraphim_analyses/TreeSimulations_",i,"_CS_env_distances.txt"), overwrite=T)
 				file.remove(paste0("TreeSimulations_",i,"_CS_env_distances.txt"))
			}
		unlink("CS_rasters", recursive=T)
		qValues = matrix(nrow=100, ncol=4); colnames(qValues) = c("LC_1","CS_1","LC_2","CS_2")
		for (i in 1:100)
			{
				tab = read.table(paste("Seraphim_analyses/TreeSimulations_",i,"_LC_env_distances.txt",sep=""), header=T)
				durations = tab[,"dispersal_times"]; distances = tab[,2]; durations4 = 4*durations; squareDistances = distances^2
				LM1 = lm(as.formula(paste("durations ~ distances", sep=""))); R2_null_1 = summary(LM1)$r.squared
				LM2 = lm(as.formula(paste("durations4 ~ squareDistances", sep=""))); R2_null_2 = summary(LM2)$r.squared
				durations = tab[,"dispersal_times"]; distances = tab[,3]; durations4 = 4*durations; squareDistances = distances^2				
				LM1 = lm(as.formula(paste("durations ~ distances", sep=""))); R2_env_1 = summary(LM1)$r.squared
				LM2 = lm(as.formula(paste("durations4 ~ squareDistances", sep=""))); R2_env_2 = summary(LM2)$r.squared
				qValues[i,"LC_1"] = R2_env_1 - R2_null_1; qValues[i,"LC_2"] = R2_env_2 - R2_null_2				
				tab = read.table(paste("Seraphim_analyses/TreeSimulations_",i,"_CS_env_distances.txt",sep=""), header=T)
				durations = tab[,"dispersal_times"]; distances = tab[,2]; durations4 = 4*durations; squareDistances = distances^2
				LM1 = lm(as.formula(paste("durations ~ distances", sep=""))); R2_null_1 = summary(LM1)$r.squared
				LM2 = lm(as.formula(paste("durations4 ~ squareDistances", sep=""))); R2_null_2 = summary(LM2)$r.squared
				durations = tab[,"dispersal_times"]; distances = tab[,3]; durations4 = 4*durations; squareDistances = distances^2				
				LM1 = lm(as.formula(paste("durations ~ distances", sep=""))); R2_env_1 = summary(LM1)$r.squared
				LM2 = lm(as.formula(paste("durations4 ~ squareDistances", sep=""))); R2_env_2 = summary(LM2)$r.squared
				qValues[i,"CS_1"] = R2_env_1 - R2_null_1; qValues[i,"CS_2"] = R2_env_2 - R2_null_2
			}
		cat(gsub("All_RRW_simulations/","",simulationDirectory)," - LC: mean Q1 = ",round(mean(qValues[,"LC_1"]),3),", p(Q1>0) = ",sum(qValues[,"LC_1"]>0)/100,sep="")
		cat("; CS: mean Q1 = ",round(mean(qValues[,"CS_1"]),3),", p(Q1>0) = ",sum(qValues[,"CS_1"]>0)/100,"\n",sep="")
		cat(gsub("All_RRW_simulations/","",simulationDirectory)," - LC: mean Q2 = ",round(mean(qValues[,"LC_2"]),3),", p(Q2>0) = ",sum(qValues[,"LC_2"]>0)/100,sep="")
		cat("; CS: mean Q2 = ",round(mean(qValues[,"CS_2"]),3),", p(Q2>0) = ",sum(qValues[,"CS_2"]>0)/100,"\n",sep="")
		# RABV_US1_4cells_T - LC: mean Q1 = 0.148, p(Q1>0) = 0.99; CS: mean Q1 = 0.048, p(Q1>0) = 0.78
		# RABV_US1_4cells_F - LC: mean Q1 = 0.134, p(Q1>0) = 0.98; CS: mean Q1 = 0.025, p(Q1>0) = 0.65
		# RABV_US1_4cells_T - LC: mean Q2 = 0.184, p(Q2>0) = 0.99; CS: mean Q2 = 0.070, p(Q2>0) = 0.82
		# RABV_US1_4cells_F - LC: mean Q2 = 0.158, p(Q2>0) = 0.96; CS: mean Q2 = 0.041, p(Q2>0) = 0.78
		setwd(wd)
	}

# 2. Simulations conducted along a MCC tree whose branch durations are impacted by environmental distances

simulationDirectory0 = "All_MCC_simulations"; nberOfSimulations = 10; tab = read.csv("RABV_US1_MCC.csv", head=T)

envVariable = raster("Elevation_16_k10.asc"); resistance = TRUE; avgResistance = resistance; fourCells = TRUE; OS="Unix"
tab1 = tab; tab2 = tab; tab3 = tab; tab4 = tab; tab5 = tab; tab6 = tab; tabs = list(); extension = ".asc"; nberOfCores_CS = 1
simulationDirectories = c("RABV_US1_LC_s_1","RABV_US1_LC_s_2","RABV_US1_LC_s_3","RABV_US1_CS_s_1","RABV_US1_CS_s_2","RABV_US1_CS_s_3")
for (i in 1:length(simulationDirectories)) simulationDirectories[i] = paste0(simulationDirectory0,"/",simulationDirectories[i])

	# 2.1. Modification of the MCC branch lengths according to different scenarios and path models

if (file.exists("All_MCC_simulations/RABV_US1_original/")) file.copy("RABV_US1_MCC.csv", "All_MCC_simulations/RABV_US1_original/RABV_US1_simu.csv")
points = tab[,c("endLon","endLat")]; hull = chull(points); hull = c(hull,hull[1]); p = Polygon(points[hull,])
ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps)); hullRaster = envVariable; hullRaster[] = hullRaster[]+1 # hull raster not used
hullRaster = raster::mask(crop(hullRaster, sps, snap="out"), sps, snap="out"); nullRaster = hullRaster; nullRaster[!is.na(nullRaster[])] = 1
envVariable_mod = envVariable; envVariable_mod[] = envVariable_mod[]+1; nullRaster = envVariable_mod; nullRaster[!is.na(nullRaster[])] = 1
trRast = transition(nullRaster, function(x) 1/mean(x), directions=4); trRastCorr = geoCorrection(trRast, type="c", multpl=F, scl=T)	
d = diag(costDistance(trRastCorr, as.matrix(tab[,c("startLon","startLat")]), as.matrix(tab[,c("endLon","endLat")]))); tab1[,"length"] = d
if (grepl("RABV_US1_on_DCs",simulationDirectory0)) tab1[,"length"] = (d^2)/4
trRast = transition(envVariable_mod, function(x) 1/mean(x), directions=4); trRastCorr = geoCorrection(trRast, type="c", multpl=F, scl=T)	
d = diag(costDistance(trRastCorr, as.matrix(tab[,c("startLon","startLat")]), as.matrix(tab[,c("endLon","endLat")]))); tab3[,"length"] = d
if (grepl("RABV_US1_on_DCs",simulationDirectory0)) tab3[,"length"] = (d^2)/4
outputName = "scenario_1"; envVariableName = paste("CS_rasters/NullRaster_",outputName,"_cs",extension,sep="")
dir.create(file.path(getwd(), "CS_rasters"), showWarnings=F); writeRaster(nullRaster, envVariableName, overwrite=T)
fromCoor = as.matrix(tab[,c("startLon","startLat")]); toCoor = as.matrix(tab[,c("endLon","endLat")]); prefix = outputName; ID = 1
d = circuitScape1(nullRaster, envVariableName, resistance, avgResistance, fourCells, fromCoor, toCoor, OS, prefix, ID, nberOfCores_CS); tab4[,"length"] = d
if (grepl("RABV_US1_on_DCs",simulationDirectory0)) tab4[,"length"] = (d^2)/4
outputName = "scenario_3"; envVariableName = paste("CS_rasters/EnvVariable_",outputName,"_cs",extension,sep="")
dir.create(file.path(getwd(), "CS_rasters"), showWarnings=F); writeRaster(envVariable_mod, envVariableName, overwrite=T)
fromCoor = as.matrix(tab[,c("startLon","startLat")]); toCoor = as.matrix(tab[,c("endLon","endLat")]); prefix = outputName; ID = 1
d = circuitScape1(envVariable_mod, envVariableName, resistance, avgResistance, fourCells, fromCoor, toCoor, OS, prefix, ID, nberOfCores_CS); tab6[,"length"] = d
if (grepl("RABV_US1_on_DCs",simulationDirectory0)) tab6[,"length"] = (d^2)/4
tabs[[1]] = tab1; tabs[[2]] = tab2; tabs[[3]] = tab3; tabs[[4]] = tab4; tabs[[5]] = tab5; tabs[[6]] = tab6
for (i in c(1,3,4,6))
	{
		tabs[[i]] = tabs[[i]][c("node1","node2","length","startLon","startLat","endLon","endLat","startYear","endYear")]
		newStartYear = rep(NA, dim(tabs[[i]])[1]); newEndYear = rep(NA, dim(tabs[[i]])[1])
		for (j in 1:dim(tabs[[i]])[1])
			{
				node2 = tabs[[i]][j,"node1"]; d2 = tabs[[i]][j,"length"]; d = d2; ancestralNode1 = FALSE
				while (ancestralNode1 == FALSE)
					{
						node1 = c(); ancestralNode2 = FALSE
						for	(k in 1:dim(tabs[[i]])[1])
							{
								if (node2 == tabs[[i]][k,"node2"])
									{
										d2 = d2 + tabs[[i]][k,"length"]
										node1 = tabs[[i]][k,"node1"]
										node2 = node1; ancestralNode2 = TRUE
									}
							}
						if (ancestralNode2 == FALSE) ancestralNode1 = TRUE
					}
				d1 = d2-d; newStartYear[j] = d1; newEndYear[j] = d2
			}
		Min0 = min(tab[,"startYear"]); Max0 = max(tab[,"endYear"]); L0 = Max0-Min0
		Min = min(newStartYear); Max = max(newEndYear); L = Max-Min
		for (j in 1:dim(tabs[[i]])[1])
			{
				 l = newEndYear[j]-newStartYear[j]
				 tabs[[i]][j,"length"] = (l/L)*L0
				 newStartYear[j] = ((newStartYear[j]/L)*L0) + Min0
				 newEndYear[j] = ((newEndYear[j]/L)*L0) + Min0
			}
		tabs[[i]][,"startYear"] = newStartYear; tabs[[i]][,"endYear"] = newEndYear
		dir.create(simulationDirectories[i], showWarnings=F)
		write.csv(tabs[[i]], paste0(simulationDirectories[i],"/RABV_US1_simu.csv"), row.names=F, quote=F)
	}
for (i in 1:dim(tabs[[2]])[1]) tabs[[2]][i,"length"] = mean(c(tabs[[1]][i,"length"],tabs[[2]][i,"length"]))
for (i in 1:dim(tabs[[5]])[1]) tabs[[5]][i,"length"] = mean(c(tabs[[4]][i,"length"],tabs[[6]][i,"length"]))
write.csv(tabs[[2]], paste0(simulationDirectories[2],"/RABV_US1_simu.csv"), row.names=F, quote=F)
write.csv(tabs[[5]], paste0(simulationDirectories[5],"/RABV_US1_simu.csv"), row.names=F, quote=F)

for (i in 1:length(simulationDirectories)) tabs[[i]] = read.csv(paste0(simulationDirectories[i],"/RABV_US1_simu.csv"))
for (i in 1:length(tabs))
	{
		dir.create(paste0(simulationDirectories[i],"/Selected_simulations"), showWarnings=F)
		dir.create(paste0(simulationDirectories[i],"/All_BEAST_analyses"), showWarnings=F)
		dir.create(paste0(simulationDirectories[i],"/Seraphim_analyses"), showWarnings=F)
		nberOfExtractionFiles = 1; localTreesDirectory = paste0(simulationDirectories[i],"/Selected_simulations")
		write.csv(tabs[[i]], paste0(localTreesDirectory,"/TreeSimulations_1.csv"), quote=F, row.names=F)
		envVariables = list(envVariable); envVariables[[1]][] = envVariables[[1]][]+1; resistances = c(TRUE); avgResistances = c(TRUE)
		fourCells = TRUE; nberOfRandomisations = 0; randomProcedure = 3; minimumConvexHull = TRUE; showingPlots = FALSE; nberOfCores = 1; unix = "OS"
		pathModel = 2; outputName = paste0("TreeSimulations_1_LC"); simulations = TRUE; randomisations = FALSE; juliaCSImplementation = FALSE
		spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
					  randomProcedure, outputName, showingPlots, nberOfCores, OS, juliaCSImplementation, simulations, randomisations, minimumConvexHull)
		file.copy(paste0("TreeSimulations_1_LC_env_distances.txt"), paste0(simulationDirectories[i],"/Seraphim_analyses/TreeSimulations_1_LC_env_distances.txt"), overwrite=T)
		pathModel = 3; outputName = paste0("TreeSimulations_1_CS"); simulations = TRUE; randomisations = FALSE; juliaCSImplementation = FALSE
		spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
					  randomProcedure, outputName, showingPlots, nberOfCores, OS, juliaCSImplementation, simulations, randomisations, minimumConvexHull)
		file.copy(paste0("TreeSimulations_1_CS_env_distances.txt"), paste0(simulationDirectories[i],"/Seraphim_analyses/TreeSimulations_1_CS_env_distances.txt"), overwrite=T)
		file.remove(paste0("TreeSimulations_1_LC_env_distances.txt")); file.remove(paste0("TreeSimulations_1_CS_env_distances.txt"))
	}
for (i in 1:length(tabs))
	{
		tab = read.table(paste0(simulationDirectories[i],"/Seraphim_analyses/TreeSimulations_1_LC_env_distances.txt"), header=T)
		durations = tab[,"dispersal_times"]; distances = tab[,2]; durations4 = 4*durations; squareDistances = distances^2
		LM1 = lm(as.formula(paste("durations ~ distances", sep=""))); R2_null_1 = summary(LM1)$r.squared
		LM2 = lm(as.formula(paste("durations4 ~ squareDistances", sep=""))); R2_null_2 = summary(LM2)$r.squared
		durations = tab[,"dispersal_times"]; distances = tab[,3]; durations4 = 4*durations; squareDistances = distances^2				
		LM1 = lm(as.formula(paste("durations ~ distances", sep=""))); R2_env_1 = summary(LM1)$r.squared
		LM2 = lm(as.formula(paste("durations4 ~ squareDistances", sep=""))); R2_env_2 = summary(LM2)$r.squared
		Q_LC_1 = round(R2_env_1-R2_null_1, 3); Q_LC_2 = round(R2_env_2-R2_null_2, 3)
		tab = read.table(paste0(simulationDirectories[i],"/Seraphim_analyses/TreeSimulations_1_CS_env_distances.txt"), header=T)
		durations = tab[,"dispersal_times"]; distances = tab[,2]; durations4 = 4*durations; squareDistances = distances^2
		LM1 = lm(as.formula(paste("durations ~ distances", sep=""))); R2_null_1 = summary(LM1)$r.squared
		LM2 = lm(as.formula(paste("durations4 ~ squareDistances", sep=""))); R2_null_2 = summary(LM2)$r.squared
		durations = tab[,"dispersal_times"]; distances = tab[,3]; durations4 = 4*durations; squareDistances = distances^2				
		LM1 = lm(as.formula(paste("durations ~ distances", sep=""))); R2_env_1 = summary(LM1)$r.squared
		LM2 = lm(as.formula(paste("durations4 ~ squareDistances", sep=""))); R2_env_2 = summary(LM2)$r.squared
		Q_CS_1 = round(R2_env_1-R2_null_1, 3); Q_CS_2 = round(R2_env_2-R2_null_2, 3)
		cat("Q_LC_1 = ",Q_LC_1,", Q_LC_2 = ",Q_LC_2,", Q_CS_1 = ",Q_CS_1,", Q_CS_2 = ",Q_CS_2,"\n",sep="")
		# Based on LDV, LC, scenario 1: Q_LC_1 = -0.289, Q_LC_2 = -0.074, Q_CS_1 = -0.157, Q_CS_2 = -0.129
		# Based on LDV, LC, scenario 2: Q_LC_1 =  0.065, Q_LC_2 =  0.035, Q_CS_1 =  0.080, Q_CS_2 =  0.031
		# Based on LDV, LC, scenario 3: Q_LC_1 =  0.329, Q_LC_2 =  0.517, Q_CS_1 =  0.322, Q_CS_2 =  0.335
		# Based on LDV, CS, scenario 1: Q_LC_1 = -0.165, Q_LC_2 =  0.019, Q_CS_1 = -0.356, Q_CS_2 = -0.255
		# Based on LDV, CS, scenario 2: Q_LC_1 =  0.065, Q_LC_2 =  0.035, Q_CS_1 =  0.080, Q_CS_2 =  0.031
		# Based on LDV, CS, scenario 3: Q_LC_1 =  0.333, Q_LC_2 =  0.223, Q_CS_1 =  0.591, Q_CS_2 =  0.362
		# Based on DCs, LC, scenario 1: Q_LC_1 = -0.403, Q_LC_2 = -0.459, Q_CS_1 = -0.115, Q_CS_2 = -0.181
		# Based on DCs, LC, scenario 2: Q_LC_1 =  0.065, Q_LC_2 =  0.035, Q_CS_1 =  0.080, Q_CS_2 =  0.031
		# Based on DCs, LC, scenario 3: Q_LC_1 =  0.237, Q_LC_2 =  0.553, Q_CS_1 =  0.172, Q_CS_2 =  0.231
		# Based on DCs, CS, scenario 1: Q_LC_1 = -0.286, Q_LC_2 = -0.112, Q_CS_1 = -0.399, Q_CS_2 = -0.392
		# Based on DCs, CS, scenario 2: Q_LC_1 =  0.065, Q_LC_2 =  0.035, Q_CS_1 =  0.080, Q_CS_2 =  0.031
		# Based on DCs, CS, scenario 3: Q_LC_1 =  0.386, Q_LC_2 =  0.315, Q_CS_1 =  0.573, Q_CS_2 =  0.412
	}	# Note: differences among environmental distances are expected due to the use of a convex hull mask

	# 2.2. Generation of Nexus tree files based on the modified extraction files of the MCC tree

for (i in 1:length(simulationDirectories))
	{
		tab = read.csv(paste0(simulationDirectories[i],"/RABV_US1_simu.csv"), head=T)
		branches1 = c(); nodes1 = c(); tipNodes = c()
		for (j in 1:dim(tab)[1])
			{
				tipNode = TRUE
				for (k in 1:dim(tab)[1])
					{
						if (tab[j,"node2"] == tab[k,"node1"]) tipNode = FALSE
					}
				if (tipNode == TRUE)
					{
						branches1 = c(branches1, paste(tab[j,"node2"],":",tab[j,"length"],sep=""))
						nodes1 =  c(nodes1, tab[j,"node2"])
						tipNodes = rbind(tipNodes, c(tab[j,"node2"], tab[j,"newStartYear"]))
					}
			}
		coalescenceEvents = length(branches1)-1
		while (length(branches1) != 1)
			{
				branches2 = branches1; nodes2 = nodes1
				nodesToRemove = c()
				for (j in 2:length(branches1))
					{
						nodeA = nodes1[j]
						for (k in 1:(j-1))
							{
								nodeB = nodes1[k]
								if (tab[(tab[,"node2"]==nodeA),"node1"] == tab[(tab[,"node2"]==nodeB),"node1"])
									{
										coalescenceEvents = coalescenceEvents-1
										nodesToRemove = c(nodesToRemove, k, j)
										node2 = tab[(tab[,"node2"]==tab[(tab[,"node2"]==nodeA),"node1"]),"node2"]
										if (coalescenceEvents > 0)
											{
												branches2 = c(branches2, paste("(",branches1[k],",",branches1[j],")",":",tab[(tab[,"node2"]==node2),"length"],sep=""))
											}	else	{
												branches2 = c(branches2, paste("(",branches1[k],",",branches1[j],")",";",tab[(tab[,"node2"]==node2),"length"],sep=""))
											}
										nodes2 = c(nodes2, node2) # print(paste(coalescenceEvents,nodeB,nodeA,node2,sep=" "))
									}
							}
					}
				branches2 = branches2[-nodesToRemove]
				nodes2 = nodes2[-nodesToRemove]
				branches1 = branches2; nodes1 = nodes2
			}
		sink(file=paste0(simulationDirectories[i],"/RABV_US1_simu_TO_BE_RE-EXPORTED_IN_NEXUS_BY_FIGTREE.tre"))
		cat("#NEXUS"); cat("\n\n")
		cat("Begin taxa;"); cat("\n")
		cat("\tDimensions ntax=47;"); cat("\n")
		cat("\tTaxlabels"); cat("\n")
		for (j in 1:length(tipNodes))
			{
				a = tab[tab[,"node2"]==j,]
				cat(paste("\t\t'US1_",round(a[1,"endYear"],6),"_",round(a[1,"endLon"],6),"_",round(a[1,"endLat"],6),"'\n",sep=""))
			}
		cat("\t\t;"); cat("\n")
		cat("End;"); cat("\n\n")
		cat("Begin trees;"); cat("\n")
		cat("\tTranslate"); cat("\n")
		for (j in 1:length(tipNodes))
			{
				a = tab[tab[,"node2"]==j,]
				cat(paste("\t\t ",j," 'US1_",round(a[1,"endYear"],6),"_",round(a[1,"endLon"],6),"_",round(a[1,"endLat"],6),"',\n",sep=""))
			}
		cat("\t\t;"); cat("\n")	
		cat(paste("tree TREE1 = [&R] ",branches1,"\n",sep=""))
		cat("End;"); cat("\n")
		sink(NULL)
	}
for (i in 1:length(simulationDirectories))
	{
		pdf(paste0(simulationDirectories[i],"/RABV_US1_simu.pdf"), width=8.0, height=4.0) # dev.new(width=10, height=4.0)
		par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(0.3,0,0,0.3), mar=c(1.5,0.5,0.5,0), lwd=0.4, col="gray30")
		tab = read.csv(paste0(simulationDirectories[i],"/RABV_US1_simu.csv"), head=T)
		tre = read.nexus(paste0(simulationDirectories[i],"/RABV_US1_simu.tre"))
		startingYear = 1972; samplingWindow = cbind(1982.2,2004.7) # should work with the rescaling step above
		startingYear = min(tab[,"startYear"]); samplingWindow = cbind(min(tab[,"startYear"]),max(tab[,"endYear"]))
		for (j in 1:dim(tab)[1])
			{
				if (tab[j,"node1"]%in%tab[,"node2"])
					{
						index = which(tab[,"node2"]==tab[j,"node1"])
						if (tab[index,"endLon"] != tab[j,"startLon"]) print(c(i,j))
					}
			}
		col_start = colorRampPalette(brewer.pal(9,"PuBu"))(101)[1]
		cols3 = colorRampPalette(brewer.pal(9,"PuBu"))(101)[(((nodeHeights(tre)[,2])/(samplingWindow[2]-startingYear))*100)+1]
		plot(tre, show.tip.label=F, edge.width=0.4, edge.col="gray30")
		for (j in 1:dim(tre$edge)[1])
			{
				if (j == 1)
					{
						nodelabels(node=tre$edge[j,1], pch=16, cex=0.7, col=col_start)
						nodelabels(node=tre$edge[j,1], pch=1, cex=0.7, col="gray30", lwd=0.25)
					}
				nodelabels(node=tre$edge[j,2], pch=16, cex=0.7, col=cols3[j])
				nodelabels(node=tre$edge[j,2], pch=1, cex=0.7, col="gray30", lwd=0.25)
			}
		axis(side=1, lwd.tick=0, cex.axis=0.65, mgp=c(0,0.1,0), lwd=0.4, tck=0, col.tick="gray30", col.axis="gray30", col="gray30", at=c(0,335), labels=rep(NA,2))
		axis(side=1, lwd.tick=0.4, cex.axis=0.65, mgp=c(0,0.1,0), lwd=0, tck=-0.016, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(-2,38,5), labels=seq(1970,2010,5))
		plot(extent(envVariable), axes=F, ann=F, col="white")
		plot(envVariable, add=T, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(101)[1:101], colNA="gray90", legend=F, useRaster=T, asp=1, axes=F, box=F)
		ext = extent(envVariable); rect(ext@xmin, ext@ymin, ext@xmax, ext@ymax, lwd=0.4, border="gray30", col=NA)
		startingYear = min(tab[,"startYear"]); col_start = colorRampPalette(brewer.pal(9,"PuBu"))(101)[1]
		cols2 = colorRampPalette(brewer.pal(9,"PuBu"))(101)[(((tab[,"endYear"]-startingYear)/(samplingWindow[2]-startingYear))*100)+1]
		for (j in 1:dim(tab)[1])
			{	
				curvedarrow(cbind(tab[j,"startLon"],tab[j,"startLat"]), cbind(tab[j,"endLon"],tab[j,"endLat"]), arr.length=0,
							arr.width=0, lwd=0.3, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
			}
		for (j in dim(tab)[1]:1)
			{
				if (j == 1)
					{
						points(cbind(tab[j,"startLon"],tab[j,"startLat"]), pch=16, cex=0.7, col=col_start)
						points(cbind(tab[j,"startLon"],tab[j,"startLat"]), pch=1, cex=0.7, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
					}
				points(cbind(tab[j,"endLon"],tab[j,"endLat"]), pch=16, cex=0.7, col=cols2[j])
				points(cbind(tab[j,"endLon"],tab[j,"endLat"]), pch=1, cex=0.7, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
			}
		plot(envVariable, legend.only=T, add=T, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(101)[1:101], legend.width=0.5, legend.shrink=0.3,
			 smallplot=c(0.060,0.965,0.079,0.094), legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
			 axis.args=list(cex.axis=0.65, lwd=0, lwd.tick=0.4, tck=-1.0, col.tick="gray30", col.axis="gray30", line=0, mgp=c(0,0.1,0), at=seq(1,10,1)))
		dev.off()
	}

	# 2.3. Simulation of genomic sequences with piBUSS based on the modified MCC trees

for (i in 1:length(simulationDirectories))
	{
		for (j in 1:nberOfSimulations)
			{
				# piBUSS simulations (require local BEAGLE compilation): 12000 bp, HKY+G (default settings), lognormal clock (ucld.mean = 0.0005, ucld.stdev = 0.00025)
				system(paste0("BEAGLEPATH=\"/usr/local/lib/\" && java -Djava.library.path=\"$BEAGLEPATH\" -jar piBUSS_version_1_4/pibuss_v1.4.jar -treeFile ",simulationDirectories[i],"/RABV_US1_simu.tre -from 1 -to 12000 -every 1 -branchSubstitutionModel HKY -HKYsubstitutionParameterValues 1.0 -siteRateModel gammaSiteRateModel -gammaSiteRateModelParameterValues 4.0 0.5 0.0 -clockRateModel lognormalRelaxedClock -lognormalRelaxedClockParameterValues 0.0005 0.00025 0 -baseFrequencies nucleotideFrequencies -nucleotideFrequencyParameterValues 0.25 0.25 0.25 0.25 : ",simulationDirectories[i],"/All_BEAST_analyses/",gsub(paste0(simulationDirectory0,"/"),"",gsub("CS_s","CS_scenario",gsub("LC_s","LC_scenario",gsub("original","original_tree",simulationDirectories[i])))),"_simulation_",j,".fas"))
			}
	}

# 3. Conducting all the MDS (multi-dimensional scaling) and cartogram transformations for the simulations

simulationDirectories = list(); envVariable = raster("Elevation_16_k10.asc"); nberOfSimulations = 50
simulationDirectory = "All_RRW_simulations"; simulationDirectories[[1]] = simulationDirectory

simulationDirectory0 = "All_MCC_simulations"; nberOfSimulations = 10; simulationDirectories = list(); envVariable = raster("Elevation_16_k10.asc")
simulationDirectories[[1]] = "RABV_US1_LC_s_1"; simulationDirectories[[2]] = "RABV_US1_LC_s_2"; simulationDirectories[[3]] = "RABV_US1_LC_s_3"
simulationDirectories[[4]] = "RABV_US1_CS_s_1"; simulationDirectories[[5]] = "RABV_US1_CS_s_2"; simulationDirectories[[6]] = "RABV_US1_CS_s_3"
for (i in 1:length(simulationDirectories)) simulationDirectories[[i]] = paste0(simulationDirectory0,"/",simulationDirectories[[i]])

fastaIDtransformation = function(fas, tab)
	{
		c = 0; new = fas
		for (i in 1:length(fas))
			{
				if (grepl(">",fas[i]))
					{
						c = c+1; new[i] = paste(unlist(strsplit(fas[i],"_"))[1],unlist(strsplit(fas[i],"_"))[2],tab[c,3],tab[c,4],sep="_")
					}
			}
		return(new)
	}
if (!file.exists("Cartogram_transf.rds"))
	{
		envVariable_mod = envVariable; envVariable_mod[] = envVariable_mod[] + 1; envVariable_mod = rasterToPolygons(envVariable_mod)
		names(envVariable_mod) = "envVariable"; envVariable_mod_sf = st_transform(st_as_sf(envVariable_mod), 3857)
		cartogram = cartogram_cont(envVariable_mod_sf, "envVariable", itermax=5); cartogram = st_transform(cartogram, crs(envVariable_mod))
		cartogram = sf:::as_Spatial(cartogram); cartograms = list(); cartograms[[1]] = cartogram; saveRDS(cartogram, "Cartogram_transf.rds")
	}	else	{
		cartogram = readRDS("Cartogram_transf.rds"); cartograms = list(); cartograms[[1]] = cartogram
	}
for (h in 1:length(simulationDirectories))
	{
		simulationDirectory = simulationDirectories[[h]]; dir.create(paste0(simulationDirectory,"/All_BEAST_analyses"), showWarnings=F)
		resistances = c(TRUE); avgResistances = c(TRUE); fourCells = TRUE; OS = "Unix"; showingPlots = FALSE
		for (i in 1:nberOfSimulations)
			{
				if (grepl("RRW",simulationDirectory))
					{							
						tab = read.csv(paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",i,".csv"), head=T)
						tab = tab[which(!tab[,"node2"]%in%tab[,"node1"]),c("node2","endYear","endLat","endLon")]
						colnames(tab) = c("tipNodeID","collection_date","latitude","longitude"); row.names(tab) = c()
						txt = tab[,c("collection_date","latitude","longitude")]; row.names(txt) = tab[,"tipNodeID"]; tmp = tab[,c("latitude","longitude")]
						write.table(txt, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_original_coordinates.txt"), sep="\t", quote=F)
						pathModel = 2; input = tab; outputName = paste0("TreeSimulations_",i); envVariables = list(envVariable); envVariables[[1]][] = envVariables[[1]][] + 1
						mdsTransformation(input, envVariables, pathModel, resistances, avgResistances, fourCells, outputName, OS)
						tmp1 = read.table(paste0("TreeSimulations_",i,"_MDS_nullRaster_LC_R.txt")); tmp2 = read.table(paste0("TreeSimulations_",i,"_MDS_Elevation_16_k10_LC_R.txt"))
						tab1 = tab; tab1[,c("latitude","longitude")] = tmp1; tab2 = tab; tab2[,c("latitude","longitude")] = tmp2		
						txt1 = tab1[,c("collection_date","latitude","longitude")]; row.names(txt1) = tab1[,"tipNodeID"]
						write.table(txt1, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_LC_null_raster.txt"), sep="\t", quote=F)
						txt2 = tab2[,c("collection_date","latitude","longitude")]; row.names(txt2) = tab2[,"tipNodeID"]
						write.table(txt2, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_LC_env_raster.txt"), sep="\t", quote=F)
						file.remove(paste0("TreeSimulations_",i,"_MDS_nullRaster_LC_R.txt")); file.remove(paste0("TreeSimulations_",i,"_MDS_Elevation_16_k10_LC_R.txt"))
						pathModel = 3; input = tab; outputName = paste0("TreeSimulations_",i); envVariables = list(envVariable); envVariables[[1]][] = envVariables[[1]][] + 1
						mdsTransformation(input, envVariables, pathModel, resistances, avgResistances, fourCells, outputName, OS)
						tmp1 = read.table(paste0("TreeSimulations_",i,"_MDS_nullRaster_CS_R.txt")); tmp2 = read.table(paste0("TreeSimulations_",i,"_MDS_Elevation_16_k10_CS_R.txt"))
						tab1 = tab; tab1[,c("latitude","longitude")] = tmp1; tab2 = tab; tab2[,c("latitude","longitude")] = tmp2		
						txt1 = tab1[,c("collection_date","latitude","longitude")]; row.names(txt1) = tab1[,"tipNodeID"]
						write.table(txt1, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_CS_null_raster.txt"), sep="\t", quote=F)
						txt2 = tab2[,c("collection_date","latitude","longitude")]; row.names(txt2) = tab2[,"tipNodeID"]
						write.table(txt2, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_CS_env_raster.txt"), sep="\t", quote=F)
						file.remove(paste0("TreeSimulations_",i,"_MDS_nullRaster_CS_R.txt")); file.remove(paste0("TreeSimulations_",i,"_MDS_Elevation_16_k10_CS_R.txt"))
						input = tab; envVariables = list(envVariable); envVariables[[1]][] = envVariables[[1]][] + 1; outputName = paste0("TreeSimulations_",i)
						cartogramTransformation(input, envVariables, resistances, outputName, cartograms)
						tmp1 = read.table(paste0("TreeSimulations_",i,"_cartogram_Elevation_16_k10_R.txt"))
						tab1 = tab; tab1[,c("latitude","longitude")] = tmp1; txt1 = tab1[,c("collection_date","latitude","longitude")]; row.names(txt1) = tab1[,"tipNodeID"]
						write.table(txt1, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_cartogram_transform.txt"), sep="\t", quote=F)
						file.remove(paste0("TreeSimulations_",i,"_cartogram_Elevation_16_k10_R.txt")); unlink("CS_rasters")
					}
				if (grepl("MCC",simulationDirectory))
					{
						fastaPrefix = gsub(paste0(simulationDirectory0,"/"),"",gsub("_s_","_scenario_",gsub("original","original_tree",simulationDirectory)))
						fas = scan(paste0(simulationDirectory,"/All_BEAST_analyses/",fastaPrefix,"_simulation_",i,".fas"), what="", sep="\n", quiet=T)
						seqIDs = gsub(">","",fas[which(grepl(">",fas))]); tab = matrix(nrow=length(seqIDs), ncol=4); colnames(tab) = c("tipNodeID","collection_date","latitude","longitude")
						for (j in 1:dim(tab)[1]) tab[j,] = c(seqIDs[j], as.numeric(unlist(strsplit(seqIDs[j],"_"))[c(2,4,3)]))
						new = fastaIDtransformation(fas, tab); write(new, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_original_coordinates.fas"))
						pathModel = 2; input = as.data.frame(tab); outputName = paste0("TreeSimulations_",i); envVariables = list(envVariable); envVariables[[1]][] = envVariables[[1]][] + 1
						mdsTransformation(input, envVariables, pathModel, resistances, avgResistances, fourCells, outputName, OS)
						tmp1 = read.table(paste0("TreeSimulations_",i,"_MDS_nullRaster_LC_R.txt")); tmp2 = read.table(paste0("TreeSimulations_",i,"_MDS_Elevation_16_k10_LC_R.txt"))
						tab1 = as.data.frame(tab); tab1[,c("latitude","longitude")] = tmp1; tab2 = as.data.frame(tab); tab2[,c("latitude","longitude")] = tmp2		
						new = fastaIDtransformation(fas, tab1); write(new, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_LC_null_raster.fas"))
						new = fastaIDtransformation(fas, tab2); write(new, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_LC_env_raster.fas"))
						file.remove(paste0("TreeSimulations_",i,"_MDS_nullRaster_LC_R.txt")); file.remove(paste0("TreeSimulations_",i,"_MDS_Elevation_16_k10_LC_R.txt"))
						pathModel = 3; input = as.data.frame(tab); outputName = paste0("TreeSimulations_",i); envVariables = list(envVariable); envVariables[[1]][] = envVariables[[1]][] + 1
						mdsTransformation(input, envVariables, pathModel, resistances, avgResistances, fourCells, outputName, OS)
						tmp1 = read.table(paste0("TreeSimulations_",i,"_MDS_nullRaster_CS_R.txt")); tmp2 = read.table(paste0("TreeSimulations_",i,"_MDS_Elevation_16_k10_CS_R.txt"))
						tab1 = as.data.frame(tab); tab1[,c("latitude","longitude")] = tmp1; tab2 = as.data.frame(tab); tab2[,c("latitude","longitude")] = tmp2		
						new = fastaIDtransformation(fas, tab1); write(new, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_CS_null_raster.fas"))
						new = fastaIDtransformation(fas, tab2); write(new, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_CS_env_raster.fas"))
						file.remove(paste0("TreeSimulations_",i,"_MDS_nullRaster_CS_R.txt")); file.remove(paste0("TreeSimulations_",i,"_MDS_Elevation_16_k10_CS_R.txt"))
						input = as.data.frame(tab); envVariables = list(envVariable); envVariables[[1]][] = envVariables[[1]][] + 1; outputName = paste0("TreeSimulations_",i)
						cartogramTransformation(input, envVariables, resistances, outputName, cartograms)
						tmp1 = read.table(paste0("TreeSimulations_",i,"_cartogram_Elevation_16_k10_R.txt")); tab1 = as.data.frame(tab); tab1[,c("latitude","longitude")] = tmp1
						new = fastaIDtransformation(fas, tab1); write(new, paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_cartogram_transform.fas"))
						file.remove(paste0("TreeSimulations_",i,"_cartogram_Elevation_16_k10_R.txt")); unlink("CS_rasters")
					}					
			}
	}

# 4. Preparing the XML files for the different continuous phylogeographic analyses to conduct with BEAST

simulationDirectories = list(); simulationDirectories[[1]] = "All_RRW_simulations"; nberOfSimulations = 50

simulationDirectory0 = "All_MCC_simulations"; nberOfSimulations = 10; simulationDirectories = list(); envVariable = raster("Elevation_16_k10.asc")
simulationDirectories[[1]] = "RABV_US1_LC_s_1"; simulationDirectories[[2]] = "RABV_US1_LC_s_2"; simulationDirectories[[3]] = "RABV_US1_LC_s_3"
simulationDirectories[[4]] = "RABV_US1_CS_s_1"; simulationDirectories[[5]] = "RABV_US1_CS_s_2"; simulationDirectories[[6]] = "RABV_US1_CS_s_3"
for (i in 1:length(simulationDirectories)) simulationDirectories[[i]] = paste0(simulationDirectory0,"/",simulationDirectories[[i]])

analyses = c("original_coordinates","MDS_LC_null_raster","MDS_LC_env_raster","MDS_CS_null_raster","MDS_CS_env_raster","cartogram_transform")
for (h in 1:length(simulationDirectories))
	{
		simulationDirectory = simulationDirectories[[h]]; jitter = 0.01
		for (i in 1:length(analyses))
			{
				for (j in 1:nberOfSimulations)
					{
						if (grepl("RRW",simulationDirectory))
							{
								templateFile = "Template1_BEAST.xml"; prefix = paste0("TreeSimulations_",j,"_",analyses[i])
								newXMLFile = paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",j,"_",analyses[i],".xml")
								metadata = read.table(paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",j,"_",analyses[i],".txt"))
								treeFile = paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",j,".tre")
								generateXMLfiles1(templateFile, newXMLFile, prefix, treeFile, metadata, jitter)
							}
						if (grepl("MCC",simulationDirectory))
							{
								templateFile = "Template2_BEAST.xml"; prefix = paste0("TreeSimulations_",j,"_",analyses[i])
								newXMLFile = paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",j,"_",analyses[i],".xml")
								fastaPrefix = gsub("All_MCC_simulations/","",gsub("_s_","_scenario_",gsub("original","original_tree",simulationDirectory)))
								fastaFile = paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",j,"_",analyses[i],".fas")
								generateXMLfiles2(templateFile, newXMLFile, prefix, fastaFile, jitter)
							}
					}
			}
	}

# 5. Performing all the different landscape phylogeographic analyses (post hoc and prior-informed analyses)

simulationDirectories = list("All_RRW_simulations"); resistance = TRUE; envVariable = raster("Elevation_16_k10.asc"); nberOfSimulations = 50; nberOfExtractionFiles = 100

simulationDirectory0 = "All_MCC_simulations"; nberOfSimulations = 10; simulationDirectories = list(); envVariable = raster("Elevation_16_k10.asc")
simulationDirectories[[1]] = "RABV_US1_LC_s_1"; simulationDirectories[[2]] = "RABV_US1_LC_s_2"; simulationDirectories[[3]] = "RABV_US1_LC_s_3"
simulationDirectories[[4]] = "RABV_US1_CS_s_1"; simulationDirectories[[5]] = "RABV_US1_CS_s_2"; simulationDirectories[[6]] = "RABV_US1_CS_s_3"
for (i in 1:length(simulationDirectories)) simulationDirectories[[i]] = paste0(simulationDirectory0,"/",simulationDirectories[[i]])

	# 5.1. All tree extractions

for (g in 1:length(simulationDirectories))
	{
		simulationDirectory = simulationDirectories[[g]]
		analyses = c("original_coordinates","MDS_LC_null_raster","MDS_LC_env_raster","MDS_CS_null_raster","MDS_CS_env_raster","cartogram_transform")
		for (h in 1:length(analyses))
			{
				for (i in 1:nberOfSimulations)
					{
						localTreesDirectory = paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_",analyses[h],"_ext"); dir.create(localTreesDirectory, showWarnings=F)		
						allTrees = scan(file=paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_",analyses[h],".trees"), what="", sep="\n", quiet=T)
						if (grepl("RRW_simulations",simulationDirectory))
							{
								mostRecentSamplingDatum = max(read.csv(paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",i,".csv"), head=T)[,"endYear"])
							}
						if (grepl("MCC_simulations",simulationDirectory))
							{
								mostRecentSamplingDatum = max(read.csv(paste0(simulationDirectory,"/RABV_US1_simu.csv"), head=T)[,"endYear"])
							}
						burnIn = 51; randomSampling = FALSE; nberOfTreesToSample = 100; coordinateAttributeName = "location"; nberOfCores = 10
						# treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores) # PROBLEM
						trees = readAnnotatedNexus(paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_",analyses[h],".trees"))
						indices = sample((burnIn+1):length(trees), 100, replace=F)
						for (j in 1:length(indices))
							{
								tab = postTreeExtractions(trees[[indices[j]]], mostRecentSamplingDatum)
								write.csv(tab, paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
							}
						if (analyses[h] == "original_coordinates")
							{
								if (grepl("RRW_simulations",simulationDirectory))
									{
										file.copy(paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",i,".pdf"),
												  paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_RRW_tree_simulation.pdf"), overwrite=T)
										tree = read.tree(paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",i,".tre"))
									}
								if (grepl("MCC_simulations",simulationDirectory))
									{
										file.copy(paste0(simulationDirectory,"/RABV_US1_simu.pdf"),
												  paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_MCC_tree_simulation.pdf"), overwrite=T)
										tree = trees[[indices[1]]]
									}
								pdf(paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_1st_tree_extraction.pdf"), width=8.0, height=4.0)
								par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(0.3,0,0,0.3), mar=c(1.5,0.5,0.5,0), lwd=0.4, col="gray30")
								tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_1.csv"), head=T); startingYear = min(tab[,"startYear"])
								for (j in 1:dim(tab)[1])
									{
										if (tab[j,"node1"]%in%tab[,"node2"])
											{
												index = which(tab[,"node2"]==tab[j,"node1"])
												if (tab[index,"endLon"] != tab[j,"startLon"]) print(c(i,j))
											}
									}
								col_start = colorRampPalette(brewer.pal(9,"PuBu"))(101)[1]; samplingWindow = cbind(1982.2,2004.7)
								cols3 = colorRampPalette(brewer.pal(9,"PuBu"))(101)[(((nodeHeights(tree)[,2])/(samplingWindow[2]-startingYear))*100)+1]
								plot(tree, show.tip.label=F, edge.width=0.4, edge.col="gray30")
								for (j in 1:dim(tree$edge)[1])
									{
										if (j == 1)
											{
												nodelabels(node=tree$edge[j,1], pch=16, cex=0.7, col=col_start)
												nodelabels(node=tree$edge[j,1], pch=1, cex=0.7, col="gray30", lwd=0.25)
											}
										nodelabels(node=tree$edge[j,2], pch=16, cex=0.7, col=cols3[j])
										nodelabels(node=tree$edge[j,2], pch=1, cex=0.7, col="gray30", lwd=0.25)
									}
								axis(side=1,lwd.tick=0,cex.axis=0.65,mgp=c(0,0.1,0),lwd=0.4,tck=0,col.tick="gray30",col.axis="gray30",col="gray30",at=c(0,335),labels=rep(NA,2))
								axis(side=1,lwd.tick=0.4,cex.axis=0.65,mgp=c(0,0.1,0),lwd=0,tck=-0.016,col.tick="gray30",col.axis="gray30",col="gray30",at=seq(-2,38,5), labels=seq(1970,2010,5))
								plot(extent(envVariable), axes=F, ann=F, col="white")
								plot(envVariable, add=T, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(101)[1:101], colNA="gray90", legend=F, useRaster=T, asp=1, axes=F, box=F)
								ext = extent(envVariable); rect(ext@xmin, ext@ymin, ext@xmax, ext@ymax, lwd=0.4, border="gray30", col=NA)
								col_start = colorRampPalette(brewer.pal(9,"PuBu"))(101)[1]
								cols2 = colorRampPalette(brewer.pal(9,"PuBu"))(101)[(((tab[,"endYear"]-startingYear)/(samplingWindow[2]-startingYear))*100)+1]
								for (j in 1:dim(tab)[1])
									{	
										curvedarrow(cbind(tab[j,"startLon"],tab[j,"startLat"]), cbind(tab[j,"endLon"],tab[j,"endLat"]), arr.length=0,
													arr.width=0, lwd=0.3, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
									}
								for (j in dim(tab)[1]:1)
									{
										if (j == 1)
											{
												points(cbind(tab[j,"startLon"],tab[j,"startLat"]), pch=16, cex=0.7, col=col_start)
												points(cbind(tab[j,"startLon"],tab[j,"startLat"]), pch=1, cex=0.7, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
											}
										points(cbind(tab[j,"endLon"],tab[j,"endLat"]), pch=16, cex=0.7, col=cols2[j])
										points(cbind(tab[j,"endLon"],tab[j,"endLat"]), pch=1, cex=0.7, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
									}
								plot(envVariable, legend.only=T, add=T, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(101)[1:101], legend.width=0.5, legend.shrink=0.3,
									 smallplot=c(0.060,0.965,0.079,0.094), legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
									 axis.args=list(cex.axis=0.65, lwd=0, lwd.tick=0.4, tck=-1.0, col.tick="gray30", col.axis="gray30", line=0, mgp=c(0,0.1,0), at=seq(1,10,1)))
								dev.off()
							}
					}
			}
	}

	# 5.2. Post hoc analyses

		# 5.2.1. Isolation-by-resistance analyses

for (g in 1:length(simulationDirectories))
	{
		simulationDirectory = simulationDirectories[[g]]
		for (i in 1:nberOfSimulations)
			{
				localTreesDirectory = paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_original_coordinates_ext")
				directory = paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_all_seraphim_analyses"); dir.create(directory, showWarnings=F)
				envVariables = list(envVariable); envVariables[[1]][] = envVariables[[1]][]+1; resistances = list(resistance); avgResistances = list(resistance)
				fourCells = TRUE; nberOfRandomisations = 1; randomProcedure = 3; showingPlots = FALSE; nberOfCores = 10; OS = "Unix"; juliaCSImplementation = FALSE
				pathModel = 2; outputName = paste0("TreeSimulations_",i,"_LC"); simulations = FALSE; randomisations = FALSE
				isolationByResistance(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances, fourCells,
									  nberOfRandomisations, randomProcedure, outputName, showingPlots, nberOfCores, OS, juliaCSImplementation)
				file.copy(paste0("TreeSimulations_",i,"_LC_linear_regression_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_seraphim_LR1.txt"), overwrite=T)
				file.copy(paste0("TreeSimulations_",i,"_LC_randomisation_BF_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_seraphim_BF1.txt"), overwrite=T)
				file.remove(paste0("TreeSimulations_",i,"_LC_linear_regression_results.txt")); file.remove(paste0("TreeSimulations_",i,"_LC_randomisation_BF_results.txt"))
				pathModel = 3; outputName = paste0("TreeSimulations_",i,"_CS"); simulations = FALSE; randomisations = FALSE
				isolationByResistance(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances, fourCells,
									  nberOfRandomisations, randomProcedure, outputName, showingPlots, nberOfCores, OS, juliaCSImplementation)
				file.copy(paste0("TreeSimulations_",i,"_CS_linear_regression_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_seraphim_LR1.txt"), overwrite=T)
				file.copy(paste0("TreeSimulations_",i,"_CS_randomisation_BF_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_seraphim_BF1.txt"), overwrite=T)
				file.remove(paste0("TreeSimulations_",i,"_CS_linear_regression_results.txt")); file.remove(paste0("TreeSimulations_",i,"_CS_randomisation_BF_results.txt"))
				unlink("CS_rasters")
			}
	}
for (g in 1:length(simulationDirectories))
	{
		simulationDirectory = simulationDirectories[[g]]; nberOfTips = rep(NA, nberOfSimulations)
		results_1 = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_1) = c("LC_beta_env","LC_R2_env","LC_Q","LC_pQ","LC_BF")
		results_2 = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_2) = c("CS_beta_env","CS_R2_env","CS_Q","CS_pQ","CS_BF")
		for (i in 1:nberOfSimulations)
			{
				directory = paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_all_seraphim_analyses")
				if (grepl("RRW_simulations",simulationDirectory))
					{
						tab = read.csv(paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",i,".csv"), head=T)
					}
				if (grepl("MCC_simulations",simulationDirectory))
					{
						tab = read.csv(paste0(simulationDirectory,"/RABV_US1_simu.csv"), head=T)
					}
				nberOfTips[i] = length(which(!tab[,"node2"]%in%tab[,"node1"]))
				directory = paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_all_seraphim_analyses")
				tab1 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_seraphim_LR1.txt"), head=T)
				tab2 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_seraphim_BF1.txt"), head=T)
				vS = tab1[,"LR_coefficients_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_1[i,"LC_beta_env"] = paste0(round(mean(vS),3)) # ," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR_R2_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_1[i,"LC_R2_env"] = paste0(round(mean(vS),3)) # ," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR_Q_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_1[i,"LC_Q"] = paste0(round(mean(vS),3)) # ," [",hds[1],", ",hds[2],"]")
				results_1[i,"LC_pQ"] = round(sum(vS>0)/length(vS),3)
				if (as.numeric(results_1[i,"LC_pQ"]) > 0.9)
					{
						results_1[i,"LC_BF"] = round(mean(tab2["Elevation_16_k10_R","BFs_Q_LR_randomisation_1"]),1)
					}	else	{
						results_1[i,"LC_BF"] = "-"
					}
				tab1 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_seraphim_LR1.txt"), head=T)
				tab2 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_seraphim_BF1.txt"), head=T)
				vS = tab1[,"LR_coefficients_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_2[i,"CS_beta_env"] = paste0(round(mean(vS),3)) # ," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR_R2_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_2[i,"CS_R2_env"] = paste0(round(mean(vS),3)) # ," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR_Q_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_2[i,"CS_Q"] = paste0(round(mean(vS),3)) # ," [",hds[1],", ",hds[2],"]")
				results_2[i,"CS_pQ"] = round(sum(vS>0)/length(vS),3)
				if (as.numeric(results_2[i,"CS_pQ"]) > 0.9)
					{
						results_2[i,"CS_BF"] = round(mean(tab2["Elevation_16_k10_R","BFs_Q_LR_randomisation_1"]),1)
					}	else	{
						results_2[i,"CS_BF"] = "-"
					}
			}
		write.table(results_1, paste0(simulationDirectory,"/Posthoc_IBR_LC.csv"), row.names=F, quote=F, sep=";")
		write.table(results_2, paste0(simulationDirectory,"/Posthoc_IBR_CS.csv"), row.names=F, quote=F, sep=";")
		results_1 = cleanResultsTable(read.csv(paste0(simulationDirectory,"/Posthoc_IBR_LC.csv"), head=T, sep=";"))
		write.table(results_1, paste0(simulationDirectory,"/Posthoc_IBR_LC.txt"), row.names=F, quote=F, sep="\t")
		results_2 = cleanResultsTable(read.csv(paste0(simulationDirectory,"/Posthoc_IBR_CS.csv"), head=T, sep=";"))
		write.table(results_2, paste0(simulationDirectory,"/Posthoc_IBR_CS.txt"), row.names=F, quote=F, sep="\t")
		vS = as.numeric(results_1[,"LC_Q"]); ci = round(quantile(vS,c(0.025,0.975)),3); BFs = as.numeric(gsub(">","",results_1[,"LC_BF"])); BFs = BFs[!is.na(BFs)]
		cat("LC: Q = ",round(mean(vS),3)," [",ci[1],", ",ci[2],"], p(Q>0) = ",sum(vS>0),"/",length(vS),", BF>3: ",sum(BFs>=3),"/",length(BFs),", BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
		vS = as.numeric(results_2[,"CS_Q"]); ci = round(quantile(vS,c(0.025,0.975)),3); BFs = as.numeric(gsub(">","",results_2[,"CS_BF"])); BFs = BFs[!is.na(BFs)]
		cat("CS: Q = ",round(mean(vS),3)," [",ci[1],", ",ci[2],"], p(Q>0) = ",sum(vS>0),"/",length(vS),", BF>3: ",sum(BFs>=3),"/",length(BFs),", BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
	}

		# 5.2.2. Lineage diffusion velocity analyses

for (g in 1:length(simulationDirectories))
	{
		simulationDirectory = simulationDirectories[[g]]
		for (i in 1:nberOfSimulations)
			{
				localTreesDirectory = paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_original_coordinates_ext")
				directory = paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_all_seraphim_analyses"); dir.create(directory, showWarnings=F)
				envVariables = list(envVariable); envVariables[[1]][] = envVariables[[1]][]+1; resistances = list(resistance); avgResistances = list(resistance)
				fourCells = TRUE; nberOfRandomisations = 1; randomProcedure = 3; juliaCSImplementation = FALSE; showingPlots = FALSE; nberOfCores = 10; OS = "Unix"
				pathModel = 2; outputName = paste0("TreeSimulations_",i,"_LC"); minimumConvexHull = TRUE; simulations = FALSE; randomisations = FALSE
				spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
							  randomProcedure, outputName, showingPlots, nberOfCores, OS, juliaCSImplementation, simulations, randomisations, minimumConvexHull)
				file.copy(paste0("TreeSimulations_",i,"_LC_linear_regression_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_seraphim_LR2.txt"), overwrite=T)
				file.copy(paste0("TreeSimulations_",i,"_LC_randomisation_BF_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_seraphim_BF2.txt"), overwrite=T)
				file.remove(paste0("TreeSimulations_",i,"_LC_linear_regression_results.txt")); file.remove(paste0("TreeSimulations_",i,"_LC_randomisation_BF_results.txt"))
				pathModel = 3; outputName = paste0("TreeSimulations_",i,"_CS"); minimumConvexHull = TRUE; simulations = FALSE; randomisations = FALSE
				spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
							  randomProcedure, outputName, showingPlots, nberOfCores, OS, juliaCSImplementation, simulations, randomisations, minimumConvexHull)
				file.copy(paste0("TreeSimulations_",i,"_CS_linear_regression_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_seraphim_LR2.txt"), overwrite=T)
				file.copy(paste0("TreeSimulations_",i,"_CS_randomisation_BF_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_seraphim_BF2.txt"), overwrite=T)
				file.remove(paste0("TreeSimulations_",i,"_CS_linear_regression_results.txt")); file.remove(paste0("TreeSimulations_",i,"_CS_randomisation_BF_results.txt"))
				pathModel = 2; outputName = paste0("TreeSimulations_",i,"_LC"); minimumConvexHull = FALSE; simulations = FALSE; randomisations = FALSE
				spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
							  randomProcedure, outputName, showingPlots, nberOfCores, OS, juliaCSImplementation, simulations, randomisations, minimumConvexHull)
				file.copy(paste0("TreeSimulations_",i,"_LC_linear_regression_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_woConvexH_LR2.txt"), overwrite=T)
				file.copy(paste0("TreeSimulations_",i,"_LC_randomisation_BF_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_woConvexH_BF2.txt"), overwrite=T)
				file.remove(paste0("TreeSimulations_",i,"_LC_linear_regression_results.txt")); file.remove(paste0("TreeSimulations_",i,"_LC_randomisation_BF_results.txt"))
				pathModel = 3; outputName = paste0("TreeSimulations_",i,"_CS"); minimumConvexHull = FALSE; simulations = FALSE; randomisations = FALSE
				spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
							  randomProcedure, outputName, showingPlots, nberOfCores, OS, juliaCSImplementation, simulations, randomisations, minimumConvexHull)
				file.copy(paste0("TreeSimulations_",i,"_CS_linear_regression_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_woConvexH_LR2.txt"), overwrite=T)
				file.copy(paste0("TreeSimulations_",i,"_CS_randomisation_BF_results.txt"), paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_woConvexH_BF2.txt"), overwrite=T)
				file.remove(paste0("TreeSimulations_",i,"_CS_linear_regression_results.txt")); file.remove(paste0("TreeSimulations_",i,"_CS_randomisation_BF_results.txt"))
				unlink("CS_rasters")
			}
	}
analysesWithoutConvexHull = FALSE
for (g in 1:length(simulationDirectories))
	{
		simulationDirectory = simulationDirectories[[g]]; nberOfTips = rep(NA, nberOfSimulations)
		results_1 = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_1) = c("LC_LR1_beta_env","LC_LR1_R2_env","LC_LR1_Q","LC_LR1_pQ","LC_LR1_BF")
		results_2 = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_2) = c("LC_LR2_beta_env","LC_LR2_R2_env","LC_LR2_Q","LC_LR2_pQ","LC_LR2_BF")
		results_3 = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_3) = c("CS_LR1_beta_env","CS_LR1_R2_env","CS_LR1_Q","CS_LR1_pQ","CS_LR1_BF")
		results_4 = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_4) = c("CS_LR2_beta_env","CS_LR2_R2_env","CS_LR2_Q","CS_LR2_pQ","CS_LR2_BF")
		for (i in 1:nberOfSimulations)
			{
				if (grepl("RRW_simulations",simulationDirectory))
					{
						tab = read.csv(paste0(simulationDirectory,"/Selected_simulations/TreeSimulations_",i,".csv"), head=T)
					}
				if (grepl("MCC_simulations",simulationDirectory))
					{
						tab = read.csv(paste0(simulationDirectory,"/RABV_US1_simu.csv"), head=T)
					}
				nberOfTips[i] = length(which(!tab[,"node2"]%in%tab[,"node1"]))
				directory = paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_all_seraphim_analyses")
				if (analysesWithoutConvexHull == FALSE)
					{
						tab1 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_seraphim_LR2.txt"), head=T)
						tab2 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_seraphim_BF2.txt"), head=T)
					}	else	{
						tab1 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_woConvexH_LR2.txt"), head=T)
						tab2 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_LC_woConvexH_BF2.txt"), head=T)						
					}
				vS = tab1[,"LR1_coefficients_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_1[i,"LC_LR1_beta_env"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR1_R2_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_1[i,"LC_LR1_R2_env"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR1_Q_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_1[i,"LC_LR1_Q"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				results_1[i,"LC_LR1_pQ"] = round(sum(vS>0)/length(vS),3)
				if (as.numeric(results_1[i,"LC_LR1_pQ"]) > 0.9)
					{
						results_1[i,"LC_LR1_BF"] = round(mean(tab2["Elevation_16_k10_R","BFs_Q_LR1_randomisation_1"]),1)
					}	else	{
						results_1[i,"LC_LR1_BF"] = "-"
					}
				vS = tab1[,"LR2_coefficients_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_2[i,"LC_LR2_beta_env"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR2_R2_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_2[i,"LC_LR2_R2_env"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR2_Q_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_2[i,"LC_LR2_Q"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				results_2[i,"LC_LR2_pQ"] = round(sum(vS>0)/length(vS),3)
				if (as.numeric(results_2[i,"LC_LR2_pQ"]) > 0.9)
					{
						results_2[i,"LC_LR2_BF"] = round(mean(tab2["Elevation_16_k10_R","BFs_Q_LR2_randomisation_1"]),1)
					}	else	{
						results_2[i,"LC_LR2_BF"] = "-"
					}
				if (analysesWithoutConvexHull == FALSE)
					{
						tab1 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_seraphim_LR2.txt"), head=T)
						tab2 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_seraphim_BF2.txt"), head=T)
					}	else	{
						tab1 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_woConvexH_LR2.txt"), head=T)
						tab2 = read.table(paste0(directory,"/TreeSimulations_",i,"_original_coordinates_CS_woConvexH_BF2.txt"), head=T)						
					}
				vS = tab1[,"LR1_coefficients_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_3[i,"CS_LR1_beta_env"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR1_R2_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_3[i,"CS_LR1_R2_env"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR1_Q_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_3[i,"CS_LR1_Q"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				results_3[i,"CS_LR1_pQ"] = round(sum(vS>0)/length(vS),3)
				if (as.numeric(results_3[i,"CS_LR1_pQ"]) > 0.9)
					{
						results_3[i,"CS_LR1_BF"] = round(mean(tab2["Elevation_16_k10_R","BFs_Q_LR1_randomisation_1"]),1)
					}	else	{
						results_3[i,"CS_LR1_BF"] = "-"
					}
				vS = tab1[,"LR2_coefficients_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_4[i,"CS_LR2_beta_env"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR2_R2_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_4[i,"CS_LR2_R2_env"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				vS = tab1[,"LR2_Q_Elevation_16_k10_R"]; hds = round(HDInterval::hdi(vS)[1:2],3)
				results_4[i,"CS_LR2_Q"] = paste0(round(mean(vS),3)," [",hds[1],", ",hds[2],"]")
				results_4[i,"CS_LR2_pQ"] = round(sum(vS>0)/length(vS),3)
				if (as.numeric(results_4[i,"CS_LR2_pQ"]) > 0.9)
					{
						results_4[i,"CS_LR2_BF"] = round(mean(tab2["Elevation_16_k10_R","BFs_Q_LR2_randomisation_1"]),1)
					}	else	{
						results_4[i,"CS_LR2_BF"] = "-"
					}
			}
		if (analysesWithoutConvexHull == FALSE)
			{
				write.table(results_1, paste0(simulationDirectory,"/Posthoc_LDV_LC.csv"), row.names=F, quote=F, sep=";")
				write.table(results_2, paste0(simulationDirectory,"/Posthoc_DCs_LC.csv"), row.names=F, quote=F, sep=";")
				write.table(results_3, paste0(simulationDirectory,"/Posthoc_LDV_CS.csv"), row.names=F, quote=F, sep=";")
				write.table(results_4, paste0(simulationDirectory,"/Posthoc_DCs_CS.csv"), row.names=F, quote=F, sep=";")
				results_1 = cleanResultsTable(read.csv(paste0(simulationDirectory,"/Posthoc_LDV_LC.csv"), head=T, sep=";"))
				write.table(results_1, paste0(simulationDirectory,"/Posthoc_LDV_LC.txt"), row.names=F, quote=F, sep="\t")
				results_2 = cleanResultsTable(read.csv(paste0(simulationDirectory,"/Posthoc_DCs_LC.csv"), head=T, sep=";"))
				write.table(results_2, paste0(simulationDirectory,"/Posthoc_DCs_LC.txt"), row.names=F, quote=F, sep="\t")
				results_3 = cleanResultsTable(read.csv(paste0(simulationDirectory,"/Posthoc_LDV_CS.csv"), head=T, sep=";"))
				write.table(results_3, paste0(simulationDirectory,"/Posthoc_LDV_CS.txt"), row.names=F, quote=F, sep="\t")
				results_4 = cleanResultsTable(read.csv(paste0(simulationDirectory,"/Posthoc_DCs_CS.csv"), head=T, sep=";"))
				write.table(results_4, paste0(simulationDirectory,"/Posthoc_DCs_CS.txt"), row.names=F, quote=F, sep="\t")
			}	else	{
				write.table(results_1, paste0(simulationDirectory,"/Posthoc_LDV_LC_woConvexH.csv"), row.names=F, quote=F, sep=";")
				write.table(results_2, paste0(simulationDirectory,"/Posthoc_DCs_LC_woConvexH.csv"), row.names=F, quote=F, sep=";")
				write.table(results_3, paste0(simulationDirectory,"/Posthoc_LDV_CS_woConvexH.csv"), row.names=F, quote=F, sep=";")
				write.table(results_4, paste0(simulationDirectory,"/Posthoc_DCs_CS_woConvexH.csv"), row.names=F, quote=F, sep=";")
			}
		vS1 = results_1[,"LC_LR1_Q"]; vS2 = as.numeric(results_1[,"LC_LR1_pQ"]); tmp = rep(NA, length(vS1))
		for (i in 1:length(vS1)) tmp[i] = as.numeric(unlist(strsplit(vS1[i]," \\["))[1])
		vS1 = tmp; ci1 = round(quantile(vS1,c(0.025,0.975)),3); ci2 = round(quantile(vS2,c(0.025,0.975)),2)
		BFs = as.numeric(gsub(">","",results_1[,"LC_LR1_BF"])); BFs = BFs[!is.na(BFs)]
		cat("LR1, LC: Q = ",round(mean(vS1),3)," [",ci1[1],", ",ci1[2],"], p(Q>0) = ",round(mean(vS2),2)," [",ci2[1],", ",ci2[2],"], ",sep="")
		cat("BF>3: ",sum(BFs>=3),"/",length(BFs)," BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
		vS1 = results_3[,"CS_LR1_Q"]; vS2 = as.numeric(results_3[,"CS_LR1_pQ"]); tmp = rep(NA, length(vS1))
		for (i in 1:length(vS1)) tmp[i] = as.numeric(unlist(strsplit(vS1[i]," \\["))[1])
		vS1 = tmp; ci1 = round(quantile(vS1,c(0.025,0.975)),3); ci2 = round(quantile(vS2,c(0.025,0.975)),2)
		BFs = as.numeric(gsub(">","",results_3[,"CS_LR1_BF"])); BFs = BFs[!is.na(BFs)]
		cat("LR1, CS: Q = ",round(mean(vS1),3)," [",ci1[1],", ",ci1[2],"], p(Q>0) = ",round(mean(vS2),2)," [",ci2[1],", ",ci2[2],"], ",sep="")
		cat("BF>3: ",sum(BFs>=3),"/",length(BFs)," BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
		vS1 = results_2[,"LC_LR2_Q"]; vS2 = as.numeric(results_2[,"LC_LR2_pQ"]); tmp = rep(NA, length(vS1))
		for (i in 1:length(vS1)) tmp[i] = as.numeric(unlist(strsplit(vS1[i]," \\["))[1])
		vS1 = tmp; ci1 = round(quantile(vS1,c(0.025,0.975)),3); ci2 = round(quantile(vS2,c(0.025,0.975)),2)
		BFs = as.numeric(gsub(">","",results_2[,"LC_LR2_BF"])); BFs = BFs[!is.na(BFs)]
		cat("LR2, LC: Q = ",round(mean(vS1),3)," [",ci1[1],", ",ci1[2],"], p(Q>0) = ",round(mean(vS2),2)," [",ci2[1],", ",ci2[2],"], ",sep="")
		cat("BF>3: ",sum(BFs>=3),"/",length(BFs)," BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
		vS1 = results_4[,"CS_LR2_Q"]; vS2 = as.numeric(results_4[,"CS_LR2_pQ"]); tmp = rep(NA, length(vS1))
		for (i in 1:length(vS1)) tmp[i] = as.numeric(unlist(strsplit(vS1[i]," \\["))[1])
		vS1 = tmp; ci1 = round(quantile(vS1,c(0.025,0.975)),3); ci2 = round(quantile(vS2,c(0.025,0.975)),2)
		BFs = as.numeric(gsub(">","",results_4[,"CS_LR2_BF"])); temp = BFs; BFs = BFs[!is.na(BFs)]
		cat("LR2, CS: Q = ",round(mean(vS1),3)," [",ci1[1],", ",ci1[2],"], p(Q>0) = ",round(mean(vS2),2)," [",ci2[1],", ",ci2[2],"], ",sep="")
		cat("BF>3: ",sum(BFs>=3),"/",length(BFs)," BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
		# plot(nberOfTips[which(!is.na(temp))], BFs); print(nberOfTips[which(is.na(temp))])
	}

	# 5.3. Prior-informed analyses

for (g in 1:length(simulationDirectories))
	{
		simulationDirectory = simulationDirectories[[g]]
		results_1a = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_1a) = c("LC_R2_null","LC_R2_env","LC_Q","LC_pQ","LC_BF")
		results_2a = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_2a) = c("CS_R2_null","CS_R2_env","CS_Q","CS_pQ","CS_BF")
		results_3a = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_3a) = c("CGT_R2_null","CGT_R2_env","CGT_Q","CGT_pQ","CGT_BF")
		results_1b = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_1b) = c("LC_R2_null","LC_R2_env","LC_Q","LC_pQ","LC_BF")
		results_2b = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_2b) = c("CS_R2_null","CS_R2_env","CS_Q","CS_pQ","CS_BF")
		results_3b = matrix(nrow=nberOfSimulations, ncol=5); colnames(results_3b) = c("CGT_R2_null","CGT_R2_env","CGT_Q","CGT_pQ","CGT_BF")
		for (i in 1:nberOfSimulations)
			{
				R2_mdsLCNul = read.table(paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_LC_null_raster.log"), header=T)[52:501,"location.squareddistTime4.Rsquared"]
				R2_mdsLCEnv = read.table(paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_LC_env_raster.log"), header=T)[52:501,"location.squareddistTime4.Rsquared"]
				hds = round(HDInterval::hdi(R2_mdsLCNul),3); meanV = mean(R2_mdsLCNul); results_1a[i,"LC_R2_null"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				hds = round(HDInterval::hdi(R2_mdsLCEnv),3); meanV = mean(R2_mdsLCEnv); results_1a[i,"LC_R2_env"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				vS = R2_mdsLCEnv-R2_mdsLCNul; hds = round(HDInterval::hdi(vS)[1:2],3); meanV = mean(vS); results_1a[i,"LC_Q"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				p = sum(R2_mdsLCEnv>R2_mdsLCNul)/length(R2_mdsLCNul); results_1a[i,"LC_pQ"] = round(p,2); results_1a[i,"LC_BF"] = round(p/(1-p),1)			
				R2_mdsCSNul = read.table(paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_CS_null_raster.log"), header=T)[52:501,"location.squareddistTime4.Rsquared"]
				R2_mdsCSEnv = read.table(paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_MDS_CS_env_raster.log"), header=T)[52:501,"location.squareddistTime4.Rsquared"]
				hds = round(HDInterval::hdi(R2_mdsCSNul),3); meanV = mean(R2_mdsCSNul); results_2a[i,"CS_R2_null"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				hds = round(HDInterval::hdi(R2_mdsCSEnv),3); meanV = mean(R2_mdsCSEnv); results_2a[i,"CS_R2_env"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				vS = R2_mdsCSEnv-R2_mdsCSNul; hds = round(HDInterval::hdi(vS)[1:2],3); meanV = mean(vS); results_2a[i,"CS_Q"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				p = sum(R2_mdsCSEnv>R2_mdsCSNul)/length(R2_mdsCSNul); results_2a[i,"CS_pQ"] = round(p,2); results_2a[i,"CS_BF"] = round(p/(1-p),1)
				R2_cartoNul = read.table(paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_original_coordinates.log"), header=T)[52:501,"location.squareddistTime4.Rsquared"]
				R2_cartoEnv = read.table(paste0(simulationDirectory,"/All_BEAST_analyses/TreeSimulations_",i,"_cartogram_transform.log"), header=T)[52:501,"location.squareddistTime4.Rsquared"]
				hds = round(HDInterval::hdi(R2_cartoNul),3); meanV = mean(R2_cartoNul); results_3a[i,"CGT_R2_null"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				hds = round(HDInterval::hdi(R2_cartoEnv),3); meanV = mean(R2_cartoEnv); results_3a[i,"CGT_R2_env"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				vS = R2_cartoEnv-R2_cartoNul; hds = round(HDInterval::hdi(vS)[1:2],3); meanV = mean(vS); results_3a[i,"CGT_Q"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				p = sum(R2_cartoEnv>R2_cartoNul)/length(R2_cartoNul); results_3a[i,"CGT_pQ"] = round(p,2); results_3a[i,"CGT_BF"] = round(p/(1-p),1)
				results = matrix(nrow=nberOfExtractionFiles, ncol=6); colnames(results) = c("MDS_LC_R2_null","MDS_LC_R2_env","MDS_CS_R2_null","MDS_CS_R2_env","CGT_R2_null","CGT_R2_env")
				for (j in 1:nberOfExtractionFiles)
					{
						tab = read.csv(paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_MDS_LC_null_raster_ext/TreeExtractions_",j,".csv"))
						x1 = tab[,"startLon"]; y1 = tab[,"startLat"]; x2 = tab[,"endLon"]; y2 = tab[,"endLat"]
						durations4 = 4*tab[,"length"]; distances = sqrt(((x2-x1)^2)+((y2-y1)^2)); squareDistances = distances^2 # --> Euclidean distances !!
			 			LM = lm(durations4 ~ squareDistances); results[j,"MDS_LC_R2_null"] = summary(LM)$r.squared
						tab = read.csv(paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_MDS_LC_env_raster_ext/TreeExtractions_",j,".csv"))
						x1 = tab[,"startLon"]; y1 = tab[,"startLat"]; durations4 = 4*tab[,"length"]; x2 = tab[,"endLon"]; y2 = tab[,"endLat"]
						durations4 = 4*tab[,"length"]; distances = sqrt(((x2-x1)^2)+((y2-y1)^2)); squareDistances = distances^2 # --> Euclidean distances !!
			 			LM = lm(durations4 ~ squareDistances); results[j,"MDS_LC_R2_env"] = summary(LM)$r.squared
						tab = read.csv(paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_MDS_CS_null_raster_ext/TreeExtractions_",j,".csv"))
						x1 = tab[,"startLon"]; y1 = tab[,"startLat"]; x2 = tab[,"endLon"]; y2 = tab[,"endLat"]
						durations4 = 4*tab[,"length"]; distances = sqrt(((x2-x1)^2)+((y2-y1)^2)); squareDistances = distances^2 # --> Euclidean distances !!
			 			LM = lm(durations4 ~ squareDistances); results[j,"MDS_CS_R2_null"] = summary(LM)$r.squared
						tab = read.csv(paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_MDS_CS_env_raster_ext/TreeExtractions_",j,".csv"))
						x1 = tab[,"startLon"]; y1 = tab[,"startLat"]; x2 = tab[,"endLon"]; y2 = tab[,"endLat"]
						durations4 = 4*tab[,"length"]; distances = sqrt(((x2-x1)^2)+((y2-y1)^2)); squareDistances = distances^2 # --> Euclidean distances !!
			 			LM = lm(durations4 ~ squareDistances); results[j,"MDS_CS_R2_env"] = summary(LM)$r.squared
						tab = read.csv(paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_original_coordinates_ext/TreeExtractions_",j,".csv"))
						x1 = cbind(tab[,"startLon"],tab[,"startLat"]); x2 = cbind(tab[,"endLon"],tab[,"endLat"])
			 			durations4 = 4*tab[,"length"]; distances = diag(rdist.earth(x1, x2, miles=F, R=NULL)); squareDistances = distances^2 # --> great-circle distances !!
			 			LM = lm(durations4 ~ squareDistances); results[j,"CGT_R2_null"] = summary(LM)$r.squared
						tab = read.csv(paste0(simulationDirectory,"/Seraphim_analyses/TreeSimulations_",i,"_cartogram_transform_ext/TreeExtractions_",j,".csv"))
						x1 = tab[,"startLon"]; y1 = tab[,"startLat"]; x2 = tab[,"endLon"]; y2 = tab[,"endLat"]
						durations4 = 4*tab[,"length"]; distances = sqrt(((x2-x1)^2)+((y2-y1)^2)); squareDistances = distances^2 # --> Euclidean distances !!
			 			LM = lm(durations4 ~ squareDistances); results[j,"CGT_R2_env"] = summary(LM)$r.squared
					}
				hds = round(HDInterval::hdi(results[,"MDS_LC_R2_null"])[1:2],3); meanV = mean(results[,"MDS_LC_R2_null"])
				results_1b[i,"LC_R2_null"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				hds = round(HDInterval::hdi(results[,"MDS_LC_R2_env"])[1:2],3); meanV = mean(results[,"MDS_LC_R2_env"])
				results_1b[i,"LC_R2_env"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				vS = results[,"MDS_LC_R2_env"]-results[,"MDS_LC_R2_null"]; hds = round(HDInterval::hdi(vS)[1:2],3); meanV = mean(vS)
				results_1b[i,"LC_Q"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				p = sum(results[,"MDS_LC_R2_env"]>results[,"MDS_LC_R2_null"])/dim(results)[1]
				results_1b[i,"LC_pQ"] = p; results_1b[i,"LC_BF"] = round(p/(1-p),1)
				hds = round(HDInterval::hdi(results[,"MDS_CS_R2_null"])[1:2],3); meanV = mean(results[,"MDS_CS_R2_null"])
				results_2b[i,"CS_R2_null"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				hds = round(HDInterval::hdi(results[,"MDS_CS_R2_env"])[1:2],3); meanV = mean(results[,"MDS_CS_R2_env"])
				results_2b[i,"CS_R2_env"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				vS = results[,"MDS_CS_R2_env"]-results[,"MDS_CS_R2_null"]; hds = round(HDInterval::hdi(vS)[1:2],3); meanV = mean(vS)
				results_2b[i,"CS_Q"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				p = sum(results[,"MDS_CS_R2_env"]>results[,"MDS_CS_R2_null"])/dim(results)[1]
				results_2b[i,"CS_pQ"] = p; results_2b[i,"CS_BF"] = round(p/(1-p),1)
				hds = round(HDInterval::hdi(results[,"CGT_R2_null"])[1:2],3); meanV = mean(results[,"CGT_R2_null"])
				results_3b[i,"CGT_R2_null"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				hds = round(HDInterval::hdi(results[,"CGT_R2_env"])[1:2],3); meanV = mean(results[,"CGT_R2_env"])
				results_3b[i,"CGT_R2_env"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				vS = results[,"CGT_R2_env"]-results[,"CGT_R2_null"]; hds = round(HDInterval::hdi(vS)[1:2],3); meanV = mean(vS)
				results_3b[i,"CGT_Q"] = paste0(round(meanV,3)," [",hds[1],", ",hds[2],"]")
				p = sum(results[,"CGT_R2_env"]>results[,"CGT_R2_null"])/dim(results)[1]
				results_3b[i,"CGT_pQ"] = p; results_3b[i,"CGT_BF"] = round(p/(1-p),1)
			}
		write.table(results_1a, paste0(simulationDirectory,"/MDS_transf_LC1.csv"), row.names=F, quote=F, sep=";")
		write.table(results_2a, paste0(simulationDirectory,"/MDS_transf_CS1.csv"), row.names=F, quote=F, sep=";")
		write.table(results_3a, paste0(simulationDirectory,"/Cartog_transfo1.csv"), row.names=F, quote=F, sep=";")
		results_1a = cleanResultsTable(read.csv(paste0(simulationDirectory,"/MDS_transf_LC1.csv"), head=T, sep=";"))
		write.table(results_1a, paste0(simulationDirectory,"/MDS_transf_LC1.txt"), row.names=F, quote=F, sep="\t")
		results_2a = cleanResultsTable(read.csv(paste0(simulationDirectory,"/MDS_transf_CS1.csv"), head=T, sep=";"))
		write.table(results_2a, paste0(simulationDirectory,"/MDS_transf_CS1.txt"), row.names=F, quote=F, sep="\t")
		results_3a = cleanResultsTable(read.csv(paste0(simulationDirectory,"/Cartog_transfo1.csv"), head=T, sep=";"))
		write.table(results_3a, paste0(simulationDirectory,"/Cartog_transfo1.txt"), row.names=F, quote=F, sep="\t")
		vS1 = results_1a[,"LC_Q"]; vS2 = as.numeric(results_1a[,"LC_pQ"]); tmp = rep(NA, length(vS1))
		for (i in 1:length(vS1)) tmp[i] = as.numeric(unlist(strsplit(vS1[i]," \\["))[1])
		vS1 = tmp; ci1 = round(quantile(vS1,c(0.025,0.975)),3); ci2 = round(quantile(vS2,c(0.025,0.975)),2)
		BFs = as.numeric(gsub(">","",results_1a[,"LC_BF"])); BFs = BFs[!is.na(BFs)]
		cat("MDS, LC-1: Q = ",round(mean(vS1),3)," [",ci1[1],", ",ci1[2],"], p(Q>0) = ",round(mean(vS2),2)," [",ci2[1],", ",ci2[2],"], ",sep="")
		cat("BF>3: ",sum(BFs>=3),"/",length(BFs),", BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
		vS1 = results_2a[,"CS_Q"]; vS2 = as.numeric(results_2a[,"CS_pQ"]); tmp = rep(NA, length(vS1))
		for (i in 1:length(vS1)) tmp[i] = as.numeric(unlist(strsplit(vS1[i]," \\["))[1])
		vS1 = tmp; ci1 = round(quantile(vS1,c(0.025,0.975)),3); ci2 = round(quantile(vS2,c(0.025,0.975)),2)
		BFs = as.numeric(gsub(">","",results_2a[,"CS_BF"])); BFs = BFs[!is.na(BFs)]
		cat("MDS, CS-1: Q = ",round(mean(vS1),3)," [",ci1[1],", ",ci1[2],"], p(Q>0) = ",round(mean(vS2),2)," [",ci2[1],", ",ci2[2],"], ",sep="")
		cat("BF>3: ",sum(BFs>=3),"/",length(BFs),", BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
		vS1 = results_3a[,"CGT_Q"]; vS2 = as.numeric(results_3a[,"CGT_pQ"]); tmp = rep(NA, length(vS1))
		for (i in 1:length(vS1)) tmp[i] = as.numeric(unlist(strsplit(vS1[i]," \\["))[1])
		vS1 = tmp; ci1 = round(quantile(vS1,c(0.025,0.975)),3); ci2 = round(quantile(vS2,c(0.025,0.975)),2)
		BFs = as.numeric(gsub(">","",results_3a[,"CGT_BF"])); BFs = BFs[!is.na(BFs)]
		cat("CGT-1: Q = ",round(mean(vS1),3)," [",ci1[1],", ",ci1[2],"], p(Q>0) = ",round(mean(vS2),2)," [",ci2[1],", ",ci2[2],"], ",sep="")
		cat("BF>3: ",sum(BFs>=3),"/",length(BFs),", BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
		write.table(results_1b, paste0(simulationDirectory,"/MDS_transf_LC2.csv"), row.names=F, quote=F, sep=";")
		write.table(results_2b, paste0(simulationDirectory,"/MDS_transf_CS2.csv"), row.names=F, quote=F, sep=";")
		write.table(results_3b, paste0(simulationDirectory,"/Cartog_transfo2.csv"), row.names=F, quote=F, sep=";")
		results_1b = cleanResultsTable(read.csv(paste0(simulationDirectory,"/MDS_transf_LC2.csv"), head=T, sep=";"))
		write.table(results_1b, paste0(simulationDirectory,"/MDS_transf_LC2.txt"), row.names=F, quote=F, sep="\t")
		results_2b = cleanResultsTable(read.csv(paste0(simulationDirectory,"/MDS_transf_CS2.csv"), head=T, sep=";"))
		write.table(results_2b, paste0(simulationDirectory,"/MDS_transf_CS2.txt"), row.names=F, quote=F, sep="\t")
		results_3b = cleanResultsTable(read.csv(paste0(simulationDirectory,"/Cartog_transfo2.csv"), head=T, sep=";"))
		write.table(results_3b, paste0(simulationDirectory,"/Cartog_transfo2.txt"), row.names=F, quote=F, sep="\t")
		vS1 = results_1b[,"LC_Q"]; vS2 = as.numeric(results_1b[,"LC_pQ"]); tmp = rep(NA, length(vS1))
		for (i in 1:length(vS1)) tmp[i] = as.numeric(unlist(strsplit(vS1[i]," \\["))[1])
		vS1 = tmp; ci1 = round(quantile(vS1,c(0.025,0.975)),3); ci2 = round(quantile(vS2,c(0.025,0.975)),2)
		BFs = as.numeric(gsub(">","",results_1b[,"LC_BF"])); BFs = BFs[!is.na(BFs)]
		cat("MDS, LC-2: Q = ",round(mean(vS1),3)," [",ci1[1],", ",ci1[2],"], p(Q>0) = ",round(mean(vS2),2)," [",ci2[1],", ",ci2[2],"], ",sep="")
		cat("BF>3: ",sum(BFs>=3),"/",length(BFs),", BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
		vS1 = results_2b[,"CS_Q"]; vS2 = as.numeric(results_2b[,"CS_pQ"]); tmp = rep(NA, length(vS1))
		for (i in 1:length(vS1)) tmp[i] = as.numeric(unlist(strsplit(vS1[i]," \\["))[1])
		vS1 = tmp; ci1 = round(quantile(vS1,c(0.025,0.975)),3); ci2 = round(quantile(vS2,c(0.025,0.975)),2)
		BFs = as.numeric(gsub(">","",results_2b[,"CS_BF"])); BFs = BFs[!is.na(BFs)]
		cat("MDS, CS-2: Q = ",round(mean(vS1),3)," [",ci1[1],", ",ci1[2],"], p(Q>0) = ",round(mean(vS2),2)," [",ci2[1],", ",ci2[2],"], ",sep="")
		cat("BF>3: ",sum(BFs>=3),"/",length(BFs),", BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
		vS1 = results_3b[,"CGT_Q"]; vS2 = as.numeric(results_3b[,"CGT_pQ"]); tmp = rep(NA, length(vS1))
		for (i in 1:length(vS1)) tmp[i] = as.numeric(unlist(strsplit(vS1[i]," \\["))[1])
		vS1 = tmp; ci1 = round(quantile(vS1,c(0.025,0.975)),3); ci2 = round(quantile(vS2,c(0.025,0.975)),2)
		BFs = as.numeric(gsub(">","",results_3b[,"CGT_BF"])); BFs = BFs[!is.na(BFs)]
		cat("CGT-2: Q = ",round(mean(vS1),3)," [",ci1[1],", ",ci1[2],"], p(Q>0) = ",round(mean(vS2),2)," [",ci2[1],", ",ci2[2],"], ",sep="")
		cat("BF>3: ",sum(BFs>=3),"/",length(BFs),", BF>20: ",sum(BFs>=20),"/",length(BFs),"\n",sep="")
	}	

