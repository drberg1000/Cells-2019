#!/usr/bin/env Rscript
#Make sure the "tiff" library is installed:
#  apt-get install libtiff5-dev
#  Rscript - <<< "install.packages('tiff'  ,,'http://www.rforge.net/'    )"
# This is loaded for image saving... May not be needed
#  Rscript - <<< "install.packages('sourcefile',repos=NULL, type='source'))"
#  Server seems to be having problems as the following command can't find the directory
#  Rscript - <<< "install.packages('deSolve'  ,,'http://desolve.r-forge.r-project.org'    )"


inColor = TRUE
plotDims  = 1080

popFile = "cat 2DTess_10000.out"   #set to "" to run simulation

if( inColor == TRUE ) {
	Colors = c("black",             #Empty
		   "blue",              #Normal
		   "red",               #Cancer
		   "green")             #Infect
	          #"white"              #NoNode
} else {
	Colors = c(rgb(0.00,0.00,0.00), #Empty 
		   rgb(0.33,0.33,0.33), #Normal
		   rgb(0.67,0.67,0.67), #Cancer
		   rgb(1.00,1.00,1.00)) #Infect
}

########################################
# Length of Simulation
########################################
time = 10000
populStep = 10
spaceStep = 50
# Used if plots are requested but simulation isn't run.
spaceSteps = c( seq(0, 250, 50) )
########################################
########################################


########################################
# DIMENSIONALITY
########################################
nDims = 3
if ( nDims == 2 ) {
	Slices = (1)
} else if ( nDims == 3 ) {
	zCompression=1
	Slices = c(0,25,50)
} else {
	stop()
}

########################################
########################################


########################################
# For Tessellation
########################################
########################################
Network_Type = "GRID"
if( Network_Type == "TESS" ){
	numNodes = 1000000
	base = paste(" 2d_vor_1000000.")  # Keep space at the beginning of the filename
} else if (Network_Type == "GRID" ) {
	width  = 100
	height = 100
	depth  = 100
} else {
	stop()
}
########################################
########################################


#######################################
# For All Models & Network Types
########################################
########################################
Infection_Type = "CENTER"
#Rates for 3-Population Equilibrium
Normal_Growth = 0.5
Cancer_Growth = 1.0
Infect_Growth = 1.2
Normal_Death  = 0.2
Cancer_Death  = 0.1
Infect_Death  = 0.1
#Enter Initial Population Counts Here
Normal = 0.9
Cancer = 0.09
Infect = 0.01
Total  = Normal + Cancer + Infect
########################################
########################################


if( Network_Type == "GRID" ){
	numNodes = width*height*depth
	print( paste("Create", width, height, depth, "Network? (Y/N)") )
	ans <- scan("stdin", character(), n=1)
	if( ans == 'Y' || ans == 'y' ){
		command = paste("../../Networks/gen_grid", width, height, depth)
		print(command)
		system(command)
	}
	if( depth == 1 ){
		base = paste0(" ", width, "x", height, "x1_nowrap.")
	}
	else if ( depth > 1 ){
		base = paste0(" ", width, "x", height, "x", depth, "_nowrap.")
	}
	else{
		print(paste("Invalid depth:", depth))
		stop()
	}

}

  datFILE = paste0( base, "dat")
  netFILE = paste0( base, "net")
coordFILE = paste0( base, "cor")


Cancer_Count = (Cancer+Infect)*numNodes
Percent_Infected = Infect / (Cancer + Infect) * 100



basecommand = paste("./Cells",
	" --NormalGrowth=",    Normal_Growth, " --NormalDeath=",     Normal_Death,
	" --CancerGrowth=",    Cancer_Growth, " --CancerDeath=",     Cancer_Death,
	" --ViralIntGrowth=",  Infect_Growth, " --ViralIntDeath=",   Infect_Death,
	" --ViralExtGrowth=",  Infect_Growth, " --ViralExtDeath=",   Infect_Death,
	" --RandomSeed=",      42,
#	" --Equalize",
#	" --CancerRadius=",    NUM,
	" --CancerCount=",     Cancer_Count,
	" --TimeInfect=",      0,
	" --PercentInfected=", Percent_Infected,
	" --InfectionType=",   Infection_Type,
#	" --PercentInterior=", NUM,
	" --TimeMax=",         time,
#	" --Grid",
#	" --Summary",
#	" --TumorSizes",
#	" --Watch",
#	" --Files",
#	" --File=", BASENAME,
  datFILE, netFILE, coordFILE, 
  sep="" )

print("Generate Simulator Spatial Plot & Files? (Y/N)")
ans <- scan("stdin", character(), n=1)
if( ans == 'Y' || ans == 'y' )
{
	print("Run Simulation? (Y/N)")
	ans <- scan("stdin", character(), n=1)
	if( ans == 'Y' || ans == 'y' )
	{
		spaceSteps = seq(0, time, spaceStep)
		command= paste(basecommand, " --OutputInterval=",  spaceStep, " --File", sep="")
		print( command )
		system( command )
	}
	print("Include Legend? (Y/N)")
	ans <- scan("stdin", character(), n=1)
	if( ans == 'Y' || ans == 'y' ) {
		IncludeLegend=TRUE
	} else { 
		IncludeLegend=FALSE
	}

	for( It in spaceSteps ){
		filename = paste("Lattice_State_", It, ".0.csv", sep="")
		data<- read.table(filename, sep=" ", nrows=numNodes)

		if( nDims == 3 ) {
			data<-data[,c(1,2,3,5)]
			if( Network_Type == "TESS" ){
				data[,3] = data[ ,3]/zCompression
				data[,1:3] <- floor( data[,1:3] )
			}
			names(data)<-c("X","Y","Z","Type")
			if( Network_Type == "TESS" ) {
				data <- data[ with(data[,], order(X, Y, Z)), ]
			}

		} else {
			data<-data[,c(1,2,4)]
			if( Network_Type == "TESS" ){
				data[,c(1,2)] = data[ , c(1,2)]/10
			}
			names(data)<-c("X","Y","Type")
			if( Network_Type == "TESS" ) {
				data <- data[ with(data[,], order(X, Y)), ]
			}
		}
		data$Color               = Colors[1]
		data$Color[data$Type==1] = Colors[2]
		data$Color[data$Type==2] = Colors[3]
		data$Color[data$Type==3] = Colors[4]



		for( zIt in Slices ) {
			if( nDims == 3 ){
				if( inColor == TRUE ) {
					slicefile = sprintf("3D%04.2f_S%d.png",    It/1000, zIt)
				} else {
					slicefile = sprintf("3D%04.2f_BW_S%d.png", It/1000, zIt)
				}
			} else {
				if( inColor == TRUE ) {
					slicefile = sprintf("2D%04.2f.png",    It/1000)
				} else {
					slicefile = sprintf("2D%04.2f_BW.png", It/1000)
				}
			}

                        if( nDims == 3 ){
				SliceX <- data[data$Z==zIt,1] 
				SliceY <- data[data$Z==zIt,2] 
				SliceT <- data[data$Z==zIt,5] 
			} else { #( nDims == 2 )
				SliceX <- data[,1] 
				SliceY <- data[,2] 
				SliceT <- data[,4] 
			}
			png(filename=slicefile, width=plotDims, height=plotDims, units="px" )
			plot(    SliceX, 
			         SliceY,
		             col=SliceT, 
			      bg=SliceT,

			      cex=2,
			      pch=22,
			      main=paste("Time =", It, "Slice =", zIt), 
			      xlab="", 
			      ylab="" )

			if( IncludeLegend ){
				legend("topright", c("Empty", "Normal", "Cancer", "Infected"), 
				       col = Colors,
				       bg="white",
				       fill= Colors)
			}
			dev.off()
		}
			
	}
}



print("Generate Simulator Population Plot? (Y/N)")
ans <- scan("stdin", character(), n=1)
if( ans == 'Y' || ans == 'y' ) {
	print("Include Legend? (Y/N)")
	ans <- scan("stdin", character(), n=1)
	if( ans == 'Y' || ans == 'y' ) {
		IncludeLegend=TRUE
	} else { 
		IncludeLegend=FALSE
	}

	print("Running")
	command= paste(basecommand, " --OutputInterval=",  populStep, " --Stats | tee 3D_Pop_", time, ".out", sep="")
	if( popFile != "" ){
		command = popFile  
	}
	print(command)

	X <- system( command, intern=TRUE )
#        X <- lapply(X, substring, c(8, 26, 59, 91), c(17, 43, 76, 120 )) #Doesn't work w/ 1000x1000
        X <- lapply(X, substring, c(8, 27, 60, 92), c(17, 43, 76, 120 ))
	X <- lapply(X, as.numeric)
	X <- do.call(rbind, X)

	png(filename=paste("2DSim_Pop_", time, ".png", sep=""), width=plotDims, height=plotDims, units="px")
	matplot(X[,1], X[,2:4], 
		type = "l",
	       	xlab = "Time", 
		ylab = "Percent of Population", 
	        ylim=c(0,0.8),
		main = "Simulator", 
		col = c("blue", "red", "green"), 
		lwd = 3)

	if ( IncludeLegend ) {
		legend("topright", c("Normal", "Cancer", "Infected"), 
		       col = c("blue", "red", "green"), 
		       bg="white",
		       lty = 1:3)
	}

	dev.off()
}




print("MeanField Plot? (Y/N)")
ans <- scan("stdin", character(), n=1)
if( ans == 'Y' || ans == 'y' )
{
	print("Include Legend? (Y/N)")
	ans <- scan("stdin", character(), n=1)
	if( ans == 'Y' || ans == 'y' ) {
		IncludeLegend=TRUE
	} else { 
		IncludeLegend=FALSE
	}
	library(deSolve)
	model <- function (time, state, pars) {
		with(as.list(c(state, pars)), {
		     ###    Normal/Cancer Grow # Death # Viral Infection
		     dUn <- Ln*Un*(1-Un-Uc-Uv) - Dn*Un
		     dUc <- Lc*Uc*(1-Un-Uc-Uv) - Dc*Uc - Lv*Uc*Uv
		     dUi <-                    - Dv*Uv + Lv*Uc*Uv
		     return(list(c(dUn, dUc, dUi)))
		})
	}

	pars <- c( Ln=Normal_Growth, Lc=Cancer_Growth, Lv=Infect_Growth,
		   Dn=Normal_Death,  Dc=Cancer_Death,  Dv=Infect_Death)
	yini <- c( Un=Normal/Total,  Uc=Cancer/Total,  Uv=Infect/Total) 

	times <- seq(0, time, by = populStep)
	out <- ode(func = model, y = yini, parms = pars, times = times)

	# TODO Change filename as appropriate
	png(filename=paste("2DMF_Pop_", time, ".png", sep=""), width=plotDims, height=plotDims, units="px")
	matplot(out[,"time"], out[,2:4], type = "l", xlab = "Time", ylab = "Percent of Population", 
		main = "MeanField", col = c("blue", "red", "green"),lwd = 3)
	if( IncludeLegend == TRUE ) {
		legend("topright", c("Normal", "Cancer", "Infected"), col = c("blue", "red", "green"), lty = 1:3)
	}
	dev.off()

}





