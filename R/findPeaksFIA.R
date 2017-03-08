#Functions to detect peaks on a single FIA-HRMS acquisition.
#The entry is an xcmsRaw object, and the ouptput is a data matrix.
#This is done by band detection and then Continuous Wavelets transform to find the value.

#Package used : pracma, xcms

#XCMS functions used : rawEIC, xcmsRaw


### Function to detect peaks in an FIA acquisition.
AntiSymTriangle <- function(x) {
	if (-1 <= x & x < 0) {
		return(-x)
	}
	if (-2 <= x & x < (-1)) {
		return(x + 2)
	}
	if (0 <= x & x < 1) {
		return(-x)
	}
	if (1 <= x & x < 2) {
		return(x - 2)
	}
	return(0)
	
}


calcAngle <- function(x, y, p1, p2, p3) {
	#cat("in")
	#cat("x",x,"\ny",y,"p",p1,"]",p2,"]",p3,"\n")
	x <- x / max(x)
	y <- y / max(y)
	a <- sqrt((x[p2] - x[p1]) ^ 2 + (y[p2] - y[p1]) ^ 2)
	b <- sqrt((x[p2] - x[p3]) ^ 2 + (y[p2] - y[p3]) ^ 2)
	c <- sqrt((x[p1] - x[p3]) ^ 2 + (y[p1] - y[p3]) ^ 2)
	val_cos_1 <- (a ^ 2 + b ^ 2 - c ^ 2) / (2 * a * b)
	valcos <- suppressWarnings(acos(val_cos_1))
	if (is.nan(valcos))
		return(pi)
	valcos
}


#Taken from the pracma package.
trapezArea <- function (x, y)
{
	if (missing(y)) {
		if (length(x) == 0)
			return(0)
		y <- x
		x <- seq(along = x)
	}
	if (length(x) == 0 && length(y) == 0)
		return(0)
	if (!(is.numeric(x) || is.complex(x)) || !(is.numeric(y) ||
											   is.complex(y)))
		stop("Arguments 'x' and 'y' must be real or complex vectors.")
	m <- length(x)
	if (length(y) != m)
		stop("Arguments 'x', 'y' must be vectors of the same length.")
	if (m <= 1)
		return(0)
	xp <- c(x, x[m:1])
	yp <- c(numeric(m), y[m:1])
	n <- 2 * m
	p1 <- sum(xp[1:(n - 1)] * yp[2:n]) + xp[n] * yp[1]
	p2 <- sum(xp[2:n] * yp[1:(n - 1)]) + xp[1] * yp[n]
	return(0.5 * (p1 - p2))
}

checkIsoValues <- function(x) {
	if (length(x) < 3) {
		return(all(x == 0))
	}
	##Binarizing the data
	shift_right <- x[-1] == 0
	shift_left <- x[-length(x)] == 0
	if (all(shift_right | shift_left))
		return(TRUE)
	return(FALSE)
}


fuseRange <- function(r1, r2, extend = TRUE) {
	if (extend)
		return(c(min(r1[1], r2[1]), max(r1[2], r2[2])))
	return(c(max(r1[1], r2[1]), min(r1[2], r2[2])))
}


#' Detect peaks in an FIA acquisition.
#'
#' Detect the peak corresponding to compounds present in the sample
#' in a Flow Injection Analysis (FIA) acquisition. The item provided
#' must be an xcmsRaw object.
#'
#' @export
#' @param xraw An xcmsRaw object as returned by \code{\link[xcms]{xcmsRaw}}.
#' @param ppm The authorized deviation between scans in ppm, this parameter
#' will also be used to fuse the bands if there are close enough.
#' @param dmz The minimum absolute value of the deviation between scans, to take into account
#' the higher diviations at low masses.
#' @param es A noise estimation object as returned by \link{estimateNoiseMS}, or NULL
#' if the parameter noise if only an threshold is supposed to be used.
#' @param  solvar Should the signal corresponding to solvent be kept ?
#' Only their maximum intensiyt will be calculated.
#' @param  solint How are the intensity of signal with bot solvent and sample
#' be treated in the injection zone region, the area of the rectangle
#' with peak-width and solvent intensity is considered.
#' @param fullInteg If no solvent is detected before the injection,
#' should the signals be integrated on the ufll chromatogram.
#' \itemize{
#'     \item poly. Half of this area is kept.
#'     \item substract. The area is removed substracted.
#'     \item remove. The area is conserved in the final value.
#' }
#' @param graphical A boolean indicating if the detected area shall be plotted.
#' @param SNT NULL by default a relative intensity of signal/intensity of solvent threshold,
#' used only if es is equal to NULL.
#' @param f method to design the filter, "TIC" means that the peak of the TIC is
#' used as a filter. "regression" means that the signal is regressed form the most
#' intense band as an Exponential modified gaussian.
#' @param pvalthresh The threshold used in p-value to discard signal, only used if
#' a noise model is provided.
#' @param scanmin The first scan to consider.
#' @param scanmax The last scan to consider.
#' @param fracPoint A filter on the number of point found in a band. 0.3 by default to allow matrix
#' effect.
#' @param sizeMin The minimum number of point considered for a band to be considred for solvent filtration.
#' @param bandlist An optional bandlist to be passed.
#' @param ... more arguments to be passed to the \link{determiningInjectionZone} function.
#' @return A numeric matrix with the following column
#' \itemize{
#'     \item mzmin the minimum value of the mass traces in the m/z dimension.
#'     \item mzmax the maximum value of the mass traces in the m/z dimension.
#'     \item scanMin the first scan on which the signal is detected.
#'     \item scanMax the last scan on which the signal is detected.
#'     \item areaIntensity the integrated area of the signal.
#'     \item maxIntensity the maximum intensity of the signal.
#'     \item solventIntensity the intensity of the solvent, 0 means that no significant
#'     solvent was detected.
#'     \item corPeak An idicator of matrix effect, if it's close to 1, the compound
#'     does not suffer from heavy matrix effect, if it is inferior to 0.5, the compound
#'     suffer from heavy matrix effect.
#'     \item timeShifted An indicator of the shifting of th epeak in the time direction.
#'     \item signalOverSolventRatio The ratio of the signal max intensity on the solvent max intensity.
#'     }
#' @aliases findFIASignal findPeaks
#' @examples
#' if(require(plasFIA)){
#'   #Getting the path of a file.
#'   path_raw<-list.files(system.file(package="plasFIA","mzML"),full.names=TRUE)[2]
#'
#'   #Opening the file with xcms
#'   xraw<-xcmsRaw(path_raw)
#'
#'   ppm<-2
#'
#'   #getting the filtered signals without noise model which is not recommended.
#'   tsignal<-findFIASignal(xraw,ppm=ppm,SNT=3)
#'
#'   #Getting the noise model un the plasSet object.
#'   data(plasSet)
#'   es<-attr(plasSet,"noiseEstimation")
#'
#'   #Getting the signal with a noise model.
#'   tsignal<-findFIASignal(xraw,ppm=2,es=es,pvalthresh=0.005)
#'   head(tsignal)
#' }
#' @useDynLib proFIA
findFIASignal <-
	function(xraw,
			 ppm,
			 es = NULL,
			 solvar = c("throw", "keep"),
			 solint = c("poly", "substract", "add"),
			 fullInteg = FALSE,
			 dmz =
			 	0.0005,
			 graphical = FALSE,
			 SNT = NULL,
			 f = c("regression", "TIC"),
			 pvalthresh = NULL,
			 scanmin = 1,
			 scanmax = length(xraw@scantime),
			 fracPoint = 0.3,
			 sizeMin = NULL,
			 bandlist = NULL,
			 ...) {
		if(class(xraw)=="character"){
			xraw <- xcmsRaw(xraw)
		}else if(class(xraw)!="xcmsRaw"){
			stop("Character giving path to a file or xcmsRaw object required for xraw argument.")
		}
		solint <- match.arg(solint)
		solvar <- match.arg(solvar)
		f <- match.arg(f)
		n <- length(xraw@scantime)
		sizepeak <-
			determiningInjectionZone(xraw, scanmin = scanmin, scanmax = scanmax, ...)
		sizepeak <- sizepeak + scanmin - 1
		if (length(sizepeak) < 3) {
			warning(paste(
				"No injection peak has been detected in the acquisition",
				basename(xraw@filepath)
			))
		}
		psip <- 1
		headPeaks <-
			c(
				"mzmin",
				"mzmax",
				"mz",
				"scanMin",
				"scanMax",
				"areaIntensity",
				"maxIntensity",
				"solventIntensity",
				"corSampPeak",
				"timeShifted",
				"signalOverSolventRatio"
			)
		
		###Which kind of quality control will be used.
		QC <- "Reduced"
		if (is.null(es) & is.null(SNT)) {
			stop("No noise model and no SNT provided.")
		} else if (is.null(es) & is.numeric(SNT)) {
			QC <- "Snt"
			message("No noise model provided, a signal/solvant threshold will be used.")
		} else if (!is.null(es@estimation)) {
			if (is.null(pvalthresh)) {
				message(
					"A noise model is provided but no threshold was specified, the threshold is therefore set to 0.01 by default"
				)
				pvalthresh <- 0.01
			} else{
				message("A noise model and a pvalue threshold are provided, SNT will be ignored.")
			}
			QC <- "Nes"
		} else{
			stop("Wrong type of noise or SNT, check the parameters.")
		}
		if (QC == "Nes") {
			headPeaks <- c(headPeaks, "signalOverSolventPvalue")
		}
		modelPeak <- NULL
		if(is.null(sizeMin)){
			sizeMin <- floor((sizepeak[3] -sizepeak[1]) * 0.5)
		}
		
		###Bands are calculated.
		if(is.null(bandlist)){
		bandlist <- findBandsFIA(
			xraw,
			ppm = ppm,
			sizeMin = sizeMin,
			dmz =
				dmz,
			beginning = sizepeak[1],
			end=sizepeak[2],
			firstScan = scanmin,
			lastScan = scanmax,
			fracMin = fracPoint
		)
		}
		pmmin = sizepeak[1]
		pmmax = sizepeak[2]
		if (f == "TIC") {
			TICv <- rawEIC(xraw, mzrange = range(xraw@env$mz))$intensity
			
			### A model peak is given by the detected TIC peak.
			modelPeak <- TICv[sizepeak[1]:sizepeak[2]]
			modelPeak <- modelPeak - modelPeak[1]
			modelPeak <- smooth.BWF(modelPeak, freq = 0.2)
			modelPeak <- modelPeak / max(modelPeak)
		} else if (f == "regression") {
			modelPeak <- getInjectionPeak(
				xraw,
				bandlist,
				sec = 2,
				iquant = 0.95,
				gpeak = sizepeak,
				graphical = graphical
			)
			#cat("klkl");print(modelPeak);cat("kllkkl")
			###Beginning and end correspond to the position where th epeak int is superior at 5%
			mpPmax <- which.max(modelPeak)
			pmmin <-
				which(modelPeak[1:mpPmax] > 0.05 * modelPeak[mpPmax])
			pmmin <- pmmin[1]
			pmmax <-
				which(modelPeak[(mpPmax + 1):length(modelPeak)] > 0.05 *
					  	modelPeak[mpPmax]) + mpPmax - 1
			pmmax <- pmmax[length(pmmax)]
			modelPeak <- modelPeak[pmmin:pmmax]
			modelPeak <- modelPeak / max(modelPeak)
			sizepeak[2] <- max(sizepeak[2], pmmax)
			
		}
		
		###defining the injection peak model.
		xseq <- seq(-2, 2, length = 1000)
		yTriangl <- sapply(xseq, AntiSymTriangle, simplify = TRUE)
		sizeInc <- floor((sizepeak[3] - sizepeak[1]))
		seqIndex <- floor(seq(0, 1000, length = sizeInc))
		modelInjPeak <- yTriangl[seqIndex]
		if (sizeInc < 5)
			modelInjPeak <- yTriangl[c(0, 250, 750, 1000)]
		nInj <- length(modelInjPeak)
		
		pmaxInjPeak <- which.max(modelInjPeak)
		pleftInjPeaks <- floor(pmaxInjPeak / 2)
		prightInjPeaks <- pmaxInjPeak + floor(pmaxInjPeak / 2)
		
		##Constant which will be used on all iteration.
		nf <- length(modelPeak)
		ninj <- length(modelInjPeak)
		countSol <- 0
		adjustZone <- floor((sizepeak[3] - sizepeak[1]) / 6)
		
		#Important positions on the iflter.
		pmaxModelPeak <- which.max(modelPeak)
		
		#Position ot start looking for the limit.
		pleftModelPeaks <- floor(pmaxModelPeak / 2)
		prightModelPeaks <-
			pmaxModelPeak + floor((nf - pmaxModelPeak) * 0.75)
		
		
		###list which will stock the result.
		ResList <- list()
		if(graphical){
			plot.new()
			title(main=paste("Peak and filters vizualisation : ",getRawName(xraw@filepath)))
			legend("center",legend = c("Raw EIC","Smoothed Eic","Integration limit","Injection peak filter coef","Triangular filter Coef (when used)"),
				   col=c("black","darkgreen","red","orange","purple"),lwd=2,lty=c(1,1,2,1,1))
		}
		
		message("Band filtering: ", appendLF = FALSE)
		memmess <- 1
		for (i in 1:nrow(bandlist)) {
			bshifted <- 0
			vact <- floor(i / nrow(bandlist) * 100) %/% 10
			if (vact != memmess) {
				memmess <- vact
				message(memmess * 10, " ", appendLF = FALSE)
			}
			vEic <-
				rawEIC(xraw, mzrange = c(bandlist[i, "mzmin"], bandlist[i, "mzmax"]))
			#Extending the sequence to at least two times the size of the
			#injection filter
			# if (length(vEic$intensity) / length(modelPeak) < 2) {
			# 	nlength <- nextn(2 * length(modelPeak), 2)
			# 	valrep <-
			# 		mean(vEic$intensity[(length(vEic$intensity) - 3):length(vEic$intensity)])
			# 	vEic$intensity <- c(vEic$intensity,
			# 						rep(valrep,
			# 							times = nlength - length(vEic$intensity)))
			# }
			
			
			
			###Determining if there is solvant in the signal.
			#seqSol = vEic$intensity[1:sizepeak[1]]
			###CHecking if there is only isolated values, 0 means no reliably detectable solvant.
			Bsol <-
				ifelse(checkIsoValues(vEic$intensity[1:sizepeak[1]]), FALSE, TRUE)
			valSol <- NULL
			if (Bsol) {
				valOkSol <- which(vEic$intensity[1:sizepeak[1]] != 0)
				valSol <- mean(vEic$intensity[valOkSol])
			} else{
				valSol <- 0
			}
			##Calculating the correlation.
			mEic <- max(vEic$intensity)
			pos0 <-
				which(vEic$intensity[sizepeak[1]:sizepeak[2]] == 0) + sizepeak[1] -
				1
			valCor <- NULL
			cinter <- pmmin:pmmax
			if (length(pos0) == 0) {
				valCor <- cor(modelPeak, vEic$intensity[cinter])
			} else{
				if ((sizepeak[2] - sizepeak[1] + 1 - length(pos0)) > 1) {
					valCor <-
						suppressWarnings(cor(modelPeak[-pos0], vEic$intensity[cinter][-pos0]))
				} else{
					valCor <- 0
				}
			}
			###Problematic case in really noisy peaks.
			if (is.na(valCor)) {
				valCor <- 0
			}
			
			##Calculating the convolution
			convolvedSeqF <-
				convolve(vEic$intensity, modelPeak, type = "filter")
			convolvedSeqInj <-
				convolve(vEic$intensity, modelInjPeak, type = "filter")
			smoothedSeq <- smooth.BWF(vEic$intensity, freq = 0.2)
			
			###Getting the local maximum
			#pMf is the position of the maximum on the convolved sequence.
			scminF <- max(1, scanmin - floor(nf / 2))
			scmaxF <- max(1, scanmax - floor(nf / 2))
			
			pMf <- which.max(convolvedSeqF[scminF:scmaxF]) + scminF - 1
			
			
			##PosMaw is the position of the maximum of the filter on the EIC.
			posMax <- pMf + pmaxModelPeak
			
			###Peaklim store the found limit of the peak.
			peaklim <- NULL
			
			###SecodMaw is the position of the second maximum on th epeak
			SecondMax <- NULL
			FirstMax <- FALSE
			PeakFound <- FALSE
			extendedFilter <- FALSE
			###If the maximum detected is in the injection peak.
			if (posMax < sizepeak[2] && posMax > sizepeak[1]) {
				###First estimation of the peak limit. peaklim will store the peak limit all along the process.
				#peaklim <- c(posMax - (floor(nf / 2)), posMax + (floor(nf / 2)))
				##TEST1
				peaklim <-
					c(pMf + pleftModelPeaks, pMf + prightModelPeaks)
				
				###Locating the limit of the peak found
				pl <-
					findPeaksLimits(smoothedSeq, peaklim[1] + 1, peaklim[2] - 1)
				
				#Once the limit of the peak is located, the peak limit is found.
				#peaklim[2] <- pl[2]
				peaklim <- fuseRange(peaklim, pl)
				###The filter coeff goes down to find if there is a second peak.
				# 				posRawMax <- which.max(convolvedSeqF[max((pMf - adjustZone), 1):(pMf +
				#                                                                                        adjustZone)]) + max(((pMf - adjustZone) - 1),1)
				#We find the limit on the convolved sequence.
				flim <-
					findPeaksLimits(convolvedSeqF, pMf - 1, pMf + 1)
				second_filter <- FALSE
				###Case where the max is the right of the injection peak.
				if (posMax < sizepeak[2] &
					posMax > sizepeak[3]) {
					###We check the second part of the interval.
					# SecondInter <- c(sizepeak[1] - floor(nf / 2),
					#                 min(flim[1] - floor(nf /
					#                                         2), pl[1]))
					SecondInter <-
						c(sizepeak[1] - floor(nf / 2), flim[1])
					
					#Necessary as flim may return -1 if out of sequence.
					SecondInter[1] = max(SecondInter[1], 0)
					
					#Checking hat there is an inter to consider.
					if (SecondInter[2] - SecondInter[1] > 1) {
						Secondpmf <-
							which.max(convolvedSeqF[SecondInter[1]:SecondInter[2]]) +
							SecondInter[1] - 1
						
						#The second maximum on the filter.
						SecondMax <- Secondpmf + pmaxModelPeak
						
						###CHecking that the second max is in the left part of the injection peak.
						if (SecondMax > sizepeak[1] &
							SecondMax < sizepeak[3]) {
							###The limit of the second peak is localised in a decent windows.
							pl <- findPeaksLimits(smoothedSeq,
												  SecondMax - 1,
												  SecondMax + 1)
							peaklim[1] <- pl[1]
							extendedFilter <- TRUE
						} else {
							### Checking hte injection filter in the good direction.
							SecondInter <-
								c(sizepeak[1] - floor(ninj / 2),
								  pl[1] - floor(ninj / 2))
							
							##Necessary because both filter have different size.
							SecondInter[1] <- max(SecondInter[1], 0)
							if ((SecondInter[2] - SecondInter[1]) > 0) {
								###Checking if we are not too close from the beginning.
								SecondMax <-
									which.max(convolvedSeqInj[SecondInter[1]:SecondInter[2]]) +
									SecondInter[1] - 1 + pmaxInjPeak
								Secondpmf <-
									SecondMax + pmaxInjPeak
								second_filter <- TRUE
								if (Secondpmf > sizepeak[1] &
									Secondpmf < pl[1]) {
									rawSecondPmf <-
										which.max(smoothedSeq[(SecondMax - adjustZone):(SecondMax +
																							adjustZone)]) + (SecondMax - adjustZone) - 1
									pl <-
										findPeaksLimits(smoothedSeq,
														rawSecondPmf - 2,
														rawSecondPmf + 2)
									extendedFilter <- TRUE
									peaklim[1] <- pl[1]
								}
							}
						}
					}
				}
				if (posMax >= sizepeak[1] &
					###Case where the detected mas is at the left of the peak.
					posMax < sizepeak[3]) {
					#In this case the injection peak is always the limit
					peaklim[1] <- sizepeak[1]
					
					SecondInter <- c(flim[2], sizepeak[2] - floor(nf / 2))
					
					#Checking that there is a second maximum in th second direction.
					SecondpMf <-
						which.max(convolvedSeqF[SecondInter[1]:SecondInter[2]]) +
						SecondInter[1] - 1
					SecondMax <- SecondpMf + pmaxModelPeak
					
					###CHecking that the second max is in the left part of the injection peak.
					if (SecondMax > sizepeak[3] &
						SecondMax <= sizepeak[2]) {
						###The limit of the second peak is localised in a decent windows.
						rawSecondPmf <-
							which.max(smoothedSeq[(SecondMax - adjustZone):(SecondMax +
																				adjustZone)]) + (SecondMax - adjustZone) - 1
						pl <- findPeaksLimits(smoothedSeq,
											  rawSecondPmf - 1,
											  rawSecondPmf + 1)
						if(pl[2]>peaklim[2]){
							extendedFilter <- TRUE
							peaklim[2] <- pl[2]
						}
					}
				}
			}
			if (posMax > sizepeak[2]) {
				##Checking that there is some retention in  the colmun.
				##If it the case the peak is wider than usual, so wider than the
				##modelPeak.
				#peaklim = floor(nf / 2)+c(pMf - (floor(nf / 2)), pMf + (floor(nf / 2)))
				###Locating the limit of the peak found
				bshifted <- 1
				pl <-
					findPeaksLimits(convolvedSeqF, pMf - 3, pMf + 3)
				peaklim <- pl + floor(nf / 2)
				#cat("mumu",posMax,"    ")
				if (diff(peaklim) < length(modelPeak))
					peaklim = NULL
			}
			
			###Quality control of the found peak.
			NoPeak <- is.null(peaklim)
			
			####Refinement of the right peak limit.
			if (!NoPeak) {
				if (!Bsol & fullInteg) {
					peaklim[2] <- length(xraw@scantime)
				} else{
					if (peaklim[2] < sizepeak[2]) {
						msol <- mean(vEic$intensity[sizepeak[2]:length(vEic$intensity)])
						sdsol <-
							sd(vEic$intensity[sizepeak[2]:length(vEic$intensity)])
						
						p_candidate <- peaklim[2]
						
						vtreshint <-msol + 2 * sdsol
						while (smoothedSeq[p_candidate]>vtreshint&smoothedSeq[p_candidate+1]>vtreshint ) {
							###We try to put
							p_candidate <- p_candidate + 1
							if(p_candidate == length(xraw@scantime)) break
							
						}
						###Three candidates are considered, the found point, the sipzekea, and orignal value.
						vcandidates <- c(peaklim[2],sizepeak[2],p_candidate)
						#cat("pl2",peaklim[2],"sp2",sizepeak[2],"pc",p_candidate,"\n")
						psup <- which(vEic$intensity[vcandidates]<vEic$intensity[posMax])
						vcandidates <- vcandidates[psup]
						vval <- sapply(vcandidates,calcAngle,
										x=xraw@scantime,y=vEic$intensity,p1=posMax,p3=length(xraw@scantime))
						###Possible adding a term to considered the area in the peak and the are outside.
						if(length(vval)!=0){
						peaklim[2] <- vcandidates[which.min(vval)]
						}
					}
					
				}
				###If there is no solvent we try to adjust the righ limit ot a 0.
				# if(!Bsol){
				# 	pp0 <- which(vEic$intensity[peaklim[2]:length(xraw@scantime)]==0)
				# 	while(vEic$intensity[peaklim[2]]!=0&peaklim[2]<length(vEic$intensity)){
				# 		peaklim[2] <- peaklim[2]+1
				# 	}
				# }
			}
			
			#Refinement of the left peak limit.
			if (!NoPeak) {
				peaklim[1] = max(peaklim[1], sizepeak[1])
				if (vEic$intensity[peaklim[1]] > valSol * 1.5)
					peaklim[1] <- sizepeak[1]
				if (peaklim[2] < peaklim[1]) {
					NoPeak <- TRUE
				}
			}
			pval <- NULL
			SNTval <- max(smoothedSeq)/ valSol
			if (NoPeak) {
				pval <- 1
				SNTval <- 0
			} else if (valSol == 0) {
				pval <- 0
				SNTval <- Inf
				###Checking that the final peak detected is not shifted.
				# if (peaklim[2] > sizepeak[2] &
				# 	abs(posMax - sizepeak[2]) < abs(posMax - sizepeak[3]) &
				# 	abs(peaklim[1] - sizepeak[3]) < abs(peaklim[1] - sizepeak[1])) {
				if ((mean(peaklim[1:2])>mean(sizepeak[1:2]))&(!extendedFilter)) {
					bshifted <- 1
				}
			} else if (QC == "Nes") {
				p1 <- max(1,peaklim[1]-1)
				p2 <- min(length(xraw@scantime),peaklim[2]+1)				
				pval <- calcPvalue(
					es,
					vEic$intensity[peaklim[1]:peaklim[2]],
					seq(smoothedSeq[p1], smoothedSeq[p2],
						length = peaklim[2] - peaklim[1] + 1)
				)
				if (pval > pvalthresh)
					next
				if (SNTval < 1.5)
					next
				###Checking that the final peak detected is not shifted.
				# if (peaklim[2] > sizepeak[2] &
				# 	abs(posMax - sizepeak[2]) < abs(posMax - sizepeak[3]) &
				# 	abs(peaklim[1] - sizepeak[3]) < abs(peaklim[1] - sizepeak[1])) {
				# 	bshifted <- 1
				# }
				if ((mean(peaklim)>mean(sizepeak[1:2]))&(!extendedFilter)) {
					bshifted <- 1
				}
				
				
			}else if(QC=="Snt"){
				if (SNTval < SNT)
					next
			}
			if (bshifted) {
				valCor <- NA_real_
			}
			###Intensities are calculated
			iArea <- 0
			solArea <- 0
			areaEIC <- vEic$intensity
			tempRes <- NULL
			if (NoPeak & solvar == "throw") {
				next
			} else{
				if(NoPeak){
					tempRes <- c(
						bandlist[i, "mzmin"],
						bandlist[i, "mzmax"],
						bandlist[i, "mz"],
						NA_real_,
						NA_real_,
						0,
						valSol,
						valSol,
						ifelse(bshifted,NA_real_,valCor),
						1,
						0
					)
				}else{
				if (peaklim[2] > length(vEic$scan)) {
					peaklim[2] <- length(vEic$scan)
				}
				if (peaklim[1] > length(vEic$scan))
					next
				iArea <-
					trapezArea(xraw@scantime[peaklim[1]:peaklim[2]], areaEIC[peaklim[1]:peaklim[2]])
				solArea <- NULL
				if (solint == "poly") {
					solArea <-
						(xraw@scantime[peaklim[2]] - xraw@scantime[peaklim[1]]) * valSol /
						2
				}
				else if (solint == "substract") {
					solArea <-
						(xraw@scantime[peaklim[2]] - xraw@scantime[peaklim[1]]) * valSol
				}
				else if (solint == "add") {
					solArea <- 0
				}
				if (solArea < iArea) {
					iArea <- iArea - solArea
				}
				tempRes <- c(
					bandlist[i, "mzmin"],
					bandlist[i, "mzmax"],
					bandlist[i, "mz"],
					peaklim[1],
					peaklim[2],
					iArea,
					mEic,
					valSol,
					valCor,
					bshifted,
					SNTval
				)
				}
				###The peak is only add if the peak is not detected as solvant.
				
			}
			if (QC == "Nes") {
				tempRes <- c(tempRes, pval)
			}
			
			ResList[[psip]] <- tempRes
			psip <- psip + 1
			
			
			###Graphical output if needed.
			if (graphical) {
				if (is.null(peaklim))
					next
				ntitle <- paste("mz :",
								sprintf("%0.4f", bandlist[i, "mzmin"]),
								"-",
								sprintf("%0.4f", bandlist[i, "mzmax"]))
				plot(
					xraw@scantime,
					vEic$intensity / mEic,
					type = "n",
					col = "black",
					main =
						ntitle,
					xlab = "Time",
					ylab = "Intensity"
				)
				
				###Smoothed sequence
				lines(
					xraw@scantime,
					smoothedSeq/mEic,
					col = "darkgreen",
					lwd = 2
				)
				
				##Peak Area is delimited by verticla bands.
				segments(xraw@scantime[c(peaklim[1],peaklim[2])],c(1,1),y1 = c(0,0),
						 col="red",lty=2,lwd=2)
				
				lines(xraw@scantime, vEic$intensity / mEic, lwd = 2)
				lines(
					xraw@scantime[(floor((nf + 1) / 2)):(length(xraw@scantime) - floor(nf /
																					   	2))],
					convolvedSeqF / max(convolvedSeqF),
					col = "orange",
					lwd = 2
				)
				if(second_filter){
				lines(
				xraw@scantime[(floor((ninj + 1) / 2)):(length(xraw@scantime) - floor(ninj /
																						 	2))],
					convolvedSeqInj / max(convolvedSeqInj),
					col = "purple",
					lwd = 2
				)
				}
				
			}
		}
		message(
			"\nSignal filtering finished.\n",
			length(ResList),
			" signals detected with an injection between: ",
			sprintf("%0.1f", xraw@scantime[sizepeak[1]]),
			"-",
			sprintf("%0.1f", xraw@scantime[sizepeak[1] + length(modelPeak)])
		)
		matPeak <- do.call("rbind", ResList)
		colnames(matPeak) <- headPeaks
		return(list(
			matrix = matPeak,
			injpeak = modelPeak,
			injscan = sizepeak[1]
		))
	}



OverlapSplit <- function(x, nsplit = 1, overlap = 2) {
	vsize <- length(x)
	nperdf <- ceiling((vsize + overlap * nsplit) / (nsplit + 1))
	start <-
		seq(1, nsplit * (nperdf - overlap) + 1, by = nperdf - overlap)
	sapply(start, function(i)
		x[c(i:(i + nperdf - 1))])
}


medianFiltering <- function(x, size = 5, complete = TRUE) {
	if (size %% 2 != 1) {
		stop("feneter size should be impair")
	}
	nf = size %/% 2
	matMed <- OverlapSplit(x, length(x) - size, size - 1)
	res <- apply(matMed, 2, median)
	if (complete) {
		res <- c(rep(res[1], nf), res, rep(res[length(res)], nf))
	}
	res
}


#' First guess of the limit of the injection peak.
#'
#' Determine a first approximation of the injection peak based
#' on the angle of the injection peak and changing of variation.
#' This peak is not used directly by findFIAsignal, but will be used to initialize the regression
#' giving the injection peak. Shloud be used carefuly.
#'
#' @export
#' @param xraw An xcmsRaw object as returned by \code{\link[xcms]{xcmsRaw}}.
#' @param threshold A relative increase born to detect the limit of the injection
#' peak.
#' @param graphical should the resulting limit be plotted.
#' @param scanmin The first scan to consider.
#' @param scanmax The last scan to consider.
#' @return A triplet composed of c(left limit,right limit, maximum) of the
#' estimated injection peak.
#' @aliases determiningInjectionZone
#' @examples
#' if(require(plasFIA)){
#'   #Getting the path of a file.
#'   path_raw <- list.files(system.file(package="plasFIA","mzML"),full.names=TRUE)[2]
#'
#'   #Opening the file with xcms
#'   xraw <- xcmsRaw(path_raw)
#'
#'   #Getting a first approximation of injection peak;
#'   sp <- determiningInjectionZone(xraw)
#' }
determiningInjectionZone <-
	function(xraw,
			 threshold = 0.05,
			 graphical = FALSE,
			 scanmin = NULL,
			 scanmax = NULL) {
		#Cutting the acquisition if necessary.
		seqTIC <- rawEIC(xraw, mzrange = range(xraw@env$mz))
		if (is.null(scanmin)) {
			scanmin <- 1
		}
		if (is.null(scanmax)) {
			scanmax <- length(xraw@scantime)
		}
		
		seqTIC$scan <- seqTIC$scan[scanmin:scanmax]
		seqTIC$intensity <- seqTIC$intensity[scanmin:scanmax]
		
		##Getting the range of variation of the TIC :
		rangInt <- range(seqTIC$intensity)
		
		###We make the supposition that first scan is solvent.
		valSol <- seqTIC$intensity[1]
		
		###Checking if the peak is finished.
		valThreshEnd <- diff(rangInt) * threshold + valSol
		
		###Beginning of increase of the peak.
		pinc <- seqTIC$intensity > valThreshEnd
		
		###True for 3 conscutives scans.
		pinc <-
			which(pinc[1:(length(pinc) - 2)] &
				  	pinc[2:(length(pinc) - 1)] & pinc[3:length(pinc)])[1]
		
		sizeMed <- min(pinc + (pinc + 1) %% 2, 7)
		vMedSeq <- medianFiltering(seqTIC$intensity, size = sizeMed)
		
		
		if (pinc[1] == 1) {
			stop("impossible to process the ifle as no scan is present.")
		}
		
		pbegin <- max(1, pinc[1] - 1)
		
		pmax <- which.max(seqTIC$intensity)
		
		candidates_limits <-
			(pmax + (pmax - pbegin)):length(seqTIC$intensity)
		
		###We normalize both dimension.
		vy <- vMedSeq / max(vMedSeq)
		vx <- seq(0, 1, length = length(xraw@scantime))
		inter <- vx[2] - vx[1]
		
		####Removing the one with an angl to the bottom.
		candidates_limits <-
			candidates_limits[which(seqTIC$intensity[candidates_limits - 1] > 
										seqTIC$intensity[candidates_limits] &
										seqTIC$intensity[candidates_limits + 1]<
										seqTIC$intensity[candidates_limits])]
		
		
		##We test with the angle between the summit and the end of the signal.
		a <-
			sqrt((vy[candidates_limits] - rep(vy[pmax], length(
				candidates_limits
			))) ^ 2 +
				(rep(vx[pmax], length(
					candidates_limits
				)) - vx[candidates_limits]) ^ 2)
		b <-
			sqrt((vy[candidates_limits] - rep(vy[length(vx)], length(
				candidates_limits
			))) ^ 2 +
				(rep(vx[length(vx)], length(
					candidates_limits
				)) - vx[candidates_limits]) ^ 2)
		c <-
			sqrt(rep((vx[pmax] - vx[length(vx)]) ^ 2 + (vy[pmax] - vy[length(vx)]) ^
					 	2,
					 length(candidates_limits)))
		
		##lines value <-
		lvalues <- ((rep(vx[length(vx)], length(
			candidates_limits
		)) - vx[candidates_limits])/(vx[length(vx)] - vx[pmax]))*
			(vy[pmax] - vy[length(vx)])+vy[length(vx)]
		
		pright <- which(lvalues>=vy[candidates_limits])
		if(length(pright)==0){
			warning("No right limit may be detected, this can be caused
					by wrongly formed unjection peak.")
			p3 <- scanmax
		}else{
			
			val_cos_1 <- (a[pright] ^ 2 + b[pright] ^
						  	2 - c[pright] ^ 2) / (2 * a[pright] *
						  						  	b[pright])
			
			valcos <- acos(val_cos_1)
			
			
			#We get the peak with the highest angle.
			p3 <- candidates_limits[pright[which.min(valcos - a[pright] * b[pright])]]
			
		}
		
		res <- c(pbegin,p3,pmax)
		if (graphical) {
			title <-
				paste("Initial guess of injection  zone",
					  getRawName(xraw@filepath))
			plot(
				xraw@scantime,
				seqTIC$intensity,
				type = "l",
				xlab = "Time",
				ylab = "Total Ion Intensity",
				main = title
			)
			abline(v = xraw@scantime[res], col = "red")
		}
		
		####CHecking the percentage of area taken by the peak of proFIA.
		solventVal <- mean(seqTIC$intensity[1:res[1]])
		
		iout <-
			c(unique(1:res[1]), unique(res[2]:length(seqTIC$intensity)))
		iin <- res[1]:res[2]
		
		areaIn <- trapezArea(xraw@scantime[iin], seqTIC$intensity[iin])
		areaOut <- trapezArea(xraw@scantime[iout],
									   seqTIC$intensity[iout])
		#Substracting the solvent.
		areaIn <-
			areaIn - trapezArea(xraw@scantime[iin], rep(solventVal, length(iin)))
		areaOut <-
			areaIn - trapezArea(xraw@scantime[iout], rep(solventVal, length(iout)))
		
		percentPres <- areaIn / (areaIn + areaOut)
		
		if (percentPres <= 0.75) {
			warnings(
				paste(
					"The detected injection peak only include only ",
					round(percentPres * 100),
					" of intensity.",
					sep = ""
				)
			)
			
		}
		return(res+scanmin-1)
	}




#' Determine the limits of the injection peak in a FIA acquisition.
#'
#' Determine a first approximation of the injection peak using the
#' Douglas-Peuker Algorithm provided in the \code{rgeos} package.
#' The object provided must be an xcmsRaw object.
#'
#' @export
#' @param xraw An xcmsRaw object as returned by \code{\link[xcms]{xcmsRaw}}.
#' @param freq The degrees of smoothing used on the TIC, corresponding to
#' the cutting frequency of the blackman windowed sync filter.
#' @param graphical should the resulting peak be plotted.
#' @param smooth Should the TIC be smoothed, recommended.
#' @param extended In case of very long tailing, shloud the research be extended.
#' @param percentSol If extended is TRUE, the limiting level of solvent for
#' @param scanmin The minimum scan to consider for peak detection.
#' @param scanmax the maximum scan to consider.
#' injection peak detection.
#' @return A triplet composed of c(left limit,right limit, maximum) of the
#' estimated injection peak.
#' @aliases determiningSizePeak.Geom determiningSizePeak
#' @examples
#' if(require(plasFIA)){
#'   #Getting the path of a file.
#'   path_raw <- list.files(system.file(package="plasFIA","mzML"),full.names=TRUE)[2]
#'
#'   #Opening the file with xcms
#'   xraw <- xcmsRaw(path_raw)
#'
#'   #Getting a first approximation of injection peak;
#'   sp <- determiningSizePeak.Geom(xraw)
#' }

determiningSizePeak.Geom <-
	function(xraw,
			 scanmin = 1,
			 scanmax = length(xraw@scantime),
			 freq = 0.15,
			 graphical = FALSE,
			 smooth = TRUE,
			 extended = FALSE,
			 percentSol = NULL) {
		TIC <- rawEIC(xraw, mzrange = range(xraw@env$mz))
		oTIC <- TIC
		if (smooth) {
			TIC$intensity <- smooth.BWF(TIC$intensity, freq = freq)
		}
		n <- length(TIC$intensity)
		
		###Subsetteing the value
		TIC$intensity <- TIC$intensity[scanmin:scanmax]
		
		###Solvant is approximates as the mean of the lower decile.
		
		sTime <-
			scale(xraw@scantime[TIC$scan][scanmin:scanmax], center = FALSE)
		sInt <- scale(TIC$intensity, center = FALSE)
		maxInt <- max(sInt)
		sval <- quantile(sInt, 0.15)
		epsilon <- (maxInt - sval) * 0.5
		if ((maxInt / sval) < 2)
			warning(
				"The ratio of the solvant levels on the maximum of the TIC for",
				basename(xraw@filepath),
				" is < 2. This can be a problem in the acquistion."
			)
		posMax <- 1
		approxSize <- 2
		segment <- numeric(0)
		numiter <- 0
		while (posMax < 3 | (approxSize - posMax) < 2) {
			###Checking that the area of the peak found is sufficient from the rest of the world.
			if (numiter == 50) {
				stop("Geometric algorithm could not find a peak. Check if there is one on the TIC.")
				plotTIC(xraw)
			}
			segment <- segmentCurve(sTime, sInt, epsilon)
			approxSize <- length(segment)
			
			simplInt <- TIC$intensity[segment]
			simplTime <- xraw@scantime[segment]
			posMax <- which.max(simplInt)
			epsilon <- epsilon * 0.9
			
			numiter <- numiter + 1
		}
		peaklim <- segment[c(posMax - 1, posMax + 1, posMax)]
		###Post processing.
		Solvantlevel <- mean(oTIC$intensity[1:(peaklim[1])])
		
		###Getting all the point
		thresh <- NULL
		if (extended) {
			thresh <- (1 + percentSol / 100) * Solvantlevel
			while (TIC$intensity[peaklim[2]] > thresh &
				   peaklim[2] < n) {
				peaklim[2] <- peaklim[2] + 1
			}
		}
		if (graphical) {
			TICv <- rawEIC(xraw, mzrange = range(xraw@env$mz))[scanmin:scanmax]
			ntitle <- basename(xraw@filepath)
			ntitle <- strsplit(ntitle, ".", fixed = TRUE)[[1]][1]
			plot(
				xraw@scantime[scanmin:scanmax],
				TICv$intensity,
				xlab = "Seconds",
				ylab = "Intensity",
				main = "TIC Chroamtogram",
				type = "l"
			)
			
			xseq <- xraw@scantime[peaklim[1]:peaklim[2]]
			yseq <- TIC$intensity[peaklim[1]:peaklim[2]]
			polygon(c(xseq, rev(xseq)), c(yseq, rep(0, length(yseq))), col = "red")
			lines(
				xraw@scantime,
				TIC$intensity,
				col = "purple",
				lwd = 2,
				lty = 2
			)
			lines(simplTime, simplInt, col = "darkgreen")
			if (extended)
				abline(v = thresh, col = "darkgrey")
			legend(
				"topright",
				c("raw TIC", "smoothed TIC", "simplified TIC"),
				col = c("black", "red", "darkgreen"),
				lwd = c(1, 1, 1)
			)
		}
		message(
			"An injection peak has been spotted:",
			sprintf("%0.1f", xraw@scantime[peaklim[1] + scanmin - 1]),
			"-",
			sprintf("%0.1f", xraw@scantime[peaklim[2] + scanmin - 1]),
			"s\n"
		)
		return(peaklim)
	}

#' Fit an injection peak to an FIA acquisition.
#'
#' Determine an injection peak as an exponential modified gaussian
#' function and a second order exponential corresponding to matrix
#' effect to the most intense signals in an acquisition.
#'
#' @export
#' @param xraw An xcmsRaw object as returned by \code{\link[xcms]{xcmsRaw}}.
#' @param bandlist A list of bands which can be used. Shall stay NULL in general use.
#' bands will be determined automatically.
#' @param sec A tolerance in sec to group the signals.
#' @param iquant The maximum intensity intensity threshold under which the peaks
#' would not be used for peak determination.
#' @param gpeak An approximation of the injection peak, if NULL
#' determining sizepeak.Geom will be used.
#' @param selIndex A seleciton of index to make the regression. We recommend to let it to NULL.
#' @param graphical shald the individually fitted components be plotted.
#' @return A vector contaning the injection peak
#' @aliases getInjectionPeak
#' @examples
#' if(require(plasFIA)){
#'   #Getting the path of a file.
#'   path_raw <- list.files(system.file(package="plasFIA","mzML"),full.names=TRUE)[2]
#'
#'   #Opening the file with xcms
#'   xraw <- xcmsRaw(path_raw)
#'
#'   #Getting the injection scan
#'   gp <- determiningInjectionZone(xraw)
#'
#'   #performing band detection.
#'   tbands <- findBandsFIA(xraw,ppm = 2,sizeMin = gp[3]-gp[1],beginning=gp[1],end=gp[2])
#'
#'   #Getting the injection peak
#'   vpeak <- getInjectionPeak(xraw,bandlist=tbands,gpeak=gp)
#'   plot(vpeak,type="l")
#' }

getInjectionPeak <-
	function(xraw,
			 bandlist,
			 sec = 2,
			 iquant = 0.95,
			 gpeak = NULL,
			 selIndex = NULL,
			 graphical = FALSE) {
		if (is.null(gpeak)) {
			gpeak <- determiningInjectionZone(xraw)
		}
		
		
		###Filtering the peak to retain the peak without oslven
		matDInt <- NULL
		if (!is.null(selIndex)) {
			cat("selection")
			matDInt <-
				apply(bandlist[selIndex, , drop = FALSE], 1, function(x, xraw) {
					requireNamespace("xcms")
					a = rawEIC(xraw, mzrange = c(x["mzmin"], x["mzmax"]))
					a$intensity
				}, xraw = xraw)
			
		} else{
			ttrue <- which(bandlist[, "meanSolvent"] == 0)
			ttrue <- bandlist[ttrue,]
			
			
			
			###Hard thresh
			ht <- quantile(ttrue[, "maxIntensity"], probs = iquant)
			ttrue <- ttrue[which(ttrue[, "maxIntensity"] >= ht),]
			###Making hte matrix with the selected threshold
			matInt <- apply(ttrue, 1, function(x, xraw) {
				requireNamespace("xcms")
				a = rawEIC(xraw, mzrange = c(x["mzmin"], x["mzmax"]))
				a$intensity
			}, xraw = xraw)
			
			
			###A number of 0 equal to the size of the injection peak is provided to avoid
			###Problem fitting the beginning of the peaks
			###Making a verfication that the retianed profiles does not have solvent in it
			vok <-
				which(as.logical(apply(matInt, 2, checkIso, niso = min(gpeak[1],3))))
			if (length(vok) >= 5) {
				matInt <- matInt[, vok, drop = FALSE]
			}
			
			###For each line getting the beginning of the injection peak
			
			###MODIFIED TO HANDLE SEVIN DATASET.
			getBeginning <- function(x, size = 3) {
					a <- which((x[3:length(x)] != 0) &
							   	(x[2:(length(x) - 1)] != 0) &
							   	(x[1:(length(x) - 2)] != 0))
				return(a[1])
			}
			vb <- apply(matInt, 2, getBeginning)
			pna <- which(is.na(vb))
			if (length(pna) > 0) {
				matInt <- matInt[,-pna]
				vb <- vb[-pna]
			}
			tvb <- table(vb)
			tvbo <- sort(tvb, decreasing = TRUE)
			ctvb <- cumsum(tvbo) / sum(tvbo)
			poscut <- which(ctvb > 0.6)[1]
			pvd <- density(vb, bw = sec / 3)
			pmax <- which.max(pvd$y)
			vp <- findPeaksLimits(pvd$y, pmax - 1, pmax + 1)
			
			####Checking the correctness of the gorup found.
			retentionInColumn <- function(xraw, pok, rangeSec = sec) {
				if (diff(range(xraw@scantime[pok])) > rangeSec) {
					warning("Too much retention, a clear injection peak may not be found.")
					return(FALSE)
				}
				return(TRUE)
			}
			
			inter <- pvd$x[vp]
			tokeep <- which(vb >= inter[1] & vb <= inter[2])
			if (length(tokeep) > 20) {
				oa <- apply(matInt[, tokeep], 2, max)
				tokeep <- tokeep[order(oa, decreasing = TRUE)[1:20]]
				
			}
			
			if (!is.null(selIndex)) {
				tokeep <- selIndex
			}
			
			
			matDInt <- matInt[, tokeep]
		}
		message(paste(
			ncol(matDInt),
			"chromatograms have been used for peak shape determination."
		))
		#densval=density(vb,0.5)
		###Clustring by the time limit.
		vmax <- apply(matDInt, 2, max)
		matSInt <- t(t(matDInt) / vmax)
		n <- ncol(matSInt)
		
		
		
		# plot(xraw@scantime,
		#      matInt[, 4] / max(matInt[, 4]),
		#      col = "green",
		#      type = "l")
		
		opar <- list()
		
		####○Making the first guess of the parameters.
		#mu sigma tau then the a b and h
		inita <- 0.05
		initb <- 5
		initpar <-
			c(
				xraw@scantime[gpeak[3] - floor((gpeak[3] - gpeak[1]) / 2)],
				(xraw@scantime[gpeak[3]] - xraw@scantime[gpeak[1]]) * 0.5,
				20,
				rep(inita, n),
				rep(initb, n),
				vmax * 1.2
			)
		#matResiduals<-function(mpp,xx,mobs,type="gaussian",n){
		weigthvec <- length(c(
			rep(2, gpeak[1] - 1),
			5,
			seq(1, 2, length = (gpeak[3] - gpeak[1]) - 1),
			seq(2, 1.5, length = (gpeak[2] - gpeak[3] +
								  	1)),
			rep(1.5, length(xraw@scantime) - gpeak[2])
		))
		lower <- c(0, 0,-Inf, rep(0, 3 * n))
		
		
		nlC <- nls.lm.control(maxiter = 100, ptol = 0.001)
		tMat <- t(t(matDInt) / apply(matDInt, 2, max))
		if (graphical) {
			matplot(tMat,
					type = "l",
					xlab = "Time (s)",
					ylab = "Scaled intensity")
		}
		parestimate <- nls.lm(
			par = initpar,
			fn = matResiduals,
			observed = tMat,
			type = "gaussian",
			beginning = gpeak[1],
			xx = xraw@scantime,
			control = nlC,
			n = n,
			lower = lower,
			weigth = weigthvec
		)
		cest <- coef(parestimate)
		parv <- NULL
		parv <- cest[1:3]
		multiplier <- numeric(n)
		for (i in 1:n) {
			parm <- cest[c(3 + i, 3 + i + n, 3 + i + 2 * n)]
			tr <- persoConv(xraw@scantime, parv, type = "gaussian")
			mat_eff <- -matrix_effect(tr, parm[1], parm[2])
			fitted <- (tr + mat_eff) * parm[3]
			fitted[1:gpeak[1]] <- 0
			tr[1:gpeak[1]] <- 0
			multiplier[i] <- max(fitted)
			if (graphical) {
				plot(xraw@scantime,
					 tMat[, i],
					 ylim = c(0, max(1, tr * parm[3])),
					 type = "l")
				lines(xraw@scantime, fitted, col = "red")
				lines(xraw@scantime, (mat_eff + 1) * parm[3], col = "blue")
				lines(xraw@scantime, tr * parm[3] , col = "green")
			}
		}
		TP <- persoConv(xraw@scantime, p = parv)
		TP[1:gpeak[1]] <- 0
		tMat <- tMat %*% diag(multiplier)
		if (graphical) {
			matplot(
				tMat,
				type = "l",
				col = rainbow(
					n,
					start = 0.25,
					end = 0.6,
					alpha = 0.5
				),
				xlab = "Scan",
				ylab = "Intensity",
				main = paste(
					"Regressed injection peak for",
					getRawName(xraw@filepath)
				)
			)
			lines(TP, col = "red", lwd = 2)
			legend("topright",
				   c("regressed signal"),
				   col = "red",
				   lty = 1)
		}
		return(TP)
	}


triangleDistribution <- function(x) {
	if (x > 0 & x < 0.5) {
		return(2 * x)
	}
	if (x >= 0.5 & x < 1) {
		return(1 - 2 * (x - 0.5))
	}
	return(0)
}


###mu sigma tau.
persoConv <- function(time, p, type = c("gaussian", "triangle")) {
	type <- match.arg(type)
	mu <- p[1]
	sigma <- p[2]
	x <- NULL
	#cat(mu,sigma,tau,"\n")
	if (type == "gaussian") {
		x <- dnorm(time, mean = mu, sd = sigma)
	}
	if (type == "triangle") {
		x <- sapply((time - mu) / sigma + 0.5, triangleDistribution)
	}
	tau <- p[3]
	yexp <- exp(-time / tau) / tau
	vv <- convolve(yexp, rev(x))
	return(vv / max(vv))
}

###Term - a added to remove the intensity in 0
matrix_effect <- function(intensity, a, b) {
	tr <-
		a * exp(b * intensity) - a#+(p$c)*exp(p$d*intensity)) #+p$c*exp(p$d*intensity))/(p$a+p$c)
	tr
}


matResiduals <-
	function(mpp,
			 xx,
			 observed,
			 type = "gaussian",
			 n,
			 beginning,
			 weigth) {
		mu <- rep(mpp[1], n)
		sigma <- rep(mpp[2], n)
		tau <- rep(mpp[3], n)
		a <- mpp[4:(3 + n)]
		b <- mpp[(4 + n):(3 + 2 * n)]
		h <- mpp[(4 + 2 * n):(3 + 3 * n)]
		matpar <- matrix(c(mu, sigma, tau, a, b, h), ncol = 6)
		trmat <- apply(matpar, 1, persoConv, time = xx)
		mmareffect <-
			sapply(1:ncol(trmat), function(x, va, vb, mat) {
				a <- va[x]
				b <- vb[x]
				- matrix_effect(mat[, x], a, b)
			}, mat = trmat, va = a, vb = b)
		#We ignore the point before the solvent peak.
		
		trmatres <- t(t(trmat + mmareffect) * h)
		matdiff <- observed - trmatres
		matdiff <- matdiff * weigth
		matdiff[1:beginning,] <- 0
		return(matdiff)
	}

# Fit an injection peak to an FIA acquisition using likehood maximization
#
# Determine an injection peak as an exponential modified gaussian
# function and a second order exponential corresponding to matrix
# effect to the most intense signals in an acquisition.
#
# @param xraw An xcmsRaw object as returned by \code{\link[xcms]{xcmsRaw}}.
# @param bandlist A list of bands which can be used. Shall stay NULL in general use.
# bands will be determined automatically.
# @param noiseModel A noiseEstimation object, usually will be passed by the proFIAset
# function;
# @param sec A tolerance in sec to group the signals.
# @param iquant The maximum intensity intensity threshold under which the peaks
# would not be used for peak determination.
# @param gpeak An approximation of the injection peak, if NULL
# determining sizepeak.Geom will be used.
# @param graphical shald the individually fitted components be plotted.
# @return A vector contaning the injection peak
# @aliases getInjectionPeak.logL
# @examples
# if(require(plasFIA)){
#   #Getting the path of a file.
#   path_raw <- list.files(system.file(package="plasFIA","mzML"),full.names=TRUE)[2]
#
#   #Opening the file with xcms
#   xraw <- xcmsRaw(path_raw)
#
#   #Getting the injection scan
#   gp <- determiningInjectionZone(xraw)
#
#   #performing band detection.
#   tbands <- findBandsFIA(xraw,ppm = 2,sizeMin = gp[3]-gp[1],beginning=gp[1],end=gp[2])
#
# 	 #Loading a noise model.
# 	 data(plasSet)
# 	 nes <- plasSet@noiseEstimation
#
#   #Getting the injection peak
#   vpeak <- getInjectionPeak(xraw,bandlist=tbands,noiseModel = nes,gpeak=gp)
#   plot(vpeak,type="l")
# }

# logL <-
# 	function(mpp,
# 			 xx,
# 			 observed,
# 			 weigth,noiseModel=NULL){
# 		n <- (length(mpp)-3)/2
# 		mu <- rep(mpp[1], n)
# 		sigma <- rep(mpp[2], n)
# 		tau <- rep(mpp[3], n)
# 		a <- mpp[4:(3 + n)]
# 		b <- mpp[(4 + n):(3 + 2 * n)]
# 		h <- mpp[(4 + 2 * n):(3 + 3 * n)]
# 		matpar <- matrix(c(mu, sigma, tau, a, b, h), ncol = 6)
# 		trmat <- apply(matpar, 1, persoConv, time = xx)
# 		mmareffect <- sapply(1:ncol(trmat), function(x, va, vb, mat) {
# 			a <- va[x]
# 			b <- vb[x]
# 			- matrix_effect(mat[, x], a, b)
# 		}, mat = trmat, va = a, vb = b)
# 		trmatres <- t(t(trmat + mmareffect) * h)
# 		matdiff <- observed - trmatres
# 		matdiff <- matdiff
# 		###Calculating the noise variance.
# 		nvar <- apply(trmatres,1,noiseModel@estimation)
# 		prob <- rbind(as.numeric(nvar),as.numeric(matdiff))
# 		prob <- apply(prob,2,function(x){
# 			log10(dnorm(x[2],0,x[1]))
# 		})
# 		###Calculating the prob
#
# 		return(prob)
# 	}
#
# getInjectionPeak.logL <-
# 	function(xraw,
# 			 bandlist=NULL,
# 			 noiseModel = NULL,
# 			 sec = 2,
# 			 iquant = 0.95,
# 			 gpeak = NULL,
# 			 graphical = FALSE, maxSig = 10) {
# 		if (is.null(gpeak)) {
# 			gpeak <- determiningInjectionZone(xraw)
# 		}
# 		if(is.null(bandlist)){
# 			message("No bandList provided")
# 			bandlist=findBandsFIA(xraw,ppm=)
# 		}
#
#
# 		###Filtering the peak to retain the peak without oslven
# 		ttrue <- which(bandlist[, "meanSolvent"] == 0)
# 		ttrue <- bandlist[ttrue,]
#
#
#
# 		###Hard thresh
# 		ht <- quantile(ttrue[, "maxIntensity"], probs = iquant)
# 		ttrue <- ttrue[which(ttrue[, "maxIntensity"] >= ht),]
# 		###Making hte matrix with the selected threshold
# 		matInt <- apply(ttrue, 1, function(x, xraw) {
# 			requireNamespace("xcms")
# 			a = rawEIC(xraw, mzrange = c(x["mzmin"], x["mzmax"]))
# 			a$intensity
# 		}, xraw = xraw)
#
#
# 		###A number of 0 equal to the size of the injection peak is provided to avoid
# 		###Problem fitting the beginning of the peaks
# 		###Making a verfication that the retianed profiles does not have oslvent in it
# 		vok <- which(as.logical(apply(matInt, 2, checkIso, niso = 3)))
# 		if (length(vok) >= 5) {
# 			matInt <- matInt[, vok, drop = FALSE]
# 		}
#
# 		###For each line getting the beginning of the injection peak
# 		getBeginning <- function(x, size = 3) {
# 			a <- which((x[3:length(x)] != 0) &
# 					   	(x[2:(length(x) - 1)] != 0) &
# 					   	(x[1:(length(x) - 2)] != 0))
# 			return(a[1])
# 		}
# 		vb <- apply(matInt, 2, getBeginning)
#
# 		tvb <- table(vb)
# 		tvbo <- sort(tvb, decreasing = TRUE)
# 		ctvb <- cumsum(tvbo) / sum(tvbo)
# 		poscut <- which(ctvb > 0.6)[1]
# 		pvd <- density(vb, bw = sec / 3)
# 		pmax <- which.max(pvd$y)
# 		vp <- findPeaksLimits(pvd$y, pmax - 1, pmax + 1)
#
# 		####Checking the correctness of the gorup found.
# 		retentionInColumn <- function(xraw, pok, rangeSec = sec) {
# 			if (diff(range(xraw@scantime[pok])) > rangeSec) {
# 				warning("Too much retention, a clear injection peak may not be found.")
# 				return(FALSE)
# 			}
# 			return(TRUE)
# 		}
#
# 		inter <- pvd$x[vp]
# 		matDInt <- matInt[, which(vb >= inter[1] & vb <= inter[2])]
# 		#print(ttrue[which(vb %in% tokeep),])
# 		#return(ttrue[which(vb %in% tokeep),])
# 		message(paste(
# 			ncol(matDInt),
# 			"chromatograms have been used for peak shape determination."
# 		))
# 		#densval=density(vb,0.5)
# 		###Clustring by the time limit.
# 		vmax <- apply(matDInt, 2, max)
# 		if(ncol(matDInt)>maxSig){
# 			matDInt[,order(vmax,decreasing=TRUE)[1:maxSig]]
# 		}
# 		matSInt <- t(t(matDInt) / vmax)
# 		n <- ncol(matSInt)
#
#
# 		####○Making the first guess of the parameters.
# 		#mu sigma tau then the a b and h
# 		inita <- 0.06
# 		initb <- 0.5
# 		initmu <- xraw@scantime[gpeak[3]]
# 		initsig <- (xraw@scantime[gpeak[1]] + xraw@scantime[gpeak[2]]) / 5
# 		inittau <- 10
# 		initpar <- c(initmu,
# 					 initsig,
# 					 inittau,
# 					 rep(inita, n),
# 					 rep(initb, n))
#
# 		h <- vmax
# 		tMat <- t(t(matDInt) / apply(matDInt, 2, max))
#
# 		#Log likehood function.
# 		#TODO pass this in C code and otpimize it.
# 		logL <-
# 			function(mpp){
# 				n <- (length(mpp)-3)/2
# 				mu <- rep(mpp[1], n)
# 				sigma <- rep(mpp[2], n)
# 				tau <- rep(mpp[3], n)
# 				a <- mpp[4:(3 + n)]
# 				b <- mpp[(4 + n):(3 + 2 * n)]
# 				matpar <- matrix(c(mu, sigma, tau, a, b, h), ncol = 6)
# 				trmat <- apply(matpar, 1, persoConv, time = xraw@scantime)
# 				mmareffect <- sapply(1:ncol(trmat), function(x, va, vb, mat) {
# 					a <- va[x]
# 					b <- vb[x]
# 					- matrix_effect(mat[, x], a, b)
# 				}, mat = trmat, va = a, vb = b)
# 				trmatres <- t(t(trmat + mmareffect) * h)
# 				matdiff <- matDInt - trmatres
# 				matdiff <- matdiff
# 				###Calculating the noise variance.
# 				nvar <- apply(abs(trmatres),1,noiseModel@estimation)
# 				prob <- rbind(as.numeric(nvar),as.numeric(matdiff))
# 				prob <- apply(prob,2,function(x){
# 					log10(dnorm(x[2],0,x[1]))
# 				})
# 				#browser()
# 				###Calculating the prob
# 				return(prob)
# 			}
# 		#
# 		# LL<-function(parv){
# 		# 	#DEBUG addition of a plot.
# 		#
# 		# 	mu = parv[1]
# 		# 	sigma = parv[2]
# 		# 	tau = parv[3]
# 		# 	alla = parv[4:(3+(length(parv)-3)/2)]
# 		# 	allb = parv[(4+(length(parv)-3)/2):(3+2*(length(parv)-3)/2)]
# 		# 	logl <- numeric(length(as.numeric(matSInt)))
# 		#
# 		# 	###We remove the h "parameter we will add it at the end.
# 		# 	for(i in 1:length(alla)){
# 		# 		TP <- persoConv(xraw@scantime, p = c(mu,sigma,tau),type="gaussian")
# 		# 		mf <- -matrix_effect(TP, alla[i], allb[i])
# 		# 		obs <- (TP+mf)
# 		# 		obs <- (obs+min(obs))
# 		# 		diff <- vmax[i]*(matSInt[,i]-obs)
# 		#
# 		# 		vvar <- noiseModel@estimation(abs(diff))
# 		# 		vdnorm <- apply(rbind(diff,vvar),2,function(x){dnorm(x[1],0,sd=x[2])})
# 		#
# 		# 		logl[((i-1)*nrow(matSInt)+1):(i*nrow(matSInt))] <- (log(vdnorm))
# 		# 		##DEBUG ADDING THE PLOT ON THE SAME PARAMETER.
# 		# 		if(i==5){
# 		# 			plot(matSInt[,i])
# 		# 			lines(obs,col="red")
# 		# 			#browser()
# 		# 		}
# 		#
# 		# 	}
# 		# 	print(sum(logl))
# 		# 	return(logl)
# 		# }
#
# 		##DEBUG
# 		###Printing the LL for all the value.
# 		# tauseq <- seq(inittau/2,inittau*1.5,length=40)
# 		# sigseq <- seq(initsig/2,initsig*1.5,length=40)
# 		# res <- numeric(length(tauseq)*length(sigseq))
# 		# for(i in 1:length(sigseq)){
# 		# 	for(j in 1:length(tauseq)){
# 		# 		vpar <- c(initmu,
# 		# 				  sigseq[i],
# 		# 					tauseq[j],
# 		# 					 rep(inita,n),
# 		# 				  		rep(initb,n))
# 		#
# 		# 		res[(i-1)*length(sigseq)+j] <- LL(vpar)
# 		# 	}
# 		# }
# 		# library(lattice)
# 		#
# 		# wireframe(res ~ tauseq*sigseq, data = NULL,
# 		# 		  xlab = "tau", ylab = "sig",
# 		# 		  main = "tausig",
# 		# 		  drape = TRUE,
# 		# 		  colorkey = TRUE,
# 		# 		  screen = list(z = -60, x = -60)
# 		# )
#
# 		###END DEBUG
#
# 		ml <- maxLik( logL,start = initpar,method="BHHH", control=list(tol=1,iterlim = 200,printLevel=2))
# 		# plot(xraw@scantime,
# 		#      matInt[, 4] / max(matInt[, 4]),
# 		#      col = "green",
# 		#      type = "l")
#
# 		opar <- list()
# 		print(summary(ml))
# 		cest <- as.numeric(coef(ml))
# 		print(cest)
# 		cat(paste("cest",cest))
# 		parv <- cest[1:3]
# 		multiplier <- numeric(n)
# 		for (i in 1:n) {
# 			parm <- cest[c(3 + i, 3 + i + n)]
# 			hi <- vmax[i]
# 			cat(paste("parm",parm,"parv",parv))
# 			tr <- persoConv(xraw@scantime, parv, type = "gaussian")
# 			mat_eff <- -matrix_effect(tr, parm[1], parm[2])
# 			fitted <- (tr + mat_eff) * hi
# 			multiplier[i] <- max(fitted) / hi
# 			if (graphical) {
# 				plot(xraw@scantime, tMat[, i], type = "l")
# 				lines(xraw@scantime, fitted, col = "red")
# 				lines(xraw@scantime, mat_eff + 1, col = "blue")
# 				lines(xraw@scantime, tr , col = "green")
# 			}
# 		}
# 		TP <- persoConv(xraw@scantime, p = parv)
# 		tMat <- tMat %*% diag(multiplier)
# 		if(graphical){
# 			matplot(
# 				tMat,
# 				type = "l",
# 				col = rainbow(
# 					n,
# 					start = 0.25,
# 					end = 0.6,
# 					alpha = 0.5
# 				),
# 				xlab = "Scan",
# 				ylab = "Intensity",
# 				main = paste("Regressed injection peak for",
# 							 getRawName(xraw@filepath))
# 			)
# 			lines(TP, col = "red", lwd = 2)
# 			legend("topright",
# 				   c("regressed signal"),
# 				   col = "red",
# 				   lty = 1)
# 		}
# 		return(TP)
# 	}





matResidualsLikehood <- function(mpp,
								 xx,
								 observed,
								 type = "gaussian",
								 n,
								 beginning,
								 weigth,
								 noisefunction) {
	mu <- rep(mpp[1], n)
	sigma <- rep(mpp[2], n)
	tau <- rep(mpp[3], n)
	a <- mpp[4:(3 + n)]
	b <- mpp[(4 + n):(3 + 2 * n)]
	h <- mpp[(4 + 2 * n):(3 + 3 * n)]
	matpar <- matrix(c(mu, sigma, tau, a, b, h), ncol = 6)
	trmat <- apply(matpar, 1, persoConv, time = xx)
	mmareffect <- sapply(1:ncol(trmat), function(x, va, vb, mat) {
		a <- va[x]
		b <- vb[x]
		- matrix_effect(mat[, x], a, b)
	}, mat = trmat, va = a, vb = b)
	trmatres <- t(t(trmat + mmareffect) * h)
	
	LL <- function(diff) {
		R = sum(log10(sapply(diff, dnorm)))
	}
	
	
	
	
	matdiff <- observed - trmatres
	
	
	###Noise estimation function.
	
	matdiff <- matdiff * weigth
	return(matdiff)
}

###Wrapper for parallelism
openAndFindPeaks <- function(fname, ppm, es = es, ...) {
	requireNamespace("xcms")
	requireNamespace("proFIA")
	xraw <- xcmsRaw(fname)
	message("processing: ", fname, "\n")
	tP <- findFIASignal(xraw, ppm, es = es, ... = ...)
	
	message("\n")
	list(
		signals = tP$matrix,
		injectionPeak = tP$injpeak,
		injectionScan = tP$injscan
	)
	
}
