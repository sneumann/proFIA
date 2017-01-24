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
    ##Binarizing the data
    shift_right <- x[-1] == 0
    shift_left <- x[-length(x)] == 0
    if (all(shift_right | shift_left))
        return(TRUE)
    return(FALSE)
}


#' Detect peaks in an FIA acquisition.
#'
#' Detect the peak corresponding to compounds present in the sample
#' in a Flow Injection Analysis (FIA) acquisition. The item furnished
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
#' a noise model is furnished.
#' @param ... more arguments to be passed to the \link{determiningSizePeak} function.
#' @param scanmin The first scan to consider.
#' @param scanmax The last scan to consider.
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
#'     \item Sig_Sol The ratio of the signal max intensity on the oslvent max intensity.
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
             solint = c("poly", "substract", "keep"),
             dmz =
                 0.0005,
             graphical = FALSE,
             SNT = NULL,
             f = c("regression", "TIC"),
             pvalthresh = NULL,
    		 scanmin = 1,
    		 scanmax = length(xraw@scantime),
             ...) {
        solint <- match.arg(solint)
        solvar <- match.arg(solvar)
        f <- match.arg(f)
        n <- length(xraw@scantime)
        sizepeak <- determiningSizePeak.Geom(xraw,scanmin = scanmin, scanmax = scanmax, ...)
        sizepeak <- sizepeak + scanmin -1
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
                "corPeak",
                "shifted",
                "Sig_Sol"
            )
        
        ###Which kind of quality control will be used.
        QC <- "Reduced"
        if (is.null(es) & is.null(SNT)) {
            stop("No noise model and no SNT furnished.")
        } else if (is.null(es) & is.numeric(SNT)) {
            message("No noise model furnished, a signal/solvant threshold will be used.")
        } else if (!is.null(es@estimation)) {
            if (is.null(pvalthresh)) {
                message(
                    "A noise model is furnished but no threshold was specified, the threshold is therefore set to 0.01 by default"
                )
                pvalthresh <- 0.01
            } else{
                message("A noise model and a pvalue threshold are furnished, SNT will be ignored.")
            }
            QC <- "Nes"
        } else{
            stop("Wrong type of noise or SNT, check the parameters.")
        }
        if (QC == "Nes") {
            headPeaks <- c(headPeaks, "pvalue")
        }
        modelPeak <- NULL
        
        AllBand <- findBandsFIA(
            xraw,
            ppm = ppm,
            sizeMin = (sizepeak[3] -
                           sizepeak[1]) * 0.75,
            dmz =
                dmz,
            beginning = sizepeak[1], firstScan = scanmin,
            lastScan = scanmax
        )
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
                AllBand,
                sec = 2,
                iquant = 0.95,
                gpeak = sizepeak,
                graphical = graphical
            )
            #cat("klkl");print(modelPeak);cat("kllkkl")
            ###Beginning and end correspond to the position where th epeak int is superior at 5%
            mpPmax <- which.max(modelPeak)
            pmmin <- which(modelPeak[1:mpPmax] > 0.05 * modelPeak[mpPmax])
            pmmin <- pmmin[1]
            pmmax <- which(modelPeak[(mpPmax + 1):length(modelPeak)] > 0.05 *
                              modelPeak[mpPmax]) + mpPmax - 1
            pmmax <- pmmax[length(pmmax)]
            modelPeak <- modelPeak[pmmin:pmmax]
            modelPeak <- modelPeak / max(modelPeak)
            sizepeak[2] <- max(sizepeak[2],pmmax)
        }
        
        ###defining the injection peak model.
        xseq <- seq(-2, 2, length = 1000)
        yTriangl <- sapply(xseq, AntiSymTriangle, simplify = TRUE)
        sizeInc <- floor((sizepeak[3] - sizepeak[1]))
        seqIndex <- floor(seq(0, 1000, length = sizeInc))
        modelInjPeak <- yTriangl[seqIndex]
        if(sizeInc < 3)
        	modelInjPeak <- yTriangl[c(0,500,1000)]
        nInj <- length(modelInjPeak)
        
        
        
        
        ##Constant which will be used on all iteration.
        nf <- length(modelPeak)
        ninj <- length(modelInjPeak)
        countSol <- 0
        adjustZone <- floor((sizepeak[3] - sizepeak[1]) / 6)
        
        ###list which will stock the result.
        ResList <- list()
        message("Band filtering: ",appendLF = FALSE)
        memmess <- 1
        for (i in 1:nrow(AllBand)) {
        	bshifted <- 0
            vact <- floor(i/nrow(AllBand)*100)%/%10
            if(vact!=memmess){
                memmess <- vact
                message(memmess*10," ",appendLF = FALSE)
             }
            vEic <- rawEIC(xraw, mzrange = c(AllBand[i, "mzmin"], AllBand[i, "mzmax"]))
            #Extending the sequence to at least two times the size of the 
            #injection filter
            if(length(vEic$intensity)/length(modelPeak)<2){
                nlength <- nextn(2*length(modelPeak),2)
                valrep <- mean(vEic$intensity[(length(vEic$intensity)-3):length(vEic$intensity)])
                vEic$intensity <- c(vEic$intensity,rep(valrep,
                                                    times=nlength-length(vEic$intensity)))
            }
            
            
            
            ###Determining if there is solvant in the signal.
            #seqSol = vEic$intensity[1:sizepeak[1]]
            ###CHecking if there is only isolated values, 0 means no reliably detectable solvant.
            Bsol <- ifelse(checkIsoValues(vEic$intensity[1:sizepeak[1]]), FALSE, TRUE)
            valSol <- NULL
            if (Bsol) {
                valOkSol <- which(vEic$intensity[1:sizepeak[1]]!=0)
                valSol <- mean(vEic$intensity[valOkSol])
            } else{
                valSol <- 0
            }
            ##Calculating the correlation.
            mEic <- max(vEic$intensity)
            pos0 <- which(vEic$intensity[sizepeak[1]:sizepeak[2]] == 0) + sizepeak[1] -
                1
            valCor <- NULL
            cinter <- pmmin:pmmax
            if (length(pos0) == 0) {
                valCor <- cor(modelPeak, vEic$intensity[cinter])
            } else{
                if ((sizepeak[2] - sizepeak[1] + 1 - length(pos0)) > 1) {
                    valCor <- suppressWarnings(cor(modelPeak[-pos0], vEic$intensity[cinter][-pos0]))
                } else{
                    valCor <- 0
                }
            }
            ###Problematic case in really noisy peaks.
            if(is.na(valCor)){valCor <- 0}
            
            ##Calculating the convolution
            convolvedSeqF <- convolve(vEic$intensity, modelPeak, type = "filter")
            convolvedSeqInj <- convolve(vEic$intensity, modelInjPeak, type = "filter")
            smoothedSeq <- smooth.BWF(vEic$intensity, freq = 0.2)
            ###Getting the local maximum.
            pMf <- which.max(convolvedSeqF[scanmin:scanmax])+scanmin
            pMmin <- which.min(convolvedSeqF)
            posMax <- pMf + (floor(nf / 2))
            posMin <- pMmin + (floor(nf / 2))
            
            ###Is the local maximum to the right of the injection peak ?
            peaklim <- NULL
            SecondMax <- NULL
            FirstMax <- FALSE
            PeakFound <- FALSE
            ###If the maximum detected is in the injection peak.
            if (posMax < sizepeak[2] && posMax > sizepeak[1]) {
                ###First estimation of the peak limit. peaklim will sotre the peak limit all along the process.
                peaklim <- c(posMax - (floor(nf / 2)), posMax + (floor(nf / 2)))
                ###Locating the limit of the peak found
                pl <- findPeaksLimits(smoothedSeq, peaklim[1] + 1, peaklim[2] - 1)
                peaklim[2] <- pl[2]
                ###The filter coeff goes down to find if there is a second peak.
                #  cat("posRawMax :",convolvedSeqF[max((posMax - adjustZone), 1):(posMax +
                #  															   	adjustZone)]," sup",(posMax - adjustZone),"\n")
                posRawMax <- which.max(convolvedSeqF[max((pMf - adjustZone), 1):(pMf +
                                                                                       adjustZone)]) + max(((pMf - adjustZone) - 1),1)
                flim <- findPeaksLimits(convolvedSeqF, posRawMax - 1, posRawMax + 1)
                if (posMax < sizepeak[2] &
                    ###Case where the max is the right of the peak.
                    posMax > sizepeak[3]) {
                    SecondInter <- c(sizepeak[1] - floor(nf / 2),
                                    min(flim[1] - floor(nf /
                                                            2), pl[1]))
                    SecondInter[1] = max(SecondInter[1], 0)
                    if (SecondInter[2] - SecondInter[1] > 1) {
                        Secondpmf <- which.max(convolvedSeqF[SecondInter[1]:SecondInter[2]]) +
                            SecondInter[1] - 1
                        SecondMax <- Secondpmf + (floor(nf / 2))
                        
                        ###CHecking that the second max is in the left part of the injection peak.
                        if (SecondMax > sizepeak[1] &
                            SecondMax < sizepeak[3]) {
                            ###The limit of the second peak is localised in a decent windows.
                            pl <- findPeaksLimits(smoothedSeq,
                                                 Secondpmf - 1,
                                                 Secondpmf + 1)
                            peaklim[1] <- pl[1]
                        } else {
                            ###Injection filter is checked in the right direction
                            SecondInter <- c(sizepeak[1] - floor(ninj / 2), pl[1])
                            SecondInter[1] <- max(SecondInter[1], 0)
                            if ((SecondInter[2] - SecondInter[1]) > 0) {
                                ###Checing if we are not too close from the beginning.
                                SecondMax <- which.max(convolvedSeqInj[SecondInter[1]:SecondInter[2]]) +
                                    SecondInter[1] - 1
                                Secondpmf <- SecondMax + (floor(ninj / 4))
                                if (Secondpmf > sizepeak[1] & Secondpmf < pl[1]) {
                                    rawSecondPmf <- which.max(smoothedSeq[(SecondMax - adjustZone):(SecondMax +
                                                                  adjustZone)]) + (SecondMax - adjustZone) - 1
                                    pl <- findPeaksLimits(smoothedSeq,
                                                         rawSecondPmf - 2,
                                                         rawSecondPmf + 2)
                                    peaklim[1] <- pl[1]
                                }
                            }
                        }
                    }
                }
                if (posMax > sizepeak[1] &
                    ###Case where the detected mas is at the left of the peak.
                    posMax < sizepeak[3]) {
                    SecondInter <- c(flim[2], sizepeak[2] - floor(nf / 2))
                    SecondMax <- which.max(convolvedSeqF[SecondInter[1]:SecondInter[2]]) +
                        SecondInter[1] - 1
                    Secondpmf <- SecondMax + (floor(nf / 2))
                    
                    ###CHecking that the second max is in the left part of the injection peak.
                    if (Secondpmf > sizepeak[1] & Secondpmf < sizepeak[3]) {
                        ###The limit of the second peak is localised in a decent windows.
                        rawSecondPmf <- which.max(smoothedSeq[(SecondMax - adjustZone):(SecondMax +
                                                                                           adjustZone)]) + (SecondMax - adjustZone) - 1
                        pl <- findPeaksLimits(smoothedSeq,
                                             rawSecondPmf - 1,
                                             rawSecondPmf + 1)
                        peaklim[2] <- pl[2]
                    }
                }
            }
            if (posMax > sizepeak[2]) {
                ##Checking that there is some retention in  the colmun.
                ##If it the case the peak is wider than usual, so wider than the
                ##modelPeak.
                #peaklim = floor(nf / 2)+c(pMf - (floor(nf / 2)), pMf + (floor(nf / 2)))
                ###Locating the limit of the peak found
            	bshifted = 1
                pl <- findPeaksLimits(convolvedSeqF, pMf - 3, pMf + 3)
                peaklim <- pl+floor(nf / 2)
                #cat("mumu",posMax,"    ")
                if(diff(peaklim)<length(modelPeak)) peaklim=NULL
            }
            
            ###Resizeing the first limit
            ###Quality control of the found peak.
            NoPeak <- is.null(peaklim)
            if (!NoPeak) {
                peaklim[1] = max(peaklim[1], sizepeak[1])
            }
            pval <- NULL
            SNTval <- smoothedSeq[posMax] / valSol
            if (NoPeak) {
                pval <- 1
                SNTval <- 0
            } else if (valSol == 0) {
                pval <- 0
                SNTval <- Inf
                ###Checking that the final peak detected is not shifted.
                if(peaklim[2]>sizepeak[2] & abs(posMax-sizepeak[2]) < abs(posMax-sizepeak[3]) &
                   abs(peaklim[1]-sizepeak[3])<abs(peaklim[1]-sizepeak[1])){
                	bshifted <- 1
                }
            } else if (QC == "Nes") {
                pval <- calcPvalue(
                    es,
                    vEic$intensity[peaklim[1]:peaklim[2]],
                    seq(smoothedSeq[peaklim[1]], smoothedSeq[peaklim[2]],
                        length = peaklim[2] - peaklim[1] + 1)
                )
                if (pval > pvalthresh)
                    next
                if (SNTval < 2)
                    next
                ###Checking that the final peak detected is not shifted.
                if(peaklim[2]>sizepeak[2] & abs(posMax-sizepeak[2]) < abs(posMax-sizepeak[3]) &
                   abs(peaklim[1]-sizepeak[3])<abs(peaklim[1]-sizepeak[1])){
                	bshifted <- 1
                }
                
                
            }
            if(bshifted){
            	valCor <- -2
            }
            ###Intensities are calculated
            iArea <- 0
            solArea <- 0
            areaEIC <- vEic$intensity
            tempRes <- NULL
            if (NoPeak & solvar == "throw") {
                next
            } else{
                if(peaklim[2]>length(vEic$scan)){
                    peaklim[2] <- length(vEic$scan)
                }
                if(peaklim[1]>length(vEic$scan)) next
                iArea <- trapezArea(xraw@scantime[peaklim[1]:peaklim[2]], areaEIC[peaklim[1]:peaklim[2]])
                solArea <- NULL
                if (solint == "poly") {
                    solArea <- (xraw@scantime[peaklim[2]] - xraw@scantime[peaklim[1]]) * valSol /
                        2
                }
                else if (solint == "substract") {
                    solArea <- (xraw@scantime[peaklim[2]] - xraw@scantime[peaklim[1]]) * valSol
                }
                else if (solint == "add") {
                    solArea <- 0
                }
                if (solArea < iArea) {
                    iArea <- iArea - solArea
                }
                tempRes <- c(
                    AllBand[i, "mzmin"],
                    AllBand[i, "mzmax"],
                    AllBand[i, "mz"],
                    peaklim[1],
                    peaklim[2],
                    iArea,
                    mEic,
                    valSol,
                    valCor,
                    bshifted,
                    SNTval
                )
                
                ###The peak is only add if the peak is not detected as solvant.
                
            }
            if (QC == "Nes") {
                tempRes <- c(tempRes, pval)
            }
            
            ResList[[psip]] <- tempRes
            psip <- psip + 1
            
            ###Eventually affine the peak limit.
            if (graphical) {
                if (is.null(peaklim))
                    next
                ntitle <- paste("mz :",
                               sprintf("%0.4f", AllBand[i, "mzmin"]),
                               "-",
                               sprintf("%0.4f", AllBand[i, "mzmax"]))
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
                yseq <- areaEIC[peaklim[1]:peaklim[2]] / mEic
                yseq <- c(yseq, rep(0, peaklim[2] -
                                       peaklim[1] + 1))
                lines(
                    xraw@scantime[sizepeak[1]:(sizepeak[2])],
                    modelPeak / max(modelPeak),
                    col =
                        "darkgrey",
                    lwd = 2
                )
                lines(xraw@scantime, vEic$intensity / mEic, lwd = 2)
                lines(
                    xraw@scantime[(floor((nf + 1) / 2)):(length(xraw@scantime) - floor(nf /
                                                                                           2))],
                    convolvedSeqF / max(convolvedSeqF),
                    col = "orange",
                    lwd = 2
                )
                lines(
                    xraw@scantime[(floor((ninj + 1) / 2)):(length(xraw@scantime) - floor(ninj /
                                                                                             2))],
                    convolvedSeqInj / max(convolvedSeqInj),
                    col = "purple",
                    lwd = 2
                )
                
            }
        }
        message(
            "\nSignal filtering finished.\n",
            length(ResList),
            " signals detected with an injection between: ",
            sprintf("%0.1f", xraw@scantime[sizepeak[1]]),
            "-",
            sprintf("%0.1f", xraw@scantime[sizepeak[1]+length(modelPeak)])
        )
        matPeak <- do.call("rbind", ResList)
        colnames(matPeak) <- headPeaks
        return(list(
            matrix = matPeak,
            injpeak = modelPeak,
            injscan = pmmin
        ))
    }

#' Determine the limits of the injection peak in a FIA acquisition.
#'
#' Determine a first approximation of the injection peak using the
#' Douglas-Peuker Algorithm provided in the \code{rgeos} package.
#' The object furnished must be an xcmsRaw object.
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
    function(xraw, scanmin=1,
    		 scanmax=length(xraw@scantime),
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
        
        sTime <- scale(xraw@scantime[TIC$scan][scanmin:scanmax], center = FALSE)
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
            segment <- segmentCurve(sTime,sInt,epsilon)
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
                type="l"
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
        message("An injection peak has been spotted:",
            sprintf("%0.1f",xraw@scantime[peaklim[1]+scanmin-1]),
            "-",
            sprintf("%0.1f",xraw@scantime[peaklim[2]+scanmin-1]),
            "s\n")
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
#'   gp <- determiningSizePeak.Geom(xraw)
#' 
#'   #performing band detection.
#'   tbands <- findBandsFIA(xraw,ppm = 2,sizeMin = gp[3]-gp[1],beginning=gp[1])
#' 
#'   #Getting the injection peak
#'   vpeak <- getInjectionPeak(xraw,bandlist=tbands,gpeak=gp)
#'   plot(vpeak,type="l")
#' }

getInjectionPeak <-
    function(xraw,
             bandlist=NULL,
             sec = 2,
             iquant = 0.95,
             gpeak = NULL,
             graphical = FALSE) {
        if (is.null(gpeak)) {
            gpeak <- determiningSizePeak.Geom(xraw)
        }
        
        ###Filtering the peak to retain the peak without oslven
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
        
        
        ###A number of 0 equal to the size of the injection peak is furnished to avoid
        ###Problem fitting the beginning of the peaks
        ###Making a verfication that the retianed profiles does not have oslvent in it
        vok <- which(as.logical(apply(matInt, 2, checkIso, niso = 3)))
        if (length(vok) >= 5) {
            matInt <- matInt[, vok, drop = FALSE]
        }
        
        ###For each line getting the beginning of the injection peak
        getBeginning <- function(x, size = 3) {
            a <- which((x[3:length(x)] != 0) &
                          (x[2:(length(x) - 1)] != 0) &
                          (x[1:(length(x) - 2)] != 0))
            return(a[1])
        }
        vb <- apply(matInt, 2, getBeginning)

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
        matDInt <- matInt[, which(vb >= inter[1] & vb <= inter[2])]
        #print(ttrue[which(vb %in% tokeep),])
        #return(ttrue[which(vb %in% tokeep),])
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
        inita <- 0.06
        initb <- 0.5
        initpar <- c(xraw@scantime[gpeak[3]-floor((gpeak[3]-gpeak[1])/2)],
        			 (gpeak[3] - gpeak[1]) / 2,
        			 20,
        			 rep(inita, n),
        			 rep(initb, n),
        			 vmax)
        #matResiduals<-function(mpp,xx,mobs,type="gaussian",n){
        weigthvec <- length(c(
        	rep(2, gpeak[1]-1),3,
        	seq(1, 0.75, length = (gpeak[3] - gpeak[1])-1),
        	seq(1, 0.75, length = (gpeak[2] - gpeak[3] +
        						   	1)),
        	rep(1, length(xraw@scantime) - gpeak[2])
        ))
        lower <- c(rep(-Inf,3),rep(0,n),rep(-Inf,2*n))
        
        
        nlC <- nls.lm.control(maxiter = 100, ptol = 0.001)
        tMat <- t(t(matDInt) / apply(matDInt, 2, max))
        if(graphical){
        matplot(tMat, type = "l",xlab = "Time (s)",ylab="Scaled intensity")
        }
        parestimate <- nls.lm(
            par = initpar,
            fn = matResiduals,
            observed = tMat,
            type = "gaussian",
            beginning = inter[1],
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
            multiplier[i] <- max(fitted)
            if (graphical) {
                plot(xraw@scantime, tMat[, i],ylim=c(0,max(1,tr*parm[3])), type = "l")
                lines(xraw@scantime, fitted*tr, col = "red")
                lines(xraw@scantime, (mat_eff + 1)*parm[3], col = "blue")
                lines(xraw@scantime, tr*parm[3] , col = "green")
            }
        }
        TP <- persoConv(xraw@scantime, p = parv)
        tMat <- tMat %*% diag(multiplier)
        if(graphical){
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
            main = paste("Regressed injection peak for",
                         getRawName(xraw@filepath))
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
    tr <- a * exp(b * intensity)#+(p$c)*exp(p$d*intensity)) #+p$c*exp(p$d*intensity))/(p$a+p$c)
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
        mmareffect <- sapply(1:ncol(trmat), function(x, va, vb, mat) {
            a <- va[x]
            b <- vb[x]
            - matrix_effect(mat[, x], a, b)
        }, mat = trmat, va = a, vb = b)
        trmatres <- t(t(trmat + mmareffect) * h)
        matdiff <- observed - trmatres
        matdiff <- matdiff * weigth
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
#   gp <- determiningSizePeak.Geom(xraw)
# 
#   #performing band detection.
#   tbands <- findBandsFIA(xraw,ppm = 2,sizeMin = gp[3]-gp[1],beginning=gp[1])
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
# 			gpeak <- determiningSizePeak.Geom(xraw)
# 		}
# 		if(is.null(bandlist)){
# 			message("No bandList furnished")
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
# 		###A number of 0 equal to the size of the injection peak is furnished to avoid
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
			 weigth,noisefunction) {
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
		
		LL<-function(diff){
			R = sum(log10(sapply(diff,dnorm)))
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
