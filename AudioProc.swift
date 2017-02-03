//
//  AudioProc.swift
//
//  Created by Michael Pilkington on 10/17/16.
//	No restrictions on usage
//

import Foundation
import AudioKit
import Surge

let PI : Float = 3.14159265359 // GLOBAL
let NUMBER_OF_CHROMA : Int = 5

/**: Main Audio Processing function
 
 Frameworks: AudioKit, Surge.
 
 - seealso:
 
 The mathmatical principles that I have used in AudioProc are derived from a webpage [1] maintained by Dan Ellis <dpwe@ee.columbia.edu> which hosts matlab files and a link to his article [2] in New Music Research (2007) which gives alot of background to this area. For the resampling algorithm, see Ron Nicholson's quick and dirty FIR resampling algorithm here [3] from [4].
 
 1. Online Source : [Music Audio Tempo Estimation and Beat Tracking](http://labrosa.ee.columbia.edu/projects/beattrack)
 
 2. D. Ellis (2007) Beat Tracking by Dynamic Programming, J. New Music Research, Special Issue on Beat and Tempo Extraction, vol. 36 no. 1, March 2007, pp. 51-60. (10pp) DOI: 10.1080/09298210701653344
 
 3. Online Source : [Ron's Digital Signal Processing Page](http://www.nicholson.com/rhn/dsp.html)
 
 4. Online Source : [CCRMA's Resampling page](https://ccrma.stanford.edu/~jos/resample/resample.html)
 
 */
class AudioProcess {
    
    //! CONSTANTS
    let newSampleRate       : Float = 8000.00
    let windowWidth         : Int = 256
    let windowHop           : Int = 32
    let numberOfMelBins     : Int = 40
    let nyquistFrequency    : Float = 4000.00
    let tightness           : Float = 400
    let oesr                : Float = 8000/32
    
    //! CLASS VARIABLES
    var file                : AKAudioFile?
    var numChannels         : Int?
    var m_delegate          : [AudioProcessMonitorDelegate]?
    var output              : outputData?
    
    //! STRUCTURES
    struct tempoData {
        var tempoSlow   : Double?
        var tempoFast   : Double?
        var pratio      : Double?
    }
    
    struct beatData {
        var beatsInSeconds : [Float]?
        var maximas : [Int]?
    }
    
    struct featureData {
        var featsInSeconds : [Float]?
        var chroma         : [Int]?
    }
    
    struct outputData {
        var tempo               : tempoData?
        var hash                : Int32
        var trackName           : NSString?
        var beatsInSeconds      : [Float]?
        var featuresInSeconds   : [Float]?
        var chroma              : [Int]?
    }
    
    /* Dummy init, TODO: Why do I need this?? */
    init() {}
    
    /* Start AudioProcess with only filename */
    init(fileName: String) {
        file = try? AKAudioFile(readFileName: fileName, baseDir: .Resources )
        if(file == nil) { print("AudioProcess: Could not open file: \(fileName)") }
        else {startProcess(file!)}
    }
    /* Start AudioProcess with filename + single delegate */
    init(fileName: String, delegate: AudioProcessMonitorDelegate?) {
        m_delegate = []
        m_delegate?.append(delegate!)
        file = try? AKAudioFile(readFileName: fileName, baseDir: .Resources )
        if(file == nil) { print("AudioProcess: Could not open file: \(fileName)") }
        else {startProcess(file!)}
    }
    /* Start AudioProcess with filename + group of delegates */
    init(fileName: String, delegate: [AudioProcessMonitorDelegate]) {
        m_delegate = delegate
        file = try? AKAudioFile(readFileName: fileName, baseDir: .Resources )
        if(file == nil) { print("AudioProcess: Could not open file: \(fileName)") }
        else {startProcess(file!)}
    }
    /* Start AudioProcess with URL + single delegate */
    init(url: NSURL, delegate: AudioProcessMonitorDelegate?) {
        m_delegate = []
        m_delegate?.append(delegate!)
        file = try? AKAudioFile(forReading: url)
        if(file == nil) {
        }
        else {startProcess(file!)}
    }
    
    /* Process "main" gets here if file is open and everything is well */
    func startProcess(file: AKAudioFile) {
        
        //! Part 1: onset strength envelope, O(t)
        progressPercChanged(0)
        progressTextChanged("Resampling..")
        let fileData = file.arraysOfFloats[0]
        let fileSampleRate  : Float = Float(file.sampleRate)
        let data = resampleToNewSampleRate(fileData, sampleRate: fileSampleRate)
        
        progressTextChanged("Spectrum...")
        let bins = packDataIntoMelBins(data)
        
        progressTextChanged("Calculus...")
        let env = differentiateAndSumAllMelBins(bins)
        
        //! Part 2: Global tempo estimation which calculates the inter-beat interval, τp
        progressTextChanged("Wizardry...")
        let tempo = performGlobalTempoEstimation(env)
        
        //! Part 3: Beat Tracking
        let beat = peformBeatTrackingAlgorithm(tempo, onsetEnvelope: env)

        progressTextChanged("Finished")
        progressPercChanged(0)

		//! YOUR METHOD FOR SAVING THE DATA HERE
       
    }
    
    /**
     This function resamples the audio file to newSampleRate given an AudioKit AudioFile [3][4]
     */
    func resampleToNewSampleRate(origData: [Float], sampleRate: Float) -> [Float] {
        
        numChannels = 1
        let outLength   : Int = Int((origData.count / sampleRate) * newSampleRate); // Capacity of the output array
        
        var outData: [Float] = Array(count: outLength, repeatedValue: 0)
        outData.reserveCapacity(outLength)
        
        let hann = hanningWindow(windowWidth)
        
        /* Generated Sampling Poisitions */
        var x_arr : [Float] = []
        x_arr.reserveCapacity(outLength)
        for i in 0 ... outLength {
            x_arr.append(Float(i) * (sampleRate/newSampleRate))
        }
        
        /* Generated sinc stuff */
        var F : [Float] = []
        F.reserveCapacity(windowWidth)
        let scalar = (2 * PI * (nyquistFrequency/sampleRate))
        for i in (-windowWidth/2) ..< (windowWidth/2) {
            F.append(sin(Float(i) * scalar)/(Float(i) * scalar))
        }
        for i in 0 ..< windowWidth {
            if(F[i].isNaN) {
                F[i] = 1 // For some reason the peak of the sinc function ends up as NaN
            }
        }
        
        let weights = mul(F,y: hann)
        
        for samples in 0 ..< outLength {
            let x = Int(x_arr[samples])
            var sPoint = Int( x - windowWidth/2 )
            var ePoint = Int( x + windowWidth/2 )
            if(sPoint < 0) {sPoint = 0}
            if(ePoint > origData.count) {ePoint = origData.count - 1}
            let slice = origData[sPoint..<ePoint]
            outData[samples] = summul(slice, y: weights) / (sampleRate / newSampleRate) 
            if(samples % 16384 == 0) {
                progressPercChanged((Float(samples)/Float(outLength))*100)
            }
        }

        return outData
    }
    
    /*
     Perform a fourier transform on the data using a hanning window windowWidth and windowHop advance between frames. This is converted to an 'approximate auditory representation' by mapping the fft spectral bins onto 40 weighted Mel bands. Mel loosely stands for melody and is equal to the range of tones that are expected to be used most often for musical compositions.
     */
    func packDataIntoMelBins(data: [Float]) -> [[Float]] {
        
        let weights = fft2mel()
        let hann = hanningWindow(windowWidth)
        let numberOfSamples = data.count
        let minus80dBs : Float = -78.7440522
        
        var bins : [[Float]] = []
        bins.reserveCapacity(numberOfSamples - windowWidth)
        
        var ftable : [Float] = []
        ftable.reserveCapacity(windowWidth)
        
        var i = 0
        while i < (numberOfSamples - windowWidth) {

            ftable = fft(mul(data[i ..< (i  + windowWidth)], y: hann)) // Fast fourier transform on array slice
            
            var D : [Float] = []
            D.reserveCapacity(numberOfMelBins)
            for j in 0 ..< numberOfMelBins { // From this construct the db-magnitude-mel spectrogram
                var product = summul(weights[j][0..<windowWidth], y: ftable)
                product = 20*log10(max(0.0000000001,product)) // in dBs
                if(product < minus80dBs) {
                    product = minus80dBs
                }
                D.append(product)
            }
            bins.append(D)
            i = i + windowHop
            if(i % 16384 == 0) {
                progressPercChanged((Float(i)/Float(numberOfSamples))*100)
            }
        }
        return bins
    }
    
    /* First order differentiation (dy/dx) along time is calculated for each bin, giving a one dimensional 'onset strength envelope' against time that responds to proportional increase in energy summed across approximately auditory frequency bins.
     */
    func differentiateAndSumAllMelBins(array: [[Float]] ) -> [Float] {
        var decisionWaveform : [Float] = []
        decisionWaveform.reserveCapacity(array.count)
        for i in 0 ..< (array.count - 1) {
            var avg : Float = 0.0
            avg = mean(clip(sub(array[i + 1], y: array[i]), low: 0, high: FLT_MAX))
            decisionWaveform.append(avg)
            if((i%32) == 0) {
                progressPercChanged((Float(i)/Float(array.count))*100)
            }
        }
        
        /* Need to remove DC component and smooth result */
        // a(1)y(n)=b(1)x(n)+b(2)x(n−1)−a(2)y(n−1)
        // filter([1 -1], [1 -.99],mm);
        let x = decisionWaveform
        for j in 1 ..< decisionWaveform.count {
            decisionWaveform[j] = x[j] + (-1*x[j-1]) + (0.99*decisionWaveform[j-1])
        }
        
        return decisionWaveform
    }
    
    /* Global tempo estimation */
    func performGlobalTempoEstimation(array: [Float]) -> tempoData {
        
        let maxd = 60
        let maxt = 120
        let acmax : Int = Int(round(4*oesr));
        let maxcol : Int = min(Int(round(maxt*oesr)),array.count)
        let mincol : Int = max(0, Int(maxcol-round(maxd*oesr)))
        // Only use the 1st 90 sec to estimate global pd
        let xcr = xcorr(Array(array[mincol ..< maxcol]), max: acmax)
        
        // Get the latter half of the xcr
        var rawxcr = Array(xcr[acmax..<xcr.count])
        
        // Rather than strictly using the raw output from the autocorrlation the algorithm uses a gaussian distribution around 120bpm which due to human psychology is the most likely beat time for a given song. This greatly improves the performance of the algorithm!
        var xcrwin : [Float] = []
        xcrwin.reserveCapacity(acmax)
        for i in 0 ..< acmax {
            let beatsPerMilisecond : Float = (60 * oesr) / (Float(i+1) + 0.1)
            let tmean : Float = 240
            let tsd : Float = 1.0
            xcrwin.append(exp ( -0.5 * ( pow( (log(beatsPerMilisecond/tmean)/log(2)/tsd) , 2))))
        }
        rawxcr = mul(rawxcr , y: xcrwin);
        
        // let xpks = localMax(rawxcr)
        // let maximumPeak = max(rawxcr)
        
        var xcr00 = rawxcr; // rawxcr padded with 2 zeros on either end
        xcr00.append(0)
        xcr00.insert(0,atIndex: 0)
        
        // Equation 7: [2] pg.10
        // TPS2(τ) = TPS(τ) + 0.5TPS(2τ) + 0.25TPS(2τ − 1) + 0.25TPS(2τ + 1)
        let xcr2Size = Int(ceil(rawxcr.count/2))
        var xcr2 : [Float] = []
        xcr2.reserveCapacity(xcr2Size)
        for τ in 1 ..< 500 {
            let i = τ-1
            xcr2.append(xcr00[τ])
            xcr2[i] = xcr2[i] + 0.5*xcr00[2*τ]
            xcr2[i] = xcr2[i] + 0.25*xcr00[(2*τ) - 1]
            xcr2[i] = xcr2[i] + 0.25*xcr00[(2*τ) + 1]
        }
        // Equation 8: [2] pg.10
        // TPS3(τ) = TPS(τ) + 0.33TPS(3τ) + 0.33TPS(3τ − 1) + 0.33TPS(3τ + 1)
        let xcr3Size = Int(ceil(rawxcr.count/3))
        var xcr3 : [Float] = []
        xcr3.reserveCapacity(xcr3Size)
        for τ in 1 ..< 334 {
            let i = τ-1
            xcr3.append(xcr00[τ])
            xcr3[i] = xcr3[i] + 0.33*xcr00[3*τ]
            xcr3[i] = xcr3[i] + 0.33*xcr00[(3*τ) - 1]
            xcr3[i] = xcr3[i] + 0.33*xcr00[(3*τ) + 1]
        }
        
        var startpd : Int = 0
        var startpd2 : Int = 0
        
        if max(xcr2) > max(xcr3) {
            startpd = maxIndex(xcr2)
            startpd2 = startpd * 2
        }
        else {
            startpd = maxIndex(xcr3)
            startpd2 = startpd * 3
        }
        
        var tempo = tempoData (
            tempoSlow: 60/((Float(startpd))/oesr),
            tempoFast: 60/((Float(startpd2))/oesr),
            pratio: rawxcr[startpd]/(rawxcr[startpd]+rawxcr[startpd2])
        )
        
        // Reorders if it comes out the wrong way around 🙃
        if(tempo.tempoFast < tempo.tempoSlow) {
            let fast = tempo.tempoFast
            tempo.tempoFast = tempo.tempoSlow
            tempo.tempoSlow = fast
        }
        progressPercChanged(100)
        
        return tempo
    }
    
    func peformBeatTrackingAlgorithm(tempo: tempoData, onsetEnvelope: [Float]) -> beatData {
        
        progressPercChanged(0)
        var startBPM : Float = 0
        if (tempo.pratio > 0.5) {
            startBPM = Float(tempo.tempoFast!) // Fast result from the autocorrelation
        }
        else {
            startBPM = Float(tempo.tempoSlow!) // Slow result from the autocorrelation
        }
        let pd = (60*oesr)/startBPM // numbeats in onsetEnvelope step, probability there is a beat in each sampling interval
        let onsetEnv = onsetEnvelope / std(onsetEnvelope) // vector which maps change in frequency
        
        /* NOTE: Gaussian probability weighting function centered around 120bpm. 
         This is the probability of the tempo of the song taken from experimental evidence */
        var gaussTempoProb : [Float] = []
        for i in Int(-pd) ..< Int(pd) {
            gaussTempoProb.append(-0.5 * pow((Float(i) / (pd/32)), 2))
        }
        gaussTempoProb = exp(gaussTempoProb)
        
        // LocalScore is a smoothed version of the onsetEnv based on the 120bpm window
        var localScore = conv(gaussTempoProb, onsetEnv);
        let sValue = Int(round(gaussTempoProb.count/2))
        let eValue = Int(onsetEnv.count)
        localScore = Array(localScore[sValue..<eValue])
        
        // Setting up some array variables for the next part
        var backLink : [Int] = Array(count: localScore.count, repeatedValue: 0)
        var cumScore : [Float] = Array(count: localScore.count, repeatedValue: 0.0)
        var prange : [Int] = []; prange.reserveCapacity(512)
        var txwt   : [Float] = []; txwt.reserveCapacity(512)
        
        // Filling the arrays
        for i in Int(round(-2*pd)) ..< Int(-round(pd/2)) {
            prange.append(Int(i))
            txwt.append(abs(pow(log(Float(i) / -pd), 2)) * -tightness)
        }
        
        var starting = true
        /* This creates a recursive estimation of each beat interval */
        for i in 0 ..< localScore.count {
            let timeRange_start = Int(round(-2*pd) + i)
            let zpad : Int = max(0, min(-timeRange_start, prange.count))
            var scorecands = txwt
            var index = 0
            for j in zpad ..< prange.count {
                scorecands[index] = scorecands[index] + cumScore[timeRange_start + j]
                index += 1
            }
            cumScore[i] = max(scorecands) + localScore[i]
            if (starting && localScore[i] < (0.01 * max(localScore))) {
                backLink[i] = -1;
            }
            else {
                backLink[i] = timeRange_start + maxIndex(scorecands)
                starting = false
            }
            progressPercChanged( Float(i) / Float(localScore.count))
        }
        
        var cumScoreHasLocalMaxima = localMax(cumScore)
        var maximum : Float = 0.0
        var minimum : Float = 0.0
        for i in 0 ..< cumScore.count - 1 {
            if(cumScoreHasLocalMaxima[i]) {
                if(cumScore[i] > maximum) {
                    maximum = cumScore[i]
                }
                if(cumScore[i] < minimum) {
                    minimum = cumScore[i]
                }
            }
        }
        // Find the best point to end the beat tracking
        var median  = minimum + ((maximum - minimum) / 2)
        var bestEndX : Int = 0
        for var i = cumScore.count - 2; i > 0; i -= 1 {
            if(cumScoreHasLocalMaxima[i] && ( cumScore[i] > 0.5*median )) {
                bestEndX = i
                i = 0;
            }
        }
        
        // This is my part for the feature extraction
        median  = min(localScore) + ((max(localScore) - min(localScore)) / 2)
        var maximas : [Int] = []
        var localScoreHasLocalMaxima = localMax(localScore)
        for i in 1 ..< (localScore.count - 2){
            if(localScoreHasLocalMaxima[i] && ( localScore[i] > 0.25*median )) {
                maximas.append(i)
            }
        }
        
        // This finally extracts the beat times
        var beatTimes : [Int] = [bestEndX]
        while(backLink[beatTimes.last!] > 0) {
            beatTimes.append(backLink[beatTimes.last!])
        }
        
        // Times it by the magic constant oesr that converts back to original sample rate
        var beatsInSeconds : [Float] = []
        for var i = beatTimes.count-1; i > 0; i -= 1 {
            beatsInSeconds.append( Float(beatTimes[i]) / oesr )
        }
        
        return beatData(beatsInSeconds: beatsInSeconds, maximas: maximas)
    }
}


