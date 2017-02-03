//
//  AudioProc_FFTToMel.swift
//
//  Created by Michael Pilkington on 10/17/16.
//	No restrictions on usage
//

import Foundation

extension AudioProcess {
    /*
     Convert frequencies in Hz to mel 'scale'.
     Adapted from matlab script : fft2melmx.m [function z = hz2mel(f,htk)]
     */
    func hzToMel(hz: Float) -> Float {
        
        let f_0    : Float = 0
        let f_sp   : Float = 200/3
        let brkfrq : Float = 1000
        let brkpt  : Float = (brkfrq - f_0)/f_sp // starting mel value for log region
        
        let logstep : Float = Float(exp(log(6.4)/27)) // magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is *almost* the ratio between 1000 Hz and the preceding linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)
        
        var mel : Float = 0
        let linputs = hz < brkpt;
        
        if (linputs == true) {
            mel = (hz - f_0)/f_sp
        }
        else {
            mel = brkpt + (log(hz/brkfrq))/log(logstep);
        }
        
        return mel
        
    }
    
    /*
     Convert values on the mel 'scale' to frequency in Hz
     Adapted from matlab script : fft2melmx.m [function z = mel2hz(f,htk)]
     */
    func melToHz(mel: Float) -> Float {
        
        let f_0    : Float = 0
        let f_sp   : Float = 200/3
        let brkfrq : Float = 1000
        let brkpt  : Float = (brkfrq - f_0)/f_sp // starting mel value for log region
        
        let logstep : Float = Float(exp(log(6.4)/27)) // magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is the ratio between 1000 Hz and the linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)
        
        var f : Float = 0
        let linputs = mel < brkpt;
        
        if (linputs == true) {
            f = f_0 + f_sp * mel
        }
        else {
            f = brkfrq*exp(log(logstep)*(mel-brkpt))
        }
        return f
    }
    /*
     Returns an array which contains the weightings of the fast fourier transform that are required to bin pack the data into bins seperated by the mel scale.
     Adapted from matlab script: fft2melmx.m [function z = fft2melmx()]
     */
    func fft2mel() -> [[Float]] {
        
        let minfrq : Float = 0
        let maxfrq : Float = newSampleRate/2
        
        // Center freqs of each FFT bin
        var fftfrqs : [Float] = []
        fftfrqs.reserveCapacity(windowWidth)
        let diff : Float = (1/Float(windowWidth)) * newSampleRate
        for i in 0 ..< windowWidth {
            fftfrqs.append(Float(i) * diff)
        }
        
        let minmel = hzToMel(minfrq);
        let maxmel = hzToMel(maxfrq);
        
        var binfreqs : [Float] = []
        binfreqs.reserveCapacity(numberOfMelBins+1)
        for i in 0 ... numberOfMelBins+1 {
            binfreqs.append(melToHz( minmel + ( ( Float(i)/(Float(numberOfMelBins)+1) ) * (maxmel - minmel) ) ))
        }
        
        var weights : [[Float]] = []
        weights.reserveCapacity(numberOfMelBins)
        
        for i in 0 ..< numberOfMelBins {
            
            var fs : [Float] = [binfreqs[i], binfreqs[i+1], binfreqs[i+2]]
            var lowSlope : Float = 0
            var highSlope : Float = 0
            var weightRow : [Float] = []
            weightRow.reserveCapacity(windowWidth)
            
            for j in 0 ..< (fftfrqs.count) {
                lowSlope = (fftfrqs[j] - fs[0]) / (fs[1] - fs[0] )
                highSlope = (fs[2] - fftfrqs[j]) /  (fs[2] - fs[1] )
                weightRow.append((2 / (fs[2] - fs[0])) * max(0, min(lowSlope, highSlope) ))
            }
            weights.append(weightRow)
        }
        return weights
    }
}
