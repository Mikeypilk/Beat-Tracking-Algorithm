//
//  AudioProc_Utils.swift
//
//  Created by Michael Pilkington on 10/17/16.
//	No restrictions on usage
//
import Foundation
import Surge
import AudioKit
import Security
import Accelerate

/* This is a quickie to speed up summing and then multiplying */
public func summul(x: ArraySlice<Float>, y: [Float]) -> Float {
    var results = [Float](count: x.count, repeatedValue: 0.0)
    x.withUnsafeBufferPointer { (xPointer: UnsafeBufferPointer) -> Void in
        vDSP_vmul(xPointer.baseAddress, 1, y, 1, &results, 1, vDSP_Length(x.count))
    }
    var product: Float = 0.0
    vDSP_sve(results, 1, &product, vDSP_Length(results.count))
    return product
}

protocol AudioProcessMonitorDelegate  {
    var progressText : String { get set}
    var percentage : Float { get set }
}

/* This creates a unique hash for the audio file */
func hashAudioFile(file : NSURL) -> Int32 {

    var hashString = ""
    let playerItem = AVPlayerItem(URL: file)
    let metadataList = playerItem.asset.commonMetadata
    for item in metadataList {
        if ((item.value as? String) != nil) {
            hashString.appendContentsOf(item.value as! String)
        }
    }
    return hashString.utf8
        .map {return $0}
        .reduce(5381) {
        ($0 << 5) &+ $0 &+ Int32($1)
    }
 
}

// Miscellaneous Utility functions
extension AudioProcess {
    
    /* Text of what is happening is returned to the view delegate */
    func progressTextChanged(text: String) {
        if(m_delegate == nil) { return }
        dispatch_async(dispatch_get_main_queue()) {
            for i in 0 ..< self.m_delegate!.count {
               self.m_delegate![i].progressText = text
            }
        }
    }
    /* Percentage of function progress returned to the view delegate */
    func progressPercChanged(perc: Float) {
        if(m_delegate == nil) { return }
        dispatch_async(dispatch_get_main_queue()) {
            for i in 0 ..< self.m_delegate!.count {
                self.m_delegate![i].percentage = perc
            }
        }
    }
    
    /* return 1 where there are local maxima in x (columnwise).
     don't include first point, maybe last point */
    func localMax(array: [Float]) -> [Bool] {
        var x : [Bool] = []
        x.reserveCapacity(array.count)
        x.append(false)
        for i in 1 ..< array.count-1 {
            if ((array[i-1] < array[i]) && (array[i+1] < array[i])) {
                x.append(true)
            }
            else {
                x.append(false)
            }
        }
        // will not include any peaks in first down slope (before goes below
        // zero for the first time) I THINK THIS IS WRONG
        /*
        for i in 1 ..< array.count-1 {
            if(array[i] > 0) {
                x[i] = false;
            }
            else {
                break
            }
        }
        */
        return x
    }
    
    /* My function to generate a hanning window of a certain width */
    func hanningWindow(width: Int) -> [Float] {
        let windowSize = Float(width)
        var window : [Float] = []
        window.reserveCapacity(width)
        for j in 1 ... width {
            window.append(0.5 - 0.5 * cos(2*PI*(Float(j)/windowSize)))
        }
        return window
    }
    
    /* returns the index of the maximum value of an array */
    func maxIndex(array: [Float]) -> Int {
        let maxVal = max(array)
        for i in 0 ..< array.count {
            if (array[i] == maxVal) {
                return i
            }
        }
        return 0
    }
    
    /* Returns the standard distribution of a floating point array */
    func std(x: [Float]) -> Float {
        var std : Float = 0.0
        var mean : Float = 0.0
        vDSP_normalize(x, 1, nil, 1, &mean, &std, vDSP_Length(x.count))
        return std
    }
    
}

