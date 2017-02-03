//
//  AudioProcQueue.swift
//
//  Created by Michael Pilkington on 10/17/16.
//	No restrictions on usage
//
import Foundation
import UIKit

// Global variable to instantiate the singleton
let APQueue = AudioProcQueue()

/** Derive this protocol to get notifications when an song was added to the audio processing queue */
protocol AudioProcQueueDelegate {
    func dataWasAddedToQueue()
}

/**: This singleton class queues up the audio processes requested by the user interface into a serial dispatch queue  */
class AudioProcQueue {
    
    var serialDispatch = dispatch_queue_create("AudioProcessingQueue", DISPATCH_QUEUE_SERIAL)
    var text = ""
    var progress : Float = 0.0
    var processingSongs : [String] = []
    var processingSongsHash : [Int32] = []
    
    var m_delegate : AudioProcQueueDelegate? = nil
    var m_monitordelegate : AudioProcessMonitorDelegate? = nil
    
    class func getInstance() -> AudioProcQueue {
        return APQueue
    }
    
    func addToQueue(nextItem: String, hash: Int32) -> Bool {
        if(processingSongsHash.contains(hash)) {
            return false
        }
        else {
            processingSongsHash.append(hash)
            processingSongs.append(nextItem)
            self.updateDelegates()
            dispatch_async(serialDispatch) {
                _ = AudioProcess(fileName: nextItem, delegate: self.m_monitordelegate)
                dispatch_async(dispatch_get_main_queue(), {
                    self.processingSongs.removeFirst()
                    self.updateDelegates()
                })
            }
            return true
        }
    }
    
    func addToQueue(nextItem: NSURL, name: String, hash: Int32) -> Bool {
        if(processingSongsHash.contains(hash)) {
            return false
        }
        else {
            processingSongsHash.append(hash)
            processingSongs.append(name)
            self.updateDelegates()
            dispatch_async(serialDispatch) {
                _ = AudioProcess(url: nextItem, delegate: self.m_monitordelegate)
                dispatch_async(dispatch_get_main_queue(), {
                    self.processingSongs.removeFirst()
                    self.updateDelegates()
                })
            }
            return true
        }
    }
    
    func songAlreadyBeingProcessed(hash: Int32) -> Bool {
        return processingSongsHash.contains(hash) ? true : false
    }
    func registerWithQueueToReceiveArrayUpdates(delegate: AudioProcQueueDelegate) {
        m_delegate = delegate
    }
    func registerWithQueueForAudioProcessMonitoring(delegate: AudioProcessMonitorDelegate) {
        m_monitordelegate = delegate
    }
    func updateDelegates() {
        if m_delegate != nil {
            m_delegate!.dataWasAddedToQueue()
        }
    }
}