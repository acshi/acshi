package main

import (
    "fmt"
    "math"
)

type autoTuneData struct {
    initialized bool
    startP, startI float64 // as yet unmodified
    startScore float64 // unmodified settings score
    newP, newI, newD float64 // current values to be used by controller
    roundNum int // the tuning trial/round we are on
    roundPart int // are we trialing a new p value (or i)
    deltaP, deltaI float64 // delta values that may be used in the current round
    deltaDeltaFactor float64 // multplies deltas each round to hone in on optimal values (< 1.0)
    scoreP, scoreI, combinedScore float64 // scores for the trials changing values of p and i
    inTrial bool
    readyForTrial bool // a trial is eligable to start if the threshold it met 
    endTrialCounts int // The number of consecutive calls the end trial requirement is met
    courseErrors []float64
    needsChange bool // awaiting a turn or change before another trial can be started
}

const (
    ROUND_BASLINE_1 = iota
    ROUND_BASLINE_2 = iota
    ROUND_PART_P_1 = iota
    ROUND_PART_P_2 = iota
    ROUND_PART_I_1 = iota
    ROUND_PART_I_2 = iota
    ROUND_COMBINED_1 = iota
    ROUND_COMBINED_2 = iota
)

// Inspiration taken from http://www.mstarlabs.com/control/self-tuning-pid.html
// We define a trial as a period of measurements that starts when the boat at a certain
// course error (it was just told to turn and then reached the designated error).
// It ends when the boat is told to turn again (larger course error).
// The tuning is done in rounds. Each round consists of three trials. The first modifying
// the p value only. The second the i value only. And then third moving both p and i according
// to the measured gradient. At the end of the round, the best of the three trials is taken.
// The next round will use smaller delta values to try and will direct the p and i values in the
// direction of downward optimization.
//
// data.deltaP and data.deltaI, the initial delta values to attempt and work done from need to be set
// before calling. The data.deltaDeltaFactor as well, which is the < 1.0 factor multipled each round to
// the deltaP and deltaI values. 
// returns true when new pid values should be read from autoTune.newP, .newI, .newD (start of new trial)
func autoTunePid(autoTune *autoTuneData, courseError, p, i float64) (newPidValues bool) {
    data := *autoTune
    newPidValues = false
    
    idGainRatio := 0.4 // ratio of d over i for output
    // if the courseError tops this value, a new trial is eligable to start when it reaches the start level
    retrialThreshold := 40.0
    // We will start taking data for a new trial when the value drops to this value...
    startTrialLevel := 40.0
    // To end the trial, the course error must be below this level this many consecutive times
    // If the retrial threshold is reached before this end trial requirement is met, the trial is thrown out.
    endTrialLevel := 5.0
    endTrialMinTimes := 50
    
    if !data.initialized {
        data.initialized = true
        data.startP = p
        data.startI = i
        data.roundNum = 1
        data.inTrial = false
        data.roundPart = ROUND_BASLINE_1
        data.courseErrors = make([]float64, 20)
        data.needsChange = true
    }
    
    if !data.inTrial && data.readyForTrial && math.Abs(courseError) <= startTrialLevel {
        data.inTrial = true
        data.readyForTrial = false
        data.needsChange = false
    } else if !data.inTrial && !data.readyForTrial && math.Abs(courseError) >= retrialThreshold {
        data.readyForTrial = true
        data.needsChange = false
        
        switch data.roundPart {
        case ROUND_BASLINE_1, ROUND_BASLINE_2:
            data.newP = data.startP
            data.newI = data.startI
        case ROUND_PART_P_1, ROUND_PART_P_2:
            data.newP = math.Max(data.startP + data.deltaP, 0.0)
            data.newI = data.startI
        case ROUND_PART_I_1, ROUND_PART_I_2:
            data.newP = data.startP
            data.newI = math.Max(data.startI + data.deltaI, 0.0)
        case ROUND_COMBINED_1, ROUND_COMBINED_2:
            if data.scoreP < data.startScore {
                data.newP = math.Max(data.startP + data.deltaP, 0.0)
            } else {
                data.newP = math.Max(data.startP - data.deltaP, 0.0)
            }
            if data.scoreI < data.startScore {
                data.newI = math.Max(data.startI + data.deltaI, 0.0)
            } else {
                data.newI = math.Max(data.startI - data.deltaI, 0.0)
            }
        }
        data.newD = data.newI * idGainRatio
        newPidValues = true
        
        fmt.Printf("Trial p: %.3f i: %.3f\n", data.newP, data.newI)
        
        data.courseErrors = data.courseErrors[:0]
    } else if data.inTrial && data.endTrialCounts > endTrialMinTimes {
        data.inTrial = false
        data.endTrialCounts = 0
        data.needsChange = true
        /*score := 0.0
        for _, val := range data.courseErrors {
            score += val * val
            fmt.Printf("%.2f ", val)
        }
        score /= float64(len(data.courseErrors))*/
        score := float64(len(data.courseErrors))
        
        fmt.Printf("Trial score: %.2f\n\n", score)
        
        switch data.roundPart {
        case ROUND_BASLINE_1:
            data.startScore = score
            data.roundPart = ROUND_BASLINE_2
        case ROUND_BASLINE_2:
            data.startScore += score
            data.roundPart = ROUND_PART_P_1
        case ROUND_PART_P_1:
            data.scoreP = score
            data.roundPart = ROUND_PART_P_2
        case ROUND_PART_P_2:
            data.scoreP += score
            data.roundPart = ROUND_PART_I_1
        case ROUND_PART_I_1:
            data.scoreI = score
            data.roundPart = ROUND_PART_I_2
        case ROUND_PART_I_2:
            data.scoreI += score
            data.roundPart = ROUND_COMBINED_1
        case ROUND_COMBINED_1:
            data.combinedScore = score
            data.roundPart = ROUND_COMBINED_2
        case ROUND_COMBINED_2:
            data.combinedScore += score
            
            // Choose the best of the options tried
            minScore := math.Min(math.Min(math.Min(data.startScore, data.scoreP), data.scoreI), data.combinedScore)
            var nextP, nextI, nextScore float64
            switch minScore {
            case data.startScore:
                nextP = data.startP
                nextI = data.startI
                nextScore = data.startScore
            case data.scoreP:
                nextP = math.Max(data.startP + data.deltaP, 0.0)
                nextI = data.startI
                nextScore = data.scoreP
            case data.scoreI:
                nextP = data.startP
                nextI = math.Max(data.startI + data.deltaI, 0.0)
                nextScore = data.scoreI
            case data.combinedScore:
                nextP = p
                nextI = i
                nextScore = data.combinedScore
            }
            
            fmt.Printf("Round settled on p: %.3f i: %.3f score: %.2f\n", nextP, nextI, nextScore)
            
            // this complete round is over
            data.roundNum++
            // Modify sign of deltas if going the wrong way
            if data.scoreP > data.startScore {
                data.deltaP = -data.deltaP
            }
            if data.scoreI > data.startScore {
                data.deltaI = -data.deltaI
            }
            // and decrease the delta values
            data.deltaP *= data.deltaDeltaFactor
            data.deltaI *= data.deltaDeltaFactor
            
            // assign new start values
            data.startP = nextP
            data.startI = nextI
            data.startScore = nextScore
            
            data.roundPart = ROUND_BASLINE_1
        }
    } else if data.inTrial {
        data.courseErrors = append(data.courseErrors, courseError)
        if math.Abs(courseError) <= endTrialLevel {
            data.endTrialCounts++
        } else {
            data.endTrialCounts = 0
        }
    }
    
    *autoTune = data
    return
}
