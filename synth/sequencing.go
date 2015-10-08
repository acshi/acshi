package synth

import "math/rand"

type SynthTimeSamples struct {
    input Synthesizer
    timeSamples int
}

func SynthTime(input Synthesizer, time float32) SynthTimeSamples {
    return SynthTimeSamples{input, int(time * samplingRate)}
}

type sequence struct {
    steps []SynthTimeSamples
    stepOn int
    stepTimeSamplesPassed int
    
    sliceRecorded int
    recordedSamples []float32
    recordedWasFinished bool
}

func Sequence(steps ...SynthTimeSamples) Synthesizer {
    return &sequence{steps: steps}
}

func (g *sequence) Synthesize (samples []float32, sliceNumber int) bool {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return g.recordedWasFinished
    }

    isFinished := true

    sampleOn := 0
    for sampleOn < len(samples) {
        if g.stepOn >= len(g.steps) {
            for i := sampleOn; i < len(samples); i++ {
                samples[i] = 0.0
            }
            break
        }

        // if we have anything to synthesize, we are not finished
        isFinished = false
        
        samplesNeeded := len(samples) - sampleOn
        stepSamplesLeft := g.steps[g.stepOn].timeSamples - g.stepTimeSamplesPassed
        
        samplesToGet := samplesNeeded
        if stepSamplesLeft < samplesToGet {
            samplesToGet = stepSamplesLeft
        }

        // Generate new slice numbers for these sub slices, as we have more than one and they must be differentiated
        g.steps[g.stepOn].input.Synthesize(samples[sampleOn : sampleOn + samplesToGet], rand.Int())
        sampleOn += samplesToGet
        stepSamplesLeft -= samplesToGet
        
        if stepSamplesLeft == 0 {
            g.stepOn++
            g.stepTimeSamplesPassed = 0
        } else {
            g.stepTimeSamplesPassed += samplesToGet
        }
    }

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
    g.recordedWasFinished = isFinished

    return isFinished
}

type transition struct {
    from Synthesizer
    to Synthesizer
    steadySamples int
    transitionSamples int
    samplesPassed int
    
    sliceRecorded int
    recordedSamples []float32
    recordedWasFinished bool
}

func lerp(from, to, progress float32) float32 {
    return from + progress * (to - from)
}

func (g *transition) Synthesize (samples []float32, sliceNumber int) bool {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return g.recordedWasFinished
    }

    g.from.Synthesize(samples, sliceNumber)

    toSamples := make([]float32, len(samples))
    g.to.Synthesize(toSamples, sliceNumber)

    // Leave from samples untouched for as long as the steady state time
    skipSamples := 0
    if g.samplesPassed < g.steadySamples {
        skipSamples = g.steadySamples - g.samplesPassed
        if skipSamples > len(samples) {
            skipSamples = len(samples)
        }
        g.samplesPassed += skipSamples
    }

    nonSteadySamplesPassed := g.samplesPassed - g.steadySamples

    currentProgress := float32(nonSteadySamplesPassed) / float32(g.transitionSamples)
    progressPerSample := 1.0 / float32(g.transitionSamples)

    samplesTransitionLeft := g.transitionSamples - nonSteadySamplesPassed

    isFinished := samplesTransitionLeft == 0

    i := skipSamples
    for ; i < len(samples) && i < samplesTransitionLeft; i++ {
        samples[i] = lerp(samples[i], toSamples[i], currentProgress)
        currentProgress += progressPerSample
    }
    g.samplesPassed += i - skipSamples

    // Copy the remaining fully transitioned, if any
    copy(samples[i:], toSamples[i:])

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
    g.recordedWasFinished = isFinished

    return isFinished
}

func Transition(from, to Synthesizer, steadyStateSeconds, transitionSeconds float32) Synthesizer {
    return &transition{from: from, to: to,
                       steadySamples: int(steadyStateSeconds * float32(samplingRate)),
                       transitionSamples: int(transitionSeconds * float32(samplingRate))}
}
