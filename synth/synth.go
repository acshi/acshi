package synth

import (
    "math"
)

type Synthesizer interface {
    // samples slice is the requested buffer to be filled
    // sliceNumber is a unique identifier so that subsequent calls with
    // the same sliceNumber are guarenteed to return the same samples
    // the return value is true if this synthesizer believes the stream of audio is finished.
    Synthesize(samples []float32, sliceNumber int) bool
}

type Constant float32

func (g Constant) Synthesize(samples []float32, sliceNumber int) bool {
    for i := 0; i < len(samples); i++ {
        samples[i] = float32(g)
    }
    return false
}

func Gain(gain float32) Constant {
    return Constant(float32(math.Pow(10, float64(gain) / 20.0)))
}

type impulse struct {
    impulseDone bool
    sliceRecorded int
    recordedSamples []float32
    recodedWasFinished bool
}

func (g *impulse) Synthesize(samples []float32, sliceNumber int) bool {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return g.recodedWasFinished
    }

    isFinished := true

    for i := 0; i < len(samples); i++ {
        samples[i] = 0
    }
    if !g.impulseDone {
        samples[0] = 1
        g.impulseDone = true
        isFinished = false
    }

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
    g.recodedWasFinished = isFinished

    return isFinished
}

func Impulse() Synthesizer {
    return new(impulse)
}

