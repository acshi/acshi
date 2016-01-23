package synth

import (
    "github.com/acshi/functional"
    "math"
)

type FilterCoefficients struct {
    As, Bs []float64
}

func convolve(xs, ys []float64) []float64 {
	zs := make([]float64, len(xs)+len(ys)-1)
	for i := range zs {
		j := 0
		if i >= len(ys) {
			j = i - (len(ys) - 1)
		}
		for ; j <= i && j < len(xs); j++ {
			zs[i] += xs[j] * ys[i-j]
		}
	}
	return zs
}

// Although b0 is not a used coefficient, space for it is still expected
func CascadeFilterStages(stage1, stage2 FilterCoefficients) (stage3 FilterCoefficients) {
    // Transfer function uses the negatives of the b coefficients
    stage1Bs := functional.Map(functional.Negate, stage1.Bs)
    stage2Bs := functional.Map(functional.Negate, stage2.Bs)

    // Transfer function has a 1 in the denominator
    stage1Bs[0], stage2Bs[0] = 1, 1

    // Polynomial multiplication of the numerator and denominator of the transfer functions is done by convolution
    stage3.As = convolve(stage1.As, stage2.As)
    stage3.Bs = convolve(stage1Bs, stage2Bs)

    functional.MapInPlace(functional.Negate, stage3.Bs)
    // just for good measure
    stage3.Bs[0] = 0

    return
}

// Although b0 is not a used coefficient, space for it is still expected
func ParallelizeFilterStages(stage1, stage2 FilterCoefficients) (stage3 FilterCoefficients) {
    // Transfer function uses the negatives of the b coefficients
    stage1Bs := functional.Map(functional.Negate, stage1.Bs)
    stage2Bs := functional.Map(functional.Negate, stage2.Bs)

    // Transfer function has a 1 in the denominator
    stage1Bs[0], stage2Bs[0] = 1, 1

    // For addition of transfer functions, numerator is stagesA1 * stagesB2 + stagesA2 * stagesB1. Denominator is the same as for cascading.
    // Polynomial multiplication of the numerator and denominator of the tranfer functions is done by convolution
    stage3.As = functional.ZipWith(functional.Add, convolve(stage1.As, stage2Bs), convolve(stage2.As, stage1Bs))
    stage3.Bs = convolve(stage1Bs, stage2Bs)

    functional.MapInPlace(functional.Negate, stage3.Bs)
    // just for good measure
    stage3.Bs[0] = 0

    return
}

func NormalizeFilterGain(stage FilterCoefficients, highpass bool) {
    sumA := 0.0
    sumB := 0.0
    for i := range stage.As {
        if highpass {
            if (i % 2) == 0 {
                sumA += stage.As[i]
                sumB += stage.Bs[i]
            } else {
                sumA -= stage.As[i]
                sumB -= stage.Bs[i]
            }
        } else {
            sumA += stage.As[i]
            sumB += stage.Bs[i]
        }
    }

    gain := sumA / (1 - sumB)

    for i := range stage.As {
        stage.As[i] /= gain
    }
}

// Chebyshev biquad (2-poles) recursive coefficients
// Adapted from The Scientist and Engineer's Guide to Digital Signal Processing, Steven W. Smith
// poleIndex = [0, poleCount)
// percentRipple in the pass band can range from 0 for a butterworth to about 0.29. Something like 0.005 is a good trade-off.
func chebyshevBiquad(freq, percentRipple float64, poleIndex, poleCount int, highpass bool) (stage FilterCoefficients) {
    // We start off by designing a low-pass filter with unity cut-off frequency

    // Location of pole on unit circle, real and imaginary parts
    // The maximally flat butterworth filter positions the poles so that
    // they form a semi-circle on the left side of the s-plane (sigma < 0)
    // The half offset keeps the poles evenly spaced and off of the sigma=0 line
    // s-plane s = sigma + i * omega = poleR + i * poleI 
    poleI, poleR := math.Sincos((float64(poleIndex) + 0.5) * math.Pi / float64(poleCount))
    poleR = -poleR

    // The chebyshev filter uses an ellipse to move all of the poles closer to the sigma=0 line
    // This causes pass-band ripple and sharpens the drop off
    // Warp coordinates from being on a circle to an ellipse
    if percentRipple != 0.0 {
        e := math.Sqrt(1 / ((1 - percentRipple) * (1 - percentRipple)) - 1)
        v := math.Asinh(1 / e) / float64(poleCount)
        k := math.Acosh(1 / e) / float64(poleCount)
        
        k = math.Cosh(k)

        poleR = poleR * math.Sinh(v) / k
        poleI = poleI * math.Cosh(v) / k
    }

    // bilinear s-domain to z-domain transformation
    t := 2 * math.Tan(0.5)
    w := 2 * math.Pi * freq
    m := poleR * poleR + poleI * poleI
    d := 4 - 4 * poleR * t + m * t * t
    x0 := t * t / d
    x1 := 2 * t * t / d
    x2 := t * t / d
    y1 := (8 - 2 * m * t * t) / d
    y2 := (-4 - 4 * poleR * t - m * t * t) / d

    // We now have the coefficients of a low-pass filter with a cutoff frequency of 1 (2 times the nyquist)...
    // We must now apply our desired frequency and convert to a high-pass filter if necessary
    // as with the bilinear tranform, these are the results of a substitution in the transfer function...

    var k float64
    if highpass {
        k = -math.Cos(w / 2 + 0.5) / math.Cos(w / 2 - 0.5)
    } else {
        k = math.Sin(0.5 - w / 2) / math.Sin(0.5 + w / 2)
    }

    d = 1 + (y1 * k - y2 * k * k)
    a0 := (x0 - x1 * k + x2 * k * k) / d
    a1 := (-2 * x0 * k + x1 + (x1 * k * k - 2 * x2 * k)) / d
    a2 := (x0 * k * k - x1 * k + x2) / d
    b1 := (2 * k + y1 + y1 * k * k - 2 * y2 * k) / d
    b2 := (-k * k - y1 * k + y2) / d
    if highpass {
        a1, b1 = -a1, -b1
    }

    // we now have the desired coefficients of our low/high pass filter with the desired cutoff frequency
    // however, the gain has not been normalized, if that is desired...
    
    stage.As = []float64{a0, a1, a2}
    stage.Bs = []float64{0, b1, b2}
    return
}

// Chebyshev recursive coefficients
// Adapted from The Scientist and Engineer's Guide to Digital Signal Processing, Steven W. Smith
func chebyshev(freq float64, highpass bool, numPoles int) (coefs []FilterCoefficients) {
    /*numPoles := 8

    freqDist := 0.25 - math.Abs(freq - 0.25)
    if freqDist >= 0.25 {
        if numPoles > 20 {
            numPoles = 20
        }
    } else if freqDist >= 0.10 {
        if numPoles > 10 {
            numPoles = 10
        }
    } else if freqDist >= 0.05 {
        if numPoles > 8 {
            numPoles = 8
        }
    } else if freqDist >= 0.02 {
        if numPoles > 4 {
            numPoles = 4
        }
    } else {
        if numPoles > 4 {
            numPoles = 4
        }
    }*/
    
    n := numPoles / 2
    for i := 0; i < n; i++ {
        stage := chebyshevBiquad(freq, 0.005, i, numPoles, highpass)
        NormalizeFilterGain(stage, highpass)
        coefs = append(coefs, stage)
        /*if i == 0 {
            coefs = []FilterCoefficients{stage}
        } else {
            coefs[0] = CascadeFilterStages(coefs[0], stage)
        }*/
    }
    //NormalizeFilterGain(coefs[0], highpass)
    return
}

type highlowPass struct {
    freq Synthesizer
    filt filter
    highpass bool
    numPoles int
    sliceRecorded int
    recordedSamples []float32
    recordedWasFinished bool
}

func (g *highlowPass) Synthesize(samples []float32, sliceNumber int) bool {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return g.recordedWasFinished
    }
    if len(samples) == 0 {
        // We haven't made anything new, so this is likely still the state
        return g.recordedWasFinished
    }

    // We sample only the very first frequency value.
    // Updating the filter coefficients for every sample would be inefficient and
    // provide unnecessary resolution
    freqIsFinished := g.freq.Synthesize(samples, sliceNumber)
    g.filt.coefs = chebyshev(float64(samples[0]) / samplingRate, g.highpass, g.numPoles)

    inputIsFinished := g.filt.Synthesize(samples, sliceNumber)

    // any one invalidates the whole
    isFinished := freqIsFinished || inputIsFinished

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
    g.recordedWasFinished = isFinished

    return isFinished
}

func LowPass(freq, input Synthesizer, numPoles int) Synthesizer {
    return &highlowPass{highpass: false, numPoles: numPoles, freq: freq, filt: filter{input: input}}
}

func HighPass(freq, input Synthesizer, numPoles int) Synthesizer {
    return &highlowPass{highpass: true, numPoles: numPoles, freq: freq, filt: filter{input: input}}
}

type bandPassStop struct {
    freq Synthesizer
    bandwidth Synthesizer
    filt filter
    bandstop bool
    numPoles int
    sliceRecorded int
    recordedSamples []float32
    recordedWasFinished bool
}

func (g *bandPassStop) Synthesize(samples []float32, sliceNumber int) bool {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return g.recordedWasFinished
    }
    if len(samples) == 0 {
        // We haven't made anything new, so this is likely still the state
        return g.recordedWasFinished
    }

    // We sample only the very first frequency value.
    // Updating the filter coefficients for every sample would be inefficient and
    // provide unnecessary resolution
    freqIsFinished := g.freq.Synthesize(samples, sliceNumber)
    centerFreq := float64(samples[0]) / samplingRate

    bandwidthIsFinished := g.bandwidth.Synthesize(samples, sliceNumber)
    halfBandwidth := float64(samples[0] * 0.5) / samplingRate

    // apparently, the decomposed lowpass and highpass filters must use a different center frequency
    // for the bandpass filter result to have the original desired one
    transformedCenterFreq := math.Sqrt(centerFreq * centerFreq + halfBandwidth * halfBandwidth)

    stageLow := chebyshev(transformedCenterFreq - halfBandwidth, !g.bandstop, g.numPoles) // highpass for bandpass
    stageHigh := chebyshev(transformedCenterFreq + halfBandwidth, g.bandstop, g.numPoles) // lowpass for bandpass

    if g.bandstop {
        for i := range stageLow {
            stageLow[i] = ParallelizeFilterStages(stageLow[i], stageHigh[i])
        }
        g.filt.coefs = stageLow
    } else {
        g.filt.coefs = append(stageLow, stageHigh...)
    }

    inputIsFinished := g.filt.Synthesize(samples, sliceNumber)

    // any one invalidates the whole
    isFinished := freqIsFinished || bandwidthIsFinished || inputIsFinished

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
    g.recordedWasFinished = isFinished

    return isFinished
}

func BandPass(freq, bandwidth, input Synthesizer, numPoles int) Synthesizer {
    return &bandPassStop{bandstop: false, numPoles: numPoles, freq: freq, bandwidth: bandwidth, filt: filter{input: input}}
}

func BandStop(freq, bandwidth, input Synthesizer, numPoles int) Synthesizer {
    halfBW := Amplify(Constant(0.5), bandwidth)
    return Mix(LowPass(Subtract(freq, halfBW), input, numPoles), HighPass(Mix(freq, halfBW), input, numPoles))
    //return &bandPassStop{bandstop: true, numPoles: numPoles, freq: freq, bandwidth: bandwidth, filt: filter{input: input}}
}

type filter struct {
    input Synthesizer
    coefs []FilterCoefficients
    lastAs, lastBs [][]float64
    sliceRecorded int
    recordedSamples []float32
    recordedWasFinished bool
}

func (g *filter) Synthesize(samples []float32, sliceNumber int) bool {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return g.recordedWasFinished
    }

    if len(g.lastAs) < len(g.coefs) {
        g.lastAs = make([][]float64, len(g.coefs))
        g.lastBs = make([][]float64, len(g.coefs))
    }

    for j, coefs := range g.coefs {
        if len(g.lastAs[j]) < len(coefs.As) {
            g.lastAs[j] = make([]float64, len(coefs.As))
            g.lastBs[j] = make([]float64, len(coefs.Bs))
        }
    }
    
    inputIsFinished := g.input.Synthesize(samples, sliceNumber)
    
    for i := 0; i < len(samples); i++ {
        a0 := float64(samples[i])
        b0 := 0.0

        // run through each stage of coefficients
        for j, coefs := range g.coefs {
            lastAs, lastBs := g.lastAs[j], g.lastBs[j]
            
            // shift the old a and b values..
            copy(lastAs[1:], lastAs[:len(lastAs) - 1])
            copy(lastBs[1:], lastBs[:len(lastBs) - 1])
        
            lastAs[0] = a0

            b0 = 0.0
            for i := 0; i < len(coefs.As) && i < len(lastAs); i++ {
                b0 += coefs.As[i] * lastAs[i]
            }
            for i := 1; i < len(coefs.Bs) && i < len(lastBs); i++ {
                b0 += coefs.Bs[i] * lastBs[i]
            }
            
            lastBs[0] = b0
            a0 = b0
        }
        samples[i] = float32(b0)
    }
    
    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
    g.recordedWasFinished = inputIsFinished

    return inputIsFinished
}

