package main

import (
    "acshi/audio"
    "acshi/functional"
    "fmt"
    "math"
    "math/rand"
    "time"
    "os"
)

const samplingRate = 44100

var (
    // The wavetables will be constructed to not aliase at frequencies below these
    // Assumed to be in ascending order
    waveTableMaxFrequencies = [...]int{220, 440, 880, 1760, 3520, 7040, 14080}
    sineWaveTable = make([][]float32, len(waveTableMaxFrequencies))
    sawtoothWaveTable = make([][]float32, len(waveTableMaxFrequencies))
    triangleWaveTable = make([][]float32, len(waveTableMaxFrequencies))
    squareWaveTable = make([][]float32, len(waveTableMaxFrequencies))
    noiseTable []float32
)

func init() {
    //samplingRate := 44100
    waveTableFrequency := 20
    waveTableSize := samplingRate / waveTableFrequency
    
    // This frequency translates to 20hz in the table
    angularFreq := 2.0 * math.Pi / float64(waveTableSize)

    // Maximum number of harmonics to avoid aliasing
    harmonicLimit := make([]int, len(waveTableMaxFrequencies))
    for i := 0; i < len(waveTableMaxFrequencies); i++ {
        harmonicLimit[i] = samplingRate / 2 / waveTableMaxFrequencies[i]
    }

    sineWaveTable[0] = make([]float32, waveTableSize)
    for i := 0; i < waveTableSize; i++ {
        sineWaveTable[0][i] = float32(math.Sin(angularFreq * float64(i)))
    }
    // All the same for sine wave...
    for i := 1; i < len(waveTableMaxFrequencies); i++ {
        sineWaveTable[i] = sineWaveTable[0]
    }

    for i := 0; i < len(harmonicLimit); i++ {
        sawtoothWaveTable[i] = make([]float32, waveTableSize)
    }
    for i := 0; i < waveTableSize; i++ {
        value := 0.0
        for j := 1; j <= harmonicLimit[0]; j++ {
            value -= 2 * math.Sin(angularFreq * float64(j) * float64(i)) / (math.Pi * float64(j))
            
            for k := 0; k < len(harmonicLimit); k++ {
                if j == harmonicLimit[k] {
                    sawtoothWaveTable[k][i] = float32(value * 0.85)
                }
            }
        }
    }

    for i := 0; i < len(harmonicLimit); i++ {
        triangleWaveTable[i] = make([]float32, waveTableSize)
    }
    for i := 0; i < waveTableSize; i++ {
        value := 0.0
        for j := 1; j <= harmonicLimit[0]; j += 2 {
            term := 8 * math.Sin(angularFreq * float64(j) * float64(i)) / (math.Pi * math.Pi * float64(j * j))
            if (j / 2) % 2 == 0 {
                value += term
            } else {
                value -= term
            }

            for k := 0; k < len(harmonicLimit); k++ {
                if j / 2 == harmonicLimit[k] / 2 {
                    triangleWaveTable[k][i] = float32(value)
                }
            }
        }
    }

    for i := 0; i < len(harmonicLimit); i++ {
        squareWaveTable[i] = make([]float32, waveTableSize)
    }
    for i := 0; i < waveTableSize; i++ {
        value := 0.0
        for j := 1; j <= harmonicLimit[0]; j += 2 {
            value += 4 * math.Sin(angularFreq * float64(j) * float64(i)) / (math.Pi * float64(j))

            for k := 0; k < len(harmonicLimit); k++ {
                if j / 2 == harmonicLimit[k] / 2 {
                    triangleWaveTable[k][i] = float32(value * 0.84)
                }
            }
        }
    }

    noiseTableSize := samplingRate * 2
    noiseTable = make([]float32, noiseTableSize)
    r := rand.New(rand.NewSource(time.Now().UnixNano()))
    for i := 0; i < noiseTableSize; i++ {
        noiseTable[i] = r.Float32() * 2 - 1
    }
    
    fmt.Println("We have inited!")
}

type Synthesizer interface {
    Synthesize(samples []float32, sliceNumber int)
}

type Constant float32

func (g Constant) Synthesize(samples []float32, sliceNumber int) {
    for i := 0; i < len(samples); i++ {
        samples[i] = float32(g)
    }
}

func Gain(gain float32) Constant {
    return Constant(float32(math.Pow(10, float64(gain) / 20.0)))
}

type impulse struct {
    impulseDone bool
    sliceRecorded int
    recordedSamples []float32
}

func (g *impulse) Synthesize(samples []float32, sliceNumber int) {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return
    }

    for i := 0; i < len(samples); i++ {
        samples[i] = 0
    }
    if !g.impulseDone {
        samples[0] = 1
        g.impulseDone = true
    }

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
}

func Impulse() Synthesizer {
    return new(impulse)
}

type amplify struct {
    inputs []Synthesizer
    sliceRecorded int
    recordedSamples []float32
}

func (g *amplify) Synthesize(samples []float32, sliceNumber int) {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return
    }

    inputSamples := make([][]float32, len(g.inputs))
    for i := 0; i < len(g.inputs); i++ {
        inputSamples[i] = make([]float32, len(samples))
        g.inputs[i].Synthesize(inputSamples[i], sliceNumber)
    }

    for i := 0; i < len(samples); i++ {
        value := float32(1.0)
        for j := 0; j < len(g.inputs); j++ {
            value *= inputSamples[j][i]
        }
        samples[i] = value
    }

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
}

func Amplify(inputs ...Synthesizer) Synthesizer {
    return &amplify{inputs: inputs}
}

type mix struct {
    inputs []Synthesizer
    sliceRecorded int
    recordedSamples []float32
}

func (g *mix) Synthesize(samples []float32, sliceNumber int) {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return
    }

    inputSamples := make([][]float32, len(g.inputs))
    for i := 0; i < len(g.inputs); i++ {
        inputSamples[i] = make([]float32, len(samples))
        g.inputs[i].Synthesize(inputSamples[i], sliceNumber)
    }

    for i := 0; i < len(samples); i++ {
        value := float32(0.0)
        for j := 0; j < len(g.inputs); j++ {
            value += inputSamples[j][i]
        }
        samples[i] = value
    }

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
}

func Mix(inputs ...Synthesizer) Synthesizer {
    return &mix{inputs: inputs}
}

type subtract struct {
    signal Synthesizer
    inputs []Synthesizer
    sliceRecorded int
    recordedSamples []float32
}

func (g *subtract) Synthesize(samples []float32, sliceNumber int) {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return
    }

    g.signal.Synthesize(samples, sliceNumber)

    inputSamples := make([][]float32, len(g.inputs))
    for i := 0; i < len(g.inputs); i++ {
        inputSamples[i] = make([]float32, len(samples))
        g.inputs[i].Synthesize(inputSamples[i], sliceNumber)
    }

    for i := 0; i < len(samples); i++ {
        for j := 0; j < len(g.inputs); j++ {
            samples[i] -= inputSamples[j][i]
        }
    }

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
}

func Subtract(signal Synthesizer, inputs ...Synthesizer) Synthesizer {
    return &subtract{signal: signal, inputs: inputs}
}

type filter struct {
    freq Synthesizer
    input Synthesizer
    a1, a2, a3, a4 float64
    b1, b2, b3, b4 float64
    sliceRecorded int
    recordedSamples []float32
}

func (g *filter) filterSynthesize(samples []float32, highPass bool, sliceNumber int) {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return
    }
    
    freqs := make([]float32, len(samples))
    g.freq.Synthesize(freqs, sliceNumber)
    g.input.Synthesize(samples, sliceNumber)
    
    for i := 0; i < len(samples); i++ {
        as, bs := chebyshev(float64(freqs[i] / 44100.0), highPass)
        a0 := float64(samples[i])
        
        b0 := as[0] * a0 +
              as[1] * g.a1 + as[2] * g.a2 + as[3] * g.a3 + as[4] * g.a4 +
              bs[1] * g.b1 + bs[2] * g.b2 + bs[3] * g.b3 + bs[4] * g.b4

        g.a1, g.a2, g.a3, g.a4 = a0, g.a1, g.a2, g.a3
        g.b1, g.b2, g.b3, g.b4 = b0, g.b1, g.b2, g.b3
        samples[i] = float32(b0)
    }
    
    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
}

type lowPass filter

func (g *lowPass) Synthesize(samples []float32, sliceNumber int) {
    (*filter)(g).filterSynthesize(samples, false, sliceNumber)
}

func LowPass(freq, input Synthesizer) Synthesizer {
    return (*lowPass)(&filter{freq: freq, input: input})
}

type highPass filter

func (g *highPass) Synthesize(samples []float32, sliceNumber int) {
    (*filter)(g).filterSynthesize(samples, true, sliceNumber)
}

func HighPass(freq, input Synthesizer) Synthesizer {
    return (*highPass)(&filter{freq: freq, input: input})
}

func BandPass(freq, bandWidth, input Synthesizer) Synthesizer {
    halfBW := Amplify(bandWidth, Constant(0.5))
    return LowPass(Mix(freq, halfBW), HighPass(Subtract(freq, halfBW), input))
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
func cascadeFilterStages(stageAs1, stageBs1, stageAs2, stageBs2 []float64) (stageAs3, stageBs3 []float64) {
    // Transfer function uses the negatives of the b coefficients
    stageBs1 = functional.Map(functional.Negate, stageBs1)
    stageBs2 = functional.Map(functional.Negate, stageBs2)

    // Transfer function has a 1 in the denominator
    stageBs1[0], stageBs2[0] = 1, 1

    // Polynomial multiplication of the numerator and denominator of the transfer functions is done by convolution
    stageAs3 = convolve(stageAs1, stageAs2)
    stageBs3 = convolve(stageBs1, stageBs2)

    functional.MapInPlace(functional.Negate, stageBs3)
    // just for good measure
    stageBs3[0] = 0

    return
}

// Although b0 is not a used coefficient, space for it is still expected
func parallelizeFilterStages(stageAs1, stageBs1, stageAs2, stageBs2 []float64) (stageAs3, stageBs3 []float64) {
    // Transfer function uses the negatives of the b coefficients
    stageBs1 = functional.Map(functional.Negate, stageBs1)
    stageBs2 = functional.Map(functional.Negate, stageBs2)

    // Transfer function has a 1 in the denominator
    stageBs1[0], stageBs2[0] = 1, 1

    // For addition of transfer functions, numerator is stagesA1 * stagesB2 + stagesA2 * stagesB1. Denominator is the same as for cascading.
    // Polynomial multiplication of the numerator and denominator of the tranfer functions is done by convolution
    stageAs3 = functional.ZipWith(functional.Add, convolve(stageAs1, stageBs2), convolve(stageAs2, stageBs1))
    stageBs3 = convolve(stageBs1, stageBs2)

    functional.MapInPlace(functional.Negate, stageBs3)
    // just for good measure
    stageBs3[0] = 0

    return
}

// Chebyshev biquad (2-poles) recursive coefficients
// Adapted from The Scientist and Engineer's Guide to Digital Signal Processing, Steven W. Smith
// poleIndex = [0, poleCount)
// percentRipple in the pass band can range from 0 for a butterworth to about 0.29. Something like 0.005 is a good trade-off.
func chebyshevBiquad(freq, percentRipple float64, poleIndex, poleCount int, highpass bool) (stageAs, stageBs []float64) {
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
    
    stageAs = []float64{a0, a1, a2}
    stageBs = []float64{0, b1, b2}
    return
}

func normalizeStageGain(stageAs, stageBs []float64, highpass bool) {
    sumA := 0.0
    sumB := 0.0
    for i := range stageAs {
        if highpass {
            if (i % 2) == 0 {
                sumA += stageAs[i]
                sumB += stageBs[i]
            } else {
                sumA -= stageAs[i]
                sumB -= stageBs[i]
            }
        } else {
            sumA += stageAs[i]
            sumB += stageBs[i]
        }
    }

    gain := sumA / (1 - sumB)

    for i := range stageAs {
        stageAs[i] /= gain
    }
}

// Chebyshev recursive coefficients
// Adapted from The Scientist and Engineer's Guide to Digital Signal Processing, Steven W. Smith
func chebyshev(freq float64, highpass bool) (as, bs []float64) {
    /*a = make([]float64, 7)
    b = make([]float64, 7)
    ta := make([]float64, 7)
    tb := make([]float64, 7)

    a[2], b[2] = 1, 1

    for p := 1; p <= 2; p++ {
        a0, a1, a2, b1, b2 := chebyshevBiquad(freq, 0.005, p, 4, highpass)

        for i := 0; i < 7; i++ {
            ta[i] = a[i]
            tb[i] = b[i]
        }

        for i := 2; i < 7; i++ {
            a[i] = a0 * ta[i] + a1 * ta[i - 1] + a2 * ta[i - 2]
            b[i] = tb[i] - b1 * tb[i - 1] - b2 * tb[i - 2]
        }
    }

    b[2] = 0
    for i := 0; i < 5; i++ {
        a[i] = a[i + 2]
        b[i] = -b[i + 2]
    }

    return a[:5], b[:5]*/

    numPoles := 4
    n := numPoles / 2
    for i := 0; i < n; i++ {
        stageAs, stageBs := chebyshevBiquad(freq, 0.005, i, numPoles, highpass)
        if i == 0 {
            as, bs = stageAs, stageBs
        } else {
            as, bs = cascadeFilterStages(as, bs, stageAs, stageBs)
        }
    }

    return
}

func Sine(freq Synthesizer) Synthesizer {
    return &WaveGen{waveTable: &sineWaveTable, input: freq}
}

func Sawtooth(freq Synthesizer) Synthesizer {
    return &WaveGen{waveTable: &sawtoothWaveTable, input: freq}
}

func Triangle(freq Synthesizer) Synthesizer {
    return &WaveGen{waveTable: &triangleWaveTable, input: freq}
}

func Square(freq Synthesizer) Synthesizer {
    return &WaveGen{waveTable: &squareWaveTable, input: freq}
}

type noise struct {
    index int
    sliceRecorded int
    recordedSamples []float32
}

func (g *noise) Synthesize(samples []float32, sliceNumber int) {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return
    }

    copied := 0
    for {
        copied += copy(samples[copied:], noiseTable[g.index:])
        g.index += copied
        if g.index >= len(noiseTable) {
            g.index = g.index % len(noiseTable)
        } else {
            break
        }
    }

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
}

func Noise() Synthesizer {
    return &noise{}
}

type WaveGen struct {
    phase float32
    waveTable *[][]float32
    input Synthesizer
    sliceRecorded int
    recordedSamples []float32
}

func (g *WaveGen) Synthesize(samples []float32, sliceNumber int) {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return
    }

    waveTable := *g.waveTable
    g.input.Synthesize(samples, sliceNumber)
    for i := 0 ; i < len(samples); i++ {
        var specificWaveTable []float32
        for j := 0; j < len(waveTableMaxFrequencies); j++ {
            if samples[i] < float32(waveTableMaxFrequencies[j]) || j == len(waveTableMaxFrequencies) - 1 {
                specificWaveTable = waveTable[j]
                break
            }
        }
        
        g.phase += samples[i] / 20.0
        if int(g.phase) >= len(specificWaveTable) {
            g.phase -= float32((int(g.phase) / len(specificWaveTable)) * len(specificWaveTable))
        } else if int(g.phase) < 0 {
            g.phase += float32((int(-g.phase) / len(specificWaveTable) + 1) * len(specificWaveTable))
        }
        samples[i] = specificWaveTable[int(g.phase)]
    }

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
}

type Vowel struct {
    Formants [4]float32
    Bandwidths [4]float32
    Gains [4]float32
}

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
}

func Sequence(steps ...SynthTimeSamples) Synthesizer {
    return &sequence{steps: steps}
}

func (g *sequence) Synthesize (samples []float32, sliceNumber int) {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return
    }

    sampleOn := 0

    for sampleOn < len(samples) {
        if g.stepOn >= len(g.steps) {
            for i := sampleOn; i < len(samples); i++ {
                samples[i] = 0
            }
            break
        }
        
        samplesNeeded := len(samples) - sampleOn
        stepSamplesLeft := g.steps[g.stepOn].timeSamples - g.stepTimeSamplesPassed
        
        samplesToGet := samplesNeeded
        if stepSamplesLeft < samplesToGet {
            samplesToGet = stepSamplesLeft
        }
        
        g.steps[g.stepOn].input.Synthesize(samples[sampleOn : sampleOn + samplesToGet], sliceNumber)
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
}

var (
    //signal = LowPass(Amplify(Mix(Triangle(Constant(1)), Constant(1)), Constant(880)), Square(Constant(220)))
    //signal = LowPass(Constant(440), Square(Constant(220)))
    //signal = Sine(Constant(220))
    rawVoice = Sawtooth(Constant(75))
    
    hbw = float32(50.0)
    
    ef1 = BandPass(Constant(275), Constant(50), rawVoice)
    ef2 = BandPass(Constant(2250), Constant(75), rawVoice)
    ef3 = BandPass(Constant(2500), Constant(100), rawVoice)
    ef4 = BandPass(Constant(3500), Constant(150), rawVoice)
    vowelEe = Amplify(Constant(1), Mix(ef1, ef2, ef3, ef4))

    if1 = BandPass(Constant(400), Constant(50), rawVoice)
    if2 = Amplify(Gain(-15), BandPass(Constant(2000), Constant(50), rawVoice))
    if3 = Amplify(Gain(-9), BandPass(Constant(2550), Constant(50), rawVoice))
    vowelIh = Mix(if1, if2, if3)

    //signal = vowelEe
    steps = Sequence(SynthTime(Sawtooth(Constant(110)), 30),
                     SynthTime(Sawtooth(Constant(440)), 30),
                     SynthTime(Sawtooth(Constant(220)), 30),
                     SynthTime(Sawtooth(Constant(880)), 30))
    signal = Amplify(Constant(0.01), steps);
    //signal = Mix(Amplify(vowelEe, Triangle(Constant(0.5))), Amplify(vowelIh, Triangle(Constant(0.25))))
    
    file, _ = os.Create("C:/Users/acsh/Documents/Go Projects/src/acshi/voice/output.raw")
)

func makeMusic(samples []int16) {
    //for j := 0; j < len(samples); j++ {
	//    samples[j] = int16(math.MaxInt16 * math.Sin(440.0 * 2.0 * math.Pi * float64(j) / 44100.0)) / 50
    //}
    rawSamples := make([]float32, len(samples))
    signal.Synthesize(rawSamples, rand.Int())
    for i := 0 ; i < len(samples); i++ {
        samples[i] = int16(math.MaxInt16 * rawSamples[i])
    }

    fileBytes := make([]byte, len(samples) * 2)
    for i := 0; i < len(samples); i++ {
        fileBytes[i * 2] = byte(samples[i] & 0xff)
        fileBytes[i * 2 + 1] = byte((int32(samples[i]) & 0xff00) >> 8)
    }
    file.Write(fileBytes)
}

func main() {
	fmt.Println("We have", audio.WaveOutGetNumDevs(), "audio output devices on this system")

	device, _ := audio.WaveOutOpen(0)


    device.Start(makeMusic)
    time.Sleep(time.Second * 120)

	device.Reset()
	device.Close()
	
	file.Close()
}
