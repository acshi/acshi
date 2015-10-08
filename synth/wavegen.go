package synth

import (
    "fmt"
    "math"
    "math/rand"
    "time"
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
    uniformNoiseTable []float32
    gaussianNoiseTable []float32
)

// From -1 to 1
func gaussianRandom(r *rand.Rand) float32 {
    // Add n numbers from r, this will give us a gaussian distribution
    // then we just need to set the variance and mean correctly
    const n = 25
    value := float32(0.0)
    for i := 0; i < n; i++ {
        value += r.Float32()
    } 
    value -= n / 2.0
    value *= 2.0 / n
    return value
}

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
                if j == harmonicLimit[k] || j + 1 == harmonicLimit[k] {
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
                if j == harmonicLimit[k] || j + 1 == harmonicLimit[k] {
                    squareWaveTable[k][i] = float32(value * 0.84)
                }
            }
        }
    }

    noiseTableSize := samplingRate * 2
    uniformNoiseTable = make([]float32, noiseTableSize)
    r := rand.New(rand.NewSource(time.Now().UnixNano()))
    for i := 0; i < noiseTableSize; i++ {
        uniformNoiseTable[i] = r.Float32() * 2 - 1
    }

    gaussianNoiseTable = make([]float32, noiseTableSize)
    for i := 0; i < noiseTableSize; i++ {
        gaussianNoiseTable[i] = gaussianRandom(r)
    }
    
    fmt.Println("We have inited!")
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
    noiseTable []float32
    sliceRecorded int
    recordedSamples []float32
}

func (g *noise) Synthesize(samples []float32, sliceNumber int) bool {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return false
    }

    copied := 0
    for {
        copied += copy(samples[copied:], g.noiseTable[g.index:])
        g.index += copied
        if g.index >= len(g.noiseTable) {
            g.index = g.index % len(g.noiseTable)
        } else {
            break
        }
    }

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber

    return false
}

func UniformNoise() Synthesizer {
    return &noise{noiseTable: uniformNoiseTable}
}

func GaussianNoise() Synthesizer {
    return &noise{noiseTable: gaussianNoiseTable}
}

type WaveGen struct {
    phase float32
    waveTable *[][]float32
    input Synthesizer
    sliceRecorded int
    recordedSamples []float32
    recordedWasFinished bool
}

func (g *WaveGen) Synthesize(samples []float32, sliceNumber int) bool {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return g.recordedWasFinished
    }

    waveTable := *g.waveTable
    inputIsFinished := g.input.Synthesize(samples, sliceNumber)
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
    g.recordedWasFinished = inputIsFinished

    return inputIsFinished
}

