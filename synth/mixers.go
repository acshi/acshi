package synth

type amplify struct {
    inputs []Synthesizer
    sliceRecorded int
    recordedSamples []float32
    recordedWasFinished bool
}

func (g *amplify) Synthesize(samples []float32, sliceNumber int) bool {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return g.recordedWasFinished
    }

    isFinished := false

    inputSamples := make([][]float32, len(g.inputs))
    for i := 0; i < len(g.inputs); i++ {
        inputSamples[i] = make([]float32, len(samples))
        partIsFinished := g.inputs[i].Synthesize(inputSamples[i], sliceNumber)
        // If a single part is finished (likely now giving all zeros)
        // then the amplified result is finished
        isFinished = isFinished || partIsFinished
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
    g.recordedWasFinished = isFinished

    return isFinished
}

func Amplify(inputs ...Synthesizer) Synthesizer {
    return &amplify{inputs: inputs}
}

type mix struct {
    inputs []Synthesizer
    sliceRecorded int
    recordedSamples []float32
    recordedWasFinished bool
}

func (g *mix) Synthesize(samples []float32, sliceNumber int) bool {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return g.recordedWasFinished
    }

    isFinished := true

    inputSamples := make([][]float32, len(g.inputs))
    for i := 0; i < len(g.inputs); i++ {
        inputSamples[i] = make([]float32, len(samples))
        partIsFinished := g.inputs[i].Synthesize(inputSamples[i], sliceNumber)
        // If there is a single part that is not finished, the result is not finished
        isFinished = isFinished && partIsFinished
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
    g.recordedWasFinished = isFinished

    return isFinished
}

func Mix(inputs ...Synthesizer) Synthesizer {
    return &mix{inputs: inputs}
}

type subtract struct {
    signal Synthesizer
    inputs []Synthesizer
    sliceRecorded int
    recordedSamples []float32
    recordedWasFinished bool
}

func (g *subtract) Synthesize(samples []float32, sliceNumber int) bool {
    if sliceNumber == g.sliceRecorded {
        copy(samples, g.recordedSamples)
        return g.recordedWasFinished
    }

    g.signal.Synthesize(samples, sliceNumber)

    isFinished := true

    inputSamples := make([][]float32, len(g.inputs))
    for i := 0; i < len(g.inputs); i++ {
        inputSamples[i] = make([]float32, len(samples))
        partIsFinished := g.inputs[i].Synthesize(inputSamples[i], sliceNumber)
        // just as with the mixer, only a single part need not be finished
        isFinished = isFinished && partIsFinished
    }

    for i := 0; i < len(samples); i++ {
        for j := 0; j < len(g.inputs); j++ {
            samples[i] -= inputSamples[j][i]
        }
    }

    g.recordedSamples = append(g.recordedSamples[:0], samples...)
    g.sliceRecorded = sliceNumber
    g.recordedWasFinished = isFinished

    return isFinished
}

func Subtract(signal Synthesizer, inputs ...Synthesizer) Synthesizer {
    return &subtract{signal: signal, inputs: inputs}
}
