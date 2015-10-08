package synth

func TdNoise() Synthesizer {
    return Amplify(Constant(2), BandPass(Constant(3500), Constant(150), GaussianNoise(), 8))
}

func SNoise() Synthesizer {
    return Amplify(Transition(Constant(0), Gain(6), 0, 0.03), BandStop(Constant(3200), Constant(1200), BandPass(Constant(7400), Constant(800), GaussianNoise(), 2), 6))
}
//BandPass(Constant(7400), Constant(800), GaussianNoise(), 2)
