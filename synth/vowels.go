package synth

type Vowel struct {
    Formants [4]float32
    Gains [4]float32
}

// f1-f3 from http://www.phon.ucl.ac.uk/home/wells/formants/table-2.htm
// f4 from http://homepages.wmich.edu/~hillenbr/Papers/HillenbrandGettyClarkWheeler.pdf
var (
    // bead
    Vowel_i = Vowel{Formants:   [4]float32{285, 2373, 3088, 3657}}

    // bid
    Vowel_I = Vowel{Formants:   [4]float32{356, 2098, 2696, 3618}}

    // bed, also written E
    Vowel_ɛ = Vowel{Formants:   [4]float32{569, 1965, 2636, 3677}}

    // bad, also written @
    Vowel_æ = Vowel{Formants:   [4]float32{748, 1746, 2460, 3624}}

    // bird, also written R
    Vowel_ɜ = Vowel{Formants:   [4]float32{474, 1379, 1710, 3334}}

    // uh... _a_bout, also written x?
    Vowel_ə = Vowel{Formants:   [4]float32{500, 1150, 1650, 3649}}

    // bud... this is different than ə? also written A
    Vowel_ʌ = Vowel{Formants:   [4]float32{677, 1083, 2340, 3557}}

    // food
    Vowel_u = Vowel{Formants:   [4]float32{309, 939, 2320, 3357}}
                    
    // good, also written U
    Vowel_ʊ = Vowel{Formants:   [4]float32{376, 950, 2440, 3400}}
    
    // born, also written c or Q
    Vowel_ɔ = Vowel{Formants:   [4]float32{599, 891, 2605, 3486}}

    // body, also written o or O
    Vowel_ɒ = Vowel{Formants:   [4]float32{449, 737, 2635, 3384}}

    // bard
    Vowel_ɑ = Vowel{Formants:   [4]float32{710, 1100, 2540, 3687}}

    // like (clear l)
    Vowel_l = Vowel{Formants:  [4]float32{300, 1225, 2950, 3500}}

    // full (dark l)
    Vowel_ɫ = Vowel{Formants:   [4]float32{450, 750, 2600, 3500}}

    // n?
    Vowel_n = Vowel{Formants:   [4]float32{150, 1307, 1900, 2638}, Gains: [4]float32{0, -6, -10, 0}}

)

func DSound(v Vowel) Vowel {
    modified := v
    modified.Formants[0] = 100
    modified.Formants[1] = 2000
    return modified
}

func GSound(v Vowel) Vowel {
    modified := v
    modified.Formants[0] = 100
    modified.Formants[1] = 1300
    return modified
}

func BSound(v Vowel) Vowel {
    modified := v
    modified.Formants[0] = 100
    modified.Formants[1] = 500
    return modified
}

func NSound(v Vowel) Vowel {
    modified := v
    modified.Gains[1] = -100
    modified.Gains[3] = -100
    return modified
}

func VowelSound(v Vowel, rawVoice Synthesizer) Synthesizer {
    f1 := BandPass(Constant(v.Formants[0]), Constant(50), rawVoice, 4)
    f2 := BandPass(Constant(v.Formants[1]), Constant(75), rawVoice, 8)
    f3 := BandPass(Constant(v.Formants[2]), Constant(100), rawVoice, 10)
    f4 := BandPass(Constant(v.Formants[3]), Constant(150), rawVoice, 10)
    if v.Gains[0] != 0.0 {
        f1 = Amplify(Gain(v.Gains[0]), f1)
    }
    if v.Gains[1] != 0.0 {
        f2 = Amplify(Gain(v.Gains[1]), f2)
    }
    if v.Gains[2] != 0.0 {
        f3 = Amplify(Gain(v.Gains[2]), f3)
    }
    return Mix(f1, f2, f3, Amplify(Gain(-18 + v.Gains[3]), f4))
}

func Diphthong(v1, v2 Vowel, steadyStateSeconds, transitionSeconds float32, rawVoice Synthesizer) Synthesizer {
    f1 := BandPass(Transition(Constant(v1.Formants[0]), Constant(v2.Formants[0]), steadyStateSeconds, transitionSeconds), Constant(50), rawVoice, 4)
    f2 := BandPass(Transition(Constant(v1.Formants[1]), Constant(v2.Formants[1]), steadyStateSeconds, transitionSeconds), Constant(75), rawVoice, 8)
    f3 := BandPass(Transition(Constant(v1.Formants[2]), Constant(v2.Formants[2]), steadyStateSeconds, transitionSeconds), Constant(100), rawVoice, 10)
    f4 := BandPass(Transition(Constant(v1.Formants[3]), Constant(v2.Formants[3]), steadyStateSeconds, transitionSeconds), Constant(150), rawVoice, 10)
    if v1.Gains[0] != 0.0 || v2.Gains[0] != 0.0 {
        f1 = Amplify(Transition(Gain(v1.Gains[0]), Gain(v2.Gains[0]), steadyStateSeconds, transitionSeconds), f1)
    }
    if v1.Gains[1] != 0.0 || v2.Gains[1] != 0.0 {
        f2 = Amplify(Transition(Gain(v1.Gains[1]), Gain(v2.Gains[1]), steadyStateSeconds, transitionSeconds), f2)
    }
    if v1.Gains[2] != 0.0 || v2.Gains[2] != 0.0 {
        f3 = Amplify(Transition(Gain(v1.Gains[2]), Gain(v2.Gains[2]), steadyStateSeconds, transitionSeconds), f3)
    }
    f4 = Amplify(Transition(Gain(v1.Gains[3] - 18), Gain(v2.Gains[3] - 18), steadyStateSeconds, transitionSeconds), f4)
    return Mix(f1, f2, f3, f4)
}


