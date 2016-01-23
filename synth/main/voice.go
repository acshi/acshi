package main

import (
    "github.com/acshi/audio"
    "github.com/acshi/synth"
    "fmt"
    "math"
    "math/rand"
    "time"
    "os"
)

var (
    rawVoice = synth.Sawtooth(synth.Constant(75))

    ef1 = synth.BandPass(synth.Constant(275), synth.Constant(50), rawVoice, 4)
    ef2 = synth.BandPass(synth.Constant(2250), synth.Constant(75), rawVoice, 8)
    ef3 = synth.BandPass(synth.Constant(2500), synth.Constant(100), rawVoice, 10)
    ef4 = synth.BandPass(synth.Constant(3500), synth.Constant(150), rawVoice, 10)
    vowelEe = synth.Amplify(synth.Constant(0.5), synth.Mix(ef1, ef2, ef3, synth.Amplify(synth.Gain(-18), ef4)))

    if1 = synth.BandPass(synth.Constant(400), synth.Constant(50), rawVoice, 4)
    if2 = synth.BandPass(synth.Constant(1900), synth.Constant(75), rawVoice, 8)
    if3 = synth.BandPass(synth.Constant(2500), synth.Constant(100), rawVoice, 10)
    if4 = synth.BandPass(synth.Constant(3500), synth.Constant(150), rawVoice, 10)
    vowelIh = synth.Amplify(synth.Constant(0.5), synth.Mix(if1, if2, if3, synth.Amplify(synth.Gain(-18), if4)))

    af1 = synth.BandPass(synth.Constant(710), synth.Constant(50), rawVoice, 4)
    af2 = synth.BandPass(synth.Constant(1100), synth.Constant(75), rawVoice, 8)
    af3 = synth.BandPass(synth.Constant(2500), synth.Constant(100), rawVoice, 10)
    af4 = synth.BandPass(synth.Constant(3500), synth.Constant(150), rawVoice, 10)
    vowelAh = synth.Amplify(synth.Constant(0.5), synth.Mix(af1, af2, af3, synth.Amplify(synth.Gain(-18), af4)))

    //signal = vowelIh
    /*steps = synth.Sequence(synth.SynthTime(synth.Sawtooth(synth.Constant(110)), 2),
                             synth.SynthTime(synth.Sawtooth(synth.Constant(440)), 2),
                             synth.SynthTime(synth.Sawtooth(synth.Constant(220)), 2),
                             synth.SynthTime(synth.Sawtooth(synth.Constant(880)), 2))
    signal = synth.Amplify(synth.Constant(0.01), steps);*/
    //signal = //synth.Sequence(synth.SynthTime(vowelAh, 0.1),
             //               synth.SynthTime(vowelEe, 0.5))

    // why
    //signal = synth.Sequence(synth.SynthTime(synth.Diphthong(synth.Vowel_u, synth.Vowel_ɑ, 0.03, 0.1, rawVoice), 0.13),
    //                        synth.SynthTime(synth.Diphthong(synth.Vowel_ɑ, synth.Vowel_i, 0, 0.5, rawVoice), 0.5))

    // rahh
    //signal = synth.Sequence(synth.SynthTime(synth.Diphthong(synth.Vowel_ɜ, synth.Vowel_ɑ, 0.05, 0.1, rawVoice), 0.6))

    // lie
    signal = synth.Sequence(synth.SynthTime(synth.Diphthong(synth.Vowel_l, synth.Vowel_ɑ, 0.05, 0.06, rawVoice), 0.11),
                            synth.SynthTime(synth.Diphthong(synth.Vowel_ɑ, synth.Vowel_i, 0, 0.4, rawVoice), 0.4))

    // all
    //signal = synth.Sequence(synth.SynthTime(synth.Diphthong(synth.Vowel_ɑ, synth.Vowel_ɫ, 0.05, 0.2, rawVoice), 0.5))

    // da?
    //signal = synth.Sequence(synth.SynthTime(synth.TdNoise(), 0.05),
    //                        synth.SynthTime(synth.Diphthong(synth.DSound(synth.Vowel_æ), synth.DSound(synth.Vowel_æ), 0.01, 0.05, synth.GaussianNoise()), 0.06),
    //                        synth.SynthTime(synth.Diphthong(synth.DSound(synth.Vowel_æ), synth.Vowel_æ, 0, 0.1, rawVoice), 0.4))
                            //synth.SynthTime(synth.VowelSound(synth.NSound(synth.Vowel_æ), rawVoice), 0.5),
                            //synth.SynthTime(synth.Diphthong(synth.NSound(synth.Vowel_æ), synth.Vowel_æ, 0.2, rawVoice), 1.0))
                            //synth.SynthTime(synth.Diphthong(synth.Vowel_æ, synth.BSound(synth.Vowel_æ), 0.05, rawVoice), 0.05))

    // sell
    //signal = synth.Sequence(synth.SynthTime(synth.SNoise(), 0.1),
    //                        synth.SynthTime(synth.VowelSound(synth.Vowel_u, rawVoice), 0.3))

    //signal = synth.Transition(synth.Constant(0), synth.Constant(0.99), 1, 2)


    // why are you all...
    //signal = synth.Sequence(synth.SynthTime(synth.Diphthong(synth.Vowel_u, synth.Vowel_ɑ, 0.03, 0.07, rawVoice), 0.1),
    //                        synth.SynthTime(synth.Diphthong(synth.Vowel_ɑ, synth.Vowel_i, 0.00, 0.5, rawVoice), 0.4),
    //                        synth.SynthTime(synth.Diphthong(synth.Vowel_ɑ, synth.Vowel_ɜ, 0.00, 0.6, rawVoice), 0.5),
    //                        synth.SynthTime(synth.Diphthong(synth.Vowel_i, synth.Vowel_u, 0.00, 0.1, rawVoice), 0.5),
    //                        synth.SynthTime(synth.Diphthong(synth.Vowel_ɑ, synth.Vowel_ɫ, 0.05, 0.2, rawVoice), 0.5))

    file, _ = os.Create("output.raw")
)

func makeMusic(samples []int16) {
    //for j := 0; j < len(samples); j++ {
	//    samples[j] = int16(math.MaxInt16 * math.Sin(440.0 * 2.0 * math.Pi * float64(j) / 44100.0)) / 50
    //}
    rawSamples := make([]float32, len(samples))
    isFinished := signal.Synthesize(rawSamples, rand.Int())
    for i := 0 ; i < len(samples); i++ {
        samples[i] = int16(math.MaxInt16 * rawSamples[i])
    }

    fileBytes := make([]byte, len(samples) * 2)
    for i := 0; i < len(samples); i++ {
        fileBytes[i * 2] = byte(samples[i] & 0xff)
        fileBytes[i * 2 + 1] = byte((int32(samples[i]) & 0xff00) >> 8)
    }
    file.Write(fileBytes)

    if isFinished {
        // give it some time to finish giving the audio driver the current stuff
        go func() {
            time.Sleep(time.Second / 10)
            os.Exit(0)
        }();
    }
}

func main() {
	fmt.Println("We have", audio.WaveOutGetNumDevs(), "audio output devices on this system")

	device, _ := audio.WaveOutOpen(0)

    device.Start(makeMusic)
    time.Sleep(time.Second * 60)

	device.Reset()
	device.Close()
	
	file.Close()
}
