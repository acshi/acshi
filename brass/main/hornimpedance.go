package main

import (
    "acshi/brass"
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/png"
	"os"
	"math"
	//"math/rand"
	//"time"
	"acshi/rosenbrock"
	"github.com/davecheney/profile"
	"runtime"
)

func writeBrassAndImpedancesImage(brassThing brass.Brass, impedances, targetFrequencies []float64, targetImpedance float64) {
    width := len(impedances)
    height := int(float64(width) / 1.618)
    // white image
    white := color.Gray{255}
    img := image.NewGray(image.Rect(0, 0, width, height))
    draw.Draw(img, img.Bounds(), &image.Uniform{white}, image.ZP, draw.Src)

    brassBounds := image.Rect(width * 2 / 10, height * 4 / 10, width * 8 / 10, height * 6 / 10)
    brass.DrawBrassOnImage(brassThing, img, brassBounds)

    brass.DrawImpedancesOnImage(impedances, targetFrequencies, targetImpedance, img)

    file, err := os.Create("impedances.png")
    if err != nil {
        fmt.Println(err)
        return
    }

    png.Encode(file, img)
    file.Close()
}

var (
    rands []float64
    nums1 []float64
    nums2 []float64
)

func funExp(x float64) {
    fmt.Println("exp", x, "for A", math.Exp(x), "for B", brass.Exp(x))
}

func main() {
    /*funExp(1)
    funExp(2)
    for i := 0; i < 10; i++ {
        funExp(rand.Float64() / 1000)
    }

    n := 8000000
    rands = make([]float64, n)
    nums1 = make([]float64, n)
    nums2 = make([]float64, n)
    
    for i := 0; i < n; i++ {
        rands[i] = rand.Float64() / 1000
    }

    startTime1 := time.Now().UnixNano()
    for i := 0; i < n; i++ {
        nums1[i] = math.Exp(rands[i])
    }
    fmt.Println("1 took", float64(time.Now().UnixNano() - startTime1) / 1000)

    startTime2 := time.Now().UnixNano()
    for i := 0; i < n; i++ {
        nums2[i] = brass.Exp(rands[i])
    }
    fmt.Println("2 took", float64(time.Now().UnixNano() - startTime2) / 1000)

    sumRelativeDifference := 0.0
    for i := 0; i < n; i++ {
        sumRelativeDifference += math.Abs(nums2[i] - nums1[i]) / math.Abs(nums2[i])
    }
    fmt.Println("average relative error", sumRelativeDifference / float64(n))
    if sumRelativeDifference != 0 {
        os.Exit(0)
    }*/

    runtime.GOMAXPROCS(runtime.NumCPU())
    fmt.Println("Attempting to use", runtime.GOMAXPROCS(0), "different cores")

    profConfig := profile.CPUProfile
    profConfig.ProfilePath = "./"
    prof := profile.Start(profConfig)
    
	/*var pipe brass.Brass
	pipe.Length = 4
	pipe.Radii = []float64{0.01335, 0.01335}
	pipe.Stepped = false*/

    // radius, length, subdivisions
    // len(seg.Coords) = subdivisions + 1s
	mouthpiece := brass.NewSegment(0.01, 0.08763, 5)
	mouthpiece.Resizable = false
	mouthpiece.Coords[0] = brass.Coord{X: 0, DeltaRadius: 0.010, Mutable: false} // Radius: 0.010
	mouthpiece.Coords[1] = brass.Coord{X: 0.008, DeltaRadius: -0.00158, Mutable: false} // Radius: 0.00842
	mouthpiece.Coords[2] = brass.Coord{X: 0.011, DeltaRadius: -0.00042, Mutable: false} // Radius: 0.008
	mouthpiece.Coords[3] = brass.Coord{X: 0.01207, DeltaRadius: -0.0005, Mutable: false} // Radius: 0.0075
	mouthpiece.Coords[4] = brass.Coord{X: 0.013, DeltaRadius: -0.00567, Mutable: false} // Radius: 0.00183
	mouthpiece.Coords[5] = brass.Coord{X: 0.08763, DeltaRadius: 0.00237, Mutable: false} // Radius: 0.0042
	body := brass.NewSegment(0.01, 1.2, 50)
	//bell := brass.NewSegment(0.01, 0.2, 20)
	//bell.Coords[9] = brass.Coord{X: 0.18, Radius: 0.05, Mutable: false}
	//bell.Coords[10] = brass.Coord{X: 0.2, Radius: 0.0635, Mutable: false}

	var trumpet brass.Brass
	trumpet.Stepped = false
	trumpet.BaseSegments = []brass.Segment{mouthpiece, body}//, bell}

    /*maxFreq := 0.238 * 344.37 / pipe.Radii[len(pipe.Radii) - 1] // original 0.58 *

    fmt.Println("maxFreq", maxFreq)

    if maxFreq > 8000 {
        maxFreq = 8000
        fmt.Println("maxFreq reduced to 8000")
    }*/

    targetFrequencies := []float64{58.27, 116.54, 233.08, 349.23, 466.16, 587.33, 698.46, 783.99, 932.33}
    targetCentsBandwidth := 10.0
    targetImpedanceValue := 40000.0

    if trumpet.Length() != 0 {
	    impedances := brass.CalculateInputImpedance(trumpet, 1, 1500, 1)
        writeBrassAndImpedancesImage(trumpet, impedances, targetFrequencies, targetImpedanceValue)
        return;
    }
    
    objectiveFunctionRunCount := 0

    objectiveFunction := func (parameters []float64) float64 {
        brass.ApplyOptimizationParameters(&trumpet, parameters)
        impedances := brass.CalculateInputImpedance(trumpet, 1, 1500, 1)

        score := 0.0
        
        // Values here count for better, lower, scores
        targetsUsed := 0
        for _, freq := range targetFrequencies {
            // 1200th root of 2 is a cent
            lowEnd := freq * math.Pow(2.0, -targetCentsBandwidth / 1200)
            highEnd := freq * math.Pow(2.0, targetCentsBandwidth / 1200)
            for i := int(lowEnd); i <= int(highEnd + 0.5); i++ {
                residual := impedances[i] - targetImpedanceValue
                score += residual * residual
                targetsUsed++
            }
        }

        // weight these targets to be equal to everything else
        nonTargets := len(impedances) - targetsUsed
        score *= float64(nonTargets) / float64(targetsUsed)

        // Everything else counts for worse, higher, scores
        // cancels out second counting done above
        for _, impedance := range impedances {
            score += impedance * impedance
        }

        // apply a small penalty for larger radii
        // because radius changes are limited, this keeps the optimizer
        // in the loop about where a larger radius it specifies is being ignored
        for _, segment := range trumpet.BaseSegments {
            currentRadius := 0.0
            for _, coords := range segment.Coords {
                currentRadius += coords.DeltaRadius
                newRadius := currentRadius + coords.SelfDeltaRadius
                score += newRadius
            }
        }
        for _, segment := range trumpet.InsertSegments {
            if !segment.Enabled {
                continue
            }
            currentRadius := 0.0
            for _, coords := range segment.Coords {
                currentRadius += coords.DeltaRadius
                newRadius := currentRadius + coords.SelfDeltaRadius
                score += newRadius
            }
        }

        // assess a penalty for being longer, to help make reasonable and short instruments!
        score *= math.Pow(trumpet.Length(), 0.1)

        // report data on progress every once in a while
        objectiveFunctionRunCount++
        if objectiveFunctionRunCount % 10 == 0 {
            fmt.Println(score)
            if objectiveFunctionRunCount % 200 == 0 {
                writeBrassAndImpedancesImage(trumpet, impedances, targetFrequencies, targetImpedanceValue)
                fmt.Println("image written")
                /*if objectiveFunctionRunCount == 10000 {
                    prof.Stop()
                    os.Exit(0)
                }*/
	        }
        }
        
        return score
    }

    parameters := brass.ExtractOptimizationParameters(trumpet)
    fmt.Println("Using", len(parameters), "parameters")
    rosenbrock.Rosenbrock(parameters, 50000, 1e-5, true, objectiveFunction)

	impedances := brass.CalculateInputImpedance(trumpet, 1, 1500, 1)
    writeBrassAndImpedancesImage(trumpet, impedances, targetFrequencies, targetImpedanceValue)

	/*rosenbrockBanana := func (x []float64) float64 {
	    a := x[1] - x[0] * x[0]
	    b := 1 - x[0]
        return 100 * a * a + b * b;
	}

	x := []float64{5, 5}

	rosenbrock.Rosenbrock(x, 5000, 1e-5, true, rosenbrockBanana)

	fmt.Println("minima with", x, "and value of", rosenbrockBanana(x))*/

	fmt.Println("Done!")

	prof.Stop()
}
