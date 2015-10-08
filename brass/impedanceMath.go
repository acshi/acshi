package brass

import (
    //"time"
	//"fmt"
	"math"
	"math/cmplx"
)

const (
	// at 22 celsius and 1 atmosphere, kg/m3
	airDensity = 1.196
	p          = airDensity
	// at 22 celsius... Pa*s, N*s/m2
	airDynamicViscosity = 1.847e-5
	n                   = airDynamicViscosity
	// at 22 celsius molar mass 28.97... m/s
	soundSpeed = 344.37
	c          = soundSpeed

    // J/(kg*K)
	specificHeatOfAirCP = 1007
	cP = specificHeatOfAirCP

    // J/(kg*K)
	isochoricSpecficHeatOfAirCV = 717.9
	cV = isochoricSpecficHeatOfAirCV

    // W/(m*K)
	thermalConductivityOfAir = 0.02587
	lm = thermalConductivityOfAir

	pi = math.Pi
)

// calculates transformation matrix for a single cone of an instrument's profile
func stageMatrix2(radiusIn, radiusOut, radiusCenter, length, radianFrequency float64) [4]complex128 {
	xi := complex(radiusIn, 0) // m
	xi1 := complex(radiusOut, 0) // m
	l := complex(length, 0) // m
	w := radianFrequency // rad / s

	centerPlaneArea := pi * radiusCenter * radiusCenter // m2
	sm := centerPlaneArea

	sphericalAreaIn := 4 * pi * radiusIn * radiusIn // m2
	si := sphericalAreaIn

	rv := math.Sqrt(p * w * sm / (n * pi)) // rad^1/2
	r0 := p * c / si // kg/(m4*s) or Pa*s, units of impedance
	k := w / c // wavenumber, rad / m

	zc := complex(r0, 0) * complex(1+0.369/rv, -0.369/rv) // Pa*s/sqrt(rad)
	g := complex(k, 0) * complex(1.045/rv, 1+1.045/rv) // sqrt(rad)/m

	gl := g * l // sqrt(rad)
	cosh_gl := cmplx.Cosh(gl) // mangled sqrt(rad)
	sinh_gl := cmplx.Sinh(gl) // mangled sqrt(rad)
	a11 := xi1 / xi * (cosh_gl - sinh_gl/(g*xi1)) // unitless (more mangled sqrt(rad))
	a12 := xi / xi1 * zc * sinh_gl // Pa*s + mangle
	a21 := 1 / zc * ((xi1/xi-1/(g*g*xi*xi))*sinh_gl + (gl*cosh_gl)/(g*g*xi*xi)) // 1/(Pa*s) + mangle
	a22 := xi / xi1 * (cosh_gl + sinh_gl/(g*xi)) // unitless + mangle

	return [4]complex128{a11, a12, a21, a22}
}

func Exp(x float64) float64 {
    const (
        a0 = 362880 * 2.75573192e-6
        a1 = 362880 * 2.75573192e-6
        a2 = 181440 * 2.75573192e-6
        a3 = 60480 * 2.75573192e-6
        a4 = 15120 * 2.75573192e-6
        a5 = 3024 * 2.75573192e-6
        a6 = 504 * 2.75573192e-6
        a7 = 72 * 2.75573192e-6
        a8 = 9 * 2.75573192e-6
    )
    return a0 + x * (a1 + x * (a2 + x * (a3 + x * (a4 + x * (a5 + x * (a6 + x * (a7 + x * a8)))))))
    //return (362880+x*(362880+x*(181440+x*(60480+x*(15120+x*(3024+x*(504+x*(72+x*(9+x))))))))) * 2.75573192e-6
}

func exp(x float64) float64 {
    const (
        a0 = 362880 * 2.75573192e-6
        a1 = 362880 * 2.75573192e-6
        a2 = 181440 * 2.75573192e-6
        a3 = 60480 * 2.75573192e-6
        a4 = 15120 * 2.75573192e-6
        a5 = 3024 * 2.75573192e-6
        a6 = 504 * 2.75573192e-6
        a7 = 72 * 2.75573192e-6
        a8 = 9 * 2.75573192e-6
    )
    return a0 + x * (a1 + x * (a2 + x * (a3 + x * (a4 + x * (a5 + x * (a6 + x * (a7 + x * a8)))))))
}

// from pkh/math/sinh.go but modified to use cheaper exp and multiplication by 0.5
func cosh(x float64) float64 {
    if x < 0 {
            x = -x
    }
    if x > 21 {
            return math.Exp(x) * 0.5
    }
    return (math.Exp(x) + math.Exp(-x)) * 0.5
}

// direct breakout of math/cmplx/sin.go
// calculate sinh and cosh
func sinhcosh(x float64) (sh, ch float64) {
    if math.Abs(x) <= 0.5 {
            return math.Sinh(x), cosh(x)
    }
    e := math.Exp(x)
    ei := 0.5 / e
    e *= 0.5
    return e - ei, e + ei
}

// combined calculation of both Cosh and Sinh
// almost directly based on math/cmplx/sin.go
func cmplxCoshSinh(x complex128) (complex128, complex128) {
    s, c := math.Sincos(imag(x))
    sh, ch := sinhcosh(real(x))
    return complex(c*ch, s*sh), complex(c*sh, s*ch)
}

var (
    pcc_inv = 1 / (p * c * c)
    c_inv = 1 / c
    y = cP / cV
    u = airDynamicViscosity
    clv_inv = c_inv * (p * c) / u
    clt_inv = c_inv * (p * c * cP) / lm
)

// calculates stage transformation matrix for a single cone of an instrument's profile
func stageMatrix(radiusIn, radiusOut, radiusCenter, length, radianFrequency float64) [5]complex128 {
    r := radiusCenter
    w := radianFrequency
    l := length

    //u := airDynamicViscosity
    //lv_inv := (p * c) / u //u / (p * c)
    //lt_inv := (p * c * cP) / lm //lm / (p * c * cP)
    //y := cP / cV
    rv_inv := 1 / (r * math.Sqrt(w * clv_inv))
    rt_inv := 1 / (r * math.Sqrt(w * clt_inv))

    //zv := complex(0, w * p) * (1 + complex(2 / rv, 0) * (1 - 1i) - complex(0, 3 / (rv * rv)))
    zv := complex(0, w * p) * (1 + complex(2 * rv_inv, -2 * rv_inv) - complex(0, 3 * (rv_inv * rv_inv)))
    //yt := complex(0, w / (p * c * c)) * (1 + complex(y - 1, 0) * (complex(math.Sqrt2 / rt, 0) * (1 - 1i) + complex(0, 1 / (rt * rt))))
    yt := complex(0, w * pcc_inv) * (1 + complex(y - 1, 0) * (complex(math.Sqrt2 * rt_inv, -math.Sqrt2 * rt_inv) + complex(0, rt_inv * rt_inv)))

    g := cmplx.Sqrt(zv * yt)
    z := cmplx.Sqrt(zv * complex128Inverse(yt))
    gl := g * complex(l, 0)
	//cosh_gl := cmplx.Cosh(gl)
	//sinh_gl := cmplx.Sinh(gl)
	cosh_gl, sinh_gl := cmplxCoshSinh(gl)
	
    a11 := cosh_gl
    a12 := z * sinh_gl
    a21 := complex128Inverse(z) * sinh_gl
    a22 := cosh_gl

    return [5]complex128{a11, a12, a21, a22, zv}
}

func complex128Inverse(x complex128) complex128 {
    r, i := real(x), imag(x)
    denom := r * r + i * i
    denom_inv := 1 / denom
    return complex(r * denom_inv, -i * denom_inv)
}

func applyStageForCone(radiusIn, radiusOut, length float64, stage [5]complex128, pressureFlow [2]complex128) [2]complex128 {
    p1 := pressureFlow[0]
    u1 := pressureFlow[1]
    zv := stage[4]
    zv_inv := complex128Inverse(zv)
    l := length
    // distance to radiusIn from the truncated cone's apex 
    x2 := complex(radiusIn * length / (radiusOut - radiusIn), 0)
    x2_inv := complex(1 / real(x2), 0)
    // distance to radiusOut...
    x1 := x2 + complex(l, 0)
    
    outMatrix := [2]complex128{p1 * x1, u1 * x1 - p1 * zv_inv}

    var newMatrix [2]complex128
	newMatrix[0] = stage[0]*outMatrix[0] + stage[1]*outMatrix[1]
	newMatrix[1] = stage[2]*outMatrix[0] + stage[3]*outMatrix[1]

	var newPressureFlow [2]complex128
	newPressureFlow[0] = newMatrix[0] * x2_inv
	newPressureFlow[1] = (newMatrix[1] + newPressureFlow[0] * zv_inv) * x2_inv

	// p2*x2 = cosh_gl * p1x1 + z*sinh_gl*u1x1 - z*sinh_gl*p1/zv
	// p2*x2 = p1 * (cosh_gl * x1 - z * sinh_gl / zv) + u1 * (z * sinh_gl * x1)
	// p2 = p1 * ((cosh_gl * x1 - z * sinh_gl / zv) / x2) + u1 * (z * sinh_gl * x1 / x2)
	
	// ... as radiusIn -> radiusOut, x1 -> x2 -> infinity
	// so x1 / x2 -> 1 and z * sinh_gl / zv / x2 -> 0
	// p2 = p1 * (cosh_gl) + u1 * (z * sinh_gl)
	// as compared to the original for a cylinder... Yay!
	// p2 = p1 * (cosh_gl) + u1 * (z * sinh_gl)

	// u2*x2 - p2/zv = p1*x1*sinh_gl/z + cosh_gl*u1*x1 - cosh_gl*p1/zv
	// u2*x2 = p1*(x1*sinh_gl/z - cosh_gl/zv) + u1*(cosh_gl*x1) + p2/zv
	// u2*x2 = p1*(x1*sinh_gl/z - cosh_gl/zv) + u1*(cosh_gl*x1) + p1 * ((cosh_gl * x1 - z * sinh_gl / zv) / x2 / zv) + u1 * (z * sinh_gl * x1 / x2 / zv)
	// u2*x2 = p1*(x1*sinh_gl/z - cosh_gl/zv + (cosh_gl * x1 - z * sinh_gl / zv) / x2 / zv) + u1*(cosh_gl*x1 + z * sinh_gl * x1 / x2 / zv)
	// u2*x2 = p1*()

	return newPressureFlow
}

func outputImpedance2At(endRadius, radianFrequency float64) complex128 {
	w := radianFrequency // rad / s

	// guessing that sl is the spherical area at the end
	sphericalAreaEnd := 4 * pi * endRadius * endRadius // m2
	sl := sphericalAreaEnd

	endPlaneArea := pi * endRadius * endRadius // m2
	sm := endPlaneArea

	si := sphericalAreaEnd // m2

	rv := math.Sqrt(p * w * sm / (n * pi)) // rad^1/2
	r0 := p * c / si // kg/(m4*s) or Pa*s, units of impedance

	al := math.Sqrt(sl / pi) // m
	zc := complex(r0, 0) * complex(1+0.369/rv, -0.369/rv) // Pa*s/sqrt(rad)
	zt := zc * complex(w*w*al*al/(4*c*c), 0.61*w*al/c) // Pa*s + rad mangle

	return zt
}

func outputImpedanceAt(endRadius, radianFrequency float64) complex128 {
    w := radianFrequency
    r := endRadius

    k := w / c
    z := k * r

    z2 := z * z
    z3 := z2 * z
    z4 := z2 * z2
    z5 := z3 * z2
    z6 := z3 * z3
        
    if z < 1.5 {
        return p * c * complex(z2 / 4 + 0.0127 * z4 + 0.082 * z4 * math.Log(z) - 0.023 * z6,
                               0.6133 * z - 0.036 * z3 + 0.034 * z3 * math.Log(z) - 0.0187 * z5)
    } else {// if z < 3.5 {
        /*if z >= 3.5 {
            fmt.Println("oy! z of", z, "for radius", r)
        }*/
        rc := math.Exp(-z) * math.Sqrt(pi * z) * (1 + 3 / 32 * (1 / z2))
        dL := 0.634 - 0.1102 * z + 0.0018 * z2 - 0.00005 * math.Pow(z, 4.9)
        return complex(0, p * c) * cmplx.Tan(complex(k * dL, 1 / 2 * math.Log(rc)))
    //} else {
    //    panic("output impedance z too large")
    }
}

func inputImpedanceAt(brass SimpleBrass, radianFrequency float64) complex128 {
	nSections := len(brass.SimpleSegments)
	endRadius := brass.SimpleSegments[nSections - 1].EndRadius
	// stepped, the end radius is ignored...
	if brass.Stepped {
        endRadius = brass.SimpleSegments[nSections - 1].StartRadius
	}
	outputImpedance := outputImpedanceAt(endRadius, radianFrequency)

	// Each stage manipulates this 2-column matrix of pressure and flow
	// whose ratio is the impedance
	pressureFlow := [2]complex128{outputImpedance, 1}

	// work from known output impedance to the front
	for i := nSections - 1; i >= 0; i-- {
	    length := brass.SimpleSegments[i].Length
	    
		radiusIn := brass.SimpleSegments[i].StartRadius
		radiusOut := brass.SimpleSegments[i].EndRadius
		radiusCenter := (radiusIn + radiusOut) / 2
		if brass.Stepped {
			radiusOut = radiusIn
			radiusCenter = radiusIn
		}
		
		var stage [5]complex128 = stageMatrix(radiusIn, radiusOut, radiusCenter, length, radianFrequency)

        if radiusIn == radiusOut {
            // calculate pressureFlow = stage * pressureFlow
            var newPressureFlow [2]complex128
		    newPressureFlow[0] = stage[0]*pressureFlow[0] + stage[1]*pressureFlow[1]
		    newPressureFlow[1] = stage[2]*pressureFlow[0] + stage[3]*pressureFlow[1]
		    pressureFlow = newPressureFlow
        } else {
            pressureFlow = applyStageForCone(radiusIn, radiusOut, length, stage, pressureFlow)
        }
	}
	return pressureFlow[0] / pressureFlow[1] // Pa*s
}

func CalculateInputImpedance(brass Brass, minHz, maxHz, deltaHz float64) []float64 {
    simpleBrass := BrassToSimpleBrass(brass)

	numFreqs := int(math.Floor(float64((maxHz - minHz) / deltaHz)))
	impedances := make([]float64, numFreqs)
	/*for i := 0; i < numFreqs; i++ {
		radianFrequency := (float64(i)*deltaHz + minHz) * 2 * pi
		inputImpedance := inputImpedanceAt(simpleBrass, radianFrequency)
		impedances[i] = float64(cmplx.Abs(inputImpedance))
	}*/

    // calculate all different frequencies in parallel
    resultChannels := make([]chan float64, numFreqs)
    for i := 0; i < numFreqs; i++ {
        i := i
        resultChannels[i] = make(chan float64, 1)
        go func() {
            radianFrequency := (float64(i)*deltaHz + minHz) * 2 * pi
		    inputImpedance := inputImpedanceAt(simpleBrass, radianFrequency)
		    resultChannels[i] <- float64(cmplx.Abs(inputImpedance))
        }()
    }

    for i := 0; i < numFreqs; i++ {
        impedances[i] = <- resultChannels[i]
    }
	
	return impedances
}

