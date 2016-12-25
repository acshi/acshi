package main

import (
	"math"
	"time"
)

// From http://www.coranac.com/2009/07/sines/
// Maximum absolute error of 0.000193
func fastSin(theta float64) float64 {
	theta /= 0.5 * math.Pi
	
	if theta < -1.0 {
		theta += 4.0
	} else if theta > 3.0 {
		theta -= 4.0
	}
	
	if theta > 1.0 && theta <= 3.0 {
		theta = 2.0 - theta
	}
	
	const A = 4.0 * (3.0 / math.Pi - 9.0 / 16.0)
	const B = 2.0 * A - 5.0 / 2.0
	const C = A - 3.0 / 2.0
	thetaSq := theta * theta
	return theta * (A - thetaSq * (B - thetaSq * C))
}

// From http://www.coranac.com/2009/07/sines/
// Maximum absolute error of 0.000193
func fastCos(theta float64) float64 {
	theta /= 0.5 * math.Pi
	theta += 1.0
	
	if theta < -1.0 {
		theta += 4.0
	} else if theta > 3.0 {
		theta -= 4.0
	}
	
	if theta > 1.0 && theta <= 3.0 {
		theta = 2.0 - theta
	}
	
	const A = 4.0 * (3.0 / math.Pi - 9.0 / 16.0)
	const B = 2.0 * A - 5.0 / 2.0
	const C = A - 3.0 / 2.0
	thetaSq := theta * theta
	return theta * (A - thetaSq * (B - thetaSq * C))
}

// Converged to minimum in 10 iterations.
// 1.40989e+01 7.56863e-01 2.55135e+01
// Final score: 7.9273e+01

// Converged to minimum in 4 iterations.
// 1.41002e+01 7.56339e-01 2.55204e+01
// Final score: 7.9270e+01

// Converged to minimum in 4 iterations.
// 1.41005e+01 7.56406e-01 2.55201e+01
// Final score: 7.9281e+01

func rotate(x, y, theta float64) (x2, y2 float64) {
	sinTheta, cosTheta := fastSin(theta), fastCos(theta)
	//sinTheta, cosTheta := math.Sincos(theta)
	x2 = x*cosTheta - y*sinTheta
	y2 = x*sinTheta + y*cosTheta
	return
}

func timeInMs() int64 {
	return int64(time.Nanosecond) * time.Now().UnixNano() / int64(time.Millisecond)
}

// linear interpolation of x in range from minX to maxX
// onto the range minVal to maxVal.
// x is coerced into the range minX to maxX if it isn't already.
func lerp(x, minX, maxX, minVal, maxVal float64) float64 {
	x = math.Min(math.Max(x, minX), maxX)
	return (x-minX)/(maxX-minX)*(maxVal-minVal) + minVal
}

// coerses the angle into the given range, using modular 360 arithmetic
// If min angle to max angle is 360 degrees, then only modular arithmetic occurs
// if min angle and max angle do not span 360 degrees, then the value is coerced
// to the closer of the two angles.
// angleToRange(0, 360, 720) = 360
// angleToRange(-40, 0, 360) = 320
// angleToRange(20, 100, 140) = 100
func angleToRange(a, minA, maxA float64) float64 {
	if a > minA && a <= maxA {
		return a
	}

	a = math.Mod(a, 360)
	modMinA := math.Mod(minA, 360)
	modMaxA := math.Mod(maxA, 360)
	if modMinA >= modMaxA {
		modMaxA += 360
	}

	if a <= modMinA && (modMinA-a > a+360-modMaxA) {
		a += 360
	} else if a > modMaxA && (a-modMaxA > modMinA-(a-360)) {
		a -= 360
	}
	a = math.Min(math.Max(a, modMinA), modMaxA)

	// Now return a to be between original min and max
	if a <= minA {
		a += math.Trunc((maxA-a)/360) * 360
	} else if a > maxA {
		a -= math.Trunc((a-minA)/360) * 360
	}
	return a
}

func normalizeAngle(a float64) float64 {
	for a <= -180 {
		a += 360
	}
	for a > 180 {
		a -= 360
	}
	return a
}

// Makes a vector from the compass direction so that north equates to -y, and east to +x.
func vecFromHeading(heading float64) vectorXyz {
	nHeading := toNativeAngle(heading)
	//sin, cos := math.Sincos(nHeading)
	sin, cos := fastSin(nHeading), fastCos(nHeading)
	return vectorXyz{cos, -sin, 0.0}
}

// Conversion to an angle in radians such that east corresponds to 0 and north to pi/2. (counterclockwise)
// Use this to get values that can be plugged into sine and cosine
func toNativeAngle(a float64) float64 {
	return -a/180*math.Pi + math.Pi/2
}

func signum(a float64) float64 {
	if a == 0.0 {
		return 0.0
	}
	return math.Abs(a) / a
}
