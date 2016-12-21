package main

import (
    "math"
    "time"
)

func rotate(x, y, theta float64) (x2, y2 float64) {
    x2 = x * math.Cos(theta) - y * math.Sin(theta)
    y2 = x * math.Sin(theta) + y * math.Cos(theta)
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
    return (x - minX) / (maxX - minX) * (maxVal - minVal) + minVal
}

// coerses the angle into the given range, using modular 360 arithmetic
// If min angle to max angle is 360 degrees, then only modular arithmetic occurs
// if min angle and max angle do not span 360 degrees, then the value is coerced
// to the closer of the two angles.
// angleToRange(0, 360, 720) = 360
// angleToRange(-40, 0, 360) = 320
// angleToRange(20, 100, 140) = 100
func angleToRange(a, minA, maxA float64) float64 {
    a = math.Mod(a, 360)
    modMinA := math.Mod(minA, 360)
    modMaxA := math.Mod(maxA, 360)
    if modMinA >= modMaxA {
        modMaxA += 360
    }

    if a < modMinA && (modMinA - a > a + 360 - modMaxA) {
        a += 360
    } else if a > modMaxA && (a - modMaxA > modMinA - (a - 360)) {
        a -= 360
    }
    a = math.Min(math.Max(a, modMinA), modMaxA)

    // Now return a to be between original min and max
    if a < minA {
        a += math.Trunc((maxA - a) / 360) * 360
    } else if a > maxA {
        a -= math.Trunc((a - minA) / 360) * 360
    }
    return a
}

func normalizeAngle(a float64) float64 {
    return angleToRange(a, -180, 180)
}

// Conversion to an angle in radians such that east corresponds to 0 and north to pi/2. (counterclockwise)
// Use this to get values that can be plugged into sine and cosine
func toNativeAngle(a float64) float64 {
    return -normalizeAngle(a) / 180 * math.Pi + math.Pi / 2
}

func signum(a float64) float64 {
    if a == 0.0 {
        return 0.0
    }
    return math.Abs(a) / a
}
