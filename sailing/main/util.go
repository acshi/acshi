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
// coerceAngleToRange(0, 360, 720) = 360
// coerceAngleToRange(-40, 0, 360) = 320
// coerceAngleToRange(20, 100, 140) = 100
func coerceAngleToRange(a, minA, maxA float64) float64 {
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

// Finds the angle between two angles in degrees, angles modulo 360
func angleBetween(a, b float64) float64 {
    a = math.Mod(a, 360)
    b = math.Mod(b, 360)
    min := math.Min(a, b)
    max := math.Max(a, b)
    if max - min > 360 {
        min += 360
    }
    if max - min > 180 {
        max, min = min + 360, max
    }
    return max - min
}

func signum(a float64) float64 {
    if a == 0.0 {
        return 0.0
    }
    return math.Abs(a) / a
}

// smallest difference between two angles
func compassDiff(head1, head2 float64) float64 {
    d := head2 - head1 // raw difference
    if d >= 0.0 { // head2 is on same compass 'rotation' as head1 and ahead of head1 in a CW sense
        if d <= 180.0 { // head1 is just a little behind head2
            return d
        } else {
            return d - 360 // head1 is very far behind head2, so easier to go CCW (must return negative!)
        }
    } else { // d < 0.0, head1 is ahead of head 2 on same 'rotation'
        if d >= -180.0 { // head1 is only a little ahead of head2
            return d 
        } else {
            return d + 360 // shorter to go CW to get there, so must be positive
        }
    }   
}