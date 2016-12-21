package main

import (
    "fmt"
    "math"
)

var _ = fmt.Print // temp

func evaluateCourse(s simData, course float64) float64 {
    s.boat.course = course
    score := 0.0

    simTime := 5.0
    for simTime > 0 {
        simTime -= s.dt
    
        performTimestep(&s)
        
        score += math.Abs(angleToRange(float64(s.boat.getVelocityDirection()) - course, -180, 180))
    }
    
    return score
}

func evaluatePid(p, i, d float64) float64 {
    s := startupSettings()
    s.proportionK, s.integralK, s.derivativeK = p, i, d
    
    score := 0.0
    for j := 0; j < 8; j++ {
        score += evaluateCourse(s, -180 + float64(j) * 22.5)
    }
    return score
}

func autoTunePid() (p, i, d float64) {
    s := startupSettings()
    p, i, d = s.proportionK, s.integralK, s.derivativeK
    
    return
}
