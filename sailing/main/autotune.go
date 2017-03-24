package main

import (
	"fmt"
	"math"
)

var _ = fmt.Print // temp

// run the simulation for 5 seconds with the given course heading
// score returned is the cumulative absolute error in course heading
// and also the cumulative rudder direction (to encourage leaving the rudder at neutral)
func evaluateCourse(s simData, course float64) float64 {
	s.boat.course = course
	score := 0.0

    // Penalty for changes in the rudder -- we want it to move as little as possible
    rudderChangePenalty := 100.0

	simTime := 40.0
	timeLeft := simTime * 0.5
	lastRudder := s.boat.rudderDirection
	for timeLeft > 0 {
		timeLeft -= s.dt

		performTimestep(&s)

		score += math.Abs(angleToRange(s.boat.getVelocityDirection()-s.boat.course, -180, 180)) * s.dt
		score += math.Abs(s.boat.rudderDirection-lastRudder) * rudderChangePenalty * s.dt
		lastRudder = s.boat.rudderDirection
		//score += math.Abs(s.boat.rudderDirection) * s.dt
		//score -= s.boat.v.Mag() * 10 * s.dt
	}

	s.boat.course = -180
	timeLeft = simTime * 0.5
	for timeLeft > 0 {
		timeLeft -= s.dt

		performTimestep(&s)

		score += math.Abs(angleToRange(s.boat.getVelocityDirection()-s.boat.course, -180, 180)) * s.dt
		score += math.Abs(s.boat.rudderDirection-lastRudder) * rudderChangePenalty * s.dt
		lastRudder = s.boat.rudderDirection
		//score += math.Abs(s.boat.rudderDirection) * s.dt
		//score -= s.boat.v.Mag() * 10 * s.dt
	}

	return score / simTime // + (s.proportionK + s.integralK + s.derivativeK) / 4.0
}

func evaluatePID(pidVals []float64) float64 {
	s := startupSettings()
	s.shouldPrint = false
	s.proportionK, s.integralK, s.derivativeK = pidVals[0], pidVals[1], pidVals[2]

	receiver := make(chan float64)

	score := 0.0
	for j := 0; j < 8; j += 2 {
		go func(j int) {
			receiver <- evaluateCourse(s, -180+float64(j)*22.5)
		}(j)
	}
	for j := 0; j < 8; j += 2 {
		score += <-receiver
	}
	return score
}

func evaluateSailPI(piVals []float64) float64 {
	s := startupSettings()
	s.shouldPrint = false
	s.sailP, s.sailI = piVals[0], piVals[1]

	receiver := make(chan float64)

	score := 0.0
	for j := 0; j < 8; j += 2 {
		go func(j int) {
			receiver <- evaluateCourse(s, -180+float64(j)*22.5)
		}(j)
	}
	for j := 0; j < 8; j += 2 {
		score += <-receiver
	}
	return score
}

func evaluatePIDSailPI(vals []float64) float64 {
	s := startupSettings()
	s.shouldPrint = false
	s.proportionK, s.integralK, s.derivativeK, s.sailP, s.sailI = vals[0], vals[1], vals[2], vals[3], vals[4]

	receiver := make(chan float64)

	score := 0.0
	for j := 0; j < 8; j += 2 {
		go func(j int) {
			receiver <- evaluateCourse(s, -180+float64(j)*22.5)
		}(j)
	}
	for j := 0; j < 8; j += 2 {
		score += <-receiver
	}

    // introduce a small penalty to incentivize low gain values.
    for _, val := range vals {
        score += math.Abs(val * 0.01);
    }

    // penalize negative values
    for _, val := range vals {
        if val < 0 {
            score -= val * 0.02;
        }
    }

	return score
}

func preconditionConjugateGradient(A [][]float64, b []float64, epsilon float64, kMax int) []float64 {
	n := len(b)
	x := make([]float64, n)
	r := make([]float64, n)
	z := make([]float64, n)
	p := make([]float64, n)
	aP := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = 1.0
		r[i] = 1.0
	}

	done := false
	k := 0

	// calculate residual
	for i := range b {
		val := b[i]
		for j := range x {
			val -= A[i][j] * x[j]
		}
		r[i] = val
	}

	rDotR := 0.0
	for i := range r {
		rDotR += r[i] * r[i]
	}

	check := math.Sqrt(rDotR)
	if check < epsilon || k+1 >= kMax {
		done = true
	}

	zDotR := 0.0

	for !done {
		// this is the preconditioning step
		// we use the diagonal of A as 'Matrix' M
		// solve M*z = r for z.
		for i := range r {
			z[i] = r[i] / A[i][i]
		}

		// on to the rest of the algorithm...
		lastZDotR := zDotR
		zDotR = 0.0
		for i := range z {
			zDotR += z[i] * r[i]
		}

		if k == 0 {
			copy(p, z)
		} else {
			delta := zDotR / lastZDotR
			for i := range z {
				p[i] = z[i] + delta*p[i]
			}
		}
		pAp := 0.0
		for i := range p {
			val := 0.0
			for j := range p {
				val = val + A[i][j]*p[j]
			}
			aP[i] = val
			pAp += p[i] * val
		}

		gamma := zDotR / pAp

		for i := range x {
			x[i] += gamma * p[i]
		}

		for i := range r {
			r[i] -= gamma * aP[i]
		}

		// check for convergence again
		rDotR = 0.0
		for i := range r {
			rDotR += r[i] * r[i]
		}

		check = math.Sqrt(rDotR)
		if check < epsilon || k+1 >= kMax {
			done = true
		}

		k++
	}
	if k >= kMax {
		fmt.Printf("Failed to converge to solution of the linear system after %d iterations.\n", kMax)
	} else {
		fmt.Printf("Converged to solution of linear system in %d iterations.\n", k)
	}
	for _, val := range x {
		fmt.Printf("%.3e ", val)
	}
	fmt.Printf("\n")

	return x
}

func calcNegGradientHessianAndFVal(f func([]float64) float64, x []float64, negGradient []float64, hessian [][]float64) float64 {
	// Since we need second derivatives, we have to use a larger value here.
	// the square root would be enough for the gradient, but so small that second derivative values would disappear
	// So we use the fourth root!
	machineC := math.Sqrt(math.Sqrt(1e-16))

	n := len(x)

	xPrime := make([]float64, n)

	fVal := f(x)
	fmt.Printf("Current score: %.5f\n", fVal)

	fNextVals := make([]float64, n)
	fPrevVals := make([]float64, n)
	for i := range fNextVals {
		theta := math.Abs(x[i])*machineC + machineC
		copy(xPrime, x)
		xPrime[i] += theta
		fNextVals[i] = f(xPrime)

		xPrime[i] -= theta * 2
		fPrevVals[i] = f(xPrime)
	}

	fmt.Printf("Gradient: ")
	for i := range negGradient {
		theta := math.Abs(x[i])*machineC + machineC
		// 2nd order centered finite difference for first derivative
		negGradient[i] = -(fNextVals[i] - fPrevVals[i]) / (2 * theta)
		fmt.Printf("%.3e ", -negGradient[i])
	}
	fmt.Printf("\n")

	for i := range hessian {
		theta1 := math.Abs(x[i])*machineC + machineC
		// The Hessian is symmetric
		for j := i; j < n; j++ {
			if i == j {
				copy(xPrime, x)
				xPrime[i] -= theta1
				fPrime2 := f(xPrime)
				// 2nd order centered finite difference for second derivative
				hessian[i][j] = (fNextVals[i] + fPrime2 - 2*fVal) / (theta1 * theta1)
			} else {
				theta2 := math.Abs(x[j])*machineC + machineC
				copy(xPrime, x)
				xPrime[i] += theta1
				xPrime[j] += theta2
				fBothNextVal := f(xPrime)
				// 2nd order forward finite difference for cross derivative
				hessian[i][j] = (fBothNextVal - fNextVals[i] - fNextVals[j] + fVal) / (theta1 * theta2)
				hessian[j][i] = hessian[i][j]
			}
		}
	}

	fmt.Printf("Hessian Matrix:\n")
	for i := range hessian {
		for j := range hessian {
			fmt.Printf("%.3e ", hessian[i][j])
		}
		fmt.Printf("\n")
	}

	return fVal
}

func newtonOptimize(f func([]float64) float64, xInitial []float64) ([]float64, float64) {
	epsilon := 1e-3
	innerEpsilon := 1e-5
	kMax := 500
	relaxation := 1.0
	relaxAttemptMax := 20
	n := len(xInitial)

	x := make([]float64, n)
	copy(x, xInitial)
	newX := make([]float64, n)
	negGradient := make([]float64, n)
	hessian := make([][]float64, n)
	for i := range hessian {
		hessian[i] = make([]float64, n)
	}

	k := 0
	check := epsilon * 10

	var fVal = 0.0
OuterLoop:
	for k < kMax && check > epsilon {
		fVal = calcNegGradientHessianAndFVal(f, x, negGradient, hessian)
		dx := preconditionConjugateGradient(hessian, negGradient, innerEpsilon, kMax)

        for _, dxVal := range dx {
            if math.IsNaN(dxVal) {
                break OuterLoop
            }
        }

		currentRelaxation := relaxation
		attemptNum := 0
		newFVal := 0.0
		for attemptNum < relaxAttemptMax {
			for i := range x {
				newX[i] = x[i] + currentRelaxation*dx[i]
			}
			newFVal = f(newX)
			if newFVal < fVal {
				break
			}
			currentRelaxation *= 0.5
			attemptNum++
		}
		copy(x, newX)
		fmt.Printf("Using relaxation: %.3e\n", currentRelaxation)

		check = math.Abs((fVal - newFVal) / fVal * 100)
		fmt.Printf("Percent change in score: %.3e\n", check)
		k++

		fmt.Printf("Iteration %d x values:\n", k)
		for _, val := range x {
			fmt.Printf("%.5f ", val)
		}
		fmt.Printf("\n")
	}

	if k >= kMax {
		fmt.Printf("Failed to converge to a minimum after %d iterations.\n", kMax)
	} else {
		fmt.Printf("Converged to minimum in %d iterations.\n", k)
	}
	for _, val := range x {
		fmt.Printf("%.5f ", val)
	}
	fmt.Printf("\n")

	fmt.Printf("Final score: %.5f\n", fVal)

	return x, fVal
}

func autotunePID() (p, i, d float64) {
	s := startupSettings()
	pidVals := []float64{s.proportionK, s.integralK, s.derivativeK}
	var bestPID []float64
	var bestScore = 0.0

	for i := 0; i < 10; i++ {
		pidVals[0] = s.proportionK + float64(i)*0.4
		pidVals[1] = s.integralK + float64(i)*0.02
		pidVals[2] = s.derivativeK + float64(i)*0.2
		newPID, newScore := newtonOptimize(evaluatePID, pidVals)
		if (newScore < bestScore || bestScore == 0.0) && newScore > 0.0 {
			bestScore = newScore
			bestPID = newPID
		}
		fmt.Println()
	}

	fmt.Printf("Best values: %.5f %.5f %.5f with score: %.5f\n", bestPID[0], bestPID[1], bestPID[2], bestScore)

	return bestPID[0], bestPID[1], bestPID[2]
}

func autotuneSailPI() (p, i float64) {
	s := startupSettings()
	piVals := []float64{s.sailP, s.sailI}
	var bestPI []float64
	var bestScore = 0.0

	for i := 0; i < 10; i++ {
		piVals[0] = s.sailP + float64(i)*0.4
		piVals[1] = s.sailI + float64(i)*0.02
		newPI, newScore := newtonOptimize(evaluateSailPI, piVals)
		if (newScore < bestScore || bestScore == 0.0) && newScore > 0.0 {
			bestScore = newScore
			bestPI = newPI
		}
		fmt.Println()
	}

	fmt.Printf("Best values: %.5f %.5f with score: %.5f\n", bestPI[0], bestPI[1], bestScore)

	return bestPI[0], bestPI[1]
}

func autotunePIDSailPI() (p, i, d, sailP, sailI float64) {
	s := startupSettings()
    vals := []float64{s.proportionK, s.integralK, s.derivativeK, s.sailP, s.sailI}
	var bestVals []float64
	var bestScore = 0.0

	for i := 0; i < 10; i++ {
        vals[0] = s.proportionK + float64(i)*0.4
		vals[1] = s.integralK + float64(i)*0.02
		vals[2] = s.derivativeK + float64(i)*0.2
		vals[3] = s.sailP + float64(i)*0.4
		vals[4] = s.sailI + float64(i)*0.02
		newVals, newScore := newtonOptimize(evaluatePIDSailPI, vals)
		if (newScore < bestScore || bestScore == 0.0) && newScore > 0.0 {
			bestScore = newScore
			bestVals = newVals
		}
		fmt.Println()
	}

	fmt.Printf("Best values: %.5f %.5f %.5f %.5f %.5f with score: %.5f\n", bestVals[0], bestVals[1], bestVals[2], bestVals[3], bestVals[4], bestScore)

	return bestVals[0], bestVals[1], bestVals[2], bestVals[3], bestVals[4]
}
