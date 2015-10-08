package rosenbrock

import (
    "fmt"
    "math"
)

func Rosenbrock(x []float64, maxIterations int, epsilon float64, verbose bool, obj func ([]float64) float64) {
    n := len(x)

    basis := make([][]float64, n)
    a := make([][]float64, n)
    for i := 0; i < n; i++ {
        basis[i] = make([]float64, n)
        a[i] = make([]float64, n)
    }

    // delta keeps track of the size of each step we make in each direction
    delta := make([]float64, n)
    // lambda keeps track of the total movement we have made since last rotating the basis matrix
    lambda := make([]float64, n)
    xBest := make([]float64, n)
    xNew := make([]float64, n)
    t := make([]float64, n)

    alpha := 0.5
    beta := 1.0 / 2.0

    numFuncEvals := 0

    // initialize the basis matrix to unity
    for i := 0; i < n; i++ {
        basis[i][i] = 1
    }

    // initialize xBest to the initial values
    copy(xBest, x)

    for i := 0; i < n; i++ {
        delta[i] = 0.1 * 0.1
    }

    // start with the initial value of the objective function to minimize
    yFirstFirst := obj(x)
    numFuncEvals++

    for {
        yBest := yFirstFirst
        for {
            yFirst := yBest
            // for each basis direction
            for i := 0; i < n; i++ {
                // Put together the new x values to try out
                for j := 0; j < n; j++ {
                    xNew[j] = xBest[j] + delta[i] * basis[i][j]
                }
                yNew := obj(xNew)
                numFuncEvals++
                if yNew < yBest {
                    lambda[i] += delta[i]
                    delta[i] *= alpha
                    yBest = yNew
                    copy(xBest, xNew)
                } else {
                    // try again in the opposite direction
                    delta[i] *= -beta
                }
            }
            
            if yBest >= yFirst {
                break
            }
        }

        // Definately restart/continue if we haven't even gotten down to minute levels of delta yet
        minDelta := epsilon + 1
        for i := 0; i < n; i++ {
            minDelta = math.Min(minDelta, math.Abs(delta[i]))
        }
        shouldRestart := minDelta > epsilon

        if yBest < yFirstFirst {
            // Definately restart if we are making progress
            minXChange := epsilon + 1
            for i := 0; i < n; i++ {
                minXChange = math.Min(minXChange, math.Abs(xBest[i] - x[i]))
            }
            shouldRestart = shouldRestart || minXChange > epsilon

            if shouldRestart {
                // rotate the basis
                for i := 0; i < n; i++ {
                    a[n-1][i] = lambda[n-1] * basis[n-1][i]
                }
                for k := n - 2; k >= 0; k-- {
                    for i := 0; i < n; i++ {
                        a[k][i] = a[k+1][i] + lambda[k] * basis[k][i]
                    }
                }
                t[n-1] = lambda[n-1] * lambda[n-1]
                for i := n - 2; i >=0; i-- {
                    t[i] = t[i+1] + lambda[i] * lambda[i]
                }
                for i := n - 1; i > 0; i-- {
                    div := math.Sqrt(t[i-1] * t[i])
                    if div != 0 {
                        for j := 0; j < n; j++ {
                            basis[i][j] = (lambda[i-1] * a[i][j] - basis[i-1][j] * t[i]) / div
                        }
                    }
                }
                div := math.Sqrt(t[0])
                for i := 0; i < n; i++ {
                    basis[0][i] = a[0][i] / div
                }

                // Reset for next run through
                copy(x, xBest)
                for i := 0; i < n; i++ {
                    lambda[i] = 0
                    delta[i] = 0.1
                }

                yFirstFirst = yBest
            }
        }

        if !shouldRestart || numFuncEvals + n > maxIterations {
            break
        }
    }
    
    if verbose {
        fmt.Printf("ROSENBROCK method for local optimization (minimization)\nnumber of evaluation of the objective function: %d\n\n", numFuncEvals);
    }
}
