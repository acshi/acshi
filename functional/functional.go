package functional

func Fold(f func(float64, float64) float64, xs []float64) float64 {
    val := xs[0]
    for i := 0; i < len(xs); i++ {
        val = f(val, xs[i])
    }
    return val
}

func Map(f func(float64) float64, xs []float64) []float64 {
    ys := make([]float64, len(xs))
    for i := 0; i < len(xs); i++ {
        ys[i] = f(xs[i])
    }
    return ys
}

func MapInPlace(f func(float64) float64, xs []float64) []float64 {
    for i := 0; i < len(xs); i++ {
        xs[i] = f(xs[i])
    }
    return xs
}

func ZipWith(f func(float64, float64) float64, xs, ys []float64) []float64 {
    n := len(xs)
    if len(ys) < n {
        n = len(ys)
    }

    zs := make([]float64, n)

    for i := 0; i < n; i++ {
        zs[i] = f(xs[i], ys[i])
    }

    return zs
}

func Negate(x float64) float64 {
    return -x
}

func Add(x, y float64) float64 {
    return x + y
}
