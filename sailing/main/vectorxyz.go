package main

import (
	"math"
)

type vectorXyz struct {
	x, y, z float64
}

func (a vectorXyz) Dot(b vectorXyz) float64 {
	return a.x*b.x + a.y*b.y + a.z*b.z
}

func (a vectorXyz) Neg() vectorXyz {
	return vectorXyz{-a.x, -a.y, -a.z}
}

func (a vectorXyz) Abs() vectorXyz {
	return vectorXyz{math.Abs(a.x), math.Abs(a.y), math.Abs(a.z)}
}

func (a vectorXyz) Add(b vectorXyz) vectorXyz {
	return vectorXyz{a.x + b.x, a.y + b.y, a.z + b.z}
}

func (a vectorXyz) Sub(b vectorXyz) vectorXyz {
	return vectorXyz{a.x - b.x, a.y - b.y, a.z - b.z}
}

func (a vectorXyz) Div(b float64) vectorXyz {
	return vectorXyz{a.x / b, a.y / b, a.z / b}
}

func (a vectorXyz) Mult(b float64) vectorXyz {
	return vectorXyz{a.x * b, a.y * b, a.z * b}
}

func (a vectorXyz) MultVec(b vectorXyz) vectorXyz {
	return vectorXyz{a.x * b.x, a.y * b.y, a.z * b.z}
}

func (a vectorXyz) Mag() float64 {
	return math.Sqrt(a.x*a.x + a.y*a.y + a.z*a.z)
}

func (a vectorXyz) vecNormal() vectorXyz {
	return a.Div(a.Mag())
}

func (a vectorXyz) ForceNormalComp(surfaceDirection vectorXyz) vectorXyz {
	surfaceNormal := surfaceDirection.TangentXy()
	return surfaceNormal.Mult(a.Dot(surfaceNormal))
}

func (a vectorXyz) TangentXy() vectorXyz {
	return vectorXyz{-a.y, a.x, 0}
}

func (a vectorXyz) ForceComponent(surfaceDirection vectorXyz) vectorXyz {
	return surfaceDirection.Mult(a.Dot(surfaceDirection))
}

// x  y  z
// Ax Ay Az
// Bx By Bz
func (a vectorXyz) Cross(b vectorXyz) vectorXyz {
	return vectorXyz{a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}
}

func (a vectorXyz) GetHeading() float64 {
	return normalizeAngle(math.Atan2(a.y, a.x)/math.Pi*180 + 90)
}
