package main

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestCanary(t *testing.T) {
	assert.True(t, true, "Canary test. Assert that true is true.")
}

func TestFastSinCos(t *testing.T) {
	assert.Equal(t, 1.0, fastSin(math.Pi / 2))
	assert.Equal(t, -1.0, fastSin(-math.Pi / 2))
	assert.Equal(t, 0.0, fastSin(0))
	assert.Equal(t, 0.0, fastSin(math.Pi))
	assert.Equal(t, 0.0, fastSin(-math.Pi))
	assert.Equal(t, 0.0, fastSin(2 * math.Pi))
	
	assert.Equal(t, 0.0, fastCos(math.Pi / 2))
	assert.Equal(t, 0.0, fastCos(-math.Pi / 2))
	assert.Equal(t, 1.0, fastCos(0))
	assert.Equal(t, -1.0, fastCos(math.Pi))
	assert.Equal(t, -1.0, fastCos(-math.Pi))
	assert.Equal(t, 1.0, fastCos(2 * math.Pi))
}

func TestNormalization(t *testing.T) {
	assert.Equal(t, 10.0, normalizeAngle(370))
	assert.Equal(t, 0.0, normalizeAngle(720))
	assert.Equal(t, 10.0, normalizeAngle(-350))
	assert.Equal(t, 0.0, normalizeAngle(-720))
	assert.Equal(t, -170.0, normalizeAngle(190))
	assert.Equal(t, 170.0, normalizeAngle(-190))
	assert.Equal(t, 180.0, normalizeAngle(-180))
	assert.Equal(t, 180.0, normalizeAngle(180))
}

func TestToNativeAngle(t *testing.T) {
	assert.InDelta(t, 0, toNativeAngle(90), 1e-8)
	assert.InDelta(t, math.Pi/2, toNativeAngle(0), 1e-8)
}

func TestVecFromHeading(t *testing.T) {
	vec := vecFromHeading(-135)
	assert.InDelta(t, -math.Sqrt(2)/2, vec.x, 1e-5)
	assert.InDelta(t, math.Sqrt(2)/2, vec.y, 1e-5)
	vec = vecFromHeading(135)
	assert.InDelta(t, math.Sqrt(2)/2, vec.x, 1e-5)
	assert.InDelta(t, math.Sqrt(2)/2, vec.y, 1e-5)
	vec = vecFromHeading(180)
	assert.InDelta(t, 0.0, vec.x, 1e-8)
	assert.InDelta(t, 1.0, vec.y, 1e-8)
	vec = vecFromHeading(-180)
	assert.InDelta(t, 0.0, vec.x, 1e-8)
	assert.InDelta(t, 1.0, vec.y, 1e-8)
}

func TestAngleToRange(t *testing.T) {
	var a float64
	a = angleToRange(100, 0, 99)
	assert.Equal(t, 99.0, a)
	a = angleToRange(100, 101, 200)
	assert.Equal(t, 101.0, a)
	a = angleToRange(0, -200, -80)
	assert.Equal(t, -80.0, a)
	a = angleToRange(200, -200, -80)
	assert.Equal(t, -160.0, a)
	a = angleToRange(560, -200, -80)
	assert.Equal(t, -160.0, a)
	a = angleToRange(560, 1000, 1200)
	assert.Equal(t, 1000.0, a)
	a = angleToRange(550, 1000, 1200)
	assert.Equal(t, 1200.0, a)
	a = angleToRange(-400, -200, -80)
	assert.Equal(t, -80.0, a)
	a = angleToRange(-500, -200, -80)
	assert.Equal(t, -140.0, a)
	a = angleToRange(-500, -560, -440)
	assert.Equal(t, -500.0, a)
	a = angleToRange(360, 0, 360)
	assert.Equal(t, 360.0, a)
	a = angleToRange(40, 0, 360)
	assert.Equal(t, 40.0, a)
	a = angleToRange(40, -360, 0)
	assert.Equal(t, -320.0, a)
}

func TestCalculateForcesTorques(t *testing.T) {
	boat := boatData{v: vectorXyz{-math.Sqrt(2) / 2.0, math.Sqrt(2) / 2.0, 0},
		yawRollHeading:  vectorXyz{0, 0, -135 / 180.0 * math.Pi},
		course:          -135,
		mainDirection:   -90,
		jibDirection:    -90,
		rudderDirection: 0}
	wind := windData{direction: 45, speed: 2.0}
	_, torques := calculateForcesTorques(boat, wind, false)
	baseZTorque := torques.z

	boat.rudderDirection = 10
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z < baseZTorque,
		"Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)

	boat.rudderDirection = -10
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z > baseZTorque,
		"Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)

	boat.rudderDirection = 0
	boat.v.x = math.Sqrt(2) / 2.0
	boat.v.y = math.Sqrt(2) / 2.0
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z < baseZTorque,
		"Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)

	boat.rudderDirection = 10
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z < baseZTorque,
		"Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)

	boat.rudderDirection = -10
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z < baseZTorque,
		"Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)

	boat.rudderDirection = 0
	boat.v.x = -math.Sqrt(2) / 2.0
	boat.v.y = -math.Sqrt(2) / 2.0
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z > baseZTorque,
		"Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)

	boat.rudderDirection = 10
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z > baseZTorque,
		"Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)

	boat.rudderDirection = -10
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z > baseZTorque,
		"Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)

	boat.rudderDirection = 0
	boat.yawRollHeading.z = -45 / 180.0 * math.Pi
	boat.v.x = -1.0
	boat.v.y = 0.0
	boat.course = -90
	_, torques = calculateForcesTorques(boat, wind, false)
	baseZTorque = torques.z

	boat.rudderDirection = 45
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z < baseZTorque,
		"Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)

	boat.rudderDirection = -45
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z > baseZTorque,
		"Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)

	boat.rudderDirection = 0
	boat.v.x = 1.0
	_, torques = calculateForcesTorques(boat, wind, false)
	baseZTorque = torques.z

	boat.rudderDirection = 45
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z > baseZTorque,
		"Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)

	boat.rudderDirection = -45
	_, torques = calculateForcesTorques(boat, wind, false)
	assert.True(t, torques.z < baseZTorque,
		"Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)
}

func TestGetHeading(t *testing.T) {
	assert.Equal(t, 0.0, vectorXyz{0, -1, 0}.GetHeading())
	assert.Equal(t, 180.0, vectorXyz{0, 1, 0}.GetHeading())
	assert.Equal(t, -90.0, vectorXyz{-1, 0, 0}.GetHeading())
	assert.Equal(t, 90.0, vectorXyz{1, 0, 0}.GetHeading())
}
