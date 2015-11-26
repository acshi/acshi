package main

import (
    "math"
    "testing"
    "github.com/stretchr/testify/assert"
)

func TestCanary(t *testing.T) {
    assert.True(t, true, "Canary test. Assert that true is true.")
}

func TestNormalization(t *testing.T) {
    assert.Equal(t, CompassDirection(10), CompassDirection(370).Normalized())
    assert.Equal(t, CompassDirection(0), CompassDirection(720).Normalized())
    assert.Equal(t, CompassDirection(10), CompassDirection(-350).Normalized())
    assert.Equal(t, CompassDirection(0), CompassDirection(-720).Normalized())
    assert.Equal(t, CompassDirection(-170), CompassDirection(190).Normalized())
    assert.Equal(t, CompassDirection(170), CompassDirection(-190).Normalized())
    assert.Equal(t, CompassDirection(180), CompassDirection(-180).Normalized())
    assert.Equal(t, CompassDirection(180), CompassDirection(180).Normalized())
}

func TestToNativeAngle(t *testing.T) {
    assert.InDelta(t, 0, CompassDirection(90).ToNativeAngle(), 1e-8)
    assert.InDelta(t, math.Pi / 2, CompassDirection(0).ToNativeAngle(), 1e-8)
}

func TestVecFromHeading(t *testing.T) {
    vec := vecFromHeading(CompassDirection(-135))
    assert.InDelta(t, -math.Sqrt(2) / 2, vec.x, 1e-8)
    assert.InDelta(t, math.Sqrt(2) / 2, vec.y, 1e-8)
    vec = vecFromHeading(CompassDirection(135))
    assert.InDelta(t, math.Sqrt(2) / 2, vec.x, 1e-8)
    assert.InDelta(t, math.Sqrt(2) / 2, vec.y, 1e-8)
    vec = vecFromHeading(CompassDirection(180))
    assert.InDelta(t, 0.0, vec.x, 1e-8)
    assert.InDelta(t, 1.0, vec.y, 1e-8)
}

func TestCoerceAngle(t *testing.T) {
    var a float64
    a = coerceAngleToRange(100, 0, 99)
    assert.Equal(t, 99.0, a)
    a = coerceAngleToRange(100, 101, 200)
    assert.Equal(t, 101.0, a)
    a = coerceAngleToRange(0, -200, -80)
    assert.Equal(t, -80.0, a)
    a = coerceAngleToRange(200, -200, -80)
    assert.Equal(t, -160.0, a)
    a = coerceAngleToRange(560, -200, -80)
    assert.Equal(t, -160.0, a)
    a = coerceAngleToRange(560, 1000, 1200)
    assert.Equal(t, 1000.0, a)
    a = coerceAngleToRange(550, 1000, 1200)
    assert.Equal(t, 1200.0, a)
    a = coerceAngleToRange(-400, -200, -80)
    assert.Equal(t, -80.0, a)
    a = coerceAngleToRange(-500, -200, -80)
    assert.Equal(t, -140.0, a)
    a = coerceAngleToRange(-500, -560, -440)
    assert.Equal(t, -500.0, a)
    a = coerceAngleToRange(360, 0, 360)
    assert.Equal(t, 0.0, a)
    a = coerceAngleToRange(40, 0, 360)
    assert.Equal(t, 40.0, a)
    a = coerceAngleToRange(40, -360, 0)
    assert.Equal(t, -320.0, a)
}

func TestAngleBetween(t *testing.T) {
    assert.Equal(t, 0.0, angleBetween(100, 100))
    assert.Equal(t, 10.0, angleBetween(90, 100))
    assert.Equal(t, 10.0, angleBetween(110, 100))
    assert.Equal(t, 10.0, angleBetween(470, 100))
    assert.Equal(t, 10.0, angleBetween(100, 470))
    assert.Equal(t, 10.0, angleBetween(470, -260))
}

func TestCalculateForcesTorques(t *testing.T) {
    boat := boatData{v: vectorXyz{-math.Sqrt(2) / 2.0, math.Sqrt(2) / 2.0, 0},
                    yawRollHeading: vectorXyz{0, 0, -135 / 180.0 * math.Pi},
                    course: CompassDirection(-135),
                    mainNormal: RelativeDirection(-90),
                    jibNormal: RelativeDirection(-90),
                    rudderDirection: RelativeDirection(0)}
    wind := windData{direction: CompassDirection(45), speed: 2.0}
    _, torques := calculateForcesTorques(boat, wind, false)
    baseZTorque := torques.z
    
    boat.rudderDirection = RelativeDirection(10)
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z < baseZTorque,
                "Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)
    
    boat.rudderDirection = RelativeDirection(-10)
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z > baseZTorque,
                "Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)
    
    boat.rudderDirection = RelativeDirection(0)
    boat.v.x = math.Sqrt(2) / 2.0
    boat.v.y = math.Sqrt(2) / 2.0
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z < baseZTorque,
                "Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)
                
    boat.rudderDirection = RelativeDirection(10)
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z < baseZTorque,
                "Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)
                
    boat.rudderDirection = RelativeDirection(-10)
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z < baseZTorque,
                "Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)
    
    boat.rudderDirection = RelativeDirection(0)
    boat.v.x = -math.Sqrt(2) / 2.0
    boat.v.y = -math.Sqrt(2) / 2.0
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z > baseZTorque,
                "Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)
                
    boat.rudderDirection = RelativeDirection(10)
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z > baseZTorque,
                "Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)
                
    boat.rudderDirection = RelativeDirection(-10)
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z > baseZTorque,
                "Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)
    
    boat.rudderDirection = RelativeDirection(0)
    boat.yawRollHeading.z = -45 / 180.0 * math.Pi
    boat.v.x = -1.0
    boat.v.y = 0.0
    boat.course = CompassDirection(-90)
    _, torques = calculateForcesTorques(boat, wind, false)
    baseZTorque = torques.z
    
    boat.rudderDirection = RelativeDirection(45)
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z < baseZTorque,
                "Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)
    
    boat.rudderDirection = RelativeDirection(-45)
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z > baseZTorque,
                "Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)
                
    boat.rudderDirection = RelativeDirection(0)
    boat.v.x = 1.0
    _, torques = calculateForcesTorques(boat, wind, false)
    baseZTorque = torques.z
    
    boat.rudderDirection = RelativeDirection(45)
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z > baseZTorque,
                "Magnitude of Z torque %.2f should be > than %.2f (torque without rudder)", torques.z, baseZTorque)
    
    boat.rudderDirection = RelativeDirection(-45)
    _, torques = calculateForcesTorques(boat, wind, false)
    assert.True(t, torques.z < baseZTorque,
                "Magnitude of Z torque %.2f should be < than %.2f (torque without rudder)", torques.z, baseZTorque)
}