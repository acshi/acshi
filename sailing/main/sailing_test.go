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

func TestCoerceAngle(t *testing.T) {
    var a, b, c float64
    a, _, _ = coerceAngleToRange(100, 0, 99)
    assert.Equal(t, 99.0, a)
    a, _, _ = coerceAngleToRange(100, 101, 200)
    assert.Equal(t, 101.0, a)
    a, _, _ = coerceAngleToRange(0, -200, -80)
    assert.Equal(t, -80.0, a)
    a, _, _ = coerceAngleToRange(200, -200, -80)
    assert.Equal(t, -160.0, a)
    a, _, _ = coerceAngleToRange(560, -200, -80)
    assert.Equal(t, -160.0, a)
    a, _, _ = coerceAngleToRange(-400, -200, -80)
    assert.Equal(t, -80.0, a)
    a, b, c = coerceAngleToRange(-500, -200, -80)
    assert.Equal(t, -140.0, a)
    assert.Equal(t, -200.0, b)
    assert.Equal(t, -80.0, c)
    a, b, c = coerceAngleToRange(-500, -560, -440)
    assert.Equal(t, -140.0, a)
    assert.Equal(t, -200.0, b)
    assert.Equal(t, -80.0, c)
}

func TestAngleBetween(t *testing.T) {
    assert.Equal(t, 0.0, angleBetween(100, 100))
    assert.Equal(t, 10.0, angleBetween(90, 100))
    assert.Equal(t, 10.0, angleBetween(110, 100))
    assert.Equal(t, 10.0, angleBetween(470, 100))
    assert.Equal(t, 10.0, angleBetween(100, 470))
    assert.Equal(t, 10.0, angleBetween(470, -260))
}