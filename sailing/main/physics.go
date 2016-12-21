package main

import (
    "math"
    "math/rand"
)

// Finds the torque applied by a force at a location (from the center of gravity)
func torqueOnPart(force, location vectorXyz) vectorXyz {
    return location.Cross(force)
}

func calculateForceCenter(baseCenter, length vectorXyz, partHeading float64, boatHeading float64) vectorXyz {
    length.x, length.y = rotate(length.x, length.y, (partHeading + boatHeading) * math.Pi / 180)
    baseCenter.x, baseCenter.y = rotate(baseCenter.x, baseCenter.y, boatHeading * math.Pi / 180)
    forceCenter := baseCenter.Add(length)
    return forceCenter
}

func calculateForcesTorques(boat boatData, wind windData, printDebug bool) (forces, torques vectorXyz) {
    boatHeading := boat.getHeading()
    
    // For the physics calculations, approximations of center of force
    mainsailRotationCenter := vectorXyz{0, -1, 9}
    mainsailLength := vectorXyz{0, 5.0, 0}
    mainsailForceCenter := calculateForceCenter(mainsailRotationCenter, mainsailLength, boat.mainDirection, boatHeading)
    
    jibRotationCenter := vectorXyz{0, -7, 5}
    jibLength := vectorXyz{0, 5.0, 0}
    jibForceCenter := calculateForceCenter(jibRotationCenter, jibLength, boat.jibDirection, boatHeading)
    
    rudderRotationCenter := vectorXyz{0, 6.5, 0}
    rudderLength := vectorXyz{0, 1, 0}
    rudderForceCenter := calculateForceCenter(rudderRotationCenter, rudderLength, boat.rudderDirection, boatHeading)
    
    keelForceCenter := vectorXyz{0, 0, -5}
    baseMassCenter := vectorXyz{0, 0, -1} // with no roll yet
    momentOfInertiaScalar := vectorXyz{0.0, 0.01, 0.1} // This will be multipled with the torques to effect angular acceleration

    // These constants (unitless) encapsulate lift/drag/density coefficients with air, water, surface areas, etc...
    mainsailConstant := 0.5 / 10
    jibConstant := 0.3 / 10
    keelConstant := 40.0 // we also give the keel credit for other lateral drag and friction
    rudderConstant := 10.0

    // Friction coefficients (dimensionless)
    axialFriction := 0.05
    angularFriction := 3.0
    forwardFriction := 0.5
    
    // strength of the noise force on the boat
    noiseConstant := 20.0 * 0

    // roughly aproximate
    gravitationalForce := vectorXyz{0, 0, 50}
    
    // negative on the wind speed because the vector is the direction wind is coming FROM
    windVector := vecFromHeading(wind.direction).Mult(-wind.speed)
    apparentWindVector := windVector.Sub(boat.v)
    
    mainsailForce := apparentWindVector.ForceNormalComp(vecFromHeading(boat.mainDirection + boatHeading)).Mult(mainsailConstant)
    jibForce := apparentWindVector.ForceNormalComp(vecFromHeading(boat.jibDirection + boatHeading)).Mult(jibConstant)
    keelForce := boat.v.Neg().ForceNormalComp(vecFromHeading(boatHeading)).Mult(keelConstant)
    axialDragForce := boat.v.Abs().MultVec(boat.v.Neg()).ForceComponent(vecFromHeading(boatHeading + 90)).Mult(axialFriction)
    forwardDragForce := boat.v.Abs().MultVec(boat.v.Neg()).Mult(forwardFriction)
    rudderForce := boat.v.Neg().ForceComponent(vecFromHeading(boat.rudderDirection + boatHeading + 90)).Mult(rudderConstant)
    noiseForce := vectorXyz{rand.Float64() * 2 - 1, rand.Float64() * 2 - 1, 0.0}.Mult(noiseConstant)
    
    mainsailTorque := momentOfInertiaScalar.MultVec(torqueOnPart(mainsailForce, mainsailForceCenter))
    jibTorque := momentOfInertiaScalar.MultVec(torqueOnPart(jibForce, jibForceCenter))
    keelTorque := momentOfInertiaScalar.MultVec(torqueOnPart(keelForce, keelForceCenter))
    rudderTorque := momentOfInertiaScalar.MultVec(torqueOnPart(rudderForce, rudderForceCenter))

    angularDragTorque := boat.w.Abs().MultVec(boat.w.Neg()).Mult(angularFriction)
    massCenter := vectorXyz{0, 0, 0}
    massCenter.x, massCenter.z = rotate(baseMassCenter.x, baseMassCenter.z, boat.yawRollHeading.y)
    gravityTorque := momentOfInertiaScalar.MultVec(torqueOnPart(gravitationalForce, massCenter))
    
    if printDebug {
        //fmt.Printf("\n\tmainFc %.3f jibGc: %.3f keelFc: %.3f rudderFc: %.3f gravityFc: %.3f ", mainsailForceCenter, jibForceCenter, keelForceCenter, rudderForceCenter, massCenter)
        //fmt.Printf("\n\tmainsailT %.3f jibT: %.3f keelT: %.3f rudderT: %.3f gravityT: %.3f angularDragT: %.3f yawRollHeading: %.3f ", mainsailTorque, jibTorque, keelTorque, rudderTorque, gravityTorque, angularDragTorque, boat.yawRollHeading)
    }
    
    forces = mainsailForce.Add(jibForce).Add(keelForce).Add(axialDragForce).Add(forwardDragForce).Add(noiseForce)
    torques = mainsailTorque.Add(jibTorque).Add(keelTorque).Add(rudderTorque).Add(gravityTorque).Add(angularDragTorque)
    return
}

func physicsUpdate(boat boatData, wind windData, dt float64, printDebug bool) (boatData, windData) {
    // Velocity verlat physics update algorithm
    forces, torques := calculateForcesTorques(boat, wind, printDebug)
    boatVHalfStep := boat.v.Add(forces.Mult(dt / 2.0))
    boatWHalfStep := boat.w.Add(torques.Mult(dt / 2.0))
    
    boat.p = boat.p.Add(boat.v.Mult(dt))
    boat.yawRollHeading = boat.yawRollHeading.Add(boat.w.Mult(dt))
    
    forcesHalfStep, torquesHalfStep := calculateForcesTorques(boat, wind, false)
    boat.v = boatVHalfStep.Add(forcesHalfStep.Mult(dt / 2.0))
    boat.w = boatWHalfStep.Add(torquesHalfStep.Mult(dt / 2.0))
    
    if printDebug {
        //fmt.Printf("veloc: %.3f percent ", boat.v.Mag() / wind.speed)
    }
    
    return boat, wind
}
