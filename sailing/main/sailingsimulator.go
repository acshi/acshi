package main

import (
    "fmt"
    //"os"
    "time"
	"image"
	"image/color"
	"image/draw"
	"image/png"
	"math"
	"math/rand"
    "net/http"
    "io"
)

const outputImageSize = 400

var (
    boatP, boatV vectorXyz // position and velocity
    boatYawRollHeading = vectorXyz{0, 0, 0}
    boatW vectorXyz // angle, angular velocity
    
    // 0 for headings means North (up)
    targetHeading = -math.Pi / 6
    windDirection = math.Pi / 6
    windSpeed = 2.0
    mainsailHeading = -math.Pi / 6
    jibHeading = -math.Pi / 6
    rudderHeading = -math.Pi
    lastPhysicsUpdateTime = int64(0);
)

type vectorXyz struct {
    x, y, z float64
}

func draw_line(start_x ,start_y ,end_x ,end_y int,col color.Gray, img *image.Gray) {
  // Bresenham's
  var cx int = start_x;
  var cy int = start_y;
 
  var dx int = end_x - cx;
  var dy int = end_y - cy;
  if dx<0 { dx = 0-dx; }
  if dy<0 { dy = 0-dy; }
 
  var sx int;
  var sy int;
  if cx < end_x { sx = 1; } else { sx = -1; }
  if cy < end_y { sy = 1; } else { sy = -1; }
  var err int = dx-dy;
 
  var n int;
  for n=0;n<1000;n++ {
    img.SetGray(cx, cy, col)
    if((cx==end_x) && (cy==end_y)) {return;}
    var e2 int = 2*err;
    if e2 > (0-dy) { err = err - dy; cx = cx + sx; }
    if e2 < dx     { err = err + dx; cy = cy + sy; }
  }
}

func sailHeadingToNative(heading float64) float64 {
    return heading - math.Pi / 2
}

func drawArrow(x, y, length int, heading float64, img *image.Gray) {
    // Convert our idea of heading (0 is north, + goes clockwise) to that of Cos and Sin's
    nativeHeading := sailHeadingToNative(heading)
    headingX := int(math.Cos(nativeHeading) * float64(length) / 2)
    headingY := int(math.Sin(nativeHeading) * float64(length) / 2)
    arrowX := int(math.Cos(nativeHeading + math.Pi / 4) * float64(length) / 3)
    arrowY := int(math.Sin(nativeHeading + math.Pi / 4) * float64(length) / 3)
    draw_line(x - headingX, y - headingY, x + headingX, y + headingY, color.Gray{0}, img)
    draw_line(x + headingX, y + headingY, x + headingX - arrowX, y + headingY - arrowY, color.Gray{0}, img)
}

func rotate(x, y, theta float64) (x2, y2 float64) {
    x2 = x * math.Cos(theta) - y * math.Sin(theta)
    y2 = x * math.Sin(theta) + y * math.Cos(theta)
    return
}

func drawBoat(x, y int, img *image.Gray) {
    // Dipicting the main sail, the jib, the rudder, and the roll and direction of the boat
    mainsailLength, jibLength, rudderLength := 30.0, 20.0, 10.0
    interspace := 8.0
    totalLength := mainsailLength + jibLength + rudderLength + interspace * 2
    // We need to make sure to use the raw heading here, because drawArrow will convert to native headings
    boatHeading := boatYawRollHeading.z

    jibX, jibY := rotate(0.0, -jibLength / 2 - interspace, boatHeading)
    drawArrow(x + int(jibX), y + int(jibY), int(jibLength), jibHeading + boatHeading, img);
    
    mainsailX, mainsailY := rotate(0.0, mainsailLength / 2, boatHeading)
    drawArrow(x + int(mainsailX), y + int(mainsailY), int(mainsailLength), mainsailHeading + boatHeading, img)
    
    rudderX, rudderY := rotate(0.0, mainsailLength + interspace + rudderLength / 2, boatHeading)
    drawArrow(x + int(rudderX), y + int(rudderY), int(rudderLength), rudderHeading + boatHeading, img)

    rollLength := boatYawRollHeading.y / (math.Pi / 2) * totalLength
    rollX, rollY := rotate(rollLength, 0, boatHeading)
    drawArrow(x + int(rollX), y + int(rollY), int(rollLength), boatHeading + math.Pi / 2, img)
}

func writeSailingImage(outputWriter io.Writer) {
    width := outputImageSize
    height := outputImageSize
    
    // make a white image
    white := color.Gray{255}
    img := image.NewGray(image.Rect(0, 0, width, height))
    draw.Draw(img, img.Bounds(), &image.Uniform{white}, image.ZP, draw.Src)

    // a grid to represent the background sea,
    // which will move as the sailboat stays in the center of the image
    // negative y motion is up
    gridResolution := 16
    startX := (int(-boatP.x * 4) % gridResolution + gridResolution) % gridResolution
    startY := (int(-boatP.y * 4) % gridResolution + gridResolution) % gridResolution
    for x := startX; x < width; x += gridResolution {
        for y := startY; y < height; y += gridResolution {
            img.SetGray(x, y, color.Gray{0});
        }
    }

    drawArrow(25, 25, 40, targetHeading, img);
    drawArrow(width - 25, 25, 40, windDirection, img);

    drawBoat(width / 2, height / 2, img);

    // file, err := os.Create("sailing.png")
    // if err != nil {
        // fmt.Println(err)
        // return
    // }

    png.Encode(outputWriter, img)
    // file.Close()
}

func timeInMs() int64 {
    return int64(time.Nanosecond) * time.Now().UnixNano() / int64(time.Millisecond)
}

func vecDot(a, b vectorXyz) float64 {
    return a.x * b.x + a.y * b.y + a.z * b.z
}

func vecNeg(a vectorXyz) vectorXyz {
    return vectorXyz{-a.x, -a.y, -a.z}
}

func vecAbs(a vectorXyz) vectorXyz {
    return vectorXyz{math.Abs(a.x), math.Abs(a.y), math.Abs(a.z)}
}

func vecAdd(a, b vectorXyz) vectorXyz {
    return vectorXyz{a.x + b.x, a.y + b.y, a.z + b.z}
}

func vecSub(a, b vectorXyz) vectorXyz {
    return vectorXyz{a.x - b.x, a.y - b.y, a.z - b.z}
}

func vecDiv(a vectorXyz, b float64) vectorXyz {
    return vectorXyz{a.x / b, a.y / b, a.z / b}
}

func vecMult(a vectorXyz, b float64) vectorXyz {
    return vectorXyz{a.x * b, a.y * b, a.z * b}
}

func vecMultVec(a, b vectorXyz) vectorXyz {
    return vectorXyz{a.x * b.x, a.y * b.y, a.z * b.z}
}

func vecFromHeading(heading float64) vectorXyz {
    nHeading := sailHeadingToNative(heading)
    return vectorXyz{math.Cos(nHeading), math.Sin(nHeading), 0.0}
}

func vecMag(a vectorXyz) float64 {
    return math.Sqrt(a.x * a.x + a.y * a.y + a.z * a.z)
}

func vecNormal(a vectorXyz) vectorXyz {
    return vecDiv(a, vecMag(a))
}

func vecForceNormalComp(force, surfaceDirection vectorXyz) vectorXyz {
    surfaceNormal := vecTangentXy(surfaceDirection)
    return vecMult(surfaceNormal, vecDot(force, surfaceNormal))
}

func vecTangentXy(a vectorXyz) vectorXyz {
    return vectorXyz{-a.y, a.x, 0}
}

func vecForceTangentComp(force, surfaceDirection vectorXyz) vectorXyz {
    return vecMult(surfaceDirection, vecDot(force, surfaceDirection))
}

// x  y  z
// Ax Ay Az
// Bx By Bz
func vecCross(a vectorXyz, b vectorXyz) vectorXyz {
    return vectorXyz{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}
}

// Finds the torque applied by a force at a location (from the center of gravity)
func torqueOnPart(force, location vectorXyz) vectorXyz {
    return vecCross(location, force)
}

func calculateForcesTorques(printDebug bool) (forces, torques vectorXyz) {
    // For the physics calculations, approximations of center of force
    mainsailForceCenter := vectorXyz{0, 15, 20}
    jibForceCenter := vectorXyz{0, -10, 15}
    keelForceCenter := vectorXyz{0, 0, -30}
    baseMassCenter := vectorXyz{0, 0, 5} // with no roll
    momentOfInertiaScalar := vectorXyz{0.0, 0.001, 0.002} // This will be multipled with the torques to effect angular acceleration

    // These constants (unitless) encapsulate lift/drag/density coefficients with air, water, surface areas, etc...
    mainsailConstant := 1.0
    jibConstant := 0.5
    keelConstant := 1.0 // we also give the keel credit for other lateral drag and friction
    rudderConstant := 10.0

    // Friction coefficients (dimensionless)
    axialFriction := 0.05
    angularFriction := 50.0

    // roughly aproximate
    gravitationalForce := vectorXyz{0, 0, -100}
    
    boatHeading := boatYawRollHeading.z

    windVector := vecMult(vecFromHeading(windDirection), windSpeed)
    apparentWindVector := vecSub(windVector, boatV)
    
    mainsailForce := vecMult(vecForceNormalComp(apparentWindVector, vecFromHeading(boatHeading + mainsailHeading)), mainsailConstant)
    jibForce := vecMult(vecForceNormalComp(apparentWindVector, vecFromHeading(boatHeading + jibHeading)), jibConstant)
    keelForce := vecMult(vecForceNormalComp(vecNeg(boatV), vecFromHeading(boatHeading)), keelConstant)
    axialDragForce := vecMult(vecForceTangentComp(vecMultVec(vecAbs(boatV), vecNeg(boatV)), vecFromHeading(boatHeading + math.Pi / 2)), axialFriction)
    rudderForce := vecMult(vecForceNormalComp(vecNeg(boatV), vecFromHeading(boatHeading + rudderHeading)), rudderConstant)
    mainsailTorque := vecMultVec(momentOfInertiaScalar, torqueOnPart(mainsailForce, mainsailForceCenter))
    jibTorque := vecMultVec(momentOfInertiaScalar, torqueOnPart(jibForce, jibForceCenter))
    keelTorque := vecMultVec(momentOfInertiaScalar, torqueOnPart(keelForce, keelForceCenter))
    
    rudderCenter := vectorXyz{0, 30, 0}
    rudderVector := vectorXyz{0, 5, 0}
    rudderVector.x, rudderVector.y = rotate(rudderVector.x, rudderVector.y, rudderHeading - math.Pi)
    rudderForceCenter := vecAdd(rudderCenter, rudderVector)
    rudderTorque := vecMultVec(momentOfInertiaScalar, torqueOnPart(rudderForce, rudderForceCenter))
    // -2 -> (|-2| * 2) = 4
    // 2 -> (|2| * -2) = -4
    angularDragTorque := vecMult(vecMultVec(vecAbs(boatW), vecNeg(boatW)), angularFriction)
    massCenter := vectorXyz{0, 0, 0}
    massCenter.x, massCenter.z = rotate(baseMassCenter.x, baseMassCenter.z, boatYawRollHeading.y)
    gravityTorque := vecMultVec(momentOfInertiaScalar, torqueOnPart(gravitationalForce, massCenter))
    
    forces = vecAdd(vecAdd(vecAdd(mainsailForce, jibForce), keelForce), axialDragForce)
    torques = vecAdd(vecAdd(vecAdd(vecAdd(vecAdd(mainsailTorque, jibTorque), keelTorque), rudderTorque), gravityTorque), angularDragTorque)
    
    if printDebug {
        // fmt.Printf("target heading: %.2f mainsail heading: %.2f jib heading: %.2f rudder heading: %.2f\n", targetHeading, mainsailHeading, jibHeading, rudderHeading)
        // fmt.Printf("apparent wind: %.2f windHeading: %.2f, windVector: %.2f\n", apparentWindVector, windDirection, windVector)
        // fmt.Printf("mainsail: %.2f jib: %.2f keel: %.2f rudder: %.2f\n", mainsailForce, jibForce, keelForce, rudderForce)
        // fmt.Printf("mainsailT: %.2f jibT: %.2f keelT: %.2f rudderT: %.2f gravityT: %.2f\n", mainsailTorque, jibTorque, keelTorque, rudderTorque, gravityTorque)
        // fmt.Printf("Pos: %.2f vel: %.2f yawRollHeading: %.2f ang vel: %.2f\n", boatP, boatV, boatYawRollHeading, boatW)
        fmt.Printf("Position: %.2f Velocity: %.2f, Force: %.2f Drag Force: %.2f\n", boatP, boatV, forces, axialDragForce)
        fmt.Printf("Angular velocity: %.2f, Torque: %.2f, Drag Torque: %.2f\n", boatW, torques, angularDragTorque)
        fmt.Printf("Rudder torque: %.2f, rudder force: %.2f, rudder force center: %.2f\n", rudderTorque, rudderForce, rudderForceCenter)
        velocityHeading := math.Atan2(boatV.y, boatV.x) + math.Pi / 2
        fmt.Printf("Boat: %.2f, Target: %.2f, Wind: %.2f, Velocity: %.2f\n", boatYawRollHeading.z, targetHeading, windDirection, velocityHeading)
        fmt.Printf("Mainsail: %.2f, Jib: %.2f, Rudder: %.2f\n\n", mainsailHeading, jibHeading, rudderHeading)
    }
    return
}

func physicsUpdate(printDebug bool) {
    currentMs := timeInMs()
    elapsedMs := currentMs - lastPhysicsUpdateTime
    lastPhysicsUpdateTime = currentMs
    elapsedSec := float64(elapsedMs) / 1000.0
    
    // wind change
    windDirection += (rand.Float64() * 2 - 1) * math.Pi / 2 * elapsedSec

    // Velocity verlat physics update algorithm
    forces, torques := calculateForcesTorques(printDebug)
    boatVHalfStep := vecAdd(boatV, vecMult(forces, elapsedSec / 2.0))
    boatWHalfStep := vecAdd(boatW, vecMult(torques, elapsedSec / 2.0))
    
    boatP = vecAdd(boatP, vecMult(boatV, elapsedSec))
    boatYawRollHeading = vecAdd(boatYawRollHeading, vecMult(boatW, elapsedSec))
    
    forcesHalfStep, torquesHalfStep := calculateForcesTorques(false)
    boatV = vecAdd(boatVHalfStep, vecMult(forcesHalfStep, elapsedSec / 2.0))
    boatW = vecAdd(boatWHalfStep, vecMult(torquesHalfStep, elapsedSec / 2.0))
    
    // % perform velocity verlat update
    // [forces, collidingParticles] = calculateForces(positions, velocities, diameter, boxL);
    // dampingForces = calculateDamping(positions, velocities, diameter, boxL, dampingB);
    // forces = forces + dampingForces;
    // velocityHalfStep = velocities + forces * dt / 2;
    // positionsStepped = positions + velocityHalfStep * dt;
    // forcesStepped = calculateForces(positionsStepped, velocityHalfStep, diameter, boxL);
    // velocityStepped = velocityHalfStep + forcesStepped * dt / 2;
}

func webHandler(w http.ResponseWriter, r *http.Request) {
    writeSailingImage(w)
}

// linear interpolation of x in range from minX to maxX
// onto the range minVal to maxVal.
// x is coerced into the range minX to maxX if it isn't already.
func lerp(x, minX, maxX, minVal, maxVal float64) float64 {
    x = math.Min(math.Max(x, minX), maxX)
    return (x - minX) / (maxX - minX) * (maxVal - minVal) + minVal
}

func main() {
    http.HandleFunc("/", webHandler)
    go http.ListenAndServe(":8080", nil)

    lastPhysicsUpdateTime = timeInMs()
    timeStepOn := 0
    dt := 0.001 // time in seconds between physics calculation updates
    controlInterval := 0.5 // time in seconds between updates to rudder/sail headings
    reportInterval := 1.0 // time in seconds between outputs of data to console

    for {
        if timeStepOn % int(controlInterval / dt) == 0 {
            // Adjust rudder to get to desired heading
            // Rudder needs to be tangential to the boat velocity to generate
            // torque needed to turn the boat
            
            // Calculate heading of the velocity
            velocityHeading := math.Atan2(boatV.y, boatV.x) + math.Pi / 2
            
            turnLeftHeading := math.Pi / 3
            turnRightHeading := -math.Pi / 3
            
            rudderHeading = lerp(targetHeading - boatYawRollHeading.z, -math.Pi / 2, math.Pi / 2,
                                 turnLeftHeading, turnRightHeading) + math.Pi + velocityHeading - boatYawRollHeading.z
            
            /*if targetHeading < boatYawRollHeading.z {
                //rudderHeading = -math.Pi + math.Pi / 3
                rudderHeading = velocityHeading - math.Pi * 2 / 3 - boatYawRollHeading.z
            } else {
                //rudderHeading = math.Pi - math.Pi / 3
                rudderHeading = velocityHeading + math.Pi * 2 / 3 - boatYawRollHeading.z
            }*/
            
            // constrain the heading values to that keep it within math.Pi/3 on either side of straight back
            rudderHeading = math.Mod(rudderHeading, 2 * math.Pi)
            rudderHeading = math.Min(math.Max(rudderHeading, math.Pi * 2 / 3), math.Pi * 4 / 3)
            
            // make sails catch full wind
            mainsailHeading = math.Mod(windDirection - math.Pi / 2 - boatYawRollHeading.z, math.Pi)
            jibHeading = mainsailHeading
            
            // but keep jib from being too closed keep it open about pi/8 radians
            jibHeading = (jibHeading / math.Abs(jibHeading)) * math.Max(math.Abs(jibHeading), 0.4)
        }
        
        //boatHeading += math.Pi / 16
        time.Sleep(time.Duration(float64(time.Second) * dt))
        
        if timeStepOn % int(reportInterval / dt) == 0 {
            //writeSailingImage()
            physicsUpdate(true)
        } else {
            physicsUpdate(false)
        }
        
        timeStepOn++
    }
}
