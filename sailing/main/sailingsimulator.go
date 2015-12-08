package main

import (
    "fmt"
    "os"
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

type boatData struct {
    p, v, yawRollHeading, w vectorXyz
    mainDirection, jibDirection, rudderDirection RelativeDirection
    course CompassDirection
}

type windData struct {
    direction CompassDirection
    speed float64
}

var (
    /*boatP, boatV vectorXyz // position and velocity
    boatYawRollHeading = vectorXyz{0, 0, -135 / 180.0 * math.Pi}
    boatW vectorXyz // angle, angular velocity
    
    // 0 for headings means North (up)
    courseDirection = CompassDirection(-135)
    // direction wind is coming FROM
    windDirection = CompassDirection(60)
    windSpeed = 2.0
    // mainsail and jib surface normals are relative to the boat's orientation
    mainsailNormal = RelativeDirection(-30)
    jibNormal = RelativeDirection(-30)
    // rudder direction is relative to the negative boat orientation
    // i.e. 0 means directly back off the end of the boat
    rudderDirection = RelativeDirection(0)*/
    //lastPhysicsUpdateTime = int64(0);
)

type vectorXyz struct {
    x, y, z float64
}

// A direction in degrees clockwise from magnetic north, so 90 corresponds to east
type CompassDirection float64

// Conversion to an angle in radians such that east corresponds to 0 and north to pi/2. (counterclockwise)
// Use this to get values that can be plugged into sine and cosine
func (compDir CompassDirection) ToNativeAngle() float64 {
    return -float64(compDir.Normalized()) / 180 * math.Pi + math.Pi / 2
}

func (compDir CompassDirection) Add(a float64) CompassDirection {
    return CompassDirection(float64(compDir) + a)
}

// Returns a normalized compass direction so that values are between in the range (-180, 180]
func (compDir CompassDirection) Normalized() CompassDirection {
    angle := float64(compDir)
    if angle <= -360 || angle > 360 {
        angle = math.Mod(float64(compDir), 360.0)
    }
    if angle <= -180 {
        return CompassDirection(angle + 360)
    } else if angle > 180 {
        return CompassDirection(angle - 360)
    }
    return CompassDirection(angle)
}

// Returns the angle between the two compass directions in the range (-180, 180]
// basically as a - b.
func (a CompassDirection) Minus(b CompassDirection) RelativeDirection {
    return RelativeDirection(float64(a) - float64(b)).Normalized()
}

// Angle in degrees relative to the boat's heading (physical CompassDirection it is pointing)
type RelativeDirection float64

func (relDir RelativeDirection) Add(compDir CompassDirection) CompassDirection {
    return CompassDirection(float64(compDir) + float64(relDir)).Normalized()
}

func (relDir RelativeDirection) Normalized() RelativeDirection {
    return RelativeDirection(CompassDirection(relDir).Normalized())
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

/*func sailHeadingToNative(heading float64) float64 {
    return heading - math.Pi / 2
}*/

// Draw an arrow centered at x, y
func drawArrow(x, y, length int, compDir CompassDirection, img *image.Gray) {
    headingX, headingY := rotate(0, float64(length) / 2, float64(compDir) * math.Pi / 180)
    arrowX, arrowY := rotate(0, float64(length) / 3, float64(compDir + 45) * math.Pi / 180)
    arrowX, arrowY = arrowX - headingX, arrowY - headingY
    draw_line(x - int(headingX), y - int(headingY), x + int(headingX), y + int(headingY), color.Gray{0}, img)
    draw_line(x - int(headingX), y - int(headingY), x + int(arrowX), y + int(arrowY), color.Gray{0}, img)
}

// Draw an arrow with the point at x, y
func drawArrowFromTop(x, y, length int, compDir CompassDirection, img *image.Gray) {
    headingX, headingY := rotate(0, float64(length), float64(compDir) * math.Pi / 180)
    arrowX, arrowY := rotate(0, float64(length) / 3, float64(compDir + 45) * math.Pi / 180)
    draw_line(x, y, x + int(headingX), y + int(headingY), color.Gray{0}, img)
    draw_line(x, y, x + int(arrowX), y + int(arrowY), color.Gray{0}, img)
}

// Draw an arrow with the base at x, y
func drawArrowFromBase(x, y, length int, compDir CompassDirection, img *image.Gray) {
    headingX, headingY := rotate(0, float64(-length), float64(compDir) * math.Pi / 180)
    arrowX, arrowY := rotate(0, float64(-length) / 3, float64(compDir + 45) * math.Pi / 180)
    arrowX, arrowY = headingX - arrowX, headingY - arrowY
    draw_line(x, y, x + int(headingX), y + int(headingY), color.Gray{0}, img)
    draw_line(x + int(headingX), y + int(headingY), x + int(arrowX), y + int(arrowY), color.Gray{0}, img)
}

func rotate(x, y, theta float64) (x2, y2 float64) {
    x2 = x * math.Cos(theta) - y * math.Sin(theta)
    y2 = x * math.Sin(theta) + y * math.Cos(theta)
    return
}

func (boat boatData) getHeading() CompassDirection {
    return CompassDirection(boat.yawRollHeading.z / math.Pi * 180).Normalized()
}

func drawBoat(boat boatData, x, y int, img *image.Gray) {
    // Dipicting the main sail, the jib, the rudder, and the roll and direction of the boat
    mainsailLength, jibLength, rudderLength := 30.0, 20.0, 16.0
    interspace := 0.0//8.0
    totalLength := mainsailLength + jibLength + rudderLength + interspace * 2
    // We need to make sure to use the raw heading here, because drawArrow will convert to native headings
    boatHeading := boat.getHeading()
    boatAngle := boat.yawRollHeading.z

    jibX, jibY := rotate(0.0, -jibLength - interspace - 5.0, boatAngle)
    drawArrowFromTop(x + int(jibX), y + int(jibY), int(jibLength), boat.jibDirection.Add(boatHeading), img);
    
    mainsailX, mainsailY := rotate(0.0, -5.0, boatAngle)
    drawArrowFromTop(x + int(mainsailX), y + int(mainsailY), int(mainsailLength), boat.mainDirection.Add(boatHeading), img)
    
    rudderX, rudderY := rotate(0.0, mainsailLength - 5.0 + interspace, boatAngle)
    drawArrowFromBase(x + int(rudderX), y + int(rudderY), int(rudderLength), boat.rudderDirection.Add(boatHeading).Add(180), img)

    rollLength := boat.yawRollHeading.y / (math.Pi / 2) * totalLength
    drawArrowFromBase(x, y, int(rollLength), boatHeading.Add(90), img)
}

func writeSailingImage(boat boatData, wind windData, outputWriter io.Writer) {
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
    startX := (int(-boat.p.x * 4) % gridResolution + gridResolution) % gridResolution
    startY := (int(-boat.p.y * 4) % gridResolution + gridResolution) % gridResolution
    for x := startX; x < width; x += gridResolution {
        for y := startY; y < height; y += gridResolution {
            img.SetGray(x, y, color.Gray{0});
        }
    }
    
    // Debugging N S W E arrows
    drawArrow(width / 2 - 15, 25, 30, CompassDirection(0), img)
    drawArrow(width / 2 - 15, height - 25, 30, CompassDirection(180), img)
    drawArrow(25, height / 2 - 15, 30, CompassDirection(-90), img)
    drawArrow(width - 25, height / 2 - 15, 30, CompassDirection(90), img)

    drawArrow(25, 25, 40, boat.course, img);
    drawArrow(width - 25, 25, 40, wind.direction, img);

    drawBoat(boat, width / 2, height / 2, img);
    
    drawArrow(width / 2 + 50, height / 2 - 50, 40, boat.getVelocityDirection(), img)
    
    // Calculate the apparent wind heading
    apparentWindDirection := getApparentWindDirection(boat, wind)
    drawArrow(width - 80, 25, 40, apparentWindDirection, img)

    png.Encode(outputWriter, img)
}

func timeInMs() int64 {
    return int64(time.Nanosecond) * time.Now().UnixNano() / int64(time.Millisecond)
}

func (a vectorXyz) Dot(b vectorXyz) float64 {
    return a.x * b.x + a.y * b.y + a.z * b.z
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

// Makes a vector from the compass direction so that north equates to -y, and east to +x.
func vecFromHeading(heading CompassDirection) vectorXyz {
    nHeading := heading.ToNativeAngle()
    return vectorXyz{math.Cos(nHeading), -math.Sin(nHeading), 0.0}
}

func (a vectorXyz) Mag() float64 {
    return math.Sqrt(a.x * a.x + a.y * a.y + a.z * a.z)
}

func (a vectorXyz) vecNormal() vectorXyz {
    return a.Div(a.Mag())
}

func (force vectorXyz) ForceNormalComp(surfaceDirection vectorXyz) vectorXyz {
    surfaceNormal := surfaceDirection.TangentXy()
    return surfaceNormal.Mult(force.Dot(surfaceNormal))
}

func (a vectorXyz) TangentXy() vectorXyz {
    return vectorXyz{-a.y, a.x, 0}
}

func (force vectorXyz) ForceComponent(surfaceDirection vectorXyz) vectorXyz {
    return surfaceDirection.Mult(force.Dot(surfaceDirection))
}

// x  y  z
// Ax Ay Az
// Bx By Bz
func (a vectorXyz) Cross(b vectorXyz) vectorXyz {
    return vectorXyz{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}
}

func (a vectorXyz) GetHeading() CompassDirection {
    return CompassDirection(math.Atan2(a.y, a.x) / math.Pi * 180 + 90).Normalized()
}

// Finds the torque applied by a force at a location (from the center of gravity)
func torqueOnPart(force, location vectorXyz) vectorXyz {
    return location.Cross(force)
}

func calculateForceCenter(baseCenter, length vectorXyz, partHeading RelativeDirection, boatHeading CompassDirection) vectorXyz {
    length.x, length.y = rotate(length.x, length.y, float64(partHeading.Add(boatHeading)) * math.Pi / 180)
    baseCenter.x, baseCenter.y = rotate(baseCenter.x, baseCenter.y, float64(boatHeading) * math.Pi / 180)
    forceCenter := baseCenter.Add(length)
    return forceCenter
}

func calculateForcesTorques(boat boatData, wind windData, printDebug bool) (forces, torques vectorXyz) {
    boatHeading := boat.getHeading()
    
    // For the physics calculations, approximations of center of force
    mainsailCenter := vectorXyz{0, -5, 20}
    mainsailLength := vectorXyz{0, 15.0 / 2, 0}
    mainsailForceCenter := calculateForceCenter(mainsailCenter, mainsailLength, boat.mainDirection, boatHeading)
    
    jibCenter := vectorXyz{0, -15, 15}
    jibLength := vectorXyz{0, 10.0 / 2, 0}
    jibForceCenter := calculateForceCenter(jibCenter, jibLength, boat.jibDirection, boatHeading)
    
    rudderCenter := vectorXyz{0, 15, 0}
    rudderLength := vectorXyz{0, 2, 0}
    rudderForceCenter := calculateForceCenter(rudderCenter, rudderLength, boat.rudderDirection, boatHeading)
    
    keelForceCenter := vectorXyz{0, 0, -30}
    baseMassCenter := vectorXyz{0, 0, 5} // with no roll
    momentOfInertiaScalar := vectorXyz{0.0, 0.01, 0.02} // This will be multipled with the torques to effect angular acceleration

    // These constants (unitless) encapsulate lift/drag/density coefficients with air, water, surface areas, etc...
    mainsailConstant := 5.0
    jibConstant := 3.0
    keelConstant := 20.0 // we also give the keel credit for other lateral drag and friction
    rudderConstant := 10.0

    // Friction coefficients (dimensionless)
    axialFriction := 0.05
    angularFriction := 1.0
    
    // strength of the noise force on the boat
    noiseConstant := 20.0 * 0

    // roughly aproximate
    gravitationalForce := vectorXyz{0, 0, -100}
    
    // negative on the wind speed because the vector is the direction wind is coming FROM
    windVector := vecFromHeading(wind.direction).Mult(-wind.speed)
    apparentWindVector := windVector.Sub(boat.v)
    
    mainsailForce := apparentWindVector.ForceNormalComp(vecFromHeading(boat.mainDirection.Add(boatHeading))).Mult(mainsailConstant)
    jibForce := apparentWindVector.ForceNormalComp(vecFromHeading(boat.jibDirection.Add(boatHeading))).Mult(jibConstant)
    keelForce := boat.v.Neg().ForceNormalComp(vecFromHeading(boatHeading)).Mult(keelConstant)
    axialDragForce := boat.v.Abs().MultVec(boat.v.Neg()).ForceComponent(vecFromHeading(boatHeading.Add(90))).Mult(axialFriction)
    rudderForce := boat.v.Neg().ForceComponent(vecFromHeading(boat.rudderDirection.Add(boatHeading.Add(90)))).Mult(rudderConstant)
    noiseForce := vectorXyz{rand.Float64() * 2 - 1, rand.Float64() * 2 - 1, 0.0}.Mult(noiseConstant)
    
    mainsailTorque := momentOfInertiaScalar.MultVec(torqueOnPart(mainsailForce, mainsailForceCenter))
    jibTorque := momentOfInertiaScalar.MultVec(torqueOnPart(jibForce, jibForceCenter))
    keelTorque := momentOfInertiaScalar.MultVec(torqueOnPart(keelForce, keelForceCenter))
    rudderTorque := momentOfInertiaScalar.MultVec(torqueOnPart(rudderForce, rudderForceCenter))
    if printDebug {
        fmt.Printf("mainsailT %.3f jibT %.3f rudderT %.3f ", mainsailTorque.z, jibTorque.z, rudderTorque.z)
    }

    angularDragTorque := boat.w.Abs().MultVec(boat.w.Neg()).Mult(angularFriction)
    massCenter := vectorXyz{0, 0, 0}
    massCenter.x, massCenter.z = rotate(baseMassCenter.x, baseMassCenter.z, boat.yawRollHeading.y)
    gravityTorque := momentOfInertiaScalar.MultVec(torqueOnPart(gravitationalForce, massCenter))
    
    forces = mainsailForce.Add(jibForce).Add(keelForce).Add(axialDragForce).Add(noiseForce)
    torques = mainsailTorque.Add(jibTorque).Add(keelTorque).Add(rudderTorque).Add(gravityTorque).Add(angularDragTorque)
    return
}

func physicsUpdate(boat boatData, wind windData, dt float64, printDebug bool) (boatData, windData) {
    //currentMs := timeInMs()
    //elapsedMs := currentMs - lastPhysicsUpdateTime
    //lastPhysicsUpdateTime = currentMs
    elapsedSec := dt//float64(elapsedMs) / 1000.0
    
    // wind change
    /*wind.direction = RelativeDirection((rand.Float64() * 2 - 1) * 360 * elapsedSec * 0 + elapsedSec * 5 * 0).Add(wind.direction)
    // we make sure that their is still possible for our boat to sail to our course direction without beating
    // wind.direction is the direction wind is FROM, so course and wind directions must differ by more than the clearance
    // Clearance can be as low as 45 degrees for some boats, but we may need to be more conservative
    courseWindAngle := float64(boat.course.Minus(wind.direction))
    ironsClearance := 45.0
    if courseWindAngle > 0 && courseWindAngle < ironsClearance {
        wind.direction = boat.course.Add(ironsClearance)
    } else if courseWindAngle < 0 && courseWindAngle > -ironsClearance {
        wind.direction = boat.course.Add(-ironsClearance)
    }*/
    
    // Velocity verlat physics update algorithm
    forces, torques := calculateForcesTorques(boat, wind, printDebug)
    boatVHalfStep := boat.v.Add(forces.Mult(elapsedSec / 2.0))
    boatWHalfStep := boat.w.Add(torques.Mult(elapsedSec / 2.0))
    
    boat.p = boat.p.Add(boat.v.Mult(elapsedSec))
    boat.yawRollHeading = boat.yawRollHeading.Add(boat.w.Mult(elapsedSec))
    
    forcesHalfStep, torquesHalfStep := calculateForcesTorques(boat, wind, false)
    boat.v = boatVHalfStep.Add(forcesHalfStep.Mult(elapsedSec / 2.0))
    boat.w = boatWHalfStep.Add(torquesHalfStep.Mult(elapsedSec / 2.0))
    
    return boat, wind
}

// linear interpolation of x in range from minX to maxX
// onto the range minVal to maxVal.
// x is coerced into the range minX to maxX if it isn't already.
func lerp(x, minX, maxX, minVal, maxVal float64) float64 {
    x = math.Min(math.Max(x, minX), maxX)
    return (x - minX) / (maxX - minX) * (maxVal - minVal) + minVal
}

// coerses the angle into the given range, using modular 360 arithmetic
func coerceAngleToRange(a, minA, maxA float64) float64 {
    a = math.Mod(a, 360)
    modMinA := math.Mod(minA, 360)
    modMaxA := math.Mod(maxA, 360)
    if modMinA >= modMaxA {
        modMaxA += 360
    }

    if a < modMinA && (modMinA - a > a + 360 - modMaxA) {
        a += 360
    } else if a > modMaxA && (a - modMaxA > modMinA - (a - 360)) {
        a -= 360
    }
    a = math.Min(math.Max(a, modMinA), modMaxA)

    // Now return a to be between original min and max
    if a < minA {
        a += math.Trunc((maxA - a) / 360) * 360
    } else if a > maxA {
        a -= math.Trunc((a - minA) / 360) * 360
    }
    return a
}

// Finds the angle between two angles in degrees, angles modulo 360
func angleBetween(a, b float64) float64 {
    a = math.Mod(a, 360)
    b = math.Mod(b, 360)
    min := math.Min(a, b)
    max := math.Max(a, b)
    if max - min > 360 {
        min += 360
    }
    if max - min > 180 {
        max, min = min + 360, max
    }
    return max - min
}

func signum(a float64) float64 {
    if a == 0.0 {
        return 0.0
    }
    return math.Abs(a) / a
}

func (boat boatData) getVelocityDirection() CompassDirection {
    return CompassDirection(math.Atan2(boat.v.y, boat.v.x) / math.Pi * 180 + 90).Normalized()
}

func getApparentWindDirection(boat boatData, wind windData) CompassDirection {
    // Calculate the apparent wind heading
    windVector := vecFromHeading(wind.direction).Mult(-wind.speed)
    apparentWindVector := windVector.Sub(boat.v)
    apparentWindDirection := float64(apparentWindVector.GetHeading()) + 180
    return CompassDirection(apparentWindDirection).Normalized()
}

func startupSettings() (boatData, windData) {
    var boat boatData
    var wind windData
    
    boat.p = vectorXyz{0, 0, 0}
    boat.v = vectorXyz{0, 0, 0}
    boat.yawRollHeading = vectorXyz{0, 0, -135 / 180.0 * math.Pi}
    boat.w = vectorXyz{0, 0, 0}
    
    boat.course = CompassDirection(-135)
    boat.mainDirection = RelativeDirection(-30)
    boat.jibDirection = RelativeDirection(-30)
    boat.rudderDirection = RelativeDirection(0)
    
    wind.direction = CompassDirection(45)
    wind.speed = 2.0
    
    return boat, wind
}

var imageBoat boatData
var imageWind windData

func main() {
    boat, wind := startupSettings()

    http.HandleFunc("/", func (w http.ResponseWriter, r *http.Request) {
        writeSailingImage(imageBoat, imageWind, w)
    })
    go http.ListenAndServe(":8080", nil)
    
    pidInitialized := false
    
    // Reset the simulation when you press enter
    go func() {
        for {
            byteBuffer := make([]byte, 1)
            bytesRead, _ := os.Stdin.Read(byteBuffer)
            if bytesRead > 0 && byteBuffer[0] == 'r' {
                boat, wind = startupSettings()
                pidInitialized = false
            } else if bytesRead > 0 {
                if byteBuffer[0] >= byte('0') && byteBuffer[0] <= byte('9') {
                    boat.course = CompassDirection(float64(byteBuffer[0] - byte('0')) * 22.5 - 180).Normalized()
                }
            }
        }
    }()

    //lastPhysicsUpdateTime = timeInMs()
    timeStepOn := 0
    dt := 0.001 // time in seconds between physics calculation updates
    controlInterval := 0.1 // time in seconds between updates to rudder/sail headings
    reportInterval := 0.5 // time in seconds between outputs of data to console

    // PID control of rudder for course correction
    previousInput := 0.0
    integralTerm := 0.0
    // control coefficients
    proportionK := 1.0
    integralK := 0.2 * controlInterval
    derivativeK := 0.5 / controlInterval
    
    for {
        if timeStepOn % int(controlInterval / dt) == 0 {
            // Calculate the limits of PID output values
            boatHeading := float64(boat.getHeading())
            velocityDirection := float64(boat.getVelocityDirection())
            
            // Use the lower-bound of the rudder as a reference point 0.
            // So an angle here of 0 would correspond to the rudder 60
            // degrees counterclockwise from straight.
            
            // The angles of the rudder most effective at effecting rotation
            // are those that are perpindicular to the velocity of the boat
            optimalTurnLeft := coerceAngleToRange((velocityDirection - 90.0) - (boatHeading - 180) + 60.0, 0, 360)
            //turningCenter := coerceAngleToRange((velocityDirection) - (boatHeading - 180) + 60.0, -180, 180)
            optimalTurnRight := coerceAngleToRange((velocityDirection + 90.0) - (boatHeading + 180) + 60.0, -360, 0)
            turningCenter := (optimalTurnLeft + optimalTurnRight) / 2.0
            
            turnLeftDirection := math.Max(coerceAngleToRange(optimalTurnLeft, 0, 120), turningCenter)
            turnLeftDirection = coerceAngleToRange(turnLeftDirection, 0, 120)
            turnRightDirection := math.Min(coerceAngleToRange(optimalTurnRight, 0, 120), turningCenter)
            turnRightDirection = coerceAngleToRange(turnRightDirection, 0, 120)
            
            // Calculate maximum and minimum PID output values currently realizable with our rudder
            // This varies over time as the boat's orientation and velocity and the wind vary
            // Note that positive PID values are negative to indicate to turn left, and positive for right
            // But the rudder direction angles are high(positive) for left and low(negative) for right
            pidLeftLimit := lerp(turnLeftDirection, 0.0, 120.0, 100, -100)
            pidRightLimit := lerp(turnRightDirection, 0.0, 120.0, 100, -100)
            
            if timeStepOn % int(reportInterval / dt) == 0 {
                fmt.Printf("boat: %.2f veloc: %.2f ", boatHeading, velocityDirection)
            }
        
            // Adjust rudder to get to desired heading
            // Rudder needs to be tangential to the boat velocity to generate
            // torque needed to turn the boat
            
            // PID control
            
            // the effect of the rudder will depend on magnitude of our velocity
            scaleFactor := 1.0// * math.Sqrt(boat.v.x * boat.v.x + boat.v.y * boat.v.y + 0 * 25 * boat.w.z * boat.w.z)
            if scaleFactor != 0.0 {
                scaleFactor = 1 / scaleFactor
            }
            scaledPropK := proportionK * scaleFactor
            scaledIntK := integralK * scaleFactor
            scaledDerivK := derivativeK * scaleFactor
            
            pidSetpoint := float64(boat.course)
            // simulate some error in our ability to read our heading
            pidInput := float64(boat.getHeading()) + (rand.Float64() * 2 - 1) * 2 * 0
            
            if !pidInitialized {
                pidInitialized = true
                previousInput = pidInput
                integralTerm = 0.0
            }
            
            courseError := float64(RelativeDirection(pidSetpoint - pidInput).Normalized())
            integralTerm += courseError * scaledIntK
            if scaledIntK != 0.0 {
                if pidRightLimit > pidLeftLimit {
                    if integralTerm >= pidRightLimit || integralTerm <= pidLeftLimit {
                        integralTerm = math.Max(math.Min(integralTerm, pidRightLimit), pidLeftLimit)
                    }
                } else {
                    if integralTerm >= pidLeftLimit || integralTerm <= pidRightLimit {
                        integralTerm = math.Max(math.Min(integralTerm, pidLeftLimit), pidRightLimit)
                    }
                }
            }
            inputDerivative := coerceAngleToRange(pidInput - previousInput, -180, 180)
            pidOutputValue := courseError * scaledPropK + integralTerm - inputDerivative * scaledDerivK
            previousInput = pidInput
            
            // conversion of that output value (from -100 to 100) to a physical rudder orientation
            
            var constrainedPidValue, newRudderDir float64
            if pidRightLimit > pidLeftLimit {
                // If we can't turn in the direction we want to, go neutral
                // Set us between the limits so rudder = 0.
                if (pidLeftLimit > 0 && pidOutputValue < 0) || (pidRightLimit < 0 && pidOutputValue > 0) {
                    constrainedPidValue = 0//(pidLeftLimit + pidRightLimit) / 2.0
                } else {
                    constrainedPidValue = math.Max(math.Min(pidOutputValue, pidRightLimit), pidLeftLimit)
                }
                newRudderDir = lerp(constrainedPidValue, -100, 100, turnLeftDirection, turnRightDirection)
                //newRudderDir = lerp(constrainedPidValue, pidLeftLimit, pidRightLimit, turnLeftDirection, turnRightDirection)
            } else {
                if (pidLeftLimit < 0 && pidOutputValue > 0) || (pidRightLimit > 0 && pidOutputValue < 0) {
                    constrainedPidValue = 0//(pidLeftLimit + pidRightLimit) / 2.0
                } else {
                    constrainedPidValue = math.Max(math.Min(pidOutputValue, pidLeftLimit), pidRightLimit)
                }
                newRudderDir = lerp(constrainedPidValue, -100, 100, turnRightDirection, turnLeftDirection)
                //newRudderDir = lerp(constrainedPidValue, pidRightLimit, pidLeftLimit, turnRightDirection, turnLeftDirection)
            }
            
            newRudderDir -= 60
            
            boat.rudderDirection = RelativeDirection(newRudderDir).Normalized();
            
            if timeStepOn % int(reportInterval / dt) == 0 {
                fmt.Printf("Course Error: %.2f OutPID: %.2f rudder: %.2f ", courseError, constrainedPidValue, boat.rudderDirection)
            }
            
            apparentWindDirection := float64(getApparentWindDirection(boat, wind))
            
            // make sails catch full wind, remember wind direction is the direction the wind is coming FROM
            mainDirection1 := coerceAngleToRange(apparentWindDirection + 90 - float64(boat.getHeading()), -180, 180)
            mainDirection2 := coerceAngleToRange(apparentWindDirection - 90 - float64(boat.getHeading()), -180, 180)
            if math.Abs(mainDirection1) < math.Abs(mainDirection2) {
                boat.mainDirection = RelativeDirection(mainDirection1)
            } else {
                boat.mainDirection = RelativeDirection(mainDirection2)
            }
            // Pull the sail fully in if it is on the side closer to the wind
            if float64(boat.mainDirection) < 0.0 && mainDirection2 < 0.0 {
                boat.mainDirection = RelativeDirection(0)
                boat.jibDirection = RelativeDirection(25)
            } else if float64(boat.mainDirection) > 0.0 && mainDirection2 < 0.0 {
                boat.mainDirection = RelativeDirection(0)
                boat.jibDirection = RelativeDirection(-25)
            } else {
                boat.jibDirection = boat.mainDirection
            }
            
            // Use the sails to help us steer to our desired course
            // Turn to the wind by pulling the mainsail in and freeing the jib
            // Turn away by pulling the jib in and freeing the mailsail
            offWindAngle := coerceAngleToRange(apparentWindDirection - float64(boat.getHeading()), -180, 180)
            turnValue := math.Abs(courseError) * 1.0
            if (courseError > 0.0 && offWindAngle > 0.0) ||
               (courseError < 0.0 && offWindAngle < 0.0) {
                // Sail more towards the wind
                mainDir := lerp(turnValue, 0, 100, float64(boat.mainDirection), 0.0)
                jibDir := lerp(turnValue, 0, 100, float64(boat.jibDirection), offWindAngle)
                boat.mainDirection = RelativeDirection(mainDir)
                boat.jibDirection = RelativeDirection(jibDir)
                if timeStepOn % int(reportInterval / dt) == 0 {
                    fmt.Printf("into main: %.2f jib: %.2f ", boat.mainDirection, boat.jibDirection)
                }
           } else {
                // Sail more away from the wind
                mainDir := lerp(turnValue, 0, 100, float64(boat.mainDirection), offWindAngle)
                jibDir := lerp(turnValue, 0, 100, float64(boat.jibDirection), 0.0)
                boat.mainDirection = RelativeDirection(mainDir)
                boat.jibDirection = RelativeDirection(jibDir)
                if timeStepOn % int(reportInterval / dt) == 0 {
                    fmt.Printf("away main: %.2f jib: %.2f ", boat.mainDirection, boat.jibDirection)
                }
           }
            
            // If the mainsail is less than 25 degrees closed, then
            // keep jib from being too closed, open at least 25 degrees
            if math.Abs(float64(boat.mainDirection)) < 25.0 && math.Abs(float64(boat.jibDirection)) < 25.0 {
                if float64(boat.jibDirection) < 0 {
                    boat.jibDirection = RelativeDirection(-25)
                } else {
                    boat.jibDirection = RelativeDirection(25)
                }
            }
            if timeStepOn % int(reportInterval / dt) == 0 {
                fmt.Printf("md1: %.2f md2: %.2f appWind: %.2f %.2f mainsail: %.2f jib: %.2f ", mainDirection1, mainDirection2, wind.direction, apparentWindDirection, float64(boat.mainDirection), float64(boat.jibDirection))
            }
        }
        
        if timeStepOn % int(reportInterval / dt) == 0 {
            boat, wind = physicsUpdate(boat, wind, dt, true)
            fmt.Printf("\n")
        } else {
            boat, wind = physicsUpdate(boat, wind, dt, false)
        }
        imageBoat, imageWind = boat, wind
        
        time.Sleep(time.Duration(float64(time.Second) * dt))
        timeStepOn++
    }
}
