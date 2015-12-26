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
    momentOfInertiaScalar := vectorXyz{0.0, 0.01, 0.1} // This will be multipled with the torques to effect angular acceleration

    // These constants (unitless) encapsulate lift/drag/density coefficients with air, water, surface areas, etc...
    mainsailConstant := 5.0 * 0.1
    jibConstant := 3.0 * 0.1
    keelConstant := 20.0 // we also give the keel credit for other lateral drag and friction
    rudderConstant := 10.0 * 0.1

    // Friction coefficients (dimensionless)
    axialFriction := 0.05
    angularFriction := 1.0
    forwardFriction := 0.5
    
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
    forwardDragForce := boat.v.Abs().MultVec(boat.v.Neg()).Mult(forwardFriction)
    rudderForce := boat.v.Neg().ForceComponent(vecFromHeading(boat.rudderDirection.Add(boatHeading.Add(90)))).Mult(rudderConstant)
    noiseForce := vectorXyz{rand.Float64() * 2 - 1, rand.Float64() * 2 - 1, 0.0}.Mult(noiseConstant)
    
    mainsailTorque := momentOfInertiaScalar.MultVec(torqueOnPart(mainsailForce, mainsailForceCenter))
    jibTorque := momentOfInertiaScalar.MultVec(torqueOnPart(jibForce, jibForceCenter))
    keelTorque := momentOfInertiaScalar.MultVec(torqueOnPart(keelForce, keelForceCenter))
    rudderTorque := momentOfInertiaScalar.MultVec(torqueOnPart(rudderForce, rudderForceCenter))
    if printDebug {
        fmt.Printf("mainsailT %.3f jibT: %.3f keelT: %.3f rudderT: %.3f ", mainsailTorque.z, jibTorque.z, keelTorque.z, rudderTorque.z)
    }

    angularDragTorque := boat.w.Abs().MultVec(boat.w.Neg()).Mult(angularFriction)
    massCenter := vectorXyz{0, 0, 0}
    massCenter.x, massCenter.z = rotate(baseMassCenter.x, baseMassCenter.z, boat.yawRollHeading.y)
    gravityTorque := momentOfInertiaScalar.MultVec(torqueOnPart(gravitationalForce, massCenter))
    
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
        fmt.Printf("veloc: %.3f percent ", boat.v.Mag() / wind.speed)
    }
    
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
    wind.speed = 4.0
    
    return boat, wind
}

type autoTuneData struct {
    initialized bool
    startP, startI float64 // as yet unmodified
    startScore float64 // unmodified settings score
    newP, newI, newD float64 // current values to be used by controller
    roundNum int // the tuning trial/round we are on
    roundPart int // are we trialing a new p value (or i)
    deltaP, deltaI float64 // delta values that may be used in the current round
    deltaDeltaFactor float64 // multplies deltas each round to hone in on optimal values (< 1.0)
    scoreP, scoreI, combinedScore float64 // scores for the trials changing values of p and i
    inTrial bool
    readyForTrial bool // a trial is eligable to start if the threshold it met 
    endTrialCounts int // The number of consecutive calls the end trial requirement is met
    courseErrors []float64
    needsChange bool // awaiting a turn or change before another trial can be started
}

const (
    ROUND_BASLINE_1 = iota
    ROUND_BASLINE_2 = iota
    ROUND_PART_P_1 = iota
    ROUND_PART_P_2 = iota
    ROUND_PART_I_1 = iota
    ROUND_PART_I_2 = iota
    ROUND_COMBINED_1 = iota
    ROUND_COMBINED_2 = iota
)

// Inspiration taken from http://www.mstarlabs.com/control/self-tuning-pid.html
// We define a trial as a period of measurements that starts when the boat at a certain
// course error (it was just told to turn and then reached the designated error).
// It ends when the boat is told to turn again (larger course error).
// The tuning is done in rounds. Each round consists of three trials. The first modifying
// the p value only. The second the i value only. And then third moving both p and i according
// to the measured gradient. At the end of the round, the best of the three trials is taken.
// The next round will use smaller delta values to try and will direct the p and i values in the
// direction of downward optimization.
//
// data.deltaP and data.deltaI, the initial delta values to attempt and work done from need to be set
// before calling. The data.deltaDeltaFactor as well, which is the < 1.0 factor multipled each round to
// the deltaP and deltaI values. 
// returns true when new pid values should be read from autoTune.newP, .newI, .newD (start of new trial)
func autoTunePid(autoTune *autoTuneData, courseError, p, i float64) (newPidValues bool) {
    data := *autoTune
    newPidValues = false
    
    idGainRatio := 0.4 // ratio of d over i for output
    // if the courseError tops this value, a new trial is eligable to start when it reaches the start level
    retrialThreshold := 40.0
    // We will start taking data for a new trial when the value drops to this value...
    startTrialLevel := 40.0
    // To end the trial, the course error must be below this level this many consecutive times
    // If the retrial threshold is reached before this end trial requirement is met, the trial is thrown out.
    endTrialLevel := 5.0
    endTrialMinTimes := 50
    
    if !data.initialized {
        data.initialized = true
        data.startP = p
        data.startI = i
        data.roundNum = 1
        data.inTrial = false
        data.roundPart = ROUND_BASLINE_1
        data.courseErrors = make([]float64, 20)
        data.needsChange = true
    }
    
    if !data.inTrial && data.readyForTrial && math.Abs(courseError) <= startTrialLevel {
        data.inTrial = true
        data.readyForTrial = false
        data.needsChange = false
    } else if !data.inTrial && !data.readyForTrial && math.Abs(courseError) >= retrialThreshold {
        data.readyForTrial = true
        data.needsChange = false
        
        switch data.roundPart {
        case ROUND_BASLINE_1, ROUND_BASLINE_2:
            data.newP = data.startP
            data.newI = data.startI
        case ROUND_PART_P_1, ROUND_PART_P_2:
            data.newP = math.Max(data.startP + data.deltaP, 0.0)
            data.newI = data.startI
        case ROUND_PART_I_1, ROUND_PART_I_2:
            data.newP = data.startP
            data.newI = math.Max(data.startI + data.deltaI, 0.0)
        case ROUND_COMBINED_1, ROUND_COMBINED_2:
            if data.scoreP < data.startScore {
                data.newP = math.Max(data.startP + data.deltaP, 0.0)
            } else {
                data.newP = math.Max(data.startP - data.deltaP, 0.0)
            }
            if data.scoreI < data.startScore {
                data.newI = math.Max(data.startI + data.deltaI, 0.0)
            } else {
                data.newI = math.Max(data.startI - data.deltaI, 0.0)
            }
        }
        data.newD = data.newI * idGainRatio
        newPidValues = true
        
        fmt.Printf("Trial p: %.3f i: %.3f\n", data.newP, data.newI)
        
        data.courseErrors = data.courseErrors[:0]
    } else if data.inTrial && data.endTrialCounts > endTrialMinTimes {
        data.inTrial = false
        data.endTrialCounts = 0
        data.needsChange = true
        /*score := 0.0
        for _, val := range data.courseErrors {
            score += val * val
            fmt.Printf("%.2f ", val)
        }
        score /= float64(len(data.courseErrors))*/
        score := float64(len(data.courseErrors))
        
        fmt.Printf("Trial score: %.2f\n\n", score)
        
        switch data.roundPart {
        case ROUND_BASLINE_1:
            data.startScore = score
            data.roundPart = ROUND_BASLINE_2
        case ROUND_BASLINE_2:
            data.startScore += score
            data.roundPart = ROUND_PART_P_1
        case ROUND_PART_P_1:
            data.scoreP = score
            data.roundPart = ROUND_PART_P_2
        case ROUND_PART_P_2:
            data.scoreP += score
            data.roundPart = ROUND_PART_I_1
        case ROUND_PART_I_1:
            data.scoreI = score
            data.roundPart = ROUND_PART_I_2
        case ROUND_PART_I_2:
            data.scoreI += score
            data.roundPart = ROUND_COMBINED_1
        case ROUND_COMBINED_1:
            data.combinedScore = score
            data.roundPart = ROUND_COMBINED_2
        case ROUND_COMBINED_2:
            data.combinedScore += score
            
            // Choose the best of the options tried
            minScore := math.Min(math.Min(math.Min(data.startScore, data.scoreP), data.scoreI), data.combinedScore)
            var nextP, nextI, nextScore float64
            switch minScore {
            case data.startScore:
                nextP = data.startP
                nextI = data.startI
                nextScore = data.startScore
            case data.scoreP:
                nextP = math.Max(data.startP + data.deltaP, 0.0)
                nextI = data.startI
                nextScore = data.scoreP
            case data.scoreI:
                nextP = data.startP
                nextI = math.Max(data.startI + data.deltaI, 0.0)
                nextScore = data.scoreI
            case data.combinedScore:
                nextP = p
                nextI = i
                nextScore = data.combinedScore
            }
            
            fmt.Printf("Round settled on p: %.3f i: %.3f score: %.2f\n", nextP, nextI, nextScore)
            
            // this complete round is over
            data.roundNum++
            // Modify sign of deltas if going the wrong way
            if data.scoreP > data.startScore {
                data.deltaP = -data.deltaP
            }
            if data.scoreI > data.startScore {
                data.deltaI = -data.deltaI
            }
            // and decrease the delta values
            data.deltaP *= data.deltaDeltaFactor
            data.deltaI *= data.deltaDeltaFactor
            
            // assign new start values
            data.startP = nextP
            data.startI = nextI
            data.startScore = nextScore
            
            data.roundPart = ROUND_BASLINE_1
        }
    } else if data.inTrial {
        data.courseErrors = append(data.courseErrors, courseError)
        if math.Abs(courseError) <= endTrialLevel {
            data.endTrialCounts++
        } else {
            data.endTrialCounts = 0
        }
    }
    
    *autoTune = data
    return
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
    
    // manage keyboard control from command line
    go func() {
        for {
            byteBuffer := make([]byte, 1)
            bytesRead, _ := os.Stdin.Read(byteBuffer)
            // Reset the simulation when you press enter
            if bytesRead > 0 && byteBuffer[0] == 'r' {
                boat, wind = startupSettings()
                pidInitialized = false
            } else if bytesRead > 0 {
                // Change the boat coarse by the numbers 0 through 9.
                if byteBuffer[0] >= byte('0') && byteBuffer[0] <= byte('9') {
                    boat.course = CompassDirection(float64(byteBuffer[0] - byte('0')) * 22.5 - 180).Normalized()
                }
            }
        }
    }()

    timeStepOn := 0
    dt := 0.001 // time in seconds between physics calculation updates
    controlInterval := 0.1 // time in seconds between updates to rudder/sail headings
    reportInterval := 2.0 // time in seconds between outputs of data to console
    
    rudderRange := 60.0 // The maximum angle the rudder can be turned in either direction off straight
    highTurnAngle := 30.0 // The angle against the velocity which will maximize turning
    
    minimumPointingAngle := 50.0 // The lowest off wind angle we will sail
    runningAngle := 165.0 // The angle at or above which we put both sails all the way out
    
    maxTurnOffset := 25.0 // The angle away from extremes that maximizes turning with the sails

    // PID control of rudder for course correction
    previousInput := 0.0
    integralTerm := 0.0
    // control coefficients
    proportionK := 1.0
    integralK := 0.4
    derivativeK := 0.8
    
    // PI control of sails for course correction
    sailITerm := 0.0
    sailP := 1.2
    sailI := 1.2
    
    /*var autoTune autoTuneData
    autoTune.deltaP = 0.1
    autoTune.deltaI = 0.001
    autoTune.deltaDeltaFactor = 0.95*/
    
    for {
        if timeStepOn % int(controlInterval / dt) == 0 {
            // Allow autotuning of the pid control constants
            /*if autoTune.needsChange {
                if float64(boat.course) == -67.5 {
                    boat.course = CompassDirection(-112.5)
                } else {
                    boat.course = CompassDirection(-67.5)
                }
            }*/
        
            // Calculate the limits of PID output values
            boatHeading := float64(boat.getHeading())
            velocityDirection := float64(boat.getVelocityDirection())
            
            // Use the lower-bound of the rudder as a reference point 0.
            // So an angle here of 0 would correspond to the rudder being rudderRange
            // degrees counterclockwise from straight.
            
            // The angles of the rudder most effective at effecting rotation
            // are those that are at highTurnAngle to the velocity of the boat
            optimalTurnLeft := coerceAngleToRange((velocityDirection - highTurnAngle) - (boatHeading - 180) + rudderRange, 0, 360)
            optimalTurnRight := coerceAngleToRange((velocityDirection + highTurnAngle) - (boatHeading + 180) + rudderRange, -360, 0)
            turningCenter := (optimalTurnLeft + optimalTurnRight) / 2.0
            
            turnLeftDirection := math.Max(coerceAngleToRange(optimalTurnLeft, 0, 2 * rudderRange), turningCenter)
            turnLeftDirection = coerceAngleToRange(turnLeftDirection, 0, 2 * rudderRange)
            turnRightDirection := math.Min(coerceAngleToRange(optimalTurnRight, 0, 2 * rudderRange), turningCenter)
            turnRightDirection = coerceAngleToRange(turnRightDirection, 0, 2 * rudderRange)
            
            // Calculate maximum and minimum PID output values currently realizable with our rudder
            // This varies over time as the boat's orientation and velocity and the wind vary
            // Note that positive PID values are negative to indicate to turn left, and positive for right
            // But the rudder direction angles are high(positive) for left and low(negative) for right
            pidLeftLimit := lerp(turnLeftDirection, 0.0, 2 * rudderRange, 100, -100)
            pidRightLimit := lerp(turnRightDirection, 0.0, 2 * rudderRange, 100, -100)
            
            if timeStepOn % int(reportInterval / dt) == 0 {
                //fmt.Printf("boat: %.2f veloc: %.2f ", boatHeading, velocityDirection)
            }
        
            // Adjust rudder to get to desired heading
            // Rudder needs to be off from the boat velocity to generate
            // torque needed to turn the boat
            
            // PID control            
            pidSetpoint := float64(boat.course)
            // simulate some error in our ability to read our heading
            pidInput := float64(boat.getHeading()) + (rand.Float64() * 2 - 1) * 2 * 0
            
            if !pidInitialized {
                pidInitialized = true
                previousInput = pidInput
                integralTerm = 0.0
            }
            
            courseError := float64(RelativeDirection(pidSetpoint - pidInput).Normalized())
            
            // Allow autotuning of the pid control constants
            /*if autoTunePid(&autoTune, courseError, proportionK, integralK) {
                proportionK = autoTune.newP
                integralK = autoTune.newI
                derivativeK = autoTune.newD
                integralTerm = 0.0
            }*/
            
            // Make the rudder integral term conditional on the sails integral
            // term already being maxed out. That gives preference to the sails
            // fixing this kind of thing.
            if math.Abs(sailITerm) == 100.0 {
                integralTerm += courseError * integralK  * controlInterval
            } else {
                integralTerm = 0.0
            }
            // Prevent the integral term from exceeding the limits of the PID output
            if integralK != 0.0 {
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
            pidOutputValue := courseError * proportionK + integralTerm - inputDerivative * derivativeK / controlInterval
            previousInput = pidInput
            
            // conversion of that output value (from -100 to 100) to a physical rudder orientation
            var constrainedPidValue, newRudderDir float64
            if pidRightLimit > pidLeftLimit {
                // If we can't turn in the direction we want to, go neutral
                // Set us between the limits so rudder = 0.
                if (pidLeftLimit > 0 && pidOutputValue < 0) || (pidRightLimit < 0 && pidOutputValue > 0) {
                    constrainedPidValue = 0
                } else {
                    constrainedPidValue = math.Max(math.Min(pidOutputValue, pidRightLimit), pidLeftLimit)
                }
                newRudderDir = lerp(constrainedPidValue, -100, 100, turnLeftDirection, turnRightDirection)
            } else {
                if (pidLeftLimit < 0 && pidOutputValue > 0) || (pidRightLimit > 0 && pidOutputValue < 0) {
                    constrainedPidValue = 0
                } else {
                    constrainedPidValue = math.Max(math.Min(pidOutputValue, pidLeftLimit), pidRightLimit)
                }
                newRudderDir = lerp(constrainedPidValue, -100, 100, turnRightDirection, turnLeftDirection)
            }
            newRudderDir -= rudderRange
            
            boat.rudderDirection = RelativeDirection(newRudderDir).Normalized();
            
            if timeStepOn % int(reportInterval / dt) == 0 {
                fmt.Printf("Course Error: %.2f pterm: %.2f iterm: %.2f dterm: %.2f OutPID: %.2f rudder: %.2f ", courseError, courseError * proportionK, integralTerm, -inputDerivative * derivativeK / controlInterval, constrainedPidValue, boat.rudderDirection)
                //fmt.Printf("Course Error: %.2f ", courseError);
            }
            
            apparentWindDirection := float64(getApparentWindDirection(boat, wind))
            
            // make sails maximize force. See points of sail for reference.
            absoluteOffWind := coerceAngleToRange(float64(wind.direction) - float64(boat.getHeading()), -180, 180)
            boat.mainDirection = RelativeDirection(lerp(math.Abs(absoluteOffWind), minimumPointingAngle, runningAngle, 0.0, 90.0))
            boat.jibDirection = boat.mainDirection
            
            // correct sign
            if absoluteOffWind < 0 {
                boat.mainDirection = -boat.mainDirection
                boat.jibDirection = -boat.jibDirection
            }
            
            // If in running angle position, jib goes opposite main
            /*if math.Abs(absoluteOffWind) >= runningAngle {
                boat.jibDirection = -boat.jibDirection
            }*/
            
            /*mainDirection1 := coerceAngleToRange(apparentWindDirection + 90 - float64(boat.getHeading()), -180, 180)
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
            }*/
            
            // Use the sails to help us steer to our desired course
            // Turn to the wind by pulling the mainsail in and freeing the jib
            // Turn away by pulling the jib in and freeing the mailsail
            offWindAngle := coerceAngleToRange(apparentWindDirection - float64(boat.getHeading()), -180, 180)
            maxTurnAngle := offWindAngle
            if maxTurnAngle > 0 {
                maxTurnAngle -= maxTurnOffset
            } else {
                maxTurnAngle += maxTurnOffset
            }
            
            // PI control
            sailPTerm := courseError * sailP
            sailITerm += courseError * sailI * controlInterval
            // Prevent value from exceeding the maximum effect value at 100
            sailITerm = math.Max(math.Min(sailITerm, 100), -100)
            // Prevent value from having the wrong effect if the course changes a lot
            // We don't want inertia in the wrong way.
            if (sailITerm > 0 && courseError < -10) || (sailITerm < 0 && courseError > 10) {
                sailITerm = 0.0
            }
            
            turnValue := math.Abs(sailPTerm + sailITerm)
            
            if timeStepOn % int(reportInterval / dt) == 0 {
                fmt.Printf("sailP: %.3f sailI: %.3f sailOut: %.3f ", sailPTerm, sailITerm, turnValue)
            }
            
            // Putting the non zero checks on courseError here help stability
            // We really do not want to oscillate between towards and wawy from the wind.
            if (courseError > 5.0 && offWindAngle > 0.0) ||
               (courseError < -5.0 && offWindAngle < 0.0) {
                // Sail more towards the wind
                mainDir := lerp(turnValue, 0, 100, float64(boat.mainDirection), maxTurnOffset)
                jibDir := lerp(turnValue, 0, 100, float64(boat.jibDirection), math.Max(math.Min(maxTurnAngle, 90 - maxTurnOffset), -90 + maxTurnOffset))
                boat.mainDirection = RelativeDirection(mainDir)
                boat.jibDirection = RelativeDirection(jibDir)
                if timeStepOn % int(reportInterval / dt) == 0 {
                    //fmt.Printf("into main: %.2f jib: %.2f ", boat.mainDirection, boat.jibDirection)
                }
           } else {
                // Sail more away from the wind
                mainDir := lerp(turnValue, 0, 100, float64(boat.mainDirection), math.Max(math.Min(maxTurnAngle, 90 - maxTurnOffset), -90 + maxTurnOffset))
                jibDir := lerp(turnValue, 0, 100, float64(boat.jibDirection), maxTurnOffset)
                boat.mainDirection = RelativeDirection(mainDir)
                boat.jibDirection = RelativeDirection(jibDir)
                if timeStepOn % int(reportInterval / dt) == 0 {
                    //fmt.Printf("away main: %.2f jib: %.2f ", boat.mainDirection, boat.jibDirection)
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
                //fmt.Printf("md1: %.2f md2: %.2f appWind: %.2f %.2f mainsail: %.2f jib: %.2f ", mainDirection1, mainDirection2, wind.direction, apparentWindDirection, float64(boat.mainDirection), float64(boat.jibDirection))
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
