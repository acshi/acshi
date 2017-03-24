package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"net/http"
	"os"
	"runtime/pprof"
	"time"
)

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")

const outputImageSize = 400

type boatData struct {
	p, v, yawRollHeading, w                      vectorXyz
	mainDirection, jibDirection, rudderDirection float64
	course                                       float64
}

type windData struct {
	direction float64
	speed     float64
}

type simData struct {
	timeStepOn int
	boat       boatData
	wind       windData

	dt, controlInterval, reportInterval        float64
	minimumPointingAngle, runningAngle         float64
	maxSailEase, maxTurnOffset                 float64
	rudderRange, highTurnAngle, minClosedAngle float64

	pidInitialized bool

	previousInput, integralTerm         float64
	proportionK, integralK, derivativeK float64
	rudderIDecay                        float64

	sailITerm, sailP, sailI float64

	lastOffWindAngle float64

	shouldPrint bool
}

func (boat boatData) getHeading() float64 {
	return normalizeAngle(boat.yawRollHeading.z / math.Pi * 180)
}

func (boat boatData) getVelocityDirection() float64 {
	return normalizeAngle(math.Atan2(boat.v.y, boat.v.x)/math.Pi*180 + 90)
}

func getApparentWindDirection(boat boatData, wind windData) float64 {
	// Calculate the apparent wind heading
	windVector := vecFromHeading(wind.direction).Mult(-wind.speed)
	apparentWindVector := windVector.Sub(boat.v)
	apparentWindDirection := apparentWindVector.GetHeading() + 180
	return normalizeAngle(apparentWindDirection)
}

func startupSettings() simData {
	var s simData

	s.boat.p = vectorXyz{0, 0, 0}
	s.boat.v = vectorXyz{0, 0, 0}
	s.boat.yawRollHeading = vectorXyz{0, 0, -135 / 180.0 * math.Pi}
	s.boat.w = vectorXyz{0, 0, 0}

	s.boat.course = -135
	s.boat.mainDirection = -30
	s.boat.jibDirection = -30
	s.boat.rudderDirection = 0

	s.wind.direction = 45
	s.wind.speed = 4.0

	s.timeStepOn = 0
	s.dt = 0.001            // time in seconds between physics calculation updates
	s.controlInterval = 0.1 // time in seconds between updates to rudder/sail headings
	s.reportInterval = 0.5  // time in seconds between outputs of data to console

	s.minimumPointingAngle = 50.0 // The lowest off wind angle we will sail
	s.runningAngle = 165.0        // The angle at or above which we put both sails all the way out
	s.maxSailEase = 90.0          // Angle from straight back to max ease of main (and for now, jib)
	s.maxTurnOffset = 5.0         // The angle away from optimal sail position used in turning
	s.rudderRange = 60.0          // The maximum angle the rudder can be turned in either direction off straight
	s.highTurnAngle = 30.0        // The angle against the velocity which will maximize turning
	s.minClosedAngle = 15.0       // The minimum angle to keep the jib open

	// PID control of rudder for course correction
	s.pidInitialized = false

	s.previousInput = 0.0
	s.integralTerm = 0.0
	// control coefficients
	s.proportionK = 4.0
	s.integralK = 0.4
	s.derivativeK = 2.0

    s.rudderIDecay = 0.95 // decay factor per second of rudder integral when sails could contribute more

	// PI control of sails for course correction
	s.sailITerm = 0.0
	s.sailP = 5.0
	s.sailI = 5.0

	if true {
		s.proportionK, s.integralK, s.derivativeK = 7.94447, 0.25778, 3.83347
        s.sailP, s.sailI = 8.90936, 0.39767
		// Best values: 5.72722 3.20107 2.12367 with score: 32.07324
		// Best values: 8.56316 -0.01802 3.26637 with score: 27.81168
		// Best values: 14.20511 10.73425 6.22322 with score: 21.01191
		// Best values: 6.84416 5.97056 2.54453 with score: 22.40385
        // Best values: 12.53007 7.61016 4.70147 13.31001 0.12787 with score: 20.72589
        // Best values: 5.94448, 0.04463, 2.83346 6.89076, 0.29762 with score: 22.79207
        // Best values: 7.94447, 0.25778, 3.83347 8.90936, 0.39767 with score: 22.61423
	}

	// Used to reduce shock from abrupt changes
	// we need to know if we might be suddenly shifting the sails from one side to the other
	s.lastOffWindAngle = 0.0

	// Whether debug info should be printed to the screen (turned off by autotune)
	s.shouldPrint = true

	return s
}

func controlUpdate(s *simData) {
	// Calculate the limits of PID output values
	boatHeading := s.boat.getHeading()
	velocityDirection := s.boat.getVelocityDirection()
	apparentWindDirection := getApparentWindDirection(s.boat, s.wind)

	if s.shouldPrint && s.timeStepOn%int(s.reportInterval/s.dt) == 0 {
		fmt.Printf("course: %.2f apwind: %.2f boat: %.2f vdir: %.2f vspeed: %.2f ", s.boat.course, apparentWindDirection, boatHeading, velocityDirection, s.boat.v.Mag())
	}

	// Use the straight back from the back of the boat as a reference point 0.
	// So a positive angle results in a left turn from the rudder, assuming forward motion

	// do not even attempt to sail into the wind
	if math.Abs(normalizeAngle(s.boat.course-apparentWindDirection)) < s.minimumPointingAngle {
		return
	}

	// The angles of the rudder most effective at effecting rotation
	// are those that are at highTurnAngle to the velocity of the boat
	// leeway is just the difference between the physical heading of the boat and the direction it is actually moving in
	leeway := normalizeAngle(velocityDirection - boatHeading)

	var optimalTurnLeft, optimalTurnRight float64
	var turnLeftDirection, turnRightDirection float64
	var pidLeftLimit, pidRightLimit float64
	// The three cases below are that 1. there are good positions for the rudder to turn
	// 2. the boat is moving backwards, so the rudder can go the opposite direction
	// 3. no matter how the rudder is positioned, the boat will be turned and we don't
	// have rudder control, so just keep the rudder neutral.

	// Calculate maximum and minimum PID output values currently realizable with our rudder
	// This varies over time as the boat's orientation and velocity and the wind vary
	// Note that PID values are negative to indicate left turning, and positive for right
	// But the rudder direction angles depend also on the leeway
	if math.Abs(leeway) <= s.rudderRange {
		optimalTurnLeft = leeway + s.highTurnAngle
		optimalTurnRight = leeway - s.highTurnAngle
		turnLeftDirection = math.Min(optimalTurnLeft, s.rudderRange)
		turnRightDirection = math.Max(optimalTurnRight, -s.rudderRange)
		pidLeftLimit = lerp(turnLeftDirection, optimalTurnRight, optimalTurnLeft, 100, -100)
		pidRightLimit = lerp(turnRightDirection, optimalTurnRight, optimalTurnLeft, 100, -100)
	} else if math.Abs(leeway) >= 180-s.rudderRange {
		oppositeTurnCenter := normalizeAngle(velocityDirection - (boatHeading + 180))
		optimalTurnLeft = oppositeTurnCenter - s.highTurnAngle
		optimalTurnRight = oppositeTurnCenter + s.highTurnAngle
		turnLeftDirection = math.Max(optimalTurnLeft, -s.rudderRange)
		turnRightDirection = math.Min(optimalTurnRight, s.rudderRange)
		pidLeftLimit = lerp(turnLeftDirection, optimalTurnLeft, optimalTurnRight, -100, 100)
		pidRightLimit = lerp(turnRightDirection, optimalTurnLeft, optimalTurnRight, -100, 100)
	} else {
		optimalTurnLeft = 0
		optimalTurnRight = 0
		turnLeftDirection = 0
		turnRightDirection = 0
		pidLeftLimit = 0
		pidRightLimit = 0
	}

	if s.shouldPrint && s.timeStepOn%int(s.reportInterval/s.dt) == 0 {
		fmt.Printf("optimalLeft: %.2f optimalRight: %.2f left1: %.2f right1: %.2f ", optimalTurnLeft, optimalTurnRight, velocityDirection+s.highTurnAngle-boatHeading+180, velocityDirection-s.highTurnAngle-boatHeading+180)
	}

	// Adjust rudder to get to desired heading
	// Rudder needs to be off from the boat velocity to generate
	// torque needed to turn the boat

	// PID control
	pidSetpoint := float64(s.boat.course)
	pidInput := float64(s.boat.getHeading())

	if !s.pidInitialized {
		s.pidInitialized = true
		s.previousInput = pidInput
		s.integralTerm = 0.0
	}

	courseError := angleToRange(pidSetpoint-pidInput, -180, 180)

	// Make adding to the rudder integral term conditional on the sails integral
	// term being maxed out. That gives preference to the sails fixing this kind of thing.
	// However, to prevent impulses, we allow changes that would decrease magnitude
	// and also cause slow decay when the sails could contribute more
	if math.Abs(s.sailITerm) >= 95.0 ||
		(s.integralTerm > 0.0 && courseError < 0.0 || s.integralTerm < 0.0 && courseError > 0.0) { // if rudderITerm would decrease in simple PID alg, allow it
		s.integralTerm += courseError * s.integralK * s.controlInterval
	} else {
		// Allow the integral to smoothly decay
		s.integralTerm *= math.Pow(s.rudderIDecay, s.controlInterval)
	}
	// Prevent the integral term from exceeding the limits of the PID output (windup)
	if s.integralTerm >= pidRightLimit || s.integralTerm <= pidLeftLimit {
		s.integralTerm = math.Max(math.Min(s.integralTerm, pidRightLimit), pidLeftLimit)
	}

	inputDerivative := angleToRange(pidInput-s.previousInput, -180, 180)
	pidOutputValue := courseError*s.proportionK + s.integralTerm - inputDerivative*s.derivativeK/s.controlInterval
	s.previousInput = pidInput

	// conversion of that output value (from -100 to 100) to a physical rudder orientation
	constrainedPIDValue := math.Max(math.Min(pidOutputValue, pidRightLimit), pidLeftLimit)
	newRudderDir := lerp(constrainedPIDValue, -100, 100, optimalTurnLeft, optimalTurnRight)

	// Smoothen out the change to the rudder
	s.boat.rudderDirection = lerp(4.0*s.controlInterval, 0, 1, s.boat.rudderDirection, newRudderDir)

	if s.shouldPrint && s.timeStepOn%int(s.reportInterval/s.dt) == 0 {
		fmt.Printf("\n\terror: %.2f pterm: %.2f iterm: %.2f dterm: %.2f Out: %.2f rudder: %.2f ", courseError, courseError*s.proportionK, s.integralTerm, -inputDerivative*s.derivativeK/s.controlInterval, constrainedPIDValue, s.boat.rudderDirection)
		//fmt.Printf("Course Error: %.2f ", courseError);
	}

	//
	// SAIL CONTROL
	//

	// Make sails maximize force. See points of sail for reference.
	// We use this as the base position for the sails from which changes are made for turning
	offWindAngle := angleToRange(apparentWindDirection-float64(s.boat.getHeading()), -180, 180)
	newMainDirection := lerp(math.Abs(offWindAngle), s.minimumPointingAngle, s.runningAngle, 0.0, 90.0) //A first (linear) approximation of where the sails should go, note right now it's always positive
	newJibDirection := newMainDirection

	// correct the sign
	if offWindAngle < 0 {
		newMainDirection = -newMainDirection
		newJibDirection = -newJibDirection
	}

	// Use the sails to help us steer to our desired course
	// Turn to the wind by pulling the mainsail in and freeing the jib
	// Turn away by pulling the jib in and freeing the mailsail

	// Make the sign of course error associate negative with away from wind, and positive into wind
	// This lets the control loop evenly progress back and forth
	sailCourseError := courseError

	if offWindAngle < 0.0 { // We're on starboard
		sailCourseError = -sailCourseError
	}

	// PI control
	sailPTerm := sailCourseError * s.sailP
	s.sailITerm += sailCourseError * s.sailI * s.controlInterval
	// Prevent value from exceeding the maximum effect at 100
	s.sailITerm = math.Max(math.Min(s.sailITerm, 100), -100)

	turnValue := sailPTerm + s.sailITerm

	if s.shouldPrint && s.timeStepOn%int(s.reportInterval/s.dt) == 0 {
		fmt.Printf("sailP: %.3f sailI: %.3f sailOut: %.3f ", sailPTerm, s.sailITerm, turnValue)
	}

	// If we are relatively close to the sails switching direction
	// We don't want that to happen besides when necessary to avoid oscillation
	if (offWindAngle < -170.0 || offWindAngle > 170.0) &&
		(offWindAngle > 0.0 && s.lastOffWindAngle < 0.0 || offWindAngle < 0.0 && s.lastOffWindAngle > 0.0) {
		offWindAngle = s.lastOffWindAngle
		newMainDirection = math.Copysign(newMainDirection, s.lastOffWindAngle)
		newJibDirection = math.Copysign(newJibDirection, s.lastOffWindAngle)
	}

	// The positions that given maximum turning torque from sails
	var maxOpenPose, maxClosePose float64
	if newMainDirection > 0 {
		maxOpenPose = math.Min(newMainDirection+s.maxTurnOffset, s.maxSailEase)
		maxClosePose = math.Max(newMainDirection-s.maxTurnOffset, 0)
	} else {
		maxOpenPose = math.Max(newMainDirection-s.maxTurnOffset, -s.maxSailEase)
		maxClosePose = math.Min(newMainDirection+s.maxTurnOffset, 0)
	}

	if turnValue > 0.0 {
		// Sail more towards the wind
		newMainDirection = lerp(turnValue, 0, 100, newMainDirection, maxClosePose)
		newJibDirection = lerp(turnValue, 0, 100, newJibDirection, maxOpenPose)
		if s.shouldPrint && s.timeStepOn%int(s.reportInterval/s.dt) == 0 {
			fmt.Printf("offwind: %.2f maxopen: %.2f _into_ main: %.2f jib: %.2f ", offWindAngle, maxOpenPose, newMainDirection, newJibDirection)
		}
	} else {
		// Sail more away from the wind
		newMainDirection = lerp(-turnValue, 0, 100, newMainDirection, maxOpenPose)
		newJibDirection = lerp(-turnValue, 0, 100, newJibDirection, maxClosePose)
		if s.shouldPrint && s.timeStepOn%int(s.reportInterval/s.dt) == 0 {
			fmt.Printf("offwind: %.2f maxopen: %.2f _away_ main: %.2f jib: %.2f ", offWindAngle, maxOpenPose, newMainDirection, newJibDirection)
		}
	}

	// To maximize flow around the sails, we keep the jib open at least some
	if math.Abs(newJibDirection) < s.minClosedAngle {
		newJibDirection = math.Copysign(s.minClosedAngle, newJibDirection)
	}

	s.lastOffWindAngle = offWindAngle
	// Smoothen out the change to the sails
	s.boat.mainDirection = lerp(4.0*s.controlInterval, 0, 1, s.boat.mainDirection, newMainDirection)
	s.boat.jibDirection = lerp(4.0*s.controlInterval, 0, 1, s.boat.jibDirection, newJibDirection)
}

func performTimestep(s *simData) {
	if s.timeStepOn%int(s.controlInterval/s.dt) == 0 {
		controlUpdate(s)
	}

	if s.shouldPrint && s.timeStepOn%int(s.reportInterval/s.dt) == 0 {
		s.boat, s.wind = physicsUpdate(s.boat, s.wind, s.dt, true)
		fmt.Printf("\n")
	} else {
		s.boat, s.wind = physicsUpdate(s.boat, s.wind, s.dt, false)
	}
	s.timeStepOn++
}

var imageBoat boatData
var imageWind windData

func main() {
	flag.Parse()
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	//newtonOptimize(func(x []float64) float64 { return 2 * math.Pow(x[0], 4) + math.Pow(x[1], 2) + x[0] * x[1] }, []float64 {2, 2})
	//autotunePID()
    //autotunePIDSailPI()
	//return

	sim := startupSettings()

	http.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {
		writeSailingImage(imageBoat, imageWind, w)
	})
	go http.ListenAndServe(":8080", nil)

	simSpeed := 1
	simSpeedWaitCount := 0

	// manage keyboard control from command line
	go func() {
		for {
			byteBuffer := make([]byte, 1)
			bytesRead, _ := os.Stdin.Read(byteBuffer)
			// Reset the simulation when you press enter
			if bytesRead == 0 {
				continue
			}
			switch byteBuffer[0] {
			case 'r':
				sim = startupSettings()
				simSpeed = 1
			case 'q':
				os.Exit(0)
			case 'f':
				simSpeed = 100
				sim.reportInterval = 0.5 * 100
			case 's':
				simSpeed = 1
				sim.reportInterval = 0.5
			default:
				// Change the boat coarse by the numbers 0 through 9.
				if byteBuffer[0] >= byte('0') && byteBuffer[0] <= byte('9') {
					sim.boat.course = normalizeAngle(float64(byteBuffer[0]-byte('0'))*22.5 - 180)
				}
			}
		}
	}()

	for {
		performTimestep(&sim)
		imageBoat, imageWind = sim.boat, sim.wind

		simSpeedWaitCount++
		if simSpeedWaitCount >= simSpeed {
			simSpeedWaitCount = 0
			time.Sleep(time.Duration(float64(time.Second) * sim.dt))
		}
	}
}
