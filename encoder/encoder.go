package main

import (
    "fmt"
    "math"
    "os"
    "image"
    "image/color"
    "image/draw"
    "image/png"
)

func rotate(x, y, theta float64) (x2, y2 float64) {
    x2 = x * math.Cos(theta) - y * math.Sin(theta)
    y2 = x * math.Sin(theta) + y * math.Cos(theta)
    return
}

// Finds the angle between two angles in degrees, a - b, angles modulo 2 Pi
func angleDiff(a, b float64) float64 {
    diff := math.Mod(a - b, 2 * math.Pi)
    if diff < -math.Pi {
        diff += 2 * math.Pi
    }
    if (diff > math.Pi) {
        diff -= 2 * math.Pi
    }
    return diff
}

type arc struct {
    center image.Point
    outerR, innerR float64
    startTheta, arcLength float64
}

func (a *arc) ColorModel() color.Model {
    return color.AlphaModel
}

func (a *arc) Bounds() image.Rectangle {
    if a.arcLength > math.Pi / 2 {
        // For larger arcs, just search the whole shape, it is harder to make a bounding rectangle
        return image.Rect(a.center.X - int(a.outerR), a.center.Y - int(a.outerR),
                          a.center.X + int(a.outerR), a.center.Y + int(a.outerR))
    } else {
        // Find min and max of all 4 end points
        innerStartX, innerStartY := math.Cos(a.startTheta) * a.innerR, math.Sin(a.startTheta) * a.innerR
        outerStartX, outerStartY := math.Cos(a.startTheta) * a.outerR, math.Sin(a.startTheta) * a.outerR
        endTheta := a.startTheta + a.arcLength
        innerEndX, innerEndY := math.Cos(endTheta) * a.innerR, math.Sin(endTheta) * a.innerR
        outerEndX, outerEndY := math.Cos(endTheta) * a.outerR, math.Sin(endTheta) * a.outerR
        return image.Rect(a.center.X + int(math.Min(math.Min(innerStartX, outerStartX), math.Min(innerEndX, outerEndX)) - 0.5),
                          a.center.Y + int(math.Min(math.Min(innerStartY, outerStartY), math.Min(innerEndY, outerEndY)) - 0.5),
                          a.center.X + int(math.Max(math.Max(innerStartX, outerStartX), math.Max(innerEndX, outerEndX)) + 0.5),
                          a.center.Y + int(math.Max(math.Max(innerStartY, outerStartY), math.Max(innerEndY, outerEndY)) + 0.5))
    }
}

func (a *arc) At(x, y int) color.Color {
    relX := float64(x - a.center.X)
    relY := float64(y - a.center.Y)
    rSquared := relX * relX + relY * relY
    theta := math.Atan2(relY, relX)
    if rSquared <= a.outerR * a.outerR && rSquared >= a.innerR * a.innerR {
        if a.arcLength >= 2 * math.Pi ||
           (angleDiff(theta, a.startTheta) >= 0 &&
           angleDiff(a.startTheta + a.arcLength, theta) >= 0) {
            //fmt.Println(x, y)
            return color.Alpha{255}
        }
    }
    return color.Alpha{0}
}

func drawArc(img *image.Gray, drawColor color.Color, centerX, centerY int, outerR, innerR, startTheta, arcLength float64) {
    colorImg := &image.Uniform{drawColor}
    a := arc{image.Point{centerX, centerY}, outerR, innerR, startTheta, arcLength}
    draw.DrawMask(img, img.Bounds(), colorImg, image.ZP, &a, image.ZP, draw.Over)
}

func main() {
    ppi := 200
    ppmm := float64(ppi) / 25.4
    
    receivingGearPitchR := 15.75 * ppmm
    receivingGearTeeth := 45
    pinionTeeth := 10
    outputDegreeResolution := 0.1
    encoderTickHeight := 1.0 * ppmm
    encoderTickEdgeOffset := 2.5 * ppmm
    
    width := int(receivingGearPitchR * 2)
    height := int(receivingGearPitchR * 2)

    // main encoder
    // make a white image
    white := color.Gray{255}
    black := color.Gray{0}
    encoderImg := image.NewGray(image.Rect(0, 0, width, height))
    draw.Draw(encoderImg, encoderImg.Bounds(), &image.Uniform{white}, image.ZP, draw.Src)
    
    // Enough ticks to give us 0.25 degree angular resolution of the output
    // There are 4 different ticks per physical encoder tick (because there are two lines of ticks)
    encoderTicks := int(math.Ceil(360.0 / outputDegreeResolution / math.Pow(float64(receivingGearTeeth) / float64(pinionTeeth), 2) / 4.0))
    fmt.Println("Encoder ticks: ", encoderTicks)
    
    encoderTickOuterR1 := receivingGearPitchR - encoderTickEdgeOffset
    encoderTickInnerR1 := encoderTickOuterR1 - encoderTickHeight
    for i := 0; i < encoderTicks; i++ {
        initialTheta := (float64(i) - 0.25) * 2.0 * math.Pi / float64(encoderTicks)
        arcLength := 2.0 * math.Pi / float64(encoderTicks) / 2.0
        drawArc(encoderImg, black, width / 2, height / 2, encoderTickOuterR1, encoderTickInnerR1, initialTheta, arcLength)
    }
    
    encoderTickOuterR2 := encoderTickOuterR1 - 2 * encoderTickHeight
    encoderTickInnerR2 := encoderTickOuterR2 - encoderTickHeight
    for i := 0; i < encoderTicks; i++ {
        initialTheta := float64(i) * 2.0 * math.Pi / float64(encoderTicks)
        arcLength := 2.0 * math.Pi / float64(encoderTicks) / 2.0
        drawArc(encoderImg, black, width / 2, height / 2, encoderTickOuterR2, encoderTickInnerR2, initialTheta, arcLength)
    }
    
    // write to file as png
    f, err := os.Create("encoder.png")
    if err != nil {
        panic(err)
    }
    defer f.Close()
    png.Encode(f, encoderImg)
    
    // single ring encoder for homing
    counterColor := color.Gray{0}
    clockwiseColor := color.Gray{64}
    centerColor := color.Gray{128}
    leftSideColor := color.Gray{192}
    rightSideColor := color.Gray{255}
    homingImg := image.NewGray(image.Rect(0, 0, width, height))
    draw.Draw(homingImg, homingImg.Bounds(), &image.Uniform{white}, image.ZP, draw.Src)
    
    drawArc(homingImg, counterColor, width / 2, height / 2, encoderTickOuterR1, encoderTickInnerR2, 0, 2 * math.Pi)
    drawArc(homingImg, clockwiseColor, width / 2, height / 2, encoderTickOuterR1, encoderTickInnerR2, 0, math.Pi / 4)
    drawArc(homingImg, clockwiseColor, width / 2, height / 2, encoderTickOuterR1, encoderTickInnerR2, -math.Pi / 2, math.Pi / 4)
    drawArc(homingImg, clockwiseColor, width / 2, height / 2, encoderTickOuterR1, encoderTickInnerR2, math.Pi / 2, math.Pi / 2)
    drawArc(homingImg, centerColor, width / 2, height / 2, encoderTickOuterR1, encoderTickInnerR2, -math.Pi / 180 * 2, math.Pi / 180 * 4)
    drawArc(homingImg, leftSideColor, width / 2, height / 2, encoderTickOuterR1, encoderTickInnerR2, -math.Pi / 180 * 92, math.Pi / 180 * 4)
    drawArc(homingImg, rightSideColor, width / 2, height / 2, encoderTickOuterR1, encoderTickInnerR2, math.Pi / 180 * 88, math.Pi / 180 * 4)
    
    //drawArc(homingImg, color.Gray{248}, width / 2, height / 2, encoderTickOuterR1, encoderTickInnerR2, 0, math.Pi)
    //drawArc(homingImg, color.Gray{140}, width / 2, height / 2, encoderTickOuterR1, encoderTickInnerR2, math.Pi, math.Pi)
    //drawArc(homingImg, color.Gray{0}, width / 2, height / 2, encoderTickOuterR1, encoderTickInnerR2, -math.Pi / 180 * 1, math.Pi / 180 * 2)
    
    f2, err2 := os.Create("homing.png")
    if err2 != nil {
        panic(err2)
    }
    defer f2.Close()
    png.Encode(f2, homingImg)
    
}