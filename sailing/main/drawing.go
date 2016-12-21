package main

import (
	"image"
	"image/color"
	"image/draw"
	"image/png"
    "io"
    "math"
)

func drawLine(startX, startY, endX, endY int, col color.Gray, img *image.Gray) {
  // Bresenham's
  var cx = startX;
  var cy = startY;
 
  var dx = endX - cx;
  var dy = endY - cy;
  if dx<0 { dx = 0-dx; }
  if dy<0 { dy = 0-dy; }
 
  var sx int;
  var sy int;
  if cx < endX { sx = 1; } else { sx = -1; }
  if cy < endY { sy = 1; } else { sy = -1; }
  var err = dx-dy;
 
  var n int;
  for n=0;n<1000;n++ {
    img.SetGray(cx, cy, col)
    if((cx==endX) && (cy==endY)) {return;}
    var e2 = 2*err;
    if e2 > (0-dy) { err = err - dy; cx = cx + sx; }
    if e2 < dx     { err = err + dx; cy = cy + sy; }
  }
}

// Draw an arrow centered at x, y
func drawArrow(x, y, length int, compDir float64, img *image.Gray) {
    headingX, headingY := rotate(0, float64(length) / 2, compDir * math.Pi / 180)
    arrowX, arrowY := rotate(0, float64(length) / 3, (compDir + 45) * math.Pi / 180)
    arrowX, arrowY = arrowX - headingX, arrowY - headingY
    drawLine(x - int(headingX), y - int(headingY), x + int(headingX), y + int(headingY), color.Gray{0}, img)
    drawLine(x - int(headingX), y - int(headingY), x + int(arrowX), y + int(arrowY), color.Gray{0}, img)
}

// Draw an arrow with the point at x, y
func drawArrowFromTop(x, y, length int, compDir float64, img *image.Gray) {
    headingX, headingY := rotate(0, float64(length), compDir * math.Pi / 180)
    arrowX, arrowY := rotate(0, float64(length) / 3, (compDir + 45) * math.Pi / 180)
    drawLine(x, y, x + int(headingX), y + int(headingY), color.Gray{0}, img)
    drawLine(x, y, x + int(arrowX), y + int(arrowY), color.Gray{0}, img)
}

// Draw an arrow with the base at x, y
func drawArrowFromBase(x, y, length int, compDir float64, img *image.Gray) {
    headingX, headingY := rotate(0, float64(-length), compDir * math.Pi / 180)
    arrowX, arrowY := rotate(0, float64(-length) / 3, (compDir + 45) * math.Pi / 180)
    arrowX, arrowY = headingX - arrowX, headingY - arrowY
    drawLine(x, y, x + int(headingX), y + int(headingY), color.Gray{0}, img)
    drawLine(x + int(headingX), y + int(headingY), x + int(arrowX), y + int(arrowY), color.Gray{0}, img)
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
    drawArrowFromTop(x + int(jibX), y + int(jibY), int(jibLength), boat.jibDirection + boatHeading, img);
    
    mainsailX, mainsailY := rotate(0.0, -5.0, boatAngle)
    drawArrowFromTop(x + int(mainsailX), y + int(mainsailY), int(mainsailLength), boat.mainDirection + boatHeading, img)
    
    rudderX, rudderY := rotate(0.0, mainsailLength - 5.0 + interspace, boatAngle)
    drawArrowFromBase(x + int(rudderX), y + int(rudderY), int(rudderLength), boat.rudderDirection + boatHeading + 180, img)

    rollLength := boat.yawRollHeading.y / (math.Pi / 2) * totalLength
    drawArrowFromBase(x, y, int(rollLength), boatHeading + 90, img)
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
    drawArrow(width / 2 - 15, 25, 30, 0, img)
    drawArrow(width / 2 - 15, height - 25, 30, 180, img)
    drawArrow(25, height / 2 - 15, 30, -90, img)
    drawArrow(width - 25, height / 2 - 15, 30, 90, img)

    drawArrow(25, 25, 40, boat.course, img);
    drawArrow(width - 25, 25, 40, wind.direction, img);

    drawBoat(boat, width / 2, height / 2, img);
    
    drawArrow(width / 2 + 50, height / 2 - 50, 40, boat.getVelocityDirection(), img)
    
    // Calculate the apparent wind heading
    apparentWindDirection := getApparentWindDirection(boat, wind)
    drawArrow(width - 80, 25, 40, apparentWindDirection, img)

    png.Encode(outputWriter, img)
}