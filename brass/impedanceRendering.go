package brass

import (
	"image"
	"image/color"
	"image/draw"
	"math"
	"acshi/functional"
	"fmt"
)



func max(xs []float64) float64 {
    return functional.Fold(math.Max, xs)
}

func min(xs []float64) float64 {
    return functional.Fold(math.Min, xs)
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

func DrawImpedancesOnImage(impedances, targetFrequencies []float64, targetImpedance float64, img *image.Gray) {
    minVal := min(impedances)
    maxVal := max(impedances)
    
    bounds := img.Bounds()
    width := bounds.Max.X - bounds.Min.X
    height := bounds.Max.Y - bounds.Min.Y

    widthPerIndex := float64(width) / float64(len(impedances))
    heightPerImpedance := float64(height) / (maxVal - minVal)
    
    gray := color.Gray{200}
    black := color.Gray{0}

    // target height for target frequencies
    targetHeight := int(float64(height) - targetImpedance * heightPerImpedance)
    draw_line(0, targetHeight, width, targetHeight, gray, img)
    draw_line(width / 2, targetHeight, width, targetHeight, gray, img)

    for _, targetFreq := range targetFrequencies {
        x := int((targetFreq - 1) * widthPerIndex)
        draw_line(x, 0, x, height - 1, gray, img)
    }

    var lastX int
    var lastY int
    for i := 0; i < len(impedances); i++ {
        x := int(float64(i) * widthPerIndex)
        y := int(float64(height) - impedances[i] * heightPerImpedance)

        if i > 0 {
            draw_line(lastX, lastY, x, y, black, img)
        }

        lastX = x
        lastY = y
    }
}

func DrawBrassOnImage(brass Brass, img *image.Gray, drawBounds image.Rectangle)  {
    width := drawBounds.Max.X - drawBounds.Min.X
    height := drawBounds.Max.Y - drawBounds.Min.Y

    yZero := height / 2 + drawBounds.Min.Y
    xZero := drawBounds.Min.X

    simpleBrass := BrassToSimpleBrass(brass)

    simpleMaxRadius := 0.0
    for _, segment := range simpleBrass.SimpleSegments {
        simpleMaxRadius = math.Max(simpleMaxRadius, segment.StartRadius)
        simpleMaxRadius = math.Max(simpleMaxRadius, segment.EndRadius)
    }

    totalLength := brass.Length()
    maxRadius := brass.MaxRadius()

    fmt.Println("Pre: max radius", maxRadius, "total length", totalLength)
    
    maxRadius = simpleMaxRadius
    fmt.Println("Post: max radius", maxRadius)

    xPerUnit := float64(width) / totalLength
    yPerUnit := float64(height / 2) / maxRadius

    black := color.Gray{0}
    lightGray := color.Gray{250}
    darkGray := color.Gray{160}
    // shade the area
    draw.Draw(img, drawBounds, &image.Uniform{lightGray}, image.ZP, draw.Src)

    // center line
    //draw_line(0, yZero, width, yZero, gray, img)

    currentX := 0.0
    lastX := 0
    lastY := 0
    for i, segment := range simpleBrass.SimpleSegments {
        startY := int(yPerUnit * segment.StartRadius)
        if i > 0 {
            draw_line(lastX + xZero, lastY + yZero, lastX + xZero, startY + yZero, black, img)
            draw_line(lastX + xZero, -lastY + yZero, lastX + xZero, -startY + yZero, black, img)
        }
        endY := int(yPerUnit * segment.EndRadius)
        currentX += xPerUnit * segment.Length
        newX := int(currentX)
        
        draw_line(lastX + xZero, startY + yZero, newX + xZero, endY + yZero, black, img)
        draw_line(lastX + xZero, -startY + yZero, newX + xZero, -endY + yZero, black, img)

        // tick marks
        draw_line(newX + xZero, endY + yZero, newX + xZero, endY + yZero - 5, darkGray, img)
        draw_line(newX + xZero, -endY + yZero, newX + xZero, -endY + yZero + 5, darkGray, img)
        
        lastX = newX
        lastY = endY
    }
}

