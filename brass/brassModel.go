package brass

import (
    "math"
    "sort"
    _ "fmt"
)

// representation of a brass instrument of a total length
// evenly divided into len(radii) - 1 sections, where radii specifies the radius
// that starts at that point. stepped specifies whether the sections are cylidrical (conical default).
// measurements in cm.
//type Brass struct {
//	Length  float64
//	Radii   []float64
//	Stepped bool
//}

type Brass struct {
    BaseSegments []Segment
	InsertSegments []Segment
	Stepped bool
}

func (brass *Brass) Length() float64 {
    length := 0.0
    for _, segment := range brass.BaseSegments {
        length += segment.Length()
    }
    for _, segment := range brass.InsertSegments {
        if segment.Enabled {
            length += segment.Length()
        }
    }
    return length
}

func (brass *Brass) MaxRadius() float64 {
    maxRadius := 0.0
    for _, segment := range brass.BaseSegments {
        maxRadius = math.Max(maxRadius, segment.MaxRadius())
    }
    for _, segment := range brass.InsertSegments {
        if segment.Enabled {
            maxRadius = math.Max(maxRadius, segment.MaxRadius())
        }
    }
    return maxRadius
}

type Coord struct {
    //Radius float64
    SelfDeltaRadius float64
    DeltaRadius float64
    X float64
    Mutable bool
}

type Segment struct {
    InsertionX float64
    Coords []Coord
    Enabled bool
    Resizable bool
}

func (segment *Segment) Length() float64 {
    // segment must be sorted first
    if len(segment.Coords) == 0 {
        return 0.0
    }
    return segment.Coords[len(segment.Coords) - 1].X - segment.Coords[0].X
}

func (segment *Segment) MaxRadius() float64 {
    maxRadius := 0.0
    currentRadius := 0.0
    for _, coords := range segment.Coords {
        currentRadius += coords.DeltaRadius
        maxRadius = math.Max(maxRadius, currentRadius + coords.SelfDeltaRadius)
        //maxRadius = math.Max(maxRadius, coords.Radius)
    }
    return maxRadius
}

func NewSegment(radius, length float64, subdivisions int) Segment {
    coords := make([]Coord, subdivisions + 1)

    deltaL := length / float64(subdivisions)
    for i := 0; i < len(coords); i++ {
        if i == 0 {
            coords[i].DeltaRadius = radius
        }
        //coords[i].Radius = radius
        coords[i].X = deltaL * float64(i)
        coords[i].Mutable = true
    }

    return Segment{0, coords, true, true}
}

type SimpleSegment struct {
    StartRadius float64
    EndRadius float64
    Length float64
}

type SimpleBrass struct {
    SimpleSegments []SimpleSegment
    Stepped bool
}

func ExtractOptimizationParameters(brass Brass) []float64 {
    parameters := make([]float64, 0, 64)

    parameters = append(parameters, brass.Length())

    for _, segment := range brass.BaseSegments {
        for _, coord := range segment.Coords {
            if coord.Mutable {
                parameters = append(parameters, coord.DeltaRadius, coord.SelfDeltaRadius)//coord.Radius)//, coord.X)
            }
        }
    }
    for _, segment := range brass.InsertSegments {
        if !segment.Enabled {
            continue
        }
        parameters = append(parameters, segment.InsertionX)
        for _, coord := range segment.Coords {
            if coord.Mutable {
                parameters = append(parameters, coord.DeltaRadius, coord.SelfDeltaRadius)//coord.Radius)//, coord.X)
            }
        }
    }

    return parameters
}

func ApplyOptimizationParameters(brass *Brass, parameters []float64) {
    i := 0

    totalLength := parameters[i]
    i++
    
    for j := range brass.BaseSegments {
        segment := &brass.BaseSegments[j]
        for k := range segment.Coords {
            coord := &segment.Coords[k]
            if coord.Mutable {
                //coord.Radius = parameters[i]
                //coord.Radius = math.Max(0.001, coord.Radius)
                coord.DeltaRadius = math.Max(0.001, parameters[i])
                coord.SelfDeltaRadius = math.Max(0.001, parameters[i + 1])
                i += 2
            }
        }
    }
    
    for j := range brass.InsertSegments {
        segment := &brass.InsertSegments[j]
        if !segment.Enabled {
            continue
        }
        segment.InsertionX = parameters[i]
        segment.InsertionX = math.Max(1e-5, segment.InsertionX)
        i++
        for k := range segment.Coords {
            coord := &segment.Coords[k]
            if coord.Mutable {
                //coord.Radius = parameters[i]
                //coord.Radius = math.Max(0.001, coord.Radius)
                coord.DeltaRadius = math.Max(0.001, parameters[i])
                coord.SelfDeltaRadius = math.Max(0.001, parameters[i + 1])
                i += 2
            }
        }
    }

    // Get a valid configuration from our changes, in case x values made coordinates switch positions
    brass.SortAll()

    // To reconstruct the desired total length, we subtract the nonresizable sections' lengths,
    // calculate the remaining as newly reconstructed, and then scale that up
    targetResizableTotalLength := totalLength
    currentResizableTotalLength := 0.0
    for _, segment := range brass.BaseSegments {
        if segment.Resizable {
            currentResizableTotalLength += segment.Length()
        } else {
            targetResizableTotalLength -= segment.Length()
        }
    }
    for _, segment := range brass.InsertSegments {
        if segment.Resizable {
            currentResizableTotalLength += segment.Length()
        } else {
            targetResizableTotalLength -= segment.Length()
        }
    }

    scaleFactor := targetResizableTotalLength / currentResizableTotalLength
    for j := range brass.BaseSegments {
        segment := &brass.BaseSegments[j]
        if segment.Resizable {
            for k := range segment.Coords {
                coord := &segment.Coords[k]
                coord.X *= scaleFactor
            }
        }
    }
    for j := range brass.InsertSegments {
        segment := &brass.InsertSegments[j]
        if segment.Resizable {
            for k := range segment.Coords {
                coord := &segment.Coords[k]
                coord.X *= scaleFactor
            }
        }
    }
    
}

// allows sorting
func (s Segment) Len() int { return len(s.Coords) }
func (s Segment) Less(i, j int) bool { return s.Coords[i].X < s.Coords[j].X }
func (s Segment) Swap(i, j int) { s.Coords[i], s.Coords[j] = s.Coords[j], s.Coords[i] }
func (b Brass) Len() int { return len(b.InsertSegments) }
func (b Brass) Less(i, j int) bool { return b.InsertSegments[i].InsertionX < b.InsertSegments[j].InsertionX }
func (b Brass) Swap(i, j int) { b.InsertSegments[i], b.InsertSegments[j] = b.InsertSegments[j], b.InsertSegments[i] }

func (b Brass) SortAll() {
    // We must sort the coordinates of each segment and the insertion segments themselves as well
    sort.Sort(b)
    for _, segment := range b.BaseSegments {
        sort.Sort(segment)
    }
    for _, segment := range b.InsertSegments {
        sort.Sort(segment)
    }
}

const maxRadiusSlope = 1

func limitedNewRadius(lastRadius, newRadius, length float64) float64 {
    // with undefined values, anything new passes
    if math.IsNaN(lastRadius + length) {
        return newRadius
    }

    radiusSlope := (newRadius - lastRadius) / length
    if math.Abs(radiusSlope) > maxRadiusSlope {
        return lastRadius + math.Copysign(maxRadiusSlope * length, radiusSlope)
        //newRadius += math.Copysign(maxDeltaRadius, deltaRadius) - deltaRadius
        //fmt.Println("would have had", nextRadius, "and now have", newRadius)
    }
    return newRadius
}

func appendSegment(newSegment Segment, segments []SimpleSegment, lastSegmentRadius float64) []SimpleSegment {
    lastRadius := lastSegmentRadius
    currentRadius := lastSegmentRadius
    segmentX := math.NaN()
    for _, coords := range newSegment.Coords {
        currentRadius += coords.DeltaRadius
        newRadius := currentRadius + coords.SelfDeltaRadius
        length := coords.X - segmentX
        if coords.Mutable {
            newLimitedRadius := limitedNewRadius(lastRadius, newRadius, length)
            // use this limiting to also limit the change transferred on from DeltaRadius to curentRadius
            currentRadius -= math.Max(0, (newRadius - newLimitedRadius) - coords.SelfDeltaRadius)
            newRadius = newLimitedRadius
        }
        if !math.IsNaN(length) {
            segments = append(segments, SimpleSegment{lastRadius, newRadius, length})
        }
        lastRadius = newRadius
        segmentX = coords.X
    }
    return segments
}

func BrassToSimpleBrass(brass Brass) SimpleBrass {
    // relies on brass being sorted

    segments := make([]SimpleSegment, 0, 64)

    // overall position through all base segments
    // this is what is used for insertion of insert segments
    overallBaseX := 0.0
    lastRadius := math.NaN()
    currentRadius := 0.0

    for _, baseSegment := range brass.BaseSegments {
        
        segmentX := math.NaN()
        for _, coords := range baseSegment.Coords {
            currentRadius += coords.DeltaRadius
            newRadius := currentRadius + coords.SelfDeltaRadius
            length := coords.X - segmentX
            if coords.Mutable {
                newLimitedRadius := limitedNewRadius(lastRadius, newRadius, length)
                // use this limiting to also limit the change transferred on from DeltaRadius to curentRadius
                currentRadius -= math.Max(0, (newRadius - newLimitedRadius) - coords.SelfDeltaRadius)
                newRadius = newLimitedRadius
            }
            if !math.IsNaN(length) {
                segments = append(segments, SimpleSegment{lastRadius, newRadius, length})
                
                // See if we have any segments to insert into the segment we just added
                for _, insertSegment := range brass.InsertSegments {
                    if !insertSegment.Enabled {
                        continue
                    }
                    distanceToInsert := insertSegment.InsertionX - overallBaseX
                    if distanceToInsert >= 0.0 && distanceToInsert <= length {
                        // shorten split segment
                        splittingSegment := &segments[len(segments) - 1]
                        splittingSegment.Length = distanceToInsert
                        overallBaseX += distanceToInsert
                        length -= distanceToInsert
                        // add insert
                        segments = appendSegment(insertSegment, segments, newRadius)
                        // add remainder of split segment
                        segments = append(segments, SimpleSegment{splittingSegment.StartRadius, splittingSegment.EndRadius, length})
                    }
                }

                overallBaseX += length
            }
            lastRadius = newRadius
            segmentX = coords.X
        }
    }
    
    return SimpleBrass{segments, brass.Stepped}
}
