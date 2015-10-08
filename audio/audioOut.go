package audio

import (
	"errors"
	"fmt"
	"syscall"
	"time"
	"unsafe"
)

type WaveOutDevice struct {
	handle   uintptr
	format   WAVEFORMATEX
	callback func([]int16)
	quit     (chan bool)
	quitted  (chan bool)
}

var (
	modWinmm                   = syscall.NewLazyDLL("winmm.dll")
	procWaveOutGetNumDevs      = modWinmm.NewProc("waveOutGetNumDevs")
	procWaveOutOpen            = modWinmm.NewProc("waveOutOpen")
	procWaveOutPrepareHeader   = modWinmm.NewProc("waveOutPrepareHeader")
	procWaveOutWrite           = modWinmm.NewProc("waveOutWrite")
	procWaveOutUnprepareHeader = modWinmm.NewProc("waveOutUnprepareHeader")
	procWaveOutReset           = modWinmm.NewProc("waveOutReset")
	procWaveOutClose           = modWinmm.NewProc("waveOutClose")

	mmErrors = map[int]string{
		MMSYSERR_ALLOCATED:    "Specified resource is already allocated.",
		MMSYSERR_BADDEVICEID:  "Specified device identifier is out of range.",
		MMSYSERR_INVALHANDLE:  "Specified device handle is invalid.",
		MMSYSERR_NODRIVER:     "No device driver is present.",
		MMSYSERR_NOMEM:        "Unable to allocate or lock memory.",
		MMSYSERR_NOTSUPPORTED: "Specified device is synchronous and does not support pausing.",
		WAVERR_BADFORMAT:      "Attempted to open with an unsupported waveform-audio format.",
		WAVERR_SYNC:           "The device is synchronous but waveOutOpen was called without using the WAVE_ALLOWSYNC flag.",
		WAVERR_STILLPLAYING:   "There are still buffers in the queue.",
	}
)

func WaveOutGetNumDevs() int {
	r, _, _ := procWaveOutGetNumDevs.Call()
	return int(r)
}

func checkMmErrors(op string, r uintptr, e error) error {
	if errno, ok := e.(syscall.Errno); ok && errno != 0 {
		return errors.New(fmt.Sprint("audio: ", op, ": ", e))
	}
	if r != MMSYSERR_NOERROR {
		return errors.New(fmt.Sprint("audio: ", op, ": ", mmErrors[int(r)]))
	}
	return nil
}

func WaveOutOpen(deviceNumber int) (WaveOutDevice, error) {
	var handle uintptr
	var format WAVEFORMATEX
	format.wFormatTag = WAVE_FORMAT_PCM
	format.nChannels = 1
	format.nSamplesPerSec = 44100
	format.wBitsPerSample = 16
	format.cbSize = 0
	format.nBlockAlign = format.nChannels * format.wBitsPerSample / 8
	format.nAvgBytesPerSec = format.nSamplesPerSec * uint32(format.nBlockAlign)

	r, _, e := procWaveOutOpen.Call(
		uintptr(unsafe.Pointer(&handle)),
		uintptr(deviceNumber),
		uintptr(unsafe.Pointer(&format)),
		uintptr(0),
		uintptr(0),
		uintptr(CALLBACK_NULL))

	return WaveOutDevice{handle: handle, format: format,
	                     quit: make(chan bool), quitted: make(chan bool)},
	                     checkMmErrors("open", r, e)
}

func (device WaveOutDevice) writeSamples(samples []int16, bufferCh (chan []int16)) error {
	var waveData WAVEHDR
	waveData.lpData = uintptr(unsafe.Pointer(&samples[0]))
	waveData.dwBufferLength = uint32(len(samples) * 2)
	r, _, e := procWaveOutPrepareHeader.Call(
		device.handle,
		uintptr(unsafe.Pointer(&waveData)),
		uintptr(unsafe.Sizeof(waveData)))

	if err := checkMmErrors("preparing samples", r, e); err != nil {
		return err
	}

	r2, _, e := procWaveOutWrite.Call(
		device.handle,
		uintptr(unsafe.Pointer(&waveData)),
		uintptr(unsafe.Sizeof(waveData)))

	go func() {
		for (waveData.dwFlags & WHDR_DONE) == 0 {
			time.Sleep(time.Second / time.Duration(len(samples)))
		}
		r, _, e := procWaveOutUnprepareHeader.Call(
			device.handle,
			uintptr(unsafe.Pointer(&waveData)),
			uintptr(unsafe.Sizeof(waveData)))
		if err := checkMmErrors("unpreparing samples", r, e); err != nil {
			fmt.Println(err)
		}
		bufferCh <- samples
	}()

	return checkMmErrors("writing samples", r2, e)
}

func (device WaveOutDevice) Start(callback func([]int16)) error {
	if device.callback != nil {
		return errors.New("audio: start: already started")
	}
	device.callback = callback
	go func() {
	    const numBuffers = 8;
		bufferCh := make(chan []int16, numBuffers)
		bufferSamples := device.format.nSamplesPerSec / numBuffers / 16
		for i := 0; i < numBuffers; i++ {
			bufferCh <- make([]int16, bufferSamples)
		}
		for {
			select {
			case <-device.quit:
			    for i := 0; i < numBuffers; i++ {
			        <-bufferCh
		        }
			    device.quitted <- true
				return
			case b := <-bufferCh:
				callback(b)
				e := device.writeSamples(b, bufferCh)
				if e != nil {
					fmt.Println(e)
				}
			default:
				time.Sleep(time.Second / time.Duration(bufferSamples) / 4)
			}
		}
	}()
	return nil
}

func (device WaveOutDevice) WriteSamples(samples []int16) {
	var waveData WAVEHDR
	waveData.lpData = uintptr(unsafe.Pointer(&samples[0]))
	waveData.dwBufferLength = uint32(len(samples) * 2)
	procWaveOutPrepareHeader.Call(
		device.handle,
		uintptr(unsafe.Pointer(&waveData)),
		uintptr(unsafe.Sizeof(waveData)))
	procWaveOutWrite.Call(
		device.handle,
		uintptr(unsafe.Pointer(&waveData)),
		uintptr(unsafe.Sizeof(waveData)))
	time.Sleep(100 * time.Millisecond)

	for (waveData.dwFlags & WHDR_DONE) == 0 {
		time.Sleep(200 * time.Millisecond)
	}

	procWaveOutUnprepareHeader.Call(
		device.handle,
		uintptr(unsafe.Pointer(&waveData)),
		uintptr(unsafe.Sizeof(waveData)))
}

func (device WaveOutDevice) Reset() error {
	device.quit <- true
	<-device.quitted
	r, _, e := procWaveOutReset.Call(device.handle)
	return checkMmErrors("reset", r, e)
}

func (device WaveOutDevice) Close() error {
	r, _, e := procWaveOutClose.Call(device.handle)
	return checkMmErrors("close", r, e)
}
