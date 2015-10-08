package audio

const (
	MMSYSERR_BASE = 0
	WAVERR_BASE   = 32

	MM_WOM_OPEN  = 0x3BB /* waveform output */
	MM_WOM_CLOSE = 0x3BC
	MM_WOM_DONE  = 0x3BD

	MM_WIM_OPEN  = 0x3BE /* waveform input */
	MM_WIM_CLOSE = 0x3BF
	MM_WIM_DATA  = 0x3C0

	/* general error return values */
	MMSYSERR_NOERROR      = 0                    /* no error */
	MMSYSERR_ERROR        = (MMSYSERR_BASE + 1)  /* unspecified error */
	MMSYSERR_BADDEVICEID  = (MMSYSERR_BASE + 2)  /* device ID out of range */
	MMSYSERR_NOTENABLED   = (MMSYSERR_BASE + 3)  /* driver failed enable */
	MMSYSERR_ALLOCATED    = (MMSYSERR_BASE + 4)  /* device already allocated */
	MMSYSERR_INVALHANDLE  = (MMSYSERR_BASE + 5)  /* device handle is invalid */
	MMSYSERR_NODRIVER     = (MMSYSERR_BASE + 6)  /* no device driver present */
	MMSYSERR_NOMEM        = (MMSYSERR_BASE + 7)  /* memory allocation error */
	MMSYSERR_NOTSUPPORTED = (MMSYSERR_BASE + 8)  /* function isn't supported */
	MMSYSERR_BADERRNUM    = (MMSYSERR_BASE + 9)  /* error value out of range */
	MMSYSERR_INVALFLAG    = (MMSYSERR_BASE + 10) /* invalid flag passed */
	MMSYSERR_INVALPARAM   = (MMSYSERR_BASE + 11) /* invalid parameter passed */
	MMSYSERR_HANDLEBUSY   = (MMSYSERR_BASE + 12) /* handle being used */

	/* flags used with waveOutOpen(), waveInOpen(), midiInOpen(), and */
	/* midiOutOpen() to specify the type of the dwCallback parameter. */
	CALLBACK_TYPEMASK = 0x00070000      /* callback type mask */
	CALLBACK_NULL     = 0x00000000      /* no callback */
	CALLBACK_WINDOW   = 0x00010000      /* dwCallback is a HWND */
	CALLBACK_TASK     = 0x00020000      /* dwCallback is a HTASK */
	CALLBACK_FUNCTION = 0x00030000      /* dwCallback is a FARPROC */
	CALLBACK_THREAD   = (CALLBACK_TASK) /* thread ID replaces 16 bit task */
	CALLBACK_EVENT    = 0x00050000      /* dwCallback is an EVENT Handle */

	/* waveform audio error return values */
	WAVERR_BADFORMAT    = (WAVERR_BASE + 0) /* unsupported wave format */
	WAVERR_STILLPLAYING = (WAVERR_BASE + 1) /* still something playing */
	WAVERR_UNPREPARED   = (WAVERR_BASE + 2) /* header not prepared */
	WAVERR_SYNC         = (WAVERR_BASE + 3) /* device is synchronous */
	WAVERR_LASTERROR    = (WAVERR_BASE + 3) /* last error in range */

	/* wave callback messages */
	WOM_OPEN  = MM_WOM_OPEN
	WOM_CLOSE = MM_WOM_CLOSE
	WOM_DONE  = MM_WOM_DONE
	WIM_OPEN  = MM_WIM_OPEN
	WIM_CLOSE = MM_WIM_CLOSE
	WIM_DATA  = MM_WIM_DATA

	/* device ID for wave device mapper */
	WAVE_MAPPER = -1

	/* flags for dwFlags parameter in waveOutOpen() and waveInOpen() */
	WAVE_FORMAT_QUERY                        = 0x0001
	WAVE_ALLOWSYNC                           = 0x0002
	WAVE_MAPPED                              = 0x0004
	WAVE_FORMAT_DIRECT                       = 0x0008
	WAVE_FORMAT_DIRECT_QUERY                 = (WAVE_FORMAT_QUERY | WAVE_FORMAT_DIRECT)
	WAVE_MAPPED_DEFAULT_COMMUNICATION_DEVICE = 0x0010

	/* flags for dwFlags field of WAVEHDR */
	WHDR_DONE      = 0x00000001 /* done bit */
	WHDR_PREPARED  = 0x00000002 /* set if this header has been prepared */
	WHDR_BEGINLOOP = 0x00000004 /* loop start block */
	WHDR_ENDLOOP   = 0x00000008 /* loop end block */
	WHDR_INQUEUE   = 0x00000010 /* reserved for driver */

	/* flags for dwSupport field of WAVEOUTCAPS */
	WAVECAPS_PITCH          = 0x0001 /* supports pitch control */
	WAVECAPS_PLAYBACKRATE   = 0x0002 /* supports playback rate control */
	WAVECAPS_VOLUME         = 0x0004 /* supports volume control */
	WAVECAPS_LRVOLUME       = 0x0008 /* separate left-right volume control */
	WAVECAPS_SYNC           = 0x0010
	WAVECAPS_SAMPLEACCURATE = 0x0020

	/* defines for dwFormat field of WAVEINCAPS and WAVEOUTCAPS */
	WAVE_INVALIDFORMAT = 0x00000000 /* invalid format */
	WAVE_FORMAT_1M08   = 0x00000001 /* 11.025 kHz, Mono,   8-bit  */
	WAVE_FORMAT_1S08   = 0x00000002 /* 11.025 kHz, Stereo, 8-bit  */
	WAVE_FORMAT_1M16   = 0x00000004 /* 11.025 kHz, Mono,   16-bit */
	WAVE_FORMAT_1S16   = 0x00000008 /* 11.025 kHz, Stereo, 16-bit */
	WAVE_FORMAT_2M08   = 0x00000010 /* 22.05  kHz, Mono,   8-bit  */
	WAVE_FORMAT_2S08   = 0x00000020 /* 22.05  kHz, Stereo, 8-bit  */
	WAVE_FORMAT_2M16   = 0x00000040 /* 22.05  kHz, Mono,   16-bit */
	WAVE_FORMAT_2S16   = 0x00000080 /* 22.05  kHz, Stereo, 16-bit */
	WAVE_FORMAT_4M08   = 0x00000100 /* 44.1   kHz, Mono,   8-bit  */
	WAVE_FORMAT_4S08   = 0x00000200 /* 44.1   kHz, Stereo, 8-bit  */
	WAVE_FORMAT_4M16   = 0x00000400 /* 44.1   kHz, Mono,   16-bit */
	WAVE_FORMAT_4S16   = 0x00000800 /* 44.1   kHz, Stereo, 16-bit */

	WAVE_FORMAT_44M08 = 0x00000100 /* 44.1   kHz, Mono,   8-bit  */
	WAVE_FORMAT_44S08 = 0x00000200 /* 44.1   kHz, Stereo, 8-bit  */
	WAVE_FORMAT_44M16 = 0x00000400 /* 44.1   kHz, Mono,   16-bit */
	WAVE_FORMAT_44S16 = 0x00000800 /* 44.1   kHz, Stereo, 16-bit */
	WAVE_FORMAT_48M08 = 0x00001000 /* 48     kHz, Mono,   8-bit  */
	WAVE_FORMAT_48S08 = 0x00002000 /* 48     kHz, Stereo, 8-bit  */
	WAVE_FORMAT_48M16 = 0x00004000 /* 48     kHz, Mono,   16-bit */
	WAVE_FORMAT_48S16 = 0x00008000 /* 48     kHz, Stereo, 16-bit */
	WAVE_FORMAT_96M08 = 0x00010000 /* 96     kHz, Mono,   8-bit  */
	WAVE_FORMAT_96S08 = 0x00020000 /* 96     kHz, Stereo, 8-bit  */
	WAVE_FORMAT_96M16 = 0x00040000 /* 96     kHz, Mono,   16-bit */
	WAVE_FORMAT_96S16 = 0x00080000 /* 96     kHz, Stereo, 16-bit */

	WAVE_FORMAT_PCM = 1
)
