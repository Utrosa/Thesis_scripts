s = selected("Sound")
s$ = selected$("Sound")
dur = Get total duration
print selected$("Sound")

silence_threshold = - 25
minimum_pause_duration = 0.12
minimum_dip_between_peaks = 10


if silence_threshold >= 0
	exitScript: "“Silence threshold” should be a negative number."
endif



int = Get intensity (dB)

if int <> undefined
	result = Copy: s$ + "-marksyllables"

	s = selected("Sound")
	dur1 = Get total duration
	ch = Get number of channels

	stt = Get start time
	if stt <> 0
		Scale times to: 0, dur1
	endif

	if dur1 > 5
		tmp1 = Extract part: 0, 5, "rectangular", 1, "no"
	endif

	if ch > 1
		tmp2 = Extract one channel: 1
	endif

	last_sel = selected("Sound")
	trimmed = nocheck nowarn Trim silences: 0.08, "no", 100, 0, -35, 0.1, 0.05, "no", "trimmed"
	if trimmed = last_sel or trimmed = undefined
		tmp3 = Copy: "tmp3"
		dur2 = Get total duration
	else
		tmp3 = trimmed
		dur2 = Get total duration
		stt = Get start time
		if stt <> 0
			Scale times to: 0, dur2
		endif
	endif

	tmp4 = Filter (stop Hann band): 20, 0, 20

	tmp5 = Extract part: 0, dur2, "Gaussian2", 1, "no"
	int = Get intensity (dB)

	if dur1 > 5
		removeObject: tmp1
	endif
	if ch > 1
		removeObject: tmp2
	endif
	removeObject: tmp3, tmp4, tmp5

	selectObject: s

	if int > 10
		Subtract mean
		tmp6 = Filter (stop Hann band): 0, 60, 20
		if dur1 > 0.5
			Fade in: 0, 0, 0.005, "yes"
			Fade out: 0, dur1, -0.005, "yes"
		endif
		selectObject: s
		Formula: "object[tmp6]"
		removeObject: tmp6
	endif

		#runScript: "declip.praat"
	extr = Get absolute extremum: 0, 0, "None"

		if extr = undefined
			extr = 0
		endif

		if number(fixed$(extr, 2)) > 0.99
			Scale peak: 0.99
		endif

	# endproc
	s = selected("Sound")
	s$ = selected$("Sound")
	dur = Get total duration
	intensity = noprogress To Intensity: 50, 0, "no"

	minint = Get minimum: 0, 0, "Parabolic"
	maxint = Get maximum: 0, 0, "Parabolic"
	max99int = Get quantile: 0, 0, 0.99

	threshold = max(max99int + silence_threshold, minint)
	threshold2 = maxint - max99int
	threshold3 = silence_threshold - threshold2

	tg = noprogress To TextGrid (silences): threshold3, minimum_pause_duration, 0.1, "pause", ""
	Rename: s$ + "-marksyllables"
	Set tier name: 1, "syllables"

	selectObject: intensity
	intensitytier = noprogress To IntensityTier (peaks)
	npoints = Get number of points
	peakcount = 0
	int[1] = 0
	for i to npoints
		db = Get value at index: i
		if db > threshold
			peakcount += 1
			int[peakcount] = db
			timepeaks[peakcount] = Get time from index: i
		endif
	endfor
	timepeaks[peakcount + 1] = dur

	selectObject: intensity
	validpeakcount = 0
	currenttime = timepeaks[1]
	currentint = int[1]
	for p to peakcount
		dip = Get minimum: currenttime, timepeaks[p + 1], "None"
		diffint = abs(currentint - dip)
		if diffint > minimum_dip_between_peaks
			validpeakcount += 1
			validtime[validpeakcount] = timepeaks[p]
		endif
		currenttime = timepeaks[p + 1]
		currentint = Get value at time: timepeaks[p + 1], "Cubic"
	endfor

	selectObject: result
	pitch = noprogress To Pitch (ac): 0.02, 30, 4, "no", 0.03, 0.25, 0.01, 0.35, 0.25, 450
	for i to validpeakcount
		pvalue[i] = Get value at time: validtime[i], "Hertz", "Linear"
	endfor

	selectObject: tg
	voicedcount = 0
	for i to validpeakcount
		whichinterval = Get interval at time: 1, validtime[i]
		whichlabel$ = Get label of interval: 1, whichinterval
		if pvalue[i] <> undefined
			if whichlabel$ <> "pause"
				voicedcount += 1
				voicedpeak[voicedcount] = validtime[i]
			endif
		endif
	endfor

	Insert point tier: 1, "peaks"
	for i to voicedcount
		Insert point: 1, voicedpeak[i], string$(i)
	endfor

	selectObject: intensity
	for i from 2 to voicedcount
		mintime[i] =  Get time of minimum: voicedpeak[i - 1], voicedpeak[i], "None"
	endfor

	selectObject: tg
	for i from 2 to voicedcount
		whichinterval = Get interval at time: 2, mintime[i]
		whichlabel$ = Get label of interval: 2, whichinterval
		if whichlabel$ <> "pause"
			Insert boundary: 2, mintime[i]
		endif
	endfor


	selectObject: result, tg
	removeObject: intensity, intensitytier, pitch
else
	result = Copy: s$ + "-marksyllables"

	noprogress Create TextGrid: 0, dur, "peaks syllables", "peaks"
	Set interval text: 2, 1, "pause"
	Rename: s$ + "-marksyllables"
	plusObject: result

endif

