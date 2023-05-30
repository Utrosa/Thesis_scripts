import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import logging
import tempfile
import gc
import shutil
import os

from matplotlib.gridspec import GridSpec
import parselmouth
from parselmouth import Sound
import tempfile
import textgrid

from scipy.signal import resample
from scipy.io import wavfile
from scipy.ndimage.filters import maximum_filter1d

from .detect_peaks import detect_peaks
from .utils import (
    ms_to_samples,
    samples_to_ms,
    apply_fades,
    save_samples_to_file,
    midi2freq,
    freq2midi,
    se_resample,
    band_pass_sharp,
    diff,
    to_contour,
    uds_contour_interval,
    diff_across,
    overlap_segments_score
)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

logging.getLogger('matplotlib.font_manager').disabled = True

mpl.style.use("bmh")

def recode_wav(file_path) -> object:
    with tempfile.NamedTemporaryFile() as temp_file:
        shutil.copyfile(file_path, temp_file.name)
        s = Sound(temp_file.name)
        s.save(file_path, "WAV")

def read_file_and_detect(
        filepath,
        extract_config,
        write=None,
        ramp_ms=150,
        start_cut=100,
        sample_rate=44100,

    ):
    recode_wav(filepath)

    fs, samples = wavfile.read(filepath)
    if samples.ndim == 2:
        samples = (samples[:, 0] + samples[:, 1]) / 2.0

    resampled, tt = se_resample(samples=samples, source_rate=fs, target_rate=sample_rate)
    cut_start = np.copy(resampled[int(ms_to_samples(start_cut, sample_rate)):])

    linearRamp = np.linspace(0, 1, round(1 + (ramp_ms / 1000.0 * sample_rate)))
    linearRampStart = linearRamp[:len(linearRamp) - 1]
    cut_start[:ms_to_samples(ramp_ms, sample_rate)] = cut_start[:ms_to_samples(ramp_ms, sample_rate)] * linearRamp

    filtered = band_pass_sharp(cut_start, sample_rate, min(extract_config["singing_bandpass_range"]), max(extract_config["singing_bandpass_range"]))
    max_normalized = filtered / max(abs(filtered))

    parselmouth_pitch_pre = Sound(max_normalized).to_pitch_ac(
        pitch_floor= midi2freq(min(extract_config["pitch_range_allowed"])),
        silence_threshold = extract_config["praat_silence_threshold"],
        voicing_threshold = 0.45,
        octave_cost = extract_config["praat_high_frequncy_favoring_octave_cost"],
        octave_jump_cost = extract_config["praat_octave_jump_cost"],
        voiced_unvoiced_cost = 0.14,
        pitch_ceiling = midi2freq(max(extract_config["pitch_range_allowed"]))
    ).to_array()[0]
    parselmouth_pitch_fs = Sound(max_normalized).sampling_frequency
    parselmouth_pitch = [x[0] for x in parselmouth_pitch_pre]

    resampled_f0 = resample(parselmouth_pitch, len(max_normalized))
    f0 = np.array([freq2midi(x) for x in resampled_f0])
    if write:
        wavfile.write(write, fs, max_normalized)

    good_tt = np.linspace(0, (len(max_normalized)) / (1.0 * sample_rate), len(max_normalized))

    praat_syllables= extract_syllables_via_praat(Sound(max_normalized),hz_from=min(extract_config["singing_bandpass_range_praat_syllable"]),hz_to=max(extract_config["singing_bandpass_range_praat_syllable"]))

    return max_normalized, f0, good_tt, praat_syllables

# Extract syllables using praat (based on code by Erika and Pol)
# erika's parameters were:hz_from = 2700, hz_to = 6000
def extract_syllables_via_praat(sound,hz_from = 80, hz_to = 6000):
    sound = parselmouth.praat.call(sound, "Filter (pass Hann band)", hz_from, hz_to, 100)
    with tempfile.TemporaryDirectory() as temp_dir:
        with open(os.path.dirname(os.path.abspath(__file__))+'/syllable_extract.praat', 'r') as file:
            script = ''.join(file.readlines())

            _, tg_obj = parselmouth.praat.run(sound, script)

            textgrid_path = temp_dir + "/temp.TextGrid"
            parselmouth.praat.call(tg_obj, "Save as text file", textgrid_path)
            tg = textgrid.TextGrid.fromFile(textgrid_path)

        onsets = []
        offsets = []
        for interval in tg[1]:
            if interval.mark == '':
                    onsets.append(interval.minTime*1000.0)
                    offsets.append(interval.maxTime*1000.0)
        return (onsets, offsets)


def calc_segments(audio, f0, tt, pp, sample_rate=44100):
    assert audio.ndim == 1
    assert len(audio) == len(tt)

    loudness = np.power(audio, 2)


    xs_orig = maximum_filter1d(loudness, size=ms_to_samples(pp["smoothing_env_window_ms"], fs=sample_rate)) # changed from 20 msec)
    xs_orig = np.power(xs_orig, pp["compresssion_power"])

    xs = np.copy(xs_orig)
    pos= (f0>(0.9+min(pp["pitch_range_allowed"]))) & (f0< (-0.9 +max(pp["pitch_range_allowed"]))) # find locations with not allowed pitches the 0.9 is to avoid including cieling or floor values
    xs[~pos]=0 # make all not detected pitches to 0 envelope
    #filtered = np.abs(band_pass_sharp(audio, sample_rate, 2700.0,6000.0))
    #filtered = filtered / max(filtered)
    #xs_filtered = 0*maximum_filter1d(filtered, size=ms_to_samples(5, fs=sample_rate))
    #xs= 0.5*xs+ 0.5*xs_filtered

    #xs=xs+0.1*xs_orig
    #xs[~pos]=-1 # 
    peaks = detect_peaks(xs, mpd=ms_to_samples(pp["peak_time_difference"], fs=sample_rate), mph=pp["minimum_peak_height"])

    # 1. find starts (this method needs a "fake" peak at the start of the peaks array)
    starts = []
    peaks = np.insert(peaks, 0, [0])
    for l in range(1, len(peaks)):
        idx = peaks[l]
        max_amp = xs[idx]
        max_amp_db = 20 * np.log10(max_amp)
        amp_thrsh_db = max_amp_db + pp["db_threshold"]
        amp_thrsh = np.power(10, amp_thrsh_db / 20)

        last = idx
        is_found = False
        for i in range(idx, peaks[l - 1], -1):
            if np.abs(xs[i]) > amp_thrsh:
                last = i

            ms_since_last = (1000 * abs(i - last) / sample_rate)
            if ms_since_last > pp["msec_silence"]:
                is_found = True
                break

        if last == idx:
            new_start = i
        else:
            new_start = last

        if is_found:
            starts.append(new_start)

    # 2. find ends (start scanning from closest peak after start)
    ends = []
    for l in range(len(starts)):
        idx = peaks[np.argmin(np.abs(peaks - starts[l]))]

        start_amp = xs[starts[l]]
        eps=1e-12
        start_amp_db = 20 * np.log10(start_amp +eps)
        start_thrsh_db = start_amp_db + pp["db_end_threshold_realtive_2note_start"]
        amp_thrsh1 = np.power(10, start_thrsh_db / 20)

        max_amp_db = 20 * np.log10(max_amp+eps)
        amp_thrsh_db = max_amp_db + pp["db_threshold"]
        amp_thrsh2 = np.power(10, amp_thrsh_db / 20)
        #importance =pp["max_vs_start_threshold_importance"]
        #amp_thrsh = (importance * amp_thrsh2) + ((1-importance) * amp_thrsh1)
        amp_thrsh= min(amp_thrsh2,amp_thrsh1) #  EXPERIEMNTAL
        #amp_thrsh = (0.9 * amp_thrsh2) + (0.1 * amp_thrsh1)

        last = idx
        for i in range(idx, len(audio)):
            if np.abs(xs[i]) > amp_thrsh:
                last = i

            if samples_to_ms(np.abs(i - last), fs=sample_rate) > pp["msec_silence"] and \
                    samples_to_ms(np.abs(i - starts[l]), fs=sample_rate):
                is_found = True
                break

        if last == idx:
            new_end = i
        else:
            new_end = last

        if is_found:
            ends.append(new_end)

    # Filter out overlaps by checking if any ends are overlapping.
    # Filter out crosses (not sure what the source of this bug is, but see example8.wav)
    pop_indices = []
    overlap_filter_index = 0
    for i in range(len(starts) - 1):
        curr_start = starts[i]
        curr_end = ends[i]
        next_end = ends[i+1]
        next_start = starts[i + 1]

        if curr_end > next_end:
            pop_indices.append(i)

        elif curr_start > curr_end:
            pop_indices.append(i)

    for pop_index in pop_indices:
        starts.pop(pop_index) # pop the associated start, too.
        ends.pop(pop_index)

    # fix a bug where the next start overlaps with the current end in this case trim the end to the earlier start
    for i in range(len(starts) - 1):
        curr_start = starts[i]
        curr_end = ends[i]
        next_end = ends[i + 1]
        next_start = starts[i + 1]

        if curr_end > next_start:
            starts[i+1]=next_start + samples_to_ms(pp["cut_pre"], fs=sample_rate)
            ends[i]=next_start - samples_to_ms(pp["cut_pre"], fs=sample_rate)

    peaks = peaks[1:]
    return xs, peaks, starts, ends


# tries to extend the segments based on the praat syllables
def extend_segments_based_on_praat(f0, extract_config, praat_syllables_tt,pitches):
    praat_syllables_starts, praat_syllables_ends = praat_syllables_tt

    for orig_index in range(len(pitches)):
        start=pitches[orig_index]["start_tt_extended"]
        end=pitches[orig_index]["end_tt_extended"]
        new_start=start
        new_end=end
        for praat_index, (start_praat, end_praat) in enumerate(zip(praat_syllables_starts,praat_syllables_ends)):
            # checks if segments overlap

            if overlap_segments_score(new_start, new_end,start_praat, end_praat)>0:
                # if the praatstart is earlier but not too earlier (extend ponly if within the  praat_extend_proximity_threshold_ms)
                if (start_praat<start) and abs(start_praat-start)<extract_config["praat_extend_proximity_threshold_ms"]:
                    if orig_index > 0:
                        if start_praat>pitches[orig_index-1]["end_tt"]: # don't extend if you ovelap with previous end
                            new_start=start_praat
                    else:
                        new_start=start_praat
                # the end can be extended as long as it's not overlapping to the next region
                if (end_praat>end) :
                    if orig_index<(len(pitches)-1):
                        if end_praat<pitches[orig_index+1]["start_tt"]: # don't extend if you ovelap with next start
                            new_end=end_praat
                    else:
                        new_end = end_praat

        pitches[orig_index]["start_tt_extended"]=new_start
        pitches[orig_index]["end_tt_extended"] = new_end

    return pitches


# tries to extend segmnts based on uniterupted pitch.
# tries to see if segments can be extended into the past/future as long as the pitch do not varry more than extend_pitch_threshold_semitones from the median
# this supose to combat situation where the voice detection decide to stop the syllable but there is good and stable pitch exraction, that menas that the sung tone is trimmed too shortly
def extend_segments_based_on_pitch( f0, fs, pp, pitches):

    for index in range(len(pitches)):
        start = int(fs * pitches[index]['start_tt_extended']/1000.0)
        end   = int(fs * pitches[index]['end_tt_extended'] / 1000.0)

        trimmed_f0 = f0[start:end]
        trimmed_f0 = trimmed_f0[~np.isnan(trimmed_f0)]
        very_trimmed_f0 = trimmed_f0[trimmed_f0 > 10]
        median_f0 = np.median(very_trimmed_f0)

        if index==0:
            old_end=0
        else:
            old_end=int(fs * pitches[index-1]['end_tt_extended'] / 1000.0)

        new_start=start
        for new_start in range(start, old_end, -1): # scan the past untill you find a place that deviates too much

            if (np.abs(f0[new_start]-median_f0)>pp["extend_pitch_threshold_semitones"]):
                break

        if index == (len(pitches)-1):
            old_start = len(f0)
        else:
            old_start = int(fs * pitches[index+1]['start_tt_extended'] / 1000.0)
            #old_start = starts[index + 1]

        new_end = end
        for new_end in range(end, old_start, +1):  # scan the future untill you find a place that deviates too much
            if (np.abs(f0[new_end] - median_f0) > pp["extend_pitch_threshold_semitones"]):
                break

        # update your pitches in the extended registry:
        pitches[index]["start_tt_extended"] = samples_to_ms(new_start, fs)
        pitches[index]["end_tt_extended"] = samples_to_ms(new_end, fs)

    return pitches

class PlotOptions():
    def __init__(
            self,
            display = False,
            save = False,
            path = "singing-analysis.png",
            format = None,
            dpi = 300
        ):
        """
        Plot options for singing-sing4me.
        """
        self.display = display
        self.save = save
        self.path = path
        self.format = format
        self.dpi = dpi


def get_pitches(
        audio,
        f0,
        fs,
        praat_syllables,
        target_pitches,
        plot_options,
        extract_config
    ):

     # making places with no f0 audio really really soft
    tt = np.arange(0, len(audio)) / fs
    xs, peaks, starts, ends = calc_segments(
        audio,
        f0,
        tt,
        pp=extract_config
    )


    onset_index = 0
    pitches = []

    status=[]
    for index, (start, end) in enumerate(zip(starts, ends)):
        status.append("<not defined>")


    for index, (start, end) in enumerate(zip(starts, ends)):
        status[index]="<not determined>"
        trim_start = int(start + ms_to_samples(extract_config["cut_pre"], fs=fs))
        trim_end = int(end - ms_to_samples(extract_config["cut_post"], fs=fs))
        trimmed_f0 = f0[trim_start:trim_end]
        trimmed_f0 = trimmed_f0[~np.isnan(trimmed_f0)]
        very_trimmed_f0 = trimmed_f0[trimmed_f0 > 10]

        start_tt = samples_to_ms(start, fs=fs)
        end_tt = samples_to_ms(end, fs=fs)

        extraction = {}
        percent_flactuation=100.0
        if target_pitches:
            try: # if the participant sings more tones than required.
                target = target_pitches[index]
                extraction['target_f0'] = target

            except IndexError:
                pass

        if len(trimmed_f0) > ms_to_samples(extract_config["minimal_segment_duration"], fs=fs):
            min_f0 = np.min(very_trimmed_f0)
            max_f0 = np.max(very_trimmed_f0)
            median_f0 = np.median(very_trimmed_f0)
            precentile_f0_5 = np.percentile(very_trimmed_f0,10) # NOTE experimental - I chnaged it to relax the criteria
            precentile_f0_95 = np.percentile(very_trimmed_f0,90) # NOTE experimental  - I chnaged it to relax the criteria
            divergence = np.abs(trimmed_f0 - median_f0)
            percent_flactuation = 100.0*np.sum(divergence > extract_config["allowed_pitch_flactuations_witin_one_tone"]) / len(divergence)
            median_flactuation= np.median(divergence)
            fluctuation_range_percentile=precentile_f0_95-precentile_f0_5 # difference between percentile tones

            if median_f0 < min(extract_config["pitch_range_allowed"]):
                status[index]="median f0 {:.2f} < min pitch_range_allowed".format(median_f0)
                continue

            if median_f0 > max(extract_config["pitch_range_allowed"]):
                status[index]="median f0 {:.2f} > max pitch_range_allowed".format(median_f0)
                continue

            if percent_flactuation> extract_config["percent_of_flcatuating_within_one_tone"]: # if tone flactuate too much (the percent of out of bounds flactuation is more than X% of the tone)
                status[index]="percent_flactuation {:.2f} > percent_of_flcatuating".format(percent_flactuation)
                continue
            if fluctuation_range_percentile> extract_config["allowed_pitch_flactuations_witin_one_tone"]:
                status[index]="fluctuation_range_percentile {:.2f} > allowed_pitch_flactuations".format(fluctuation_range_percentile)
                continue

            extraction['onset_index'] = onset_index
            extraction['start_tt'] = start_tt
            extraction['end_tt'] = end_tt
            extraction['start_tt_extended'] = start_tt # originaly the extended is identical
            extraction['end_tt_extended'] = end_tt # originaly the extended is identical
            extraction['start_trim'] = 1000.0 * (trim_start ) / fs
            extraction['end_trim'] =  1000.0 * (trim_end )/ fs
            extraction['median_f0'] = median_f0
            extraction['min_f0'] = min_f0
            extraction['max_f0'] = max_f0
            extraction['precentile_f0_5'] = precentile_f0_5
            extraction['precentile_f0_95'] = precentile_f0_95
            extraction['percent_flactuation'] = percent_flactuation
            extraction['median_flactuation'] = median_flactuation

            status[index]="OK"


            pitches.append(extraction)
            onset_index += 1
        else:
            status[index]="less than minimal_segment_duration"


    pitches = extend_segments_based_on_pitch(f0, fs, extract_config, pitches)

    pitches = extend_segments_based_on_praat(f0, extract_config, praat_syllables, pitches)

    if plot_options.display or plot_options.save:
        plot_pitches(
            audio=audio,
            f0=f0,
            fs=fs,
            pitches=pitches,
            praat_syllables=praat_syllables,
            targets=target_pitches,
            tt=tt,
            xs=xs,
            peaks=peaks,
            starts=starts,
            ends=ends,
            status=status,
            plot_options=plot_options,
            extract_config=extract_config,
        )

    return pitches


def extract_onsets(raw, min_duration__sec):
    # support function to extract timings from extracted raw audio
    # currently in testing state

    # get start/ end onsets
    start_onsets = [x["start_tt"] for x in raw]
    end_onsets = [x["end_tt"] for x in raw]

    # get ISI
    ISI_ms = np.diff(start_onsets)

    # note durations
    note_durations_ms = [element2 - element1 for (element2, element1) in zip(end_onsets, start_onsets)]

    # silent durations
    silence_durations_ms = [element2 - element1 for (element2, element1) in zip(ISI_ms, note_durations_ms)]

    # validate that they are the same
    ISI_ms_validate1 = [element2 + element1 for (element2, element1) in zip(note_durations_ms, silence_durations_ms)]
    ISI_ms_validate1 = [round(num, 1) for num in ISI_ms_validate1]
    ISI_ms_validate2 = [round(num, 1) for num in ISI_ms]
    assert ISI_ms_validate1 == ISI_ms_validate2

    # convert to sec
    note_durations_sec = [(i / 1000) for i in note_durations_ms]
    silence_durations_sec = [(i / 1000) for i in silence_durations_ms]
    ISI_sec = [element2 + element1 for (element2, element1) in
                         zip(note_durations_sec, silence_durations_sec)]

    # min possible allowed = min_duration__sec
    note_durations_sec_corrected = [min_duration__sec if i < min_duration__sec else i for i in
                                    note_durations_sec]
    silence_durations_sec_corrected = [min_duration__sec if i < min_duration__sec else i for i in
                                       silence_durations_sec]
    ISI_sec_corrected = [element2 + element1 for (element2, element1) in
                         zip(note_durations_sec_corrected, silence_durations_sec_corrected)]

    # add a final silence (this does not really matter as there isn't a note after this)
    silence_durations_sec_corrected.append(0.5)

    # get differences from raw ISIs to corrected ISIs
    diff_in_ISIs = []
    for i, j in zip(ISI_sec, ISI_sec_corrected):
        diff = abs(i - j)
        diff_in_ISIs.append(diff)

    # correct note and silence duration so the ISI stays the same after correction
    # this may compromise the real duration of the silence (first) or note duraiton (second)
    silence_durations_sec_matched = []
    note_durations_sec_matched = []
    diff_in_ISIs.append(0)  # append a 0 so it matches the len of the other lists
    for a, b, c in zip(silence_durations_sec_corrected, note_durations_sec_corrected, diff_in_ISIs):
        if a > min_duration__sec:
            silence = a - c
            note = b
        else:
            silence = a
            note = b - c
        silence_durations_sec_matched.append(silence)
        note_durations_sec_matched.append(note)

    return silence_durations_sec_matched, note_durations_sec_matched, ISI_sec


def plot_pitches_cogsci(
        audio,
        f0,
        fs,
        pitches,
        targets,
        tt,
        xs,
        peaks,
        starts,
        ends,
        plot_options,
        extract_config
    ):
    #audio[np.isnan(f0)]=0.000001*audio[np.isnan(f0)] # making places with no f0 audio really really soft
    xs, peaks, starts, ends = calc_segments(
        audio,
        f0,
        tt,
        pp=extract_config
    )
    if plot_options.display or plot_options.save:
        fig = plt.figure(figsize=(12, 7))
        gs = GridSpec(nrows=3, ncols=1)

        # audio_ax = fig.add_subplot(gs[0])
        # env_ax = fig.add_subplot(gs[1])
        # midi_ax = fig.add_subplot(gs[2])

        midi_ax = fig.add_subplot(gs[0])
        env_ax = fig.add_subplot(gs[1])
        audio_ax = fig.add_subplot(gs[2])

        midi_ax.grid(False)
        env_ax.grid(False)
        audio_ax.grid(False)


        ##########
        # Raw Audio Plot
        audio_ax.plot(tt, audio)
        audio_ax.axvspan(0, extract_config["silence_beginning_ms"]/1000.0, alpha=0.2, facecolor="r")
        if extract_config["silence_beginning_ms"] > 5:
            audio_ax.text(0, min(audio), "silenced", c="r", rotation=90)

        audio_ax.set_ylabel("Amplitude", fontname="Times", fontsize=12)
        audio_ax.set_xlabel("Time (s)", fontname="Times", fontsize=12)

        #########
        # Envelope Amplitude Plot
        env_ax.plot(tt,xs)
        env_ax.set_ylabel("Env. Amplitude", fontname="Times", fontsize=12)

        for peak in peaks:
             env_ax.plot(peak/fs, xs[peak], 'rx')

        for start, end in zip(starts, ends):
            env_ax.axvline(x=start/fs, c="g", linewidth=1)
            env_ax.axvline(x=end/fs, c="r", linewidth=1)
            env_ax.axvspan(start/fs, end/fs, alpha=0.1)

       
        for resp in pitches:
            env_ax.axvline(x=resp["start_trim"]/1000.0, c="k", linewidth=1)
            env_ax.axvline(x=resp["end_trim"]/1000.0, c="k", linewidth=1)
            env_ax.axvspan(resp["start_trim"]/1000.0, resp["end_trim"]/1000.0, alpha=0.051)

        #########
        # MIDI plot

        # midi_ax.plot(tt, f0)
        midi_ax.set_xlim([0,4])
        midi_ax.set_ylim(extract_config["pitch_range_allowed"])

        midi_ax.set_ylabel("MIDI", fontname="Times", fontsize=12)

        if plot_options.display:
            plt.show()

        if plot_options.save:
            fig.savefig(
                plot_options.path,
                format = plot_options.format,
                dpi = plot_options.dpi,
                transparent=True
            )
            print("Plot saved")

        if not plot_options.display:
            fig.clf()

        plt.close(fig)
        del fig
        gc.collect()


def plot_pitches(
        audio,
        f0,
        fs,
        pitches,
        praat_syllables,
        targets,
        tt,
        xs,
        peaks,
        starts,
        ends,
        status,
        plot_options,
        extract_config
    ):
    #audio[np.isnan(f0)]=0.000001*audio[np.isnan(f0)] # making places with no f0 audio really really soft
    xs, peaks, starts, ends = calc_segments(
        audio,
        f0,
        tt,
        pp=extract_config
    )
    if plot_options.display or plot_options.save:
        fig = plt.figure(figsize=(12, 7))
        gs = GridSpec(nrows=3, ncols=1)

        midi_ax = fig.add_subplot(gs[0])
        env_ax = fig.add_subplot(gs[1])
        audio_ax = fig.add_subplot(gs[2])
        # input_ax = fig.add_subplot(gs[3])


        ##########
        # MANU BUILDING
        # Synthesized Melody Plot
        # y_range = np.arange(0,71)

        # input_ax.plot(tt, f0)

        # input_ax.set_ylabel("MIDI", fontname="Times", fontsize=12)
        # input_ax.set_xlabel("Input pitches", fontname="Times", fontsize=12)

        # targets_invented = [{"onset_index": 0, "start_tt": 986.6439909297052, "end_tt": 1476.530612244898, "start_trim": 1016.6666666666666, "end_trim": 1426.5079365079366, "median_f0": 64.73153951152149, "min_f0": 64.01820521341912, "max_f0": 66.00966891663933, "precentile_f0_5": 64.23091308146273, "precentile_f0_95": 65.52699525200333, "percent_flactuation": 0.0, "median_flactuation": 0.19010747362720792}, {"onset_index": 1, "start_tt": 1726.3492063492063, "end_tt": 2347.4603174603176, "start_trim": 1756.3718820861677, "end_trim": 2297.437641723356, "median_f0": 61.441085049869294, "min_f0": 60.75770383209931, "max_f0": 62.45869700373345, "precentile_f0_5": 60.88454097252891, "precentile_f0_95": 61.84385311136434, "percent_flactuation": 0.0, "median_flactuation": 0.17825451824215577}, {"onset_index": 2, "start_tt": 2535.578231292517, "end_tt": 3037.1882086167802, "start_trim": 2565.6009070294785, "end_trim": 2987.165532879819, "median_f0": 62.37766340553658, "min_f0": 61.65028404041872, "max_f0": 64.00443277055022, "precentile_f0_5": 61.90951796472775, "precentile_f0_95": 63.18136018533443, "percent_flactuation": 0.0, "median_flactuation": 0.21131772674184646}]

        # if targets:
        # 	for response in targets_invented:
        # 		interval = [response['start_trim']/1000.0, response['end_trim']/1000.0]
        # 		input_ax.plot(interval, [response['median_f0'], response['median_f0']], 'k-', linewidth=1.5)
        # 		input_ax.plot(interval, [response['min_f0'], response['min_f0']], 'k--', linewidth=1.5)
        # 		input_ax.plot(interval, [response['max_f0'], response['max_f0']], 'k--', linewidth=1.5)

        # 		text_loc_x = response['start_trim'] / 1000.0
        # 		text_loc_y = response['precentile_f0_95']+ 0.8*(response['precentile_f0_95']-response['precentile_f0_5'])
        # 		input_ax.text(text_loc_x, text_loc_y, "{:2.1f}".format(response['median_f0']))

        # 		input_ax.axvline(x=response['start_trim']/1000.0, c='k', linewidth=1.5)
        # 		input_ax.axvline(x=response['end_trim']/1000.0, c='k', linewidth=1.5)
        # 		input_ax.axvspan(response['start_trim']/1000.0, response['end_trim']/1000.0, alpha=0.051)

        # for start, end in zip(starts, ends):
        # 	input_ax.axvline(x=start/fs, c="g", linewidth=1.5)
        # 	input_ax.axvline(x=end/fs, c="r", linewidth=1.5)
        # 	input_ax.axvspan(start/fs, end/fs, alpha=0.1)

        ##########
        # Raw Audio Plot
        audio_ax.plot(tt, audio)
        audio_ax.axvspan(0, extract_config["silence_beginning_ms"]/1000.0, alpha=0.2, facecolor="r")
        if extract_config["silence_beginning_ms"] > 5:
            audio_ax.text(0, min(audio), "silenced", c="r", rotation=90)

        audio_ax.set_ylabel("Amplitude", fontname="Times", fontsize=12)
        audio_ax.set_xlabel("Time (s)", fontname="Times", fontsize=12)


        for i in range(len(praat_syllables[0])):
            yc=0.9+0.1*i/len(praat_syllables[0])
            audio_ax.plot([  praat_syllables[0][i]/1000.0, (praat_syllables[1][i])/1000.0], [yc,yc] )

        for i in range(len(pitches)):
            start=pitches[i]["start_tt_extended"]/1000.0
            end = pitches[i]["end_tt_extended"]/1000.0
            audio_ax.axvline(x=start , c="c", linewidth=1.2)
            audio_ax.axvline(x=end , c="m", linewidth=1)
            audio_ax.axvspan(start , end , facecolor="y",alpha=0.2)
            yn=0.8+0.1*(i % 2)
            audio_ax.plot([start,end],[yn,yn],c="c")

        for start, end in zip(starts, ends):
            audio_ax.axvline(x=start / fs, c="g", linewidth=1.5)
            audio_ax.axvline(x=end / fs, c="r", linewidth=1.5)
            audio_ax.axvspan(start / fs, end / fs, alpha=0.1)

        #########
        # Envelope Amplitude Plot
        env_ax.plot(tt,xs)
        env_ax.set_ylabel("Env. Amplitude", fontname="Times", fontsize=12)

        for peak in peaks:
             env_ax.plot(peak/fs, xs[peak], 'rx')


        for resp in pitches:
            start = resp["start_tt_extended"] / 1000.0
            end = resp["end_tt_extended"] / 1000.0
            env_ax.axvline(x=start, c="c", linewidth=1.2)
            env_ax.axvline(x=end, c="m", linewidth=1)
            env_ax.axvspan(start, end, facecolor="y", alpha=0.2)

        for start, end in zip(starts, ends):
            env_ax.axvline(x=start / fs, c="g", linewidth=1.7)
            env_ax.axvline(x=end / fs, c="r", linewidth=1.7)
            env_ax.axvspan(start / fs, end / fs, alpha=0.1)

        for start, end, stat in zip(starts, ends,status):
            env_ax.text( (start/fs + end/fs)/2, 0, stat, c="r", fontsize='xx-small',rotation=90)


        for resp in pitches:
            env_ax.axvline(x=resp["start_trim"] / 1000.0, c="k", linewidth=1.5)
            env_ax.axvline(x=resp["end_trim"] / 1000.0, c="k", linewidth=1.5)
            env_ax.axvspan(resp["start_trim"] / 1000.0, resp["end_trim"] / 1000.0, alpha=0.051)

        #########
        # MIDI plot
        midi_ax.plot(tt,f0)
        midi_ax.set_ylim(extract_config["pitch_range_allowed"])
        for response in pitches:
            interval = [response['start_trim']/1000.0, response['end_trim']/1000.0]
            midi_ax.plot(interval, [response['median_f0'], response['median_f0']], 'k-', linewidth=1.5)
            midi_ax.plot(interval, [response['min_f0'], response['min_f0']], 'k--', linewidth=1.5)
            midi_ax.plot(interval, [response['max_f0'], response['max_f0']], 'k--', linewidth=1.5)

            text_loc_x = response['start_trim'] / 1000.0
            text_loc_y = response['precentile_f0_95']+ 0.8*(response['precentile_f0_95']-response['precentile_f0_5'])
            midi_ax.text(text_loc_x, text_loc_y, "{:2.1f}".format(response['median_f0']))

            midi_ax.axvline(x=response['start_trim']/1000.0, c='k', linewidth=1.5)
            midi_ax.axvline(x=response['end_trim']/1000.0, c='k', linewidth=1.5)
            midi_ax.axvspan(response['start_trim']/1000.0, response['end_trim']/1000.0, alpha=0.051)

        if targets:
            for i in range(len(targets)):
                target = targets[i]
                if len(pitches) == len(targets):
                    response = pitches[i]
                    interval = [response['start_trim'] / 1000.0, response['end_trim'] / 1000.0]
                    target_text_loc_y = target +  0.8*(response['precentile_f0_95']-response['precentile_f0_5'])
                else:
                    interval = [((max(tt)-min(tt))*i)/len(targets), ((max(tt)-min(tt))*(i+1)) / len(targets)]
                    target_text_loc_y = target + 1.5

                midi_ax.plot(interval, [target, target], "g", linewidth=1.5)
                target_text_loc_x = max(interval)
                midi_ax.text(target_text_loc_x, target_text_loc_y, "{:2.1f}".format(target), c="g")


        for k in range(len(pitches[:-1])):
            note_a = pitches[k]['median_f0']
            note_b = pitches[k+1]['median_f0']
            end_a = pitches[k]['end_trim']/1000.0
            start_b = pitches[k+1]['start_trim']/1000.0
            midi_ax.plot([end_a,start_b], [note_a,note_b], 'k--', linewidth=1.5)
            midi_ax.text(sum([end_a,start_b])/2.0, max([note_a,note_b])," {:2.1f}".format(note_b-note_a), fontsize=8)

        midi_ax.set_title("Extracted pitches", fontname="Times", fontsize=16)
        midi_ax.set_ylabel("MIDI", fontname="Times", fontsize=12)

        if plot_options.display:
            plt.show()

        if plot_options.save:
            fig.savefig(
                plot_options.path,
                format = plot_options.format,
                dpi = plot_options.dpi
            )
            print("Plot saved")

        if not plot_options.display:
            fig.clf()

        plt.close(fig)
        del fig
        gc.collect()


def analyze(
        filename,
        extract_config,
        target_pitches=None,
        plot_options=PlotOptions(),
        sample_rate=44100,
    ):
    max_normalized, resampled_f0, tt, praat_syllables = read_file_and_detect(
        filename,
        extract_config,
    )


    analysis = get_pitches(
        audio=max_normalized,
        f0=resampled_f0,
        fs=sample_rate,
        praat_syllables=praat_syllables,
        target_pitches=target_pitches,
        plot_options=plot_options,
        extract_config=extract_config
    )
    return analysis

def compute_error(raw_diffs):
    abs_differences = [abs(x) for x in raw_diffs]
    max_error = max(abs_differences, default=999)
    root_mean_squared = 999
    if abs_differences:
        mean_squared_differences = np.mean([x**2 for x in abs_differences])
        root_mean_squared = np.sqrt(mean_squared_differences)
    return max_error, root_mean_squared

def compute_stats(sung_pitches, target_pitches, sung_intervals, target_intervals):
    raw_interval_diffs = diff_across(target_intervals, sung_intervals)
    raw_pitch_diffs = diff_across(target_pitches, sung_pitches)


    target_contour = uds_contour_interval(target_intervals)
    sung_contour = uds_contour_interval(sung_intervals)

    direction_accuracy_points = 0
    for sung_contour_elem, target_contour_elem in zip(sung_contour, target_contour):
        if sung_contour_elem == target_contour_elem:
            direction_accuracy_points += 1

    direction_accuracy = (direction_accuracy_points / (len(target_pitches) - 1)) * 100

    pitch_max_error, pitch_root_mean_squared = compute_error(raw_pitch_diffs)
    interval_max_error, interval_root_mean_squared = compute_error(raw_interval_diffs)

    return {
        "num_sung_pitches": len(sung_pitches),
        "num_target_pitches": len(target_pitches),
        "raw_pitch_diffs": raw_pitch_diffs,
        "mean_pitch_diffs": float(np.mean(raw_pitch_diffs) if raw_pitch_diffs else 999),
        "max_abs_pitch_error": pitch_max_error,
        "root_mean_squared_pitch": float(pitch_root_mean_squared),
        "raw_interval_diffs": raw_interval_diffs,
        "mean_interval_diff": float(np.mean(raw_interval_diffs) if raw_interval_diffs else 999),
        "max_abs_interval_error": interval_max_error,
        "root_mean_squared_interval": float(interval_root_mean_squared),
        "direction_accuracy": direction_accuracy
    }


###################################################################################################
# Stimulus Generation
def beep(f, duration, fs, attack=5, decay=5, sample_rate=44100):
    """
    Generate a sine wave of frequency f (Hz) for duration (seconds)
    The ramp is symmetric.
    """
    raw_sine = np.sin(2 * np.pi * np.arange(fs * duration) * f / fs)
    return apply_fades(raw_sine, attack, sample_rate)

def generate_stimulus_samples(spec, fs):
    midis = spec['midis']
    durations = spec['durations']
    onsets = spec['onsets']

    sound_properties = {}

    assert len(midis) == len(durations) == len(onsets)

    f0s = [midi2freq(midi) for midi in midis]
    signals = [
        beep(f=f0, duration=duration / 1000, fs=fs)
        for f0, duration in zip(f0s, durations)
    ]

    stimulus_end_tt = max(onsets) + durations[-1] + 500
    stimulus_end_samples = ms_to_samples(stimulus_end_tt, fs=fs)

    stimulus = np.zeros(stimulus_end_samples)

    for onset, signal in zip(onsets, signals):
        onset_samples = ms_to_samples(onset, fs=fs)
        stimulus[onset_samples:(onset_samples + len(signal))] = signal

    return stimulus

def generate_sine_tones(midis, durations, onsets, sample_rate, output_file):
    assert len(midis) == len(durations)
    assert len(durations) == len(onsets)
    spec = {
        "midis": midis,
        "durations": [round(x * 1000) for x in durations],
        "onsets": [round(x * 1000) for x in onsets]
    }
    fs = round(sample_rate)
    samples = generate_stimulus_samples(spec, fs)
    save_samples_to_file(samples, output_file, fs)
