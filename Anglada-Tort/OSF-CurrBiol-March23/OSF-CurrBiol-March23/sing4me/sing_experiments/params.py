note_duration = 0.5  # ISI = note_duration + note_silence
note_silence = 0.3

# Mean for the distribution of the melody's roving reference tone.
roving_mean_oldest = dict(
    default=57.5,
    male=51.5,
    female=63.5,
)

roving_mean_old = dict(
    default=55.5,
    male=49.5,
    female=61.5,
)

roving_mean = dict(
    default=55,
    low=49,
    high=61,
)

# Expressed as semitones either side of the roving_mean.
# For example, if my roving_mean isƒtime 60, and my roving_width is 2.5,
# I will end up selecting pitches from the range [57.5, 62.5].
roving_width = dict(
    default=2.5,
)

# Maximum interval sizes for different representation schemes
max_interval_size = dict(
    reference_is_first_note=9.5,
    reference_is_previous_note=4.5,
)

# Maximum overall pitch range for the melody
max_melody_pitch_range = dict(
    default=10
)

# Duration of the silent gap between melodies in e.g. 2AFC tasks
inter_stimulus_interval = 1.25

transpose_within_trial_min_step_size = 0.0
transpose_within_trial_max_step_size = 2.5

# Specified in semitones
slider = dict(
    width=5,
    jitter=2,
    max_trial_bonus=0.05,
    bonus_threshold_semitones=1,
)

###############################################
# prams for specific paradigms
###############################################
# singing 2 intervals
singing_2intervals = dict(
    # singing task
    sing_duration=4,
    max_abs_interval_error_treshold=5.5,  # we used 5 in pilot 2-3
    max_melody_pitch_range=99,  # we used 9 in pilot 2-3
    num_int=2,
    reference_mode="pitch_mode",  # three options: first_note, previous_note, or "pitch_mode"
    max_pitch_height_seed=9.5,
    max_pitch_height=15,
    discrete=False,  # True means that everything is quantized to the 12-tone scale
    max_mean_interval_error=5.5,
    # singing sing4me
    sample_rate=44100,
    peak_time_difference=70,  # Calculated with: ms_to_samples(70, fs=STANDARD_FS),
    minimum_peak_height=0.05,
    db_threshold=-30, # used to be -22 (-30 means probably that the pitch extraction blobs are mainly used for the syllable calc_segments)
    db_end_threshold_realtive_2note_start=-15, # used to be -10 (changed )
    #max_vs_start_threshold_importance = 0.9, #I got rid of this parameters as I am calculating min as threshold
    msec_silence=30,  # was 90 (minimal silence between segments - increase to require more time between segments
    silence_beginning_ms= 50,
    extend_pitch_threshold_semitones=2.0,
    praat_extend_proximity_threshold_ms=150.0,  # (ms) max allowed extension of onset of the sung tone based on paraat - if the onset computed from paraat is deviating more than this threshold the onset would not be extended!
    #cut_pre=110,  # time ignored at start of seg (ms)
    #cut_post=70,  # time ignored at end of seg (ms)
    cut_pre=40,  # EXPERIMENTAL used to be 30 - time ignored at start of seg (ms)
    cut_post=50,  # time ignored at end of seg (ms)
#    minimal_segment_duration=35,  # minimal time in sec that segment should have (ms)
    #minimal_segment_duration=70,  # minimal time in sec that segment should have (ms)   
    minimal_segment_duration=40, # experimental  used to be 35
    pitch_range_allowed=[36, 75],
    # used to be [36,84] range of acceptable pitches - this is used in the detection algorithm to eliminate areas of no allowed pitch detected - Nori modified the highest pitch to 75 as  praat defualt i detecting up to 600 Hz
    singing_bandpass_range=[80, 6000],  # in Hz - we use this to bandpass filtering the audio
    #singing_bandpass_range_praat_syllable= [1200,7000], # in Hz - we use this to bandpass filtering the audio for syllable extraction
    singing_bandpass_range_praat_syllable= [40,8000], # in Hz - we use this to bandpass filtering the audio for syllable extraction ; this paramters used by Erika in her tpping experiments and found to be good.
    smoothing_env_window_ms=40, # EXEPERIMENTAL TODO: used to be 50 in msec. changed from 20 msec.
    compresssion_power=0.5, # EXEPERIMENTAL TODO : non linear compression - no compression = 1
    #allowed_pitch_flactuations_witin_one_tone=6,  # in semitones - max flactuations that is conisdered OK
    allowed_pitch_flactuations_witin_one_tone=8.0,  # used to be 6 
    #percent_of_flcatuating_within_one_tone=20.0,
    percent_of_flcatuating_within_one_tone=35.0, # 
    # exclude tones with more than that percent if they flactuate more than allowed_pitch_flactuations_witin_one_tone
    praat_octave_jump_cost=0.55,
    # praat default is 0.35 but we suspect that this is not good enough (too many octave doubling "jumps")
    praat_high_frequncy_favoring_octave_cost=0.03,
    # control octave_cost parms parrat dedault is 0.01  see https://www.fon.hum.uva.nl/praat/manual/Sound__To_Pitch__ac____.html
    praat_silence_threshold=0.03,  # praat defaule: silence_threshold = 0.03 this
    # progress bar
    end_melody=3.5,
    end_singing=3.5,
    end_recording=0.5,
    upload=0.5,

)


if singing_2intervals["reference_mode"] == "first_note":
    singing_2intervals["max_interval_size"] = max_interval_size["reference_is_first_note"]
elif singing_2intervals["reference_mode"] == "previous_note":
    singing_2intervals["max_interval_size"] = 6.5  # we used 4.5 in pilots 2-3
elif singing_2intervals["reference_mode"] == "pitch_mode":
    singing_2intervals["max_interval_size"] = 6.5  # we used 4.5 in pilots 2-3


# singing 4 intervals
singing_long_melodies = singing_2intervals.copy()
singing_long_melodies["num_int"] = 4
singing_long_melodies["sing_duration"] = 6


# singing 1 interval
singing_1interval = singing_2intervals.copy()
singing_1interval["max_interval_size"] = 13.5
singing_1interval["num_int"] = 1
singing_1interval["max_mean_interval_error"] = 5.5
singing_1interval["sing_duration"] = 3.5
singing_1interval["end_melody"] = 2.5
singing_1interval["end_singing"] = 3
singing_1interval["end_recording"] = 0.5
singing_1interval["upload"] = 0.5


# Nori-replication1: pleasatnes ratings
nori_replication1 = dict(
    num_int=1,
    note_duration=0.4,
    note_silence=0.1,
    reference_mode="previous_note",
    max_interval_size=8.5,
    max_melody_pitch_range=99,  # we used 9 in pilot 2-3
    discrete=False,  # True means that everything is quantized to the 12-tone scale
)
nori_replication1["melody_duration"] = (nori_replication1["note_duration"] + nori_replication1["note_silence"]) * \
                                       nori_replication1["num_int"]


probe_1interval = dict(
    num_int=1,
    # previous_note reference means that each note is expressed as an interval in semitones
    # from the previous note. So an ascending major triad would be written [4, 3],
    # because we start on the tonic, go up by 4 semitones for the third, and
    # go up by another 3 semitones for the fifth.
    #
    # first_note reference means that each note is expressed as an interval in semitones
    # from the first note in the melody. So an ascending major triad would be written [4, 7].
    reference_mode="first_note",
    max_melody_pitch_range=99,  # we used 9 in pilot 2-3
    # We randomise ('rove') the starting note of the melodies in continuous pitch space
    # so as to avoid context effects.
    # roving_mean is the mean starting note (in MIDI)
    max_interval_size=14.5,
    # roving_width is the width of the roving distribution
    # note_duration is the length of each note in the melody (seconds)
    # note_duration=params.rating_config["note_duration"],
    # note_silence is the length of the (optional) silent gap between each note in the melody
    # note_silence=params.rating_config["note_silence"],
    # post_prime_silence_duration defines how long we wait after the prime before playing the melody
    post_prime_silence_duration=1.0,
    # These transposition variables are not relevant for this rating task.
    # transpose_within_trial=True,
    # transpose_within_trial_min_step_size=0.25,
    # transpose_within_trial_max_step_size=2.5,
)
