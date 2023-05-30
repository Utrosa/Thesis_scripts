import random
import math
import numpy as np


# supporting functions for singing experiments
def sample_interval_sequence(
        n_int,
        max_interval_size,
        max_melody_pitch_range,
        discrete,
        reference_mode,
):
    """
    Generates a random interval sequence, defining the pitch content of a melody.
    We use this for example for generating random seeds for iterated singing chains.

    Parameters
    ----------
    n_int:
        Number of intervals in the melody.

    max_interval_size:
        Maximum size interval allowed in the melody; see params.py.

    max_melody_pitch_range:
        Maximum pitch range allowed in the melody; see params.py

    discrete:
        If ``True``, the intervals are constrained to take integer values.

    reference_mode:
        Can be ``"first_note"`` or ``"previous_note"``.
        When the ``reference_mode`` is ``"first_note"``, this means that each interval is expressed relative
        to the first note in the melody. For example, a major triad would be expressed as ``[4, 7]``.
        When the ``reference_mode`` is ``"previous_note"``, this means that each interval is expressed relative
        to the previous note in the melody. For example, a major triad would be expressed as ``[4, 3]``.

    Returns
    -------

    A list of intervals, each expressed as numbers.
    """
    n_try = 1e5
    counter = 0
    while counter < n_try:
        candidate = [
            sample_interval(
                max_interval_size,
                discrete,
            )
            for _ in range(n_int)
        ]
        if is_valid_interval_sequence(candidate, n_int, max_interval_size, max_melody_pitch_range, reference_mode):
            return candidate
        counter += 1
    raise RuntimeError(
        "Failed to generate valid interval sequence with the following constraints"
    )


def sample_interval(max_interval_size, discrete):
    assert max_interval_size >= 0
    if discrete:
        return random.randint(- max_interval_size, max_interval_size)
    else:
        return random.uniform(- max_interval_size, max_interval_size)


def sample_reference_pitch(roving_mean, roving_width):
    return random.uniform(roving_mean - roving_width, roving_mean + roving_width)


def sample_num_pitches(min_num, max_num):
    return int(random.uniform(min_num, max_num))


def sample_absolute_pitches(reference_pitch, max_interval2reference, num_pitches):
    target_pitches = []
    for i in range(num_pitches):
        target_pitch = sample_absolute_pitch(
            reference_pitch=reference_pitch,
            max_interval2reference=max_interval2reference
        )
        target_pitches.append(target_pitch)
    return(target_pitches)


def sample_absolute_pitch(reference_pitch, max_interval2reference):
    interval2reference = random.uniform(-max_interval2reference, max_interval2reference)
    return reference_pitch + interval2reference


def is_valid_interval_sequence(
        intervals,
        n_int,
        max_interval_size,
        max_melody_pitch_range,
        reference_mode,
    ):
    """
    Checks whether an interval sequence is valid according to the constraints
    expressed in the options.
    This can be used, for example, for making sure that sequences stay within bounds in
    serial reproduction experiments.

    Parameters
    ----------

    intervals:
        A list of intervals, with each interval expressed as a number.

    max_interval_size:
        Maximum size interval allowed in the melody; see params.py.

    max_melody_pitch_range:
        Maximum pitch range allowed in the melody; see params.py

    discrete:
        If ``True``, the intervals are constrained to take integer values.

    reference_mode:
        Can be ``"first_note"`` or ``"previous_note"``.
        When the ``reference_mode`` is ``"first_note"``, this means that each interval is expressed relative
        to the first note in the melody. For example, a major triad would be expressed as ``[4, 7]``.
        When the ``reference_mode`` is ``"previous_note"``, this means that each interval is expressed relative
        to the previous note in the melody. For example, a major triad would be expressed as ``[4, 3]``.
    """
    try:
        assert len(intervals) == n_int
        for interval in intervals:
            assert abs(interval) <= max_interval_size
        assert get_melody_pitch_range(intervals, reference_mode) <= max_melody_pitch_range
        return True
    except AssertionError:
        return False


def is_valid_pitch_range(
        reference_pitch,
        pitches,
        max_pitch_range
    ):
    """
    Checks whether a pitch sequence is valid according to the maximum pitch range.

    Parameters
    ----------

    reference_pitch:
        The reference pitch used to sample the pitches,expressed as a midi number.

    pitches:
        A list of pitches, with each pitch expressed as a midi number.

    max_pitch_range:
        Maximum pitch height allowed in the melody; see params.py.

    """
    if isinstance(reference_pitch, list):
        reference_pitch = reference_pitch.pop(0)
    try:
        for pitch in pitches:
            assert abs(reference_pitch - pitch) <= max_pitch_range
        return True
    except AssertionError:
        return False


def get_melody_pitch_range(intervals, reference_mode):
    """
    Gets the range of pitches used in an interval sequence.
    """
    example_pitches = convert_interval_sequence_to_absolute_pitches(
        intervals,
        reference_pitch=60,
        reference_mode=reference_mode
    )
    return max(example_pitches) - min(example_pitches)


def convert_interval_sequence_to_absolute_pitches(intervals, reference_pitch, reference_mode):
    """
    Takes an interval sequence and converts it to a set of absolute pitches,
    i.e. MIDI note numbers where 60 corresponds to middle C (C4) and integers correspond to semitones.

    Parameters
    ----------

    intervals:
        A list of intervals, with each interval expressed as a number.

    reference_pitch:
        The reference pitch to use, expressed as a MIDI note number. The interpretation of this pitch
        is determined by the latter ``reference_mode`` argument.

    reference_mode:
        Can be ``"first_note"`` or ``"previous_note"``.
        When the ``reference_mode`` is ``"first_note"``, this means that each interval is expressed relative
        to the first note in the melody. For example, a major triad would be expressed as ``[4, 7]``.
        When the ``reference_mode`` is ``"previous_note"``, this means that each interval is expressed relative
        to the previous note in the melody. For example, a major triad would be expressed as ``[4, 3]``.

    Returns
    -------

    A list of absolute pitches, expressed as MIDI note numbers.
    """
    if reference_mode == "first_note":
        first_note = reference_pitch
        return [first_note] + [interval + reference_pitch for interval in intervals]
    elif reference_mode == "previous_note":
        pitches = [reference_pitch]
        for interval in intervals:
            last_pitch = pitches[-1]
            new_pitch = last_pitch + interval
            pitches.append(new_pitch)
        return pitches
    else:
        raise ValueError(f"Unrecognized reference_mode: {reference_mode}.")


def convert_intervals2reference_to_absolute_pitches(intervals2refernece, reference_pitch):
    """
    Takes a sequence of intervals2reference and generates absolute pitches given a new reference pitch

    Parameters
    ----------

    intervals2refernece:
        A list of intervals relative to the refernece pitch of the melody
    reference_pitch:
        The reference pitch to generate the original melodies, expressed as a MIDI note number.

    Returns
    -------

    A list of absolute pitch heights
    """
    if isinstance(reference_pitch, list):
        reference_pitch = reference_pitch.pop(0)
    pitches = []
    for interval2reference in intervals2refernece:
        pitch = reference_pitch + interval2reference
        pitches.append(pitch)
    return pitches


def convert_absolute_pitches_to_interval_sequence(pitches, reference_mode):
    """
    Takes a sequence of absolute pitches and converts it to an interval

    Parameters
    ----------

    pitches:
        A list of absolute pitches, expressed as MIDI note numbers.

    reference_mode:
        Can be ``"first_note"`` or ``"previous_note"``.
        When the ``reference_mode`` is ``"first_note"``, this means that each interval is expressed relative
        to the first note in the melody. For example, a major triad would be expressed as ``[4, 7]``.
        When the ``reference_mode`` is ``"previous_note"``, this means that each interval is expressed relative
        to the previous note in the melody. For example, a major triad would be expressed as ``[4, 3]``.

    Returns
    -------

    A list of absolute pitches, expressed as MIDI note numbers.
    """
    intervals = []
    for i in range(1, len(pitches)):
        pitch = pitches[i]
        if reference_mode == "first_note":
            reference_pitch = pitches[0]
        elif reference_mode == "previous_note":
            reference_pitch = pitches[i - 1]
        else:
            raise ValueError(f"Invalid reference_mode: '{reference_mode}'")
        intervals.append(pitch - reference_pitch)
    return intervals


def convert_absolute_pitches_to_intervals2reference(pitches, reference_pitch):
    """
    Takes a sequence of absolute pitches and converts it to a list of intervals compared 
    against the refernece pitch

    Parameters
    ----------

    pitches:
        A list of absolute pitches, expressed as MIDI note numbers.
    reference_pitch:
        The reference pitch to generate the original melodies, expressed as a MIDI note number.

    Returns
    -------

    A list of absolute pitch heights
    """
    if isinstance(reference_pitch, list):
        reference_pitch = reference_pitch.pop(0)
    intervals2reference = []
    for pitch in pitches:
        interval2reference = pitch - reference_pitch
        intervals2reference.append(interval2reference)
    return intervals2reference


def as_native_type(x):
    if type(x).__module__ == np.__name__:
        return x.item()
    return x


def diff(x):
    return [j - i for i, j in zip(x[:-1], x[1:])]


def failing_criteria(
        sung_intervals,
        sung_pitches,
        reference_pitch,
        num_target_intervals,
        max_interval_size,
        max_melody_pitch_range,
        reference_mode,
        stats,
        max_abs_interval_error_treshold,
        max_pitch_range
):
    """ Perform basic failing criteria for singing production tasks

    Parameters
    ----------

    sung_intervals : list
        A list of the sung intervals, with each interval expressed as a number.
    num_target_intervals : float
        The number of target intervals
    max_interval_size:
        Maximum size interval allowed in the melody; see params.py.
    max_melody_pitch_range:
        Maximum pitch range allowed in the melody; see params.py
    reference_mode:
        Can be ``"first_note"`` or ``"previous_note"``.
        When the ``reference_mode`` is ``"first_note"``, this means that each interval is expressed relative
        to the first note in the melody. For example, a major triad would be expressed as ``[4, 7]``.
        When the ``reference_mode`` is ``"previous_note"``, this means that each interval is expressed relative
        to the previous note in the melody. For example, a major triad would be expressed as ``[4, 3]``.
    stats: dict
        Output of compute_stats function in singing_extract
    max_abs_interval_error_treshold : float
        maximum value allowed for max absolute interval error

    Returns
    --------

    is_failed : dict
        A dictionary with the output of the failing criteria
    """
    # check number of notes
    num_sung_pitches = stats["num_sung_pitches"]
    num_target_pitches = stats["num_target_pitches"]
    correct_num_notes = num_sung_pitches == num_target_pitches

    # check max interval error
    max_interval_error_ok = stats["max_abs_interval_error"] < max_abs_interval_error_treshold

    # check interval sequence
    if reference_mode == "pitch_mode":
        valid_sequence = is_valid_pitch_range(
            reference_pitch,
            sung_pitches,
            max_pitch_range
        )
        message_valid_sequence = "A pitch height is too large"
    else:
        valid_sequence = is_valid_interval_sequence(
            sung_intervals,
            num_target_intervals,
            max_interval_size,
            max_melody_pitch_range,
            reference_mode
            )
        message_valid_sequence = "One interval is too large"


    failed_options = [
        correct_num_notes,
        max_interval_error_ok,
        valid_sequence
    ]
    reasons = [
        (f"Wrong number of sung notes: {num_sung_pitches}  sung out of {num_target_pitches} notes in melody"),
        "Max absolute interval error is too large",
        message_valid_sequence
    ]
    if False in failed_options:
        failed = True
        index = failed_options.index(False)
        reason = reasons[index]
    else:
        failed = False
        reason = "All good"
    is_failed = {
        "failed": failed,
        "reason": reason
    }
    return is_failed


def failing_criteria_unconstrained_notes(stats, max_num_diff_pitches):
    """ Perform failing criteria for singing paradigms with unconstrained notes (this is by definition more relaxed)
    """
    # check number of notes
    num_sung_pitches = stats["num_sung_pitches"]
    num_target_pitches = stats["num_target_pitches"]

    diff_sung_target_notes = abs(num_sung_pitches - num_target_pitches)
    correct_num_notes = num_sung_pitches == num_target_pitches

    if diff_sung_target_notes > max_num_diff_pitches:
        failed = True
        reason = "Too many changes in number of notes"
    elif num_sung_pitches < 3:
        failed = True
        reason = "Too few detected notes"
    elif num_sung_pitches > 15:
        failed = True
        reason = "Too many detected notes"
    else:
        failed = False
        reason = "All good"

    is_failed = {
        "failed": failed,
        "reason": reason,
    }
    return is_failed


def feedback_generator(stats, max_abs_interval_error_treshold):
    """ Generate feedback based on performance
    stats: dict
        Output of compute_stats function in singing_extract
    max_abs_interval_error_treshold : float
        maximum value allowed for max absolute interval error

    Returns
    --------

    feedback : dict
        A dictionary with the feedback to return to participants based on their peroformance
    """
    num_sung_pitches = stats["num_sung_pitches"]
    num_target_pitches = stats["num_target_pitches"]

    too_many_sung_notes = num_target_pitches < num_sung_pitches
    too_few_sung_notes = num_target_pitches > num_sung_pitches

    failed_options = [too_many_sung_notes, too_few_sung_notes,]
    reasons = [
    f"You have sung too many notes: the melody had {num_target_pitches} notes and you sung {num_sung_pitches} notes.",
    f"You have not sung enough notes: the melody had {num_target_pitches} notes and you sung {num_sung_pitches} notes.",
    ]

    if True in failed_options:
        status = "bad"
        index = failed_options.index(True)
        message = reasons[index]
    elif  stats["max_abs_interval_error"] < max_abs_interval_error_treshold:
        status = "perfect"
        message = "Your singing accuracy is Excellent."
    elif stats["direction_accuracy"] > 75:
        status = "ok"
        message = "Your singing accuracy is ok, but you can do better."
    else:
        status = "bad"
        message = "Your singing accuracy is not good. To continue with the experiment, you will need to improve your performance by singing each note in the melody accurately."

    feedback = {
        "status": status,
        "message": message
    }
    return feedback


def get_duration(seq, default_duration=None, default_silence=None):
    """
    Gets the duration of a JSSynth sequence.
    """
    sec = 0.0
    for elt in seq:
        duration = elt["duration"]
        if duration == "default":
            if default_duration is None:
                raise RuntimeError("Need to specify default_duration")
            duration = default_duration

        silence = elt["silence"]
        if silence == "default":
            if default_silence is None:
                raise RuntimeError("Need to specify default_duration")
            silence = default_silence

        sec += duration + silence
    return sec


def midi2freq(midi_number):
    return (440 / 32) * (2**((midi_number - 9) / 12))


def freq2midi(f0):
    if f0 <= 0:
        return 0
    else:
        exp = result=12 * (np.log2(f0) - np.log2(440)) + 69
        if exp < 0:
            return 0 # get's weird log2 values otherwise...
        else:
            return exp
