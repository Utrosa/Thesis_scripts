import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from IPython.display import clear_output
from rich import print
from rich.progress import track


def build_diatonic_PDF(params, custom_intervals=None, custom_weights=None):

    """

    Build a PDF (evaluated at x=cents) representing diatonic intervals

    Args:

        params (string):
            Accepts 'major','minor', 'random', or 'custom'. Major and minor
            settings use the K-K profiles (ratings), diatonic intervals only.
            Random chooses seven intervals and seven weights at random, and
            overlaps are common. Custom uses provided intervals and weights.

        custom_intervals (numpy vector of floats, 0-12):
            Locations of the "diatonic" intervals/notes in semitones, to
            nearest cent. Values must be between 0 and 11.99. Ignored if params
            is not "custom".

        custom_weights (numpy vec of floats):
            For each interval, what is the relative expectation or prevalence?
            Must have the same number of elements as 'custom_intervals'. Ignored
            if params is not "custom".

    Returns:

        PDF (numpy array of floats):
            Summed evaluations of gaussians centered at 'means', with
            evaluation in steps of 1 cent, for 1200 cents. Vector is normalized
            across the entire length (sums to 1).

        intervals (numpy vector of floats):
            Locations of the intervals/notes in semitones

        weights (numpy vector of floats):
            For each interval, what is the relative expectation or prevalence?
            For example, uses Krumhansl-Kessler maj/minor ratings, or generates
            random values. Could use any other set of weights (e.g., Temperley
            profile). Note, only the values at positions 'intervals' are chosen.

    """

    # preconditions
    assert params in ['major','minor','random','custom'], "Invalid parameter"

    # choose the settings for intervals/weights
    if params == 'major':
        weights = [6.35,2.23,3.48,2.33,4.38,4.09,2.52,5.19,2.39,3.66,2.29,2.88]
        intervals = np.array([0,2,4,5,7,9,11])
        weights = np.array([weights[i] for i in intervals]) # diatonic intervals
    elif params == 'minor':
        weights = [6.33,2.68,3.52,5.38,2.6,3.53,2.54,4.75,3.98,2.69,3.34,3.17]
        intervals = np.array([0,2,3,5,7,8,10])
        weights = np.array([weights[i] for i in intervals]) # diatonic intervals
    elif params == 'random':
        weights = np.random.random(7)
        intervals = np.random.random(size=7) * 12
    elif params == 'custom':
        if (custom_weights is None) | (custom_intervals is None):
            raise ValueError('Must provide custom intervals and weights')
        if len(custom_weights) != len(custom_intervals):
            raise ValueError('Number of custom intervals and weights must match')
        weights = custom_weights
        intervals = custom_intervals
    else:
        raise ValueError('Must choose "major", "minor", "random", or "custom"')

    # use three octaves, to give continuity before and after octave of interest
    step_size = .01 # cents, given that integers are semitones in MIDI
    x_points = np.arange(0, 12*3, step_size)
    idxs = np.concatenate([intervals + 12*i for i in np.arange(0,3)])

    # evaluate the gaussians
    # ["magic number" alert]
    # sd is estimated from pitch sd within middle section of all notes in study
    # (i.e., mean of 40,032 note SDs = 26.15 cents)
    sd = .2615
    mod_by = len(intervals)
    PDF = np.zeros_like(x_points)
    for i,m in enumerate(idxs):
        PDF = PDF + stats.norm.pdf(x_points,
                                   loc = m,
                                   scale = sd) * weights[np.mod(i,mod_by)]

    # get rid of the buffer octaves
    PDF = PDF[1200:2400]

    # normalize points to sum to 1
    PDF /= PDF.sum()

    return PDF, intervals, weights


def MIDI_to_transpose_array(pitch_vec):

    """

    Convert MIDI pitch values to cents and transpose across the octave in
    steps of 1 cent, then change to pitch classes from C=0 (resolution 1 cent)


    Args:

        pitch_vec (numpy vector of floats):
            The pitch values in MIDI form, to the nearest cent (i.e.,
            measurement precision = 2). Example: 69.45 is A + 45 cents.

    Returns:

        pitch_array (numpy array of ints, shape 1200 x len(pitch_vec)):
            The pitch values in cents from C = 0, where maximum is 1199 (i.e.,
            B + 99 cents) and higher values are wrapped back to zero. Each row
            is a transposition of one cent. First row is the original. Will be
            used as a means of indexing, hence the int type.

    """

    # preconditions
    assert len(pitch_vec) >= 1, "Must include at least one pitch value"
    assert np.all(pitch_vec >= 0), "All pitches must be >= MIDI value 0"
    assert np.all(pitch_vec <= 127), "All pitches must be <= MIDI value 127"

    # function code
    pitch_array = np.array(
        np.mod(
          np.add(
              np.tile(pitch_vec * 100, (1200, 1)),
              np.tile(np.arange(0,1200),(len(pitch_vec),1)).transpose()
          ),
          1200 # remainder boundary
        ),
        dtype = 'int64'
    )

    # postconditions
    assert pitch_array.shape == (1200, len(pitch_vec)), "Incorrect array shape"

    return pitch_array


def nll_tonality(PDF, pitch_array):

    """

    Compute negative log likelihood of values given evaluated points PDF

    Args:

        PDF (numpy vector of floats, shape 1200 x 1):
            The evaluated points of a collection of weighted gaussian curves,
            with 1200 points evaluated per octave, normalized across points.

        pitch_array (numpy array of ints, shape 1200 x len(pitch_vec)):
            The pitch values in cents from tonic = 0, where maximum is 1199
            (e.g, if tonic C = 0, then max is B + 99 cents) and higher values
            are wrapped back to zero. Each row is a transposition of one cent.
            The first row is the original.

    Returns:

        nll_array (numpy array of floats, shape of pitch_array):
            The negative log likelihood for each pitch value, given the PDF

    """

    # preconditions
    assert np.shape(PDF) == (1200,), "Incorrect PDF shape"
    assert np.shape(pitch_array)[0] == 1200, "Incorrect pitch_array shape"

    # Calculate -log of likelihood for each pitch value
    nll_array = -np.log(PDF[pitch_array])

    # postconditions
    assert np.shape(nll_array) == np.shape(pitch_array)

    return nll_array


def optimal_shift(nll_array, pitch_durs=None):

    """

    Calculate mean of -log likelihood array, weighting pitches optional

        Args:

            nll_array (numpy array of floats):
                Negative log likelihood for each pitch value (columns) at
                transpositions of 1 cent across the octave (rows).

            pitch_durs (numpy vector of floats, length = len(nll_array[0])):
                Optional, can be used to weight each pitch column. For example,
                if pitch values are central tendencies (mean/mode), then each
                note should be weighted by its duration. Not applicable if
                using continuous pitch data.

        Returns:

            nll_means (numpy vector of floats):
                The mean (weighted) negative log likelihood at each of 1200
                transpositions.

            nll_min_idx (integer):
                The index (0-1199) of the minimum negative log likelihood. At
                this transposition in cents, the pitch values fit the PDF best.
                Taking the value of 'mean_nll' at this index returns a measure
                of how well the pitch values fit the distribution.

            nll_weighted (numpy array of floats):
                The nll_array, weighted by the relative duration of each note.
                If no durations are provided, this will match nll_array.

    """

    # setting all durations to 1 will ignore duration weighting
    if pitch_durs is None:
        pitch_durs = np.ones_like(nll_array[0])

    # preconditions
    if any(pitch_durs < 0):
        raise ValueError('Negative durations not allowed')

    # weight the nll values by relative duration (lower = better)
    nll_weighted = nll_array * pitch_durs

    # take the mean across each transposition (row), divided by the durs used
    # for weighting to standardize across different durations
    nll_means = np.sum(nll_weighted, axis = 1) / np.sum(pitch_durs)

    # find the best fit(s) across transpositions
    nll_min_idx = np.where(nll_means == nll_means.min())[0]

    return nll_means, nll_min_idx, nll_weighted

def min_idx_to_tonic(nll_min_idx):

    """

    Calculate the tonic from the transposition minimum index (0-1199), where
    0=C.

    The function optimal_shift() transposes pitches upward 1 cent per step. The
    minimum index (nll_min_idx) indicates the amount of *upward* transposition
    to align the pitches (modulus 1200) to a PDF distribution where C=0. Thus, a
    melody with a tonic two cents below B major (tonic 11, offset -0.02) would
    have an upward transposition (nll_min_idx) of 102 cents to reach C. A melody
    with a tonic 15 cents above C# (tonic 1, offset +0.15) would have an upward
    transposition of 1085 cents to reach C, etc.

    Args:

        nll_min_idx (integer):
            The index (0-1199) of the minimum negative log likelihood. At this
            transposition in cents, the pitch values fit the PDF best. Taking
            the value of 'mean_nll' at this index returns a measure of how well
            the pitch values fit the distribution. See optimal_shift().

    Returns:

        tonic (int):
            Whole number from 0 (C) to 11 (B) identifying the closest tonic.

        tonic_offset (float):
            Number from -0.49 to 0.50 expressing the offset from the closest
            tonic in cents.

    """

    # preconditions
    assert nll_min_idx >= 0, "Input index out of bounds"
    assert nll_min_idx <= 1199, "Input index out of bounds"

    # calculate tonic and tonic offset
    a_ = nll_min_idx / 100
    b_ = np.trunc(a_ + np.copysign(0.5, a_)) # avoiding "banker's rounding"
    tonic_offset = np.round(b_ - (nll_min_idx / 100), decimals = 2)
    tonic = np.mod(12 - b_, 12).astype(int)

    # postconditions
    assert tonic >= 0, "Error in tonic calculation"
    assert tonic <= 12, "Error in tonic calculation"
    assert tonic_offset >= -0.49, "Error in tonic offset calculation"
    assert tonic_offset <= 0.50, "Error in tonic offset calculation"

    return tonic, tonic_offset


def notes_in_key_major_minor(pitch_vec, tonic, tonic_offset, mode):

    """

    Calculate proportion of in-key notes for major or minor profiles.

    Args:

        pitch_vec (numpy vector of floats):
            The pitch values in MIDI form, to the nearest cent (i.e.,
            measurement precision = 2). Example: 69.45 is A + 45 cents

        tonic (int): Number from 0 (C) to 11 (B) identifying the
            closest tonic.

        tonic_offset (float):
            Number from -0.49 to 0.50 expressing the offset from the
            closest tonic in cents.

        mode (str):
            'major' or 'minor', to look up the diatonic intervals

    Returns:

        prop_in_key (float):
            Proportion of notes within 50 cents of a diatonic interval

        key (str):
            Key name, uppercase is major.

    """

    # keynames
    name_major = ['C','C#','D','D#','E','F','F#','G','G#','A','A#','B']
    name_minor = ['c','c#','d','d#','e','f','f#','g','g#','a','a#','b']

    # define scales if tonic = 0
    if mode == 'major':
        scale = np.array([0,2,4,5,7,9,11])
        key = name_major[int(tonic)]
    elif mode == 'minor':
        scale = np.array([0,2,3,5,7,8,10])
        key = name_minor[int(tonic)]

    # shift the scale based on the tonic
    scale = np.mod(scale + tonic, 12)

    # center the pitch vector, round to the nearest semitone, and restrict to
    # 0-11
    pitch_vec_c = np.mod(np.round(pitch_vec - tonic_offset), 12)

    # identify which notes are in key
    in_key_vec = [p in scale for p in pitch_vec_c]

    # calculate the proportion of in key notes
    prop_in_key = np.round(np.mean(in_key_vec), decimals = 3)

    return prop_in_key, key


def run_PDF(params, pitch_vec, pitch_durs=None, custom_intervals=None,
    custom_weights=None):

    """

    Convenience function to build a model under some PDF distribution.

    Args:

        params (string):
            Accepts 'major','minor', 'random', or 'custom'. Major and minor
            settings use the K-K profiles (ratings), diatonic intervals only.
            Random chooses seven intervals and seven weights at random, and
            overlaps are common. Custom uses provided intervals and weights.

        pitch_vec (numpy vector of floats):
            The pitch values in MIDI form, to the nearest cent (i.e.,
            measurement precision = 2). Example: 69.45 is A + 45 cents

        pitch_durs (numpy vector of floats, length = len(nll_array[0])):
            Optional, can be used to weight each pitch column. For example,
            if pitch values are central tendencies (mean/mode), then each
            note should be weighted by its duration. Not applicable if
            using continuous pitch data.

        custom_intervals (numpy vector of floats, 0-12):
            Locations of the "diatonic" intervals/notes in semitones, to
            nearest cent. Values must be between 0 and 11.99. Ignored if params
            is not "custom".

        custom_weights (numpy vec of floats):
            For each interval, what is the relative expectation or prevalence?
            Must have the same number of elements as 'custom_intervals'. Ignored
            if params is not "custom".

    Returns:

        model_out (dict):
            Data structure with the results from the model, for plotting,
            saving, etc. Includes:
                params: from args
                pitch_vec: from args
                pitch_durs: from args
                PDF (numpy vector floats): the evaluated probability density
                    function, normalized, length 1200 spanning 1 octave
                intervals (numpy vector floats): the intervals used/generated
                weights (numpy vector floats): the weights of intervals
                    used/generated
                pitch_array (numpy array): pitch_vec transposed across octave
                nll_means (numpy vector floats): negative log likelihood at
                    every shift 0-1199
                nll_min_idx (numpy array int): the index of the optimal shift,
                    value(s) between 0-1199
                nll (numpy array float): negative log likelihood at nll_min_idx,
                    if there are multiple minima, it takes the lowest index
                tonic (int): Number from 0 (C) to 11 (B) identifying the
                    closest tonic. Value is 'nan' if params='random'.
                tonic_offset (float):
                    Number from -0.49 to 0.50 expressing the offset from the
                    closest tonic in cents. Value is 'nan' if params='random'.

    """

    # build PDF
    PDF, intervals, weights = build_diatonic_PDF(params=params)

    # build transposition array
    pitch_array = MIDI_to_transpose_array(pitch_vec)

    # evaluate negative log likelihood
    nll_array = nll_tonality(PDF, pitch_array)
    nll_means, nll_min_idx, nll_weighted = optimal_shift(nll_array, pitch_durs)

    # save minimum value
    # [WARNING] if there are multiple minima, this arbitrarily takes the lowest
    # index value (all are saved with the output regardless)
    # this shouldn't be an issue with the K-S ratings (major/minor PDF) but an
    # artificial PDF with matching weights and intervals (i.e., peaks with the
    # same height and same spacing) could lead to multiple minima.
    nll = nll_means[nll_min_idx[0]]

    # for major and minor PDFs only, we can estimate the tonic and offset, and
    # then determine the name of the key and the proportion of notes in key
    if (params == 'major') | (params == 'minor'):
        tonic, tonic_offset = min_idx_to_tonic(nll_min_idx)
        prop_in_key, key = notes_in_key_major_minor(
            pitch_vec = pitch_vec,
            tonic = tonic,
            tonic_offset = tonic_offset,
            mode = params,
        )
    else:
        tonic = np.nan
        tonic_offset = np.nan
        prop_in_key = np.nan
        key = 'NA'

    # construct the output dictionary
    model_out = {
        'params': params,
        'pitch_vec': pitch_vec,
        'pitch_durs': pitch_durs,
        'PDF': PDF,
        'intervals': intervals,
        'weights': weights,
        'pitch_array': pitch_array,
        'nll_array': nll_array,
        'nll_means': nll_means,
        'nll_weighted': nll_weighted,
        'nll_min_idx': nll_min_idx,
        'nll': nll,
        'tonic': tonic,
        'tonic_offset': tonic_offset,
        'prop_in_key': prop_in_key,
        'key': key,
    }

    return model_out


def plot_model_output(model, model_id=None, dir_out_detailed=None):

    """

    Plot the results of run_PDF().

    Args:

        model (dict):
            Data structure with the results of the model from run_PDF().

        model_id (str, optional):
            The name or ID of the model, prepending each filename.

        dir_out_detailed (str, optional):
            A path of where to save the figure.

    """

    # unpack the model output
    params = model['params']
    pitch_vec = model['pitch_vec']
    pitch_durs = model['pitch_durs']
    PDF = model['PDF']
    pitch_array = model['pitch_array']
    intervals = model['intervals']
    weights = model['weights']
    nll_means = model['nll_means']
    nll_min_idx = model['nll_min_idx']
    nll = model['nll']
    tonic = model['tonic']
    tonic_offset = model['tonic_offset']

    # set plot
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312, sharex=ax1)
    ax3 = fig.add_subplot(313, sharex=ax1)
    fig.tight_layout(pad=3.0)

    # subplot 1
    ax1.title.set_text('PDF evaluated once per cent')
    ax1.plot(PDF)
    ax1.axvline(x = 0, c = 'k', linestyle = '--')
    ax1.axvline(x = 1200, c = 'k', linestyle = '--')

    # subplot 2
    ax2.title.set_text(
        'Pitch distribution before (blue) and after (red) shifting to best fit'
    )
    pitch_to_hist = np.concatenate(
        [np.repeat(p, int(pitch_durs[i]*100)) for i,p in enumerate(pitch_array[nll_min_idx][0])]
    )
    ax2.hist(
        pitch_to_hist,
        bins=np.arange(0,1201,10),
        align='left',
        color='r'
    )
    pitch_to_hist = np.concatenate(
        [np.repeat(p, int(pitch_durs[i]*100)) for i,p in enumerate(pitch_array[0])]
    )
    ax2.hist(
        pitch_to_hist,
        bins=np.arange(0,1201,10),
        color='b',
        alpha=0.15
    )

    # subplot 3
    ax3.title.set_text(
        'Weighted mean negative log likelihood at each transposition'
    )
    ax3.plot(nll_means)
    ax3.axvline(x = nll_min_idx, c = 'r', linestyle = '-')
    ax3.text(
        x = np.mod(nll_min_idx + 25, 800),
        y = 8,
        s = 'minimum at index = ' + str(nll_min_idx),
        bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 3.0)
    )

    if dir_out_detailed is not None:
        filepath = os.path.join(dir_out_detailed, model_id + '.png')
        plt.savefig(filepath, dpi=75)
        plt.close()
    else:
        plt.show()

    return


def save_model_info(model, model_id, dir_out_detailed, save_method='summary',
    plot_model=False):

    """

    Save information generated by run_PDF() in a directory specified by
    dir_out_detailed, with files prepended by model_id.

    [!!! WARNING !!!] will create dir_out_detailed if it does not exist.

    Args:

        model (dict):
            Data structure with results from run_PDF().

        model_id (str):
            The name or ID of the model, prepending each filename. For example,
            a string identifying the melody that was used.

        dir_out_detailed (str):
            Path to directory to save the model info, will create if
            nonexistent.

        save_method (str, 'summary' or 'full'):
            Method of saving.
            'summary': the default, which saves only the summary CSV file
            containing the key, mode, tonic, tonic offset, and proportion of
            in-key notes (PIKN).
            'full': saves all information about the model.

        plot_model (boolean):
            Turn on or off (default) the creation of a plot of model behaviour.
            It can slow things down considerably.

    Returns:

        summary (pandas dataframe):
            A pandas dataframe with the values for key, mode, tonic, tonic
            offset, and proportion in key notes.

    """

    # preconditions
    assert model is not None, "Provide model output from run_PDF()"
    assert model_id is not None, "Provide name of this model (str)"
    assert dir_out_detailed is not None, "Provide model output"
    assert save_method in ['summary','full'], "Invalid parameter: save_method"

    # create the output directory if it does not exist
    if not os.path.exists(dir_out_detailed):
        os.makedirs(dir_out_detailed)

    # save the data in appropriate formats
    os.chdir(dir_out_detailed)

    # the most important summary info gets returned and could be saved alone to
    # save space
    summary = pd.DataFrame({
        'key': model['key'],
        'mode': model['params'],
        'tonic': model['tonic'],
        'tonic_offset': model['tonic_offset'],
        'prop_in_key': model['prop_in_key'],
    })

    summary.to_csv('model_' + model_id + '_summary.csv', index = False)

    if save_method == 'summary':

        pass

    elif save_method == 'full':

        for item in model:

            data = model[item]

            if type(data) is np.ndarray:
                if data.ndim == 1:
                    df = pd.DataFrame(model[item], columns = [item])
                else:
                    columns = ['note_' + str(i) for i in range(data[1].size)]
                    df = pd.DataFrame(model[item], columns = columns)
            else:
                df = pd.DataFrame({item: [model[item]]})

            df.to_csv('model_' + model_id + '_' + item + '.csv', index = False)

    if plot_model is True:
        plot_model_output(model, model_id, dir_out_detailed)

    return summary


def PDF_major_minor(pitch_vec, pitch_durs, dir_out_detailed=None,
    save_method='summary', plot_model=False):

    """

    A convenience function that calculates the tonic using the best fit from
    either the major or minor settings, wrapping run_PDF().

    Args:

        pitch_vec (numpy vector of floats):
            The pitch values in MIDI form, to the nearest cent (i.e.,
            measurement precision = 2). Example: 69.45 is A + 45 cents

        pitch_durs (numpy vector of floats, length = len(nll_array[0])):
            Optional, can be used to weight each pitch column. For example, if
            pitch values are central tendencies (mean/mode), then each note
            should be weighted by its duration. Not applicable if using
            continuous pitch data.

        dir_out_detailed (str):
            Path to directory to save the model info, will create if necessary,
            see save_model_info().

        save_method (str, 'summary' or 'full'):
            Method of saving if dir_out_detailed is provided.
            'summary': the default, which saves only the summary CSV file
            containing the key, mode, tonic, tonic offset, and proportion of
            in-key notes (PIKN).
            'full': saves all information about the model.
            See save_model_info().

        plot_model (boolean):
            Save a visual summary of the model. It is slow, and could be created
            later from the model info if needed. Requires that dir_out_detailed
            has been provided. See plot_model_output().

    Returns:

        model_out (dict):
            Data structure with the results of the model from run_PDF() using
            the best fitting of either major or minor models

    """

    # run each version of the model
    model_maj = run_PDF('major', pitch_vec, pitch_durs)
    model_min = run_PDF('minor', pitch_vec, pitch_durs)

    # find the best fit between the two (lowest nll score)
    if model_maj['nll'] == model_min['nll']:
        raise ValueError('! nll major == minor, cannot determine best fit !')
    elif model_maj['nll'] < model_min['nll']:
        model_out = model_maj
    elif model_maj['nll'] > model_min['nll']:
        model_out = model_min

    # save the model
    if dir_out_detailed is not None:
        save_model_info(
            model = model_out,
            model_id = 'actual',
            dir_out_detailed = dir_out_detailed,
            save_method = save_method,
            plot_model = plot_model,
        )

    return model_out


def generate_random_pitches(pitch_vec, strategy='flat_octave'):

    """

    Generate a sequence of random pitch values matching the length of pitch_vec.

    Args:

        pitch_vec (numpy vector of floats):
            The pitch values in MIDI form, to the nearest cent (i.e.,
            measurement precision = 2). Example: 69.45 is A + 45 cents.

        strategy (string of 'flat_octave', 'flat_range', or 'gaussian_range'):
            'flat_octave': notes are drawn from across the octave (0-11.99)
            with equal probability. Default.
            'flat_range': notes are drawn from the range of the input pitches
            (min to max) with equal probability
            'gaussian_range': notes are drawn from the range of the input pitches
            (min to max) but from a Gaussian distribution centered on the M/SD of
            the input pitches.

    Returns:

        pitch_vec_rand (numpy vector of floats):
            A set of pitch values determined using one of the randomization
            strategies, matching the shape of pitch_vec.

    """

    # preconditions
    assert len(pitch_vec) >= 1, "Must include at least one pitch value"
    assert np.all(pitch_vec >= 0), "All pitches must be >= MIDI value 0"
    assert np.all(pitch_vec <= 127), "All pitches must be <= MIDI value 127"
    assert strategy in ['flat_octave','flat_range','gaussian_range'], "Invalid parameter"

    # function code to build pitch_vec_rand based on strategy
    if strategy == 'flat_octave':

        # generate random pitch values across the octave (0-12)
        pitch_vec_rand = np.round(
            np.random.rand(len(pitch_vec)) * 12,
            decimals = 2
        )

    elif strategy == 'flat_range':

        # generate random pitch values from flat dist. in performance range
        low_note = np.min(pitch_vec)
        high_note = np.max(pitch_vec)
        size_notes = len(pitch_vec)
        pitch_vec_rand = np.round(
            np.random.uniform(
                low = low_note, high = high_note, size = size_notes
            ),
            decimals = 2
        )

    elif strategy == 'gaussian_range':

        # generate random pitch values from gaussian dist. centered around
        # performance mean/sd and truncated by min/max pitches
        low_note = np.min(pitch_vec)
        high_note = np.max(pitch_vec)
        mu = np.mean(pitch_vec)
        sigma = np.std(pitch_vec)
        size_notes = len(pitch_vec)
        pitch_vec_rand = np.round(
            stats.truncnorm.rvs(
                (low_note - mu) / sigma,
                (high_note - mu) / sigma,
                loc = mu,
                scale = sigma,
                size = size_notes
            ), decimals = 2
        )

    # postconditions
    assert pitch_vec.shape == pitch_vec_rand.shape, "Incorrect array shape"

    return pitch_vec_rand


def PDF_null_notes(pitch_vec, pitch_durs, strategy=None, n=1000,
    dir_out_detailed=None, save_method='summary', plot_model=False):

    """

    A convenience function that calculates a null distribution of models using
    PDF_major_minor(), but generating random pitches for pitch_vec using
    generate_random_pitches().

    Args:

        pitch_vec (numpy vector of floats):
            The pitch values in MIDI form, to the nearest cent (i.e.,
            measurement precision = 2). Example: 69.45 is A + 45 cents

        pitch_durs (numpy vector of floats, length = len(nll_array[0])):
            Optional, can be used to weight each pitch column. For example, if
            pitch values are central tendencies (mean/mode), then each note
            should be weighted by its duration. Not applicable if using
            continuous pitch data.

        strategy (string of 'flat_octave', 'flat_range', or 'gaussian_range'):
            'flat_octave': notes are drawn from across the octave (0-11.99)
            with equal probability. Default if not specified.
            'flat_range': notes are drawn from the range of the input pitches
            (min to max) with equal probability
            'gaussian_range': notes are drawn from the range of the input pitches
            (min to max) but from a Gaussian distribution centered on the M/SD of
            the input pitches.
            See generate_random_pitches().

        n (int):
            Number of iterations, default is 1000

        dir_out_detailed (str) [OPTIONAL]:
            Path to directory to save the model info, will create if necessary

        save_method (str, 'summary' or 'full'):
            Method of saving if dir_out_detailed is provided.
            'summary': the default, which saves only the summary CSV file
            containing the key, mode, tonic, tonic offset, and proportion of
            in-key notes (PIKN).
            'full': saves all information about the model.
            See save_model_info().

        plot_model (boolean):
            Save a visual summary of the model. It is slow, and could be created
            later from the model info if needed. Requires that dir_out_detailed
            has been provided. See plot_model_output().

    Returns:

        null_summary (pandas dataframe):
            Summary of results from each iteration, including the iteration
            number, the resulting key, mode, tonic, tonic offset, and proportion
            of in-key notes.

        null_notes (pandas dataframe):
            The random notes generated by the algorithm at each iteration.

    """

    if strategy is None:
        strategy = 'flat_octave'

    # build null models
    null_summary = []
    null_pitches = []
    for i in track(range(n), description=f"{n} null models, random pitches..."):

        # generate random pitch values according to strategy
        pitch_vec_rand = generate_random_pitches(pitch_vec, strategy)

        # build model
        model_n = PDF_major_minor(pitch_vec_rand, pitch_durs)

        # record the full model results - quite slow
        if dir_out_detailed is not None:
            save_model_info(
                model = model_n,
                model_id = 'randNotes_'+str(i),
                dir_out_detailed = dir_out_detailed,
                save_method = save_method,
                plot_model = plot_model,
            )

        # append to output
        null_summary.append(
            pd.DataFrame({
                'model': f"randNotes_{i:03d}",
                'key': model_n['key'],
                'mode': model_n['params'],
                'tonic': model_n['tonic'],
                'tonic_offset': model_n['tonic_offset'],
                'prop_in_key': model_n['prop_in_key'],
            })
        )
        null_pitches.append(
            pd.DataFrame({
                f"randNotes_{i:03d}": model_n['pitch_vec'],
            })
        )

    null_summary = pd.concat(null_summary)
    null_pitches = pd.concat(null_pitches, axis=1)

    return null_summary, null_pitches


def z_score_null_notes(actual_PIKN, null_PIKN, p_max=None):

    """

    Convert the proportion of notes that are tonal in the actual performance
    into a z-score relative to a null distribution of proportions.

    Args:

        actual_PIKN (float):
            The proportion of in-key notes (PIKN) in the actual performance.

        null_PIKN (numpy vector of floats):
            The proportions of in-key notes (PIKN) from models generated under
            randomly-assigned pitch values but otherwise matching the original
            melody in number of notes and note durations.

        p_max (float, OPTIONAL):
            An optional maximum value for the proportion representing
            the maximum percentile value that can be used to calculate a z-score
            from percentiles. If None, it will be based on the number of items
            in the null distribution [1 - (1/n)/2]. It has a strong effect on
            where the z-scores are truncated.

    Returns:

        z_out (dict):
            Data structure with the results from the z-scoring, for plotting,
            saving, etc. Includes:

                z (float): z-score representing proportion of tonal notes
                relative to chance, calculated in the standard method from the
                mean and standard deviation of the null distribution. This
                method can create large z-scores, especially with longer
                sequences.

                z_pctl (float): z-score representing proportion of tonal notes
                relative to chance, calculated from the percentile rank of the
                PIKN against the null distribution. Proportions of 1.0
                (percentiles of 100.0) must be corrected to avoid infinite or
                indeterminate results. This has an effect of truncating z-scores
                beyond a certain value regardless of sequence length, which may
                be desirable, but may violate assumptions of normality.

                actual_PIKN (float): the proportion of in-key notes for the
                actual pitch values (pre- z-score).

                null_PIKN (vector): a vector of floats, each value being the
                proportion of in-key notes for one of the null models.

    """

    if p_max is not None:
        assert p_max > 0, "p_max must be greater than 0"
        assert p_max < 1.0, "p_max must be less than 1.0"

    # calculate the z-score in the standard method, (score - M) / SD
    z_PIKN = (actual_PIKN - np.mean(null_PIKN)) / np.std(null_PIKN)

    # calculate the z-score from the percentile rank of the score in its null
    # distribution, by putting that value (as proportion) through norm.ppf()
    # this generally agrees with the standard method of z-score calculation, but
    # it is limited by a correction from proportion=1.0 (inf). Here, it is
    # corrected by taking [1-(1/n)/2] where n is the number of items in the null
    # array (i.e., with n=1000, max proportion=0.9995 and max abs z=3.29)
    if p_max is None:
        n = len(null_PIKN)
        p_max = 1 - (1/n)/2
    pctl = stats.percentileofscore(null_PIKN, actual_PIKN)
    z_pctl = stats.norm.ppf(np.min([pctl/100, p_max]))

    # prepare for saving
    z_out = {
        'z': z_PIKN,
        'z_pctl': z_pctl,
        'actual_PIKN': actual_PIKN,
        'null_PIKN': null_PIKN,
        'pctl': pctl,
    }

    return z_out


def plot_PIKN_distribution(value, distribution, model_id=None,
    dir_out_detailed=None):

    """

    Visualize proportion of in-key notes value in the null distribution.

    Args:

        value (float):
            Number representing the value of interest

        distribution (numpy array of floats):
            Array of values from null distribution

        model_id (str):
            (Optional) The name or ID of the model, prepending each filename

        dir_out_detailed (string):
            (Optional) Where to save the figure

    """

    # generate a figure that visualizes dv
    dv = 'Proportion of in-key notes'
    plt.figure()
    plt.axvline(x = value,
                color = 'black',
                linestyle = ':',
                alpha = 0.5)
    plt.hist(distribution, bins=50);
    plt.xlim([-0.05,1.05])
    plt.title(dv)
    plt.xlabel('Proportion')
    plt.ylabel('Frequency')

    # save and close (optional)
    if dir_out_detailed is not None:
        filepath = os.path.join(dir_out_detailed, model_id + '.png')
        plt.savefig(filepath, dpi=3000)
        plt.close()
    else:
        plt.show()

    return