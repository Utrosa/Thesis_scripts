"""Read hums, get f0 and segment to notes."""
from itertools import groupby

import librosa
import librosa.display
import matplotlib.pyplot as plt
import mir_eval
import numpy as np
import scipy
import soundfile

# DONE: Add time to rotational maps using color of the lines
# FAIL: add curve to the lines and arrows to see both directions
# DONE: play the pitch derivation and original sound together
plt.plot()
plt.close()
# read in file
# filename = "/Users/kdkeithd/Downloads/badsinging/data/multitrial_badsinging_2022-07-22_15h02.44.618_mic_recorded/recording_mic_2022-07-22_15h04.01.229.wav"
filename = "recording_999_mic_1.wav"
y, sr = librosa.load(filename)
# get spectrogram
D = librosa.stft(y)  # STFT of y
S_db = librosa.amplitude_to_db(np.abs(D), ref=np.max)


###############################################################################
#                                    get f0                                   #
###############################################################################
# plot spectrogram
fig, ax = plt.subplots(figsize=(10, 4))
img = librosa.display.specshow(S_db, x_axis="time", y_axis="log", ax=ax, sr=sr)
fig.colorbar(img, ax=ax, format="%2.f dB")
f0, vid, vpd = librosa.pyin(
    y,
    sr=sr,
    fmin=librosa.note_to_hz("C1"),
    fmax=librosa.note_to_hz("C7"),
    resolution=0.1,
)
tf = np.linspace(0, len(y) / sr, len(f0))
ax.plot(tf, f0)
ax.set(title=filename, ylim=[75, 475])
fig.show()
fig.savefig(f"f0_{filename.replace('.wav','')}.png")

###############################################################################
#                      Identify chroma at each time point                     #
###############################################################################
# optimize filtering
chroma = librosa.feature.chroma_cqt(y=y, sr=sr, norm=None)
chroma_filter = np.minimum(
    chroma, librosa.decompose.nn_filter(chroma, aggregate=np.median, metric="cosine")
)
kernel = np.tile(np.array([-0.5, 1, -0.5])[:, None], (1, 9))
kernel /= np.abs(kernel).sum()
chroma_smooth = scipy.ndimage.filters.convolve(chroma_filter, kernel)
# threshold and plot
chroma_thresh = chroma_smooth
chroma_thresh[chroma_thresh < 0.1] = -np.inf
fig, ax = plt.subplots(nrows=2, sharex=True, figsize=(10, 6))
img = librosa.display.specshow(S_db, x_axis="time", y_axis="log", ax=ax[0], sr=sr)
fig.colorbar(img, ax=ax[0], format="%2.f dB")
tf = np.linspace(0, len(y) / sr, len(f0))
ax[0].plot(tf, f0)
ax[0].set(title=filename, ylim=[75, 475])
img = librosa.display.specshow(chroma_thresh, y_axis="chroma", x_axis="time", ax=ax[1])
cb = fig.colorbar(img, ax=ax[1])
ax[1].set(title="Chroma With Filtering and Smoothing")
index = np.argmax(chroma_thresh, axis=0).astype(float)
nanid = np.all(np.isinf(chroma_thresh), axis=0)
index[nanid] = np.nan
inList = list()
for k, g in groupby(index):
    seg = [*g]
    if len(seg) > 3:
        inList += seg
    else:
        inList += [np.nan] * len(seg)
index = np.array(inList)
ax[1].plot(tf, index, ".")
ax[1].set_xlim([0, 14])
fig.show()
fig.savefig(f"f0vChroma_{filename.replace('.wav', '')}")

###############################################################################
#   Math octaves maintaining duration and compare transcription in wavefile   #
###############################################################################
index_labels = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
hz = list()
midi_notes = list()
for k, g in groupby(zip(index, f0), key=lambda x: x[0]):
    freqs = [t[1] for t in g]
    if np.isnan(k):
        hz.append(k)
        midi_notes.append(k)
    else:
        currlab = index_labels[int(k)]
        posF = librosa.note_to_hz([f"{currlab}{i}" for i in range(1, 7)])
        avF = np.nanmean(freqs)
        octI = np.argmin(np.abs(posF - avF)) + 1
        note_name = f"{currlab}{octI}"
        freq = librosa.note_to_hz(note_name)
        midi = librosa.note_to_midi(note_name)
        validSamp = ~np.isnan(freqs)
        hz += list(freq * validSamp)
        midi_notes += [midi if vS else np.nan for vS in validSamp]
hz = np.array(hz)
hz[np.isnan(hz)] = -1
y_chroma = mir_eval.sonify.pitch_contour(tf, hz, sr, length=len(y))
rms_y = np.sqrt(np.mean(y ** 2))
rms_yC = np.sqrt(np.mean(y_chroma ** 2))
y_chroma *= 2 * rms_y / rms_yC
side_by_side = np.vstack([y, y_chroma]).T
soundfile.write(filename.replace(".wav", "_transcript.wav"), side_by_side, sr)

###############################################################################
#         Plot transition probabilites from chroma to chroma in polar         #
###############################################################################

fig, ax = plt.subplots(ncols=2, subplot_kw={"projection": "polar"})
fphase = 2 * np.pi * (f0 - 125) / 125
ax[0].plot(fphase, tf)
ax[0].set(title="raw F0")
cphase = 2 * np.pi * index / 12
ax[1].plot(cphase, tf)
ax[1].set(title="Chroma")
fig.show()

sequence = np.array(
    [[k, len([*g])] for k, g in groupby(index[~np.isnan(index)])]
).astype(int)
trans = np.zeros([12, 12])
for i in range(1, len(sequence)):
    curr, prev = sequence[i, 0], sequence[i - 1, 0]
    trans[prev, curr] += 1

med_durations = np.zeros([12, 12])
for prev in range(12):
    for curr in range(12):
        # get median time length for each
        select = np.insert(
            (sequence[:-1, 0] == prev) & (sequence[1:, 0] == curr), 0, False
        )
        med_durations[prev, curr] = np.median(sequence[select, 1])
fig, ax = plt.subplots(ncols=2, subplot_kw={"projection": "polar"})
ax[0].plot(sequence[:, 0], np.ones(len(sequence)))
ax[0].set(title="Overlapping sequence")
cmap = plt.cm.viridis
vmax = 30.0
for prev in range(12):
    pPh = 2 * np.pi * prev / 12
    for curr in range(12):
        cPh = 2 * np.pi * curr / 12
        thisTrans = trans[prev, curr]
        thisDur = med_durations[prev, curr]
        if thisTrans > 0:
            ax[1].plot(
                [pPh, cPh],
                [1, 1],
                lw=6 * thisTrans / trans.max(),
                color=cmap(np.minimum(thisDur / vmax, 1.0)),
            )
values, size = np.unique(sequence[:, 0], return_counts=True)
vPh = 2 * np.pi * values / 12
ax[1].scatter(vPh, np.ones_like(vPh), color="k", s=6 * size, zorder=100)
# plot the sequence and the transition probabilities.
ax[1].set(title="Transition probabilities")
sm = plt.cm.ScalarMappable(plt.Normalize(0, vmax * 512 / sr), cmap)
plt.tight_layout()
cax = fig.add_subplot(10, 2, 20)
cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticklocation="top")
fig.show()


###############################################################################
#            keep track of octaves to make relative pitch work well           #
###############################################################################
# for each note find the mean f0 value to best estimate absolute interval
index_labels = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
note_names = list()
for k, g in groupby(
    zip(index[~np.isnan(index)], f0[~np.isnan(index)]), key=lambda x: x[0]
):
    freq = [t[1] for t in g]
    nSamps = sum(~np.isnan(freq))
    if nSamps > 0:
        currlab = index_labels[int(k)]
        posF = librosa.note_to_hz([f"{currlab}{i}" for i in range(1, 7)])
        avF = np.nanmean(freq)
        octI = np.argmin(np.abs(posF - avF)) + 1
        note_names.append(f"{currlab}{octI}")
note_names = np.array(note_names)
midi_index = librosa.note_to_midi(note_names)
semitones = np.diff(midi_index)
fig = plt.figure()
plt.hist(semitones, bins=np.arange(semitones.min() - 0.5, semitones.max() + 0.5, 1))
fig.show()


###############################################################################
#                        plot transition probabilities                        #
###############################################################################
xy = range(np.min([-12, semitones.min()]), np.max([13, semitones.max() + 1]))
semitrans = np.zeros([len(xy)] * 2)
for i in range(1, len(semitones)):
    curr, prev = semitones[i], semitones[i - 1]
    semitrans[xy.index(prev), xy.index(curr)] += 1.0
fig, ax = plt.subplots()
semitrans[semitrans == 0] = np.nan
# semitrans /= np.nansum(semitrans, axis=1, keepdims=True)
semitrans /= np.nansum(semitrans)
img = ax.pcolormesh(xy, xy, semitrans, shading="nearest", cmap=plt.cm.cividis)
ax.set(xlabel="Current Interval", ylabel="Previous Interval")
fig.colorbar(img, ax=ax)
fig.show()
