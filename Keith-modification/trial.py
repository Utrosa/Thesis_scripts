"""Read hums, get f0 and segment to notes."""
from itertools import groupby

import librosa
import librosa.display
import matplotlib.pyplot as plt
import mir_eval
import numpy as np
import scipy
import soundfile
import pandas as pd
import os

# DONE: Add time to rotational maps using color of the lines
# FAIL: add curve to the lines and arrows to see both directions
# DONE: play the pitch derivation and original sound together
plt.plot()
plt.close()

# assign directory
directory = "C:/Users/monik/Documents/Paris/Internship/Analyses/999/Improv"
 
# iterate over files in
# that directory
filenames = []
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    # checking if it is a file
    if os.path.isfile(f):
        filenames.append(f)

# while loop might do it faster
for file in filenames:
    y, sr = librosa.load(f)
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
    #plt.show()
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
    #plt.show()
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

    ######## ----------------------- save values into csv file
         

    # dictionary of lists  
    dict = {'pitch': hz, 'note': midi_notes}  
           
    df = pd.DataFrame(dict)
        
    # saving the dataframe 
    df.to_csv('trial_999_hz_and_notes.csv') 


print("Done. Exiting!")