## IMPORT --------------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

from glob import glob

import librosa
import librosa.display
import IPython.display as ipd

from itertools import cycle

sns.set_theme(style="white", palette=None)
color_pal = plt.rcParams["axes.prop_cycle"].by_key()["color"]
color_cycle = cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])

## LOAD & DESCRIBE ------------------------------------------------------------------------------

y, sr = librosa.load('C:/Users/monik/Documents/Paris/Internship/Analyses/999/Improv/999_recording_mic_23.wav')
print(f'y: {y[:10]}')
print(f'shape y: {y.shape}')
print(f'sr: {sr}')

## SPECTROGRAM ---------------------------------------------------------------------------------

D = librosa.stft(y)
S_db = librosa.amplitude_to_db(np.abs(D), ref=np.max)
S_db.shape

# Plot the transformed audio data
fig, ax = plt.subplots(figsize=(10, 5))
img = librosa.display.specshow(S_db,
                              x_axis='time',
                              y_axis='log',
                              ax=ax)
ax.set_title('Spectogram Example', fontsize=20)
fig.colorbar(img, ax=ax, format=f'%0.2f')
plt.savefig("Spectrogram-999-Improvised-Song-23.png", format="png")
#plt.show()

## MELSPECTROGRAM -----------------------------------------------------------------------------

hop_length = 512
S = librosa.feature.melspectrogram(y=y,
                                   sr=sr,
                                   hop_length=hop_length)
S_db_mel = librosa.power_to_db(S, ref=np.max)

fig, ax = plt.subplots(figsize=(10, 5))
# Plot the mel spectogram
img = librosa.display.specshow(S_db_mel,
                              x_axis='time',
                              y_axis='mel', # y_axis can either be mel or log which affects the scale
                              hop_length=512)

ax.set_title('Mel Spectogram Example', fontsize=20)
fig.colorbar(img, ax=ax, format=f'%0.2f')
plt.tight_layout()
plt.savefig("Mel-Spectrogram-999-Improvised-Song-23.png", format="png")

## PLOT ----------------------------------------------------------------------------------------
#plt.show()