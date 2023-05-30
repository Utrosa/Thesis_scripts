''' Cyclic tempogram—A mid-level tempo representation for musicsignals '''

## IMPORT ------------------------------------------------------------------------

import matplotlib.pyplot as plt
import librosa
import numpy as np

## LOAD --------------------------------------------------------------------------

filename = "C:/Users/monik/Documents/Paris/Internship/Analyses/999/Improv/999_recording_mic_23.wav"
y, sr = librosa.load(filename)

## FFT COMPUTATION ----------------------------------------------------------------
    # Compute local onset autocorrelation

hop_length = 512
oenv = librosa.onset.onset_strength(y=y, sr=sr, hop_length=hop_length)
tempogram = librosa.feature.fourier_tempogram(onset_envelope=oenv, sr=sr,
                                              hop_length=hop_length)

    # Compute the auto-correlation tempogram, unnormalized to make comparison easier
ac_tempogram = librosa.feature.tempogram(onset_envelope=oenv, sr=sr,
                                         hop_length=hop_length, norm=None)

## PLOT --------------------------------------------------------------------------

fig, ax = plt.subplots(nrows=3, sharex=True)
ax[0].plot(librosa.times_like(oenv), oenv, label='Onset strength')
ax[0].legend(frameon=True)
ax[0].label_outer()
librosa.display.specshow(np.abs(tempogram), sr=sr, hop_length=hop_length,
                         x_axis='time', y_axis='fourier_tempo', cmap='magma',
                         ax=ax[1])
ax[1].set(title='Fourier tempogram')
ax[1].label_outer()

librosa.display.specshow(ac_tempogram, sr=sr, hop_length=hop_length,
                         x_axis='time', y_axis='tempo', cmap='magma',
                         ax=ax[2])
ax[2].set(title='Autocorrelation tempogram')

plt.tight_layout()
plt.savefig("Beat-Auto-Correlation-999-Improvised-Song-23.png", format="png")

## SHOW ----------------------------------------------------------------
plt.show()