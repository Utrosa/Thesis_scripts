
''' Dynamic programming beat tracker '''

''' Modifed example code available online: librosa beat module documentation'''

# IMPORT ----------------------------------------------------------------------------------------
import librosa
import matplotlib.pyplot as plt
import numpy as np

## LOAD ----------------------------------------------------------------------------------------
    #To analyze a subset of audio files, select the filenames with the glob function. Then loop.
    # audio_files = glob('C:/Users/monik/Documents/Paris/Internship/Analyses/999/Improv/*.wav')
    # The * selects all files in that directory. 
    # The order of items in the list of the filepaths is alphabetical. audio_files[3] will not
    # be the fourth song because of missing zeros 01, 02, 03, ... !!

filename = "C:/Users/monik/Documents/Paris/Internship/Analyses/999/Improv/999_recording_mic_23.wav"

    # Load the audio as a waveform `y`
    # Store the sampling rate as `sr`

y, sr = librosa.load(filename)

## BEAT ----------------------------------------------------------------------------------------
    # Run the beat tracker using a pre-computed onset envelope
    # Default: tempo, beat_frames = librosa.beat.beat_track(y=y, sr=sr)

onset_env = librosa.onset.onset_strength(y=y, sr=sr,
                                         aggregate=np.median)
tempo, beats = librosa.beat.beat_track(onset_envelope=onset_env,
                                       sr=sr)

print('Estimated tempo: {:.2f} beats per minute'.format(tempo))

    # Convert the frame indices of beat events into timestamps
beat_times = librosa.frames_to_time(beats, sr=sr)

    # Plot the beat events against the onset strength envelope
hop_length = 512
fig, ax = plt.subplots(nrows=2, sharex=True)
times = librosa.times_like(onset_env, sr=sr, hop_length=hop_length)


## MELSPECTROGRAM ------------------------------------------------------------------------------
M = librosa.feature.melspectrogram(y=y, 
                                    sr=sr, 
                                    hop_length=hop_length)
M_db_mel = librosa.power_to_db(M, ref=np.max) # or amplitude_to_db
librosa.display.specshow(M_db_mel,
                         x_axis='time',
                         y_axis='mel',  # or log
                         hop_length=hop_length,
                         ax=ax[0])

## PLOT ----------------------------------------------------------------------------------------

ax[0].label_outer()
ax[0].set(title='Mel spectrogram')

ax[1].set(title='Beat tracking')
ax[1].plot(times, librosa.util.normalize(onset_env),
         label='Onset strength')
ax[1].vlines(times[beats], 0, 1, alpha=0.5, color='r',
           linestyle='--', label='Beats')
ax[1].legend()

plt.savefig("Beat-tracking-999-Improvised-Song-23.png", format="png")

## SHOW ----------------------------------------------------------------------------------------
#plt.show()