import sys
import librosa

from audio_to_midi.monophonic import wave_to_midi

print("Starting...")
file_in = sys.argv[1]
file_out = sys.argv[2]
y, sr = librosa.load(filename = "recording_999_mic_1.wav", sr=None)
print("Audio file loaded!")
midi = wave_to_midi(y, sr=sr)
print("Conversion finished!")
with open (file_out, 'wb') as f:
    midi.writeFile(f)
print("Done. Exiting!")