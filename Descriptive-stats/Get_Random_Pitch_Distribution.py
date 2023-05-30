import numpy as np
import pandas as pd

# Define the length of the song in seconds
song_length = 30

# Define the number of notes in the Western scale
scale_length = 12

# Define the number of songs to generate
num_songs = 30

# Generate random distributions of notes for each song
songs = []
for i in range(num_songs):
    notes = np.random.randint(low=0, high=scale_length, size=song_length)
    songs.append(notes)

# Convert the list of songs into a Pandas DataFrame
df = pd.DataFrame(songs)

# Rename the columns of the DataFrame to indicate the time in seconds
df.columns = [f"sec_{i+1}" for i in range(song_length)]

# Display the DataFrame
print(df)
