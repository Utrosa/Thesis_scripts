# In this example, the primes are chords.
# from psynet.js_synth import Chord, Note

# intervals
PRIME_INTERVALS = {
    "major-triad": [4, 3],
    "minor-triad": [3, 4],
}

# triads
# PRIMES_TRIADS = {
#     # female
#     # "major-triad": [58, 62, 65],
#     # "minor-triad": [58, 61, 65],
#     "major-triad": [60, 64, 67],  # C major
#     "minor-triad": [60, 63, 67],  # C minor
# }

# PRIMES_TRIADS = {
#     "major-triad": [
#         Chord([58, 62, 65], duration=2.0, silence=0.0)
#     ],
#     "minor-triad": [
#         Chord([58, 61, 65], duration=2.0, silence=0.0)
#     ],
# }

# scales
# PRIMES_SCALES = {
#     "major-scale": [60, 62, 64, 65, 67, 69, 71, 72],
#     "minor-scale": [60, 62, 63, 65, 67, 68, 71, 72],
# }

#
# PRIMES_SCALES = {
#     "Major scale": [
#         Note(60, duration=1.2, silence=0.0),
#         Note(62, duration=0.6, silence=0.0),
#         Note(64, duration=0.6, silence=0.0),
#         Note(65, duration=0.6, silence=0.0),
#         Note(67, duration=0.6, silence=0.0),
#         Note(69, duration=0.6, silence=0.0),
#         Note(71, duration=0.6, silence=0.0),
#         Note(72, duration=1.2, silence=0.0),
#     ],
#     "Minor scale": [
#         Note(60, duration=1.2, silence=0.0),
#         Note(62, duration=0.6, silence=0.0),
#         Note(63, duration=0.6, silence=0.0),
#         Note(65, duration=0.6, silence=0.0),
#         Note(67, duration=0.6, silence=0.0),
#         Note(68, duration=0.6, silence=0.0),
#         Note(71, duration=0.6, silence=0.0),
#         Note(72, duration=1.2, silence=0.0),
#     ],
# }

# others
# PRIMES = {
#     prime_id: [Chord(pitches, duration=2.0, silence=0.0, timbre="prime")]
#     for prime_id, pitches in {
#         "Major triad": [58, 62, 65],
#         "Minor triad": [58, 61, 65],
#         # "Augmented triad": [58, 62, 66],
#         # "Diminished triad": [58, 61, 65]
#     }.items()
# }


# PRIMES = {
#     prime_id: [Chord(pitches, duration=2.0, silence=0.0)]
#     for prime_id, pitches in {
#         "Major triad": [60, 64, 67],
#         "Minor triad": [60, 63, 67],
#         "Augmented triad": [60, 64, 68],
#         "Diminished triad": [60, 63, 66],
#         "Phrygian triad": [60, 61, 67],
#         "Lydian triad": [60, 66, 67],
#         "Locrian triad 1": [60, 61, 66],
#         "Locrian triad 2": [60, 65, 66]
#     }.items()
# }

# PRIMES = {
#     prime_id: [
#         Note(pitch, duration=duration, silence=0.0)
#         for pitch, duration in zip(pitches, [1.0] * 3)
#     ]
#     for prime_id, pitches in {
#         "Major triad": [60, 64, 67],
#         "Minor triad": [60, 63, 67],
#         "Augmented triad": [60, 64, 68],
#         "Diminished triad": [60, 63, 66],
#         "Phrygian triad": [60, 61, 67],
#         "Lydian triad": [60, 66, 67],
#         "Locrian triad 1": [60, 61, 66],
#         "Locrian triad 2": [60, 65, 66]
#     }.items()
# }

# In this example, the primes are melodies.
# Here the primes are written out explicitly, instead of using any comprehensions.
#
# PRIMES = {
#     "Major arpeggio": [
#         Note(60, duration=1.2, silence=0.0),
#         Note(64, duration=0.6, silence=0.0),
#         Note(67, duration=0.6, silence=0.0),
#         Note(72, duration=1.2, silence=0.0),
#     ],
#     "Minor arpeggio": [
#         Note(60, duration=1.2, silence=0.0),
#         Note(63, duration=0.6, silence=0.0),
#         Note(67, duration=0.6, silence=0.0),
#         Note(72, duration=1.2, silence=0.0),
#     ],
# }
