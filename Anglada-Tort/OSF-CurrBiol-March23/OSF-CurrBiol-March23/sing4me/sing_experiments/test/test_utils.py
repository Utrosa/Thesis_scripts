from .. import melodies


def test_sample_interval_sequence():
    for i in range(100):
        x = melodies.sample_interval_sequence(
            n_int = 5,
            max_interval_size=8,
            max_melody_pitch_range=10,
            discrete=False,
            reference_mode="first_note",
        )
        assert len(x) == 5
        for interval in x:
            assert abs(interval) <= 8


def test_is_valid_interval_sequence():
    par = dict(
        n_int=2,
        max_interval_size=7,
        max_melody_pitch_range=20,
        reference_mode="previous_note",
    )
    assert melodies.is_valid_interval_sequence(intervals=[4, 7], **par)
    assert not melodies.is_valid_interval_sequence(intervals=[4, 8], **par)

    par["reference_mode"] = "first_note"
    par["max_melody_pitch_range"] = 5

    assert melodies.is_valid_interval_sequence(intervals=[4, 5], **par)
    assert melodies.is_valid_interval_sequence(intervals=[-4, -5], **par)

    par["reference_mode"] = "previous_note"

    assert not melodies.is_valid_interval_sequence(intervals=[4, 2], **par)
    assert not melodies.is_valid_interval_sequence(intervals=[-4, -2], **par)
