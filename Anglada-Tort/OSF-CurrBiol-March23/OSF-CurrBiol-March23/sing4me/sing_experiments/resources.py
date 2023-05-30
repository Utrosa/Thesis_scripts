# shared resources for melody experiments
import os
import shutil
import jsonpickle
import json
import numpy as np
from statistics import mean
from flask import Markup
from psynet.page import InfoPage
from psynet.trial.audio import AudioRecordTrial
from psynet.modular_page import (
    AudioMeterControl,
    ModularPage,
    PushButtonControl,
    Prompt,
    AudioPrompt,
    AudioRecordControl,
)
from psynet.timeline import (
    conditional,
    CodeBlock,
    join,
    Event,
    Module,
    PageMaker,
    ProgressDisplay,
    ProgressStage,
)
from psynet.trial.static import (
    StaticTrial,
    StaticTrialMaker,
    StimulusSet,
    StimulusSpec,
    StimulusVersionSpec
)
from psynet.utils import get_logger
from psynet.js_synth import JSSynth, Note
from psynet.js_synth import InstrumentTimbre, HarmonicTimbre
# sing4me
from . import params
from . import melodies
from sing4me import singing_extract as sing
config = params.singing_1interval  # 1 interval
logger = get_logger()


SYLLABLE = "TA"


def select_timbre(selected_timbre):
    # function to select timbre, thee options: piano, complex_short_ISI_long, complex_long_ISI_long, complex_short_ISI_short, complex_mid_ISI_long
    if selected_timbre == "piano":
        TIMBRE = InstrumentTimbre("piano")
        note_duration_tonejs = params.note_duration
        note_silence_tonejs = params.note_silence
    elif selected_timbre == "complex_short_ISI_long":
        note_duration_tonejs = 0.8
        note_silence_tonejs = 0
        TIMBRE = dict(
            default=HarmonicTimbre(
                attack=0.01,  # Attack phase duration in seconds
                decay=0.05,  # Decay phase duration in seconds
                sustain_amp=0.8,  # Amplitude fraction to decay to relative to max amplitude --> 0.4, 0.7
                release=0.7,  # Release phase duration in seconds
                num_harmonics=10,  # Actual number of partial harmonics to use
                roll_off=14,  # Roll-off in units of dB/octave,
            )
        )
    elif selected_timbre == "complex_mid_ISI_long":
        note_duration_tonejs = 0.8
        note_silence_tonejs = 0
        TIMBRE = dict(
            default=HarmonicTimbre(
                attack=0.01,  # Attack phase duration in seconds
                decay=0.05,  # Decay phase duration in seconds
                sustain_amp=0.8,  # Amplitude fraction to decay to relative to max amplitude --> 0.4, 0.7
                release=0.55,  # Release phase duration in seconds
                num_harmonics=10,  # Actual number of partial harmonics to use
                roll_off=14,  # Roll-off in units of dB/octave,
            )
        )
    elif selected_timbre == "complex_long_ISI_long":
        note_duration_tonejs = 0.8
        note_silence_tonejs = 0
        TIMBRE = dict(
            default=HarmonicTimbre(
                attack=0.01,  # Attack phase duration in seconds
                decay=0.05,  # Decay phase duration in seconds
                sustain_amp=0.8,  # Amplitude fraction to decay to relative to max amplitude --> 0.4, 0.7
                release=0.45,  # Release phase duration in seconds
                num_harmonics=10,  # Actual number of partial harmonics to use
                roll_off=14,  # Roll-off in units of dB/octave,
            )
        )
    elif selected_timbre == "complex_short_ISI_short":  # MacDermott match
        note_duration_tonejs = 0.5
        note_silence_tonejs = 0
        TIMBRE = dict(
            default=HarmonicTimbre(
                attack=0.01,  # Attack phase duration in seconds
                decay=0.05,  # Decay phase duration in seconds
                sustain_amp=0.8,  # Amplitude fraction to decay to relative to max amplitude --> 0.4, 0.7
                release=0.45,  # Release phase duration in seconds
                num_harmonics=10,  # Actual number of partial harmonics to use
                roll_off=14,  # Roll-off in units of dB/octave,
            )
        )
    return note_duration_tonejs, note_silence_tonejs, TIMBRE


def estimate_time_per_trial(
    # estiamte time for trials: melody and singing duration
        pitch_duration,
        num_pitches,
        time_after_singing
):
    melody_duration = pitch_duration * num_pitches
    singing_duration = melody_duration + time_after_singing
    return melody_duration, singing_duration


class ToneJSVolumeTest(Module):
    """
    This is a volume calibration test to be used when implementing experiments with tonejs.

    Parameters
    ----------
    label : string, optional
        The label for the REPPVolumeCalibration test, default: "tonejs_volume_calibration_test".

    time_estimate_per_trial : float, optional
        The time estimate in seconds per trial, default: 5.0.

    min_time_before_submitting : float, optional
        Minimum time to wait (in seconds) while the music plays and the participant cannot submit a response, default: 5.0.

    number_intervals : float, optional
        Number of intervals to generate the melody, default: 9.0.

    max_interval_size : float, optional
        Maximum interval size to between pitches, default: 8.5.

    timbre : float, optional
        Timbre to play the melody, from js_synth, default: piano

    """

    def __init__(
        self,
        label="tonejs_volume_calibration_test",
        time_estimate_per_trial: float = 5.0,
        min_time_before_submitting: float = 5.0,
        number_intervals: int = 9,
        max_interval_size: float = 8.5,
        timbre: dict = InstrumentTimbre("piano"),
        note_duration: int = 0.8,
        note_silence: int = 0.0,
    ):
        self.label = label
        self.elts = PageMaker(
            lambda participant: ModularPage(
                "tone_js_volume_test",
                JSSynth(
                    Markup(
                        """
                            <h3>Volume calibration</h3>
                            <hr>
                            Set the volume in your laptop to a level in which you can hear each note properly.
                            <hr>
                            """
                    ),
                    sequence=[
                        Note(x)
                        for x in melodies.convert_interval_sequence_to_absolute_pitches(
                            intervals=melodies.sample_interval_sequence(
                                n_int=number_intervals,
                                max_interval_size=max_interval_size,
                                max_melody_pitch_range=config["max_melody_pitch_range"],
                                discrete=config["discrete"],
                                reference_mode="first_note",
                            ),
                            reference_pitch=melodies.sample_reference_pitch(
                                params.roving_mean["default"],
                                params.roving_width["default"],
                            ),
                            reference_mode="first_note",
                        )
                    ],
                    timbre=timbre,
                    default_duration=note_duration,
                    default_silence=note_silence,
                ),
                time_estimate=time_estimate_per_trial,
                events={
                    "restartMelody": Event(
                        is_triggered_by="promptEnd",
                        delay=1.0,
                        js="psynet.trial.restart()"
                    ),
                    "submitEnable": Event(is_triggered_by="trialStart", delay=min_time_before_submitting)
                }
            ),
            time_estimate=time_estimate_per_trial,
        )
        super().__init__(self.label, self.elts)


# Singing Calibration
class SingingTestControl(AudioMeterControl):
    decay = {"display": 0.1, "high": 0.1, "low": 0.1}
    threshold = {"high": -5, "low": -22}  # this is matching the db_threshold in singing_extract
    grace = {"high": 0.2, "low": 1.5}
    warn_on_clip = False
    msg_duration = {"high": 0.25, "low": 0.25}


class SingingCalibration(Module):
    """
    This is a singing calibration test to be used when implementing singing experiments with sing4me.

    Parameters
    ----------
    label : string, optional
        The label for the REPPVolumeCalibration test, default: "singing_calibration_test".

    time_estimate_per_trial : float, optional
        The time estimate in seconds per trial, default: 5.0.

    min_time_before_submitting : float, optional
        Minimum time to wait (in seconds) while the music plays and the participant cannot submit a response, default: 5.0.

    """

    def __init__(
        self,
        label="singing_calibration_test",
        time_estimate_per_trial: float = 5.0,
        min_time_before_submitting: float = 5.0,
    ):
        self.label = label
        self.elts = ModularPage(
            "singing_calibration",
            Markup(
                """
                <h3>Singing Calibration</h3>
                <hr>
                Please speak into your microphone and check that the sound is registered
                properly. If the sound is too quiet, try moving your microphone
                closer or increasing the input volume on your computer.
                <br><br>
                <b><b>Attention:</b></b> If the sound is not registered properly,
                you won't be able to continue with the experiment.
                <br><br>
                <hr>
                """
            ),
            SingingTestControl(calibrate=False),
            events={
                "submitEnable": Event(is_triggered_by="trialStart", delay=min_time_before_submitting)
                },
            time_estimate=time_estimate_per_trial
        )
        super().__init__(self.label, self.elts)


# singing recording test
class CustomStimulusVersionSpec(StimulusVersionSpec):
    has_media = True
    media_ext = ".wav"

    @classmethod
    def generate_media(cls, definition, output_path):
        shutil.copyfile(definition["local_file"], output_path)


silent_prompt = [
    StimulusSpec(
        definition={},
        version_specs=[
            CustomStimulusVersionSpec(
                definition={
                    "local_file": os.path.join("static", "audio", "silence_1s.wav"),
                }
            )
        ],
        phase="recording_test"
    )
]

stimulus_set_silence = StimulusSet(
    "recording_test", silent_prompt, version="stimuli", s3_bucket="iterated-singing-demo"
)


class SingingRecordingTrial(AudioRecordTrial, StaticTrial):
    __mapper_args__ = {"polymorphic_identity": "singing_recording_trial"}

    time_estimate = 8

    def show_trial(self, experiment, participant):
        return ModularPage(
            "singing_recording_trial",
            AudioPrompt(
                self.definition["url_audio"],
                Markup(
                    f"""
                    <h3>Recording Test</h3>
                    Sing 2 notes to the syllable '{SYLLABLE}'.<br> 
                    <i>leave a silent gap between the notes</i>
                    """
                ),
            ),
            AudioRecordControl(
                duration=self.definition["duration_rec_sec"],
                s3_bucket="iterated-singing",
                show_meter=True,
                public_read=True,
                controls=False,
                auto_advance=False,
            ),
            time_estimate=4,
            progress_display=ProgressDisplay(
                stages=[
                    ProgressStage(self.definition["duration_rec_sec"], "Recording.. Sing two notes!", "red"),
                ],
            ),
        )

    def analyze_recording(self, audio_file: str, output_plot: str):
        raw = sing.analyze(
            audio_file,
            config,
            plot_options=sing.PlotOptions(
                save=True, path=output_plot, format="png"
            ),
        )
        raw = [
            {key: melodies.as_native_type(value) for key, value in x.items()} for x in raw
        ]
        sung_pitches = [x["median_f0"] for x in raw]

        if len(sung_pitches) >= self.definition["min_num_notes_to_detect"]:
            correct_num_notes = True
            failed = False
        else:
            correct_num_notes =False
            failed = True

        return {
            "failed": failed,
            "correct_num_notes": correct_num_notes,
            "num_sung_pitches": len(sung_pitches),
        }


class SingingRecordingTrialFeedback(SingingRecordingTrial):
    __mapper_args__ = {"polymorphic_identity": "custom_trial_practice_feedback"}

    wait_for_feedback = True

    def gives_feedback(self, experiment, participant):
        return True

    def show_feedback(self, experiment, participant):
        output_analysis =json.loads(self.details["analysis"])
        if output_analysis["num_sung_pitches"] == 2:
            return InfoPage(
                Markup(
                    f"""
                    <h3>Excellent!</h3>
                    <hr>
                    <br><br>
                    We could detect {output_analysis["num_sung_pitches"]} notes in your recording.
                    <hr>
                    <img style="width:25%" src="/static/images/happy.png"  alt="happy">
                    """
                ),
                time_estimate=5
            )
        if output_analysis["num_sung_pitches"] == 0:
            return InfoPage(
                Markup(
                    f"""
                   <h3>Bad...</h3>
                   <hr>
                   We could not detect any note in your recording.
                   <br><br>
                   Please try to do one or more of the following:
                    <ol><li>Sing each note clearly using the syllable '{SYLLABLE}'.</li>
                        <li>Make sure you computer microphone is working and you are in a quiet environment.</li>
                        <li>Leave a silent gap between the notes.</li>
                        <li>Sing each note for about 1 second.</li>
                    </ol>
                    <hr>
                   <img style="width:25%" src="/static/images/sad.png"  alt="sad">
                   """
                ),
                time_estimate=5
            )
        else:
            return InfoPage(
                Markup(
                    f"""
                    <h3>Ok, but you can do better..</h3>
                    <hr>
                    We detected {output_analysis["num_sung_pitches"]} notes in your recording, but we asked
                    you to <b><b>sing 2 notes</b></b>.
                    <br><br>
                       Please try to do one or more of the following:
                        <ol><li>Sing each note clearly using the syllable '{SYLLABLE}'.</li>
                            <li>Make sure you computer microphone is working and you are in a quiet environment.</li>
                            <li>Leave a silent gap between the notes.</li>
                            <li>Sing each note for about 1 second.</li>
                        </ol>
                        <hr>
                    <img style="width:25%" src="/static/images/neutral.png"  alt="neutral">
                    """
                ),
                time_estimate=5
            )


class SingingRecordingTest(Module):
    """
    This pre-screening test is designed to quickly determine whether participants
    are able to provide valid singing data. The task is also efficient in determining whether
    participants are following the instruction and use hardware
    and software that meet the technical requirements for singing experiments.
    To make the most out of it, the test should be used at the
    beginning of the experiment, after providing general instructions.
    By default, we start with warming trials, where we first play the recording back to the participant
    and then give them feedback. Finally, we use a third trial to exclude participants based on their performance.

    Parameters
    ----------

    label : string, optional
        The label for the test, default: "singing_recording_test".

    performance_threshold : int, optional
        The performance threshold, default: 1.

    duration_rec_sec : float, optional
        Time of the recording, default: 6 sec.

    min_num_notes_to_detect : float, optional
        Mininum number of sung notes to detect in order to pass the trial, default: 1.

    num_repeat_trials : float, optional
        Number of trials to repeat in the main trial maker, default: 0.

    """

    def __init__(
        self,
        label="singing_recording_test",
        performance_threshold: int = 1,
        duration_rec_sec: float = 6,
        min_num_notes_to_detect: int = 1,
        num_repeat_trials: int = 0,
    ):
        self.label = label
        self.elts = join(
            self.familiarization_phase(),
            self.trial_maker_feedback(
                duration_rec_sec,
                min_num_notes_to_detect,
            ),
            self.instructions_test(),
            self.trial_maker_test(
                performance_threshold,
                duration_rec_sec,
                min_num_notes_to_detect,
                num_repeat_trials,
            ),
        )
        super().__init__(self.label, self.elts)

    def familiarization_phase(self):
        return join(
            InfoPage(
                Markup(
                    f"""
                    <h3>Recording Example</h3>
                    <hr>
                    Now we will test if we can record your voice with your microphone. 
                    <br><br>
                    When ready, go to the next page and <b><b>sing 2 notes</b></b> using the syllable '{SYLLABLE}'.<br>
                    (separate each note with a silence). 
                    <hr>
                    """
                ),
                time_estimate=5,
            ),
            ModularPage(
                "singing_record_example",
                Markup(
                    f"""
                    <h3>Recording Example</h3>
                    Sing 2 notes to the syllable '{SYLLABLE}'.<br> 
                    <i>leave a silent gap between the notes</i>
                    """),
                AudioRecordControl(
                    duration=5.0,
                    s3_bucket="iterated-singing-demo",
                    show_meter=True,
                    public_read=True,
                    controls=False,
                    auto_advance=False,
                ),
                time_estimate=5,
                progress_display=ProgressDisplay(
                    stages=[
                        ProgressStage(5, "Recording.. Sing 2 notes!", "red"),
                    ],
                ),
            ),
            PageMaker(
                lambda participant: ModularPage(
                    "playback",
                    AudioPrompt(
                        participant.answer["url"],
                        Markup("""
                        <h3>Can you hear your recording?</h3>
                        <hr>
                        If you do not hear your recording, please make sure
                        to use a working microphone so we can record your voice and continue with the experiment. 
                        <hr>
                        """),
                    ),
                ),
                time_estimate=5,
            ),
            InfoPage(
                Markup(
                    f"""
                    <h3>Recording Example</h3>
                    <hr>
                    Now we will analyse your recording and provide feedback accordingly. 
                    <br><br>
                    In the next page, <b><b>sing 2 notes</b></b> using the syllable '{SYLLABLE}'.
                    <hr>
                    """
                ),
                time_estimate=5,
            ),
        )

    def trial_maker_feedback(
            self,
            duration_rec_sec: int,
            min_num_notes_to_detect: int,
    ):
        class RecordingTestTrialMakerFeedback(StaticTrialMaker):
            performance_check_type = "performance"
            performance_check_threshold = 0
            give_end_feedback_passed = False

        return RecordingTestTrialMakerFeedback(
            id_="singing_recording_feedback",
            trial_class=SingingRecordingTrialFeedback,
            phase="recording_test",
            stimulus_set = self.get_stimulus_set(duration_rec_sec, min_num_notes_to_detect),
            fail_trials_on_premature_exit = False,
            fail_trials_on_participant_performance_check = False,
            check_performance_at_end=True
        )

    def instructions_test(self):
        return InfoPage(
            Markup(
                f"""
                <h3>Test</h3>
                <hr>
                <b><b>Be careful:</b></b> This is a recording test!<br><br>
                In the next page, we will ask you again to <b><b>sing 2 notes</b></b> to the
                syllable '{SYLLABLE}'. <br><br>
                We will analyse your recording and test if we can detect your singing.
                <hr>
                """
                ),
                time_estimate=5,
        )

    def trial_maker_test(
            self,
            performance_threshold: int,
            duration_rec_sec: int,
            min_num_notes_to_detect: int,
            num_repeat_trials: int,
    ):
        class RecordingTestTrialMakerTest(StaticTrialMaker):
            performance_check_type = "performance"
            performance_check_threshold = performance_threshold
            give_end_feedback_passed = False

        return RecordingTestTrialMakerTest(
            id_="singing_recording_trialmaker",
            trial_class=SingingRecordingTrial,
            phase="recording_test",
            stimulus_set = self.get_stimulus_set(duration_rec_sec, min_num_notes_to_detect),
            target_num_participants=0,
            recruit_mode="num_participants",
            fail_trials_on_premature_exit = False,
            fail_trials_on_participant_performance_check = False,
            check_performance_at_end=True,
            num_repeat_trials=num_repeat_trials
        )

    def get_stimulus_set(self, duration_rec_sec: int, min_num_notes_to_detect: int):
        return StimulusSet(
            "silence_wav",
            [
                StimulusSpec(
                    definition={
                        "url_audio": "https://s3.amazonaws.com/repp-materials/silence_1s.wav",
                        "duration_rec_sec": duration_rec_sec,
                        "min_num_notes_to_detect": min_num_notes_to_detect,
                    },
                    phase="recording_test",
                )
            ],
        )


# singing performance test
class SingingPerformanceTrial(AudioRecordTrial, StaticTrial):
    __mapper_args__ = {"polymorphic_identity": "singing_performance_trial"}

    time_estimate = 8

    def show_trial(self, experiment, participant):
        current_trial = self.position + 1
        if "feedback" in self.trial_maker_id:
            total_num_trials = self.definition["num_trials_feedback"]
        else:
            total_num_trials = self.definition["num_trials_test"] + self.definition["num_repeat_trials"]
        show_current_trial = f'<br><br>Trial number {current_trial} out of {total_num_trials} trials.'

        return ModularPage(
            "singing_performance_trial",
            JSSynth(
                Markup(
                    f"""
                    <h3>Imitate the melody</h3>
                    This melody has two notes: <b><b>Sing each note back to the syllable '{SYLLABLE}'.</b></b><br>
                    <i>leave a silent gap between the notes</i>
                    <br><br>
                    {show_current_trial}
                    """
                ),
                sequence=[Note(x) for x in self.definition["target_pitches"]],
                timbre=jsonpickle.decode(self.definition["timbre"]),
                default_duration=self.definition["note_duration"],
                default_silence=self.definition["note_silence"],
            ),
            AudioRecordControl(
                duration=3.5,
                s3_bucket="iterated-singing",
                public_read=True,
                show_meter=False,
                controls=False,
                auto_advance=False,
            ),
            time_estimate=8,
            events={
                "promptStart": Event(is_triggered_by="trialStart"),
                "recordStart": Event(is_triggered_by="promptEnd", delay=0.25)
            },
             progress_display=ProgressDisplay(
                    stages=[
                        ProgressStage(2.5, "Listen to the melody...", "orange"),
                        ProgressStage(3, "Recording... START SINGING!", "red"),
                        ProgressStage(0.5, "Finished recording.",
                                      "green",
                                      persistent=True
                                      ),
                        ],
                    ),
                )

    def analyze_recording(self, audio_file: str, output_plot: str):
        raw = sing.analyze(
            audio_file,
            config,
            target_pitches=self.definition["target_pitches"],
            plot_options=sing.PlotOptions(
                save=True, path=output_plot, format="png"
            ),
        )
        raw = [
            {key: melodies.as_native_type(value) for key, value in x.items()} for x in raw
        ]
        sung_pitches = [x["median_f0"] for x in raw]
        sung_intervals = melodies.convert_absolute_pitches_to_interval_sequence(
            sung_pitches,
            "first_note"
        )
        target_intervals = melodies.convert_absolute_pitches_to_interval_sequence(
            self.definition["target_pitches"],
            "first_note"
        )
        stats = sing.compute_stats(
            sung_pitches,
            self.definition["target_pitches"],
            sung_intervals,
            target_intervals
        )

        # failing criteria
        correct_num_notes =  stats["num_sung_pitches"] == stats["num_target_pitches"]
        max_interval_error_ok = stats["max_abs_interval_error"] < 3
        direction_accuracy_ok = stats["direction_accuracy"] == 100

        failed_options = [
            correct_num_notes,
            max_interval_error_ok,
            direction_accuracy_ok
        ]
        reasons = [
            "Wrong number of sung notes",
            "max interval error is larger than 3",
            "direction accuyracy is wrong"
        ]
        if False in failed_options:
            failed = True
            index = failed_options.index(False)
            reason = reasons[index]
        else:
            failed = False
            reason = "All good"

        return {
            "failed": failed,
            "reason": reason,
            "target_pitches": self.definition["target_pitches"],
            "target_intervals": target_intervals,
            "sung_pitches": sung_pitches,
            "sung_intervals": sung_intervals,
            "raw": raw,
            # "stats": stats,
            "mean_pitch_diffs": stats["mean_pitch_diffs"],
            "max_abs_pitch_error": stats["max_abs_pitch_error"],
            "mean_interval_diff": stats["mean_interval_diff"],
            "max_abs_interval_error": stats["max_abs_interval_error"],
            "direction_accuracy": stats["direction_accuracy"],
        }


class SingingPerformanceFeedback(SingingPerformanceTrial):
    __mapper_args__ = {"polymorphic_identity": "singing_performance_feedback"}

    wait_for_feedback = True

    def gives_feedback(self, experiment, participant):
        return True

    def show_feedback(self, experiment, participant):
        output_analysis =json.loads(self.details["analysis"])
        num_sung_pitches = output_analysis["num_sung_pitches"]
        if num_sung_pitches == 2:
            return InfoPage(
                Markup(
                    f"""
                    <h3>Excellent!</h3>
                    <hr>
                    <br><br>
                    We detected {num_sung_pitches} notes in your recording.
                    <hr>
                    <img style="width:25%" src="/static/images/happy.png"  alt="happy">
                    """
                ),
                time_estimate=5
            )
        if num_sung_pitches == 0:
            return InfoPage(
                Markup(
                    f"""
                   <h3>Bad...</h3>
                   <hr>
                   We could not detect any note in your recording.<br><br> 
                   Please try to do one or more of the following:
                    <ol><li>Sing each note clearly using the syllable '{SYLLABLE}'.</li>
                        <li>Make sure you computer microphone is working and you are in a quiet environment.</li>
                        <li>Leave a silent gap between the notes.</li>
                        <li>Sing each note for about 1 second.</li>
                    </ol>
                    <hr>
                    <img style="width:25%" src="/static/images/sad.png"  alt="sad">
                    """
                    ),
                time_estimate=5
            )
        else:
            return InfoPage(
                Markup(
                    f"""
                    <h3>Ok, but you can do better...</h3>
                    <hr>
                    We detected {num_sung_pitches} notes in your recording, but we asked
                    you to <b><b>sing 2 notes</b></b>.<br><br> 
                    Please try to do one or more of the following:
                    <ol><li>Sing each note clearly using the syllable '{SYLLABLE}'.</li>
                        <li>Make sure you computer microphone is working and you are in a quiet environment.</li>
                        <li>Leave a silent gap between the notes.</li>
                        <li>Sing each note for about 1 second.</li>
                    </ol>
                    <hr>
                    <img style="width:25%" src="/static/images/neutral.png"  alt="neutral">
                    """
                ),
                time_estimate=5
            )

    def analyze_recording(self, audio_file: str, output_plot: str):
        raw = sing.analyze(
            audio_file,
            config,
            plot_options=sing.PlotOptions(
                save=True, path=output_plot, format="png"
            ),
        )
        raw = [
            {key: melodies.as_native_type(value) for key, value in x.items()} for x in raw
        ]
        sung_pitches = [x["median_f0"] for x in raw]

        if len(sung_pitches) >= 1:
            correct_num_notes = True
            failed = False
        else:
            correct_num_notes =False
            failed = True

        return {
            "failed": failed,
            "correct_num_notes": correct_num_notes,
            "num_sung_pitches": len(sung_pitches),
        }


class SingingPerformanceTest(Module):
    """
    This is a singing performance test to determine participants' singing accuracy and also
     exclude particpants who do not follow the instructions or do not provide valid singing data.

    Parameters
    ----------

    label : string, optional
        The label for the singing test, default: "singing_test".

    time_estimate_per_trial : float, optional
        The time estimate in seconds per trial, default: 8.0.

    performance_threshold : int, optional
        The performance threshold to exclude participants, default: 3.

    num_trials_feedback : int, optional
        The total number of trials to display in the feedback phase, default: 2.

    num_trials_test : int, optional
        The total number of trials to display in the test phase, default: 6.

    num_repeat_trials : int, optional
        Number of trials to repeat, default: 0.

    timbre : dict
        Dicitonary specifying the timbre, default: "default".

    trial_maker_id : str
        The trial maker id, default: "singing_performance_trialmaker".


    """

    def __init__(
        self,
        label="singing_perforamnce_test",
        performance_threshold: int = 5,
        num_trials_feedback: int = 2,
        num_trials_test: int = 10,
        num_repeat_trials: int = 0,
        timbre: dict = InstrumentTimbre("piano"),
        note_duration: float = 0.8,
        note_silence: float = 0.0,
        roving_width: float = 2.5,
        roving_mean_low: float = 49,
        roving_mean_high: float = 61,
        num_target_participants: int = 0,
        trial_maker_id: str = "singing_performance_trialmaker",
    ):
        self.label = label
        self.elts = join(
            self.familiarization_phase(),
            self.trial_maker_feedback(
                num_trials_feedback,
                num_trials_test,
                num_repeat_trials,
                timbre,
                note_duration,
                note_silence,
                roving_width,
                trial_maker_id
            ),
            self.instructions_test(num_trials_test, num_repeat_trials),
            self.trial_maker(
                num_trials_feedback,
                num_trials_test,
                num_repeat_trials,
                performance_threshold,
                timbre,
                note_duration,
                note_silence,
                roving_width,
                roving_mean_low,
                roving_mean_high,
                num_target_participants,
                trial_maker_id,
            ),
        )
        super().__init__(self.label, self.elts)

    def familiarization_phase(self):
        return join(
            InfoPage(
                Markup(
                    f"""
                    <h3>Recording Example</h3>
                    <hr>
                    First we will test if you can record your voice with the computer microphone. 
                    <br><br>
                    When ready, go to the next page and <b><b>sing 2 notes</b></b> using the syllable '{SYLLABLE}'.<br>
                    (separate each note with a silence). 
                    <hr>
                    """
                ),
                time_estimate=5,
            ),
            ModularPage(
                "singing_record_example",
                Markup(
                    f"""
                    <h3>Recording Example</h3>
                    Sing 2 notes to the syllable {SYLLABLE}.<br> 
                    <i>leave a silent gap between the notes</i>
                    """),
                AudioRecordControl(
                    duration=5.0,
                    s3_bucket="iterated-singing-demo",
                    show_meter=True,
                    public_read=True,
                    controls=False,
                    auto_advance=False,
                ),
                time_estimate=5,
                progress_display=ProgressDisplay(
                    stages=[
                        ProgressStage(5, "Recording.. Sing 2 notes!", "red"),
                    ],
                ),
            ),
            PageMaker(
                lambda participant: ModularPage(
                    "playback",
                    AudioPrompt(
                        participant.answer["url"],
                        Markup("""
                        <h3>Can you hear your recording?</h3>
                        <hr>
                        If you do not hear your recording, please make sure
                        to use a working microphone so we can record your voice and continue with the experiment. 
                        <hr>
                        """),
                    ),
                ),
                time_estimate=5,
            ),
            InfoPage(
                Markup(
                    f"""
                    <h3>Singing practice</h3>
                    <hr>
                    In each trial, you will hear a melody with 2 notes:<br>
                    <b><b>Your goal is to sing each note back as accurately as possible.</b></b><br>
                    <i>Important:</i> Use the syllable {SYLLABLE} to sing each note and leave a silent gap between notes.
                    <br><br>
                    We will analyse your recording and provide feedback.
                    <hr>
                    When ready, click <b><b>next</b></b> to start singing.
                    """
                ),
                time_estimate=5,
            ),
        )

    def trial_maker_feedback(
            self,
            num_trials_feedback: int,
            num_trials_test: int,
            num_repeat_trials: int,
            timbre: dict,
            note_duration: int,
            note_silence: float,
            roving_width: float,
            trial_maker_id: str,
    ):
        class SingingPerformanceFeedbackTrialMaker(StaticTrialMaker):
            performance_check_type = "performance"
            performance_check_threshold = 0.6
            give_end_feedback_passed = False

        return SingingPerformanceFeedbackTrialMaker(
            id_=f"{trial_maker_id}_feedback",
            trial_class=self.trial_class_feedback,
            phase="singing_performance",
            stimulus_set=self.get_stimulus_set(
                num_trials_test, num_repeat_trials, num_trials_feedback,
                timbre, note_duration, note_silence,
                roving_width
            ),
            max_trials_per_block=num_trials_feedback,
            fail_trials_on_premature_exit = False,
            fail_trials_on_participant_performance_check = False,
            check_performance_at_end=True
        )

    trial_class_feedback = SingingPerformanceFeedback

    def instructions_test(self, num_trials_test, num_repeat_trials):
        return InfoPage(
            Markup(
                f"""
                <h3>Singing Test</h3>
                <hr>
                Now we will test your singing performance in a total of {num_trials_test + num_repeat_trials} trials.
                <br><br> 
                Your goal is to listen to each melody and sing it back to the syllable {SYLLABLE}.
                <hr>
                When ready, click <b><b>next</b></b> to start the singing test.
                """
            ),
            time_estimate=5,
        )

    def trial_maker(
            self,
            num_trials_feedback: int,
            num_trials_test: int,
            num_repeat_trials: int,
            performance_threshold: int,
            timbre: dict,
            note_duration: int,
            note_silence: int,
            roving_width: int,
            roving_mean_low: int,
            roving_mean_high: int,
            num_target_participants: int,
            trial_maker_id: str,
    ):
        class SingingPerformanceTrialMaker(StaticTrialMaker):
            give_end_feedback_passed = True
            performance_check_type = "performance"

            def performance_check(self, experiment, participant, participant_trials):
                """Should return a tuple (score: float, passed: bool)"""
                score = 0
                # list_max_abs_pitch_diff = []
                list_max_abs_interval_error = []
                list_direction_accuracy = []
                list_sung_pitches = []

                for trial in participant_trials:
                    # list_max_abs_pitch_diff.append(trial.analysis["max_abs_pitch_error"])
                    list_max_abs_interval_error.append(trial.analysis["max_abs_interval_error"])
                    list_direction_accuracy.append(trial.analysis["direction_accuracy"])
                    list_sung_pitches.append(trial.analysis["sung_pitches"])
                    if trial.analysis["failed"] == False:
                        score += 1
                passed = score >= performance_threshold

                # store variables in particaipnts table
                participant.var.set("singing_performance", score)
                participant.var.set("list_max_abs_interval_error", list_max_abs_interval_error)
                participant.var.set("list_direction_accuracy", list_direction_accuracy)

                # determine register
                flat_list_sung_pitches = [x for xs in list_sung_pitches for x in xs]
                median_pitch = np.median(flat_list_sung_pitches)
                distance_to_low_register = abs(median_pitch - roving_mean_low)
                distance_to_high_register = abs(median_pitch - roving_mean_high)

                if distance_to_low_register < distance_to_high_register:
                    predicted_register = "low"
                elif distance_to_low_register > distance_to_high_register:
                    predicted_register = "high"
                else:
                    predicted_register = "undefined"

                participant.var.set("sung_median_pitch", median_pitch)
                participant.var.set("distance_to_low_register", distance_to_low_register)
                participant.var.set("distance_to_high_register", distance_to_high_register)
                participant.var.set("predicted_register", predicted_register)

                
                return {"score": score, "passed": passed}

            def get_end_feedback_passed_page(self, score):
                score_to_display = "NA" if score is None else f"{((score / (num_trials_test + num_repeat_trials)) * 100):.0f}"

                return InfoPage(
                    Markup(
                        f"""
                        <h4>Congratulations, you passed the singing test.</h4>
                        Your perforamnce score was <strong>{score_to_display}&#37;</strong>.
                        """
                    ),
                    time_estimate=3,
                )

        return SingingPerformanceTrialMaker(
            id_=trial_maker_id,
            trial_class=self.trial_class,
            phase="singing_performance",
            stimulus_set=self.get_stimulus_set(
                num_trials_test, num_repeat_trials, num_trials_feedback,
                timbre, note_duration, note_silence,
                roving_width
            ),
            max_trials_per_block=num_trials_test,
            allow_repeated_stimuli=True,
            check_performance_at_end=True,
            check_performance_every_trial=False,
            target_num_participants=num_target_participants,
            target_num_trials_per_stimulus=None,
            recruit_mode="num_participants",
            num_repeat_trials=num_repeat_trials
        )

    trial_class = SingingPerformanceTrial

    def get_stimulus_set(self,
                         num_trials_test: int,
                         num_repeat_trials: int,
                         num_trials_feedback: int,
                         timbre: dict,
                         note_duration: int,
                         note_silence: int,
                         roving_width: int
                         ):
        return StimulusSet(
            "intervals",
            [
                StimulusSpec(
                    definition={
                        "num_trials_test": num_trials_test,
                        "num_repeat_trials": num_repeat_trials,
                        "num_trials_feedback": num_trials_feedback,
                        "interval": interval,
                        "note_duration": note_duration,
                        "note_silence": note_silence,
                        "target_pitches": melodies.convert_interval_sequence_to_absolute_pitches(
                            intervals=[interval],
                            reference_pitch=melodies.sample_reference_pitch(
                                params.roving_mean[register],
                                roving_width
                            ),
                            reference_mode="previous_note",
                        ),
                        "timbre": jsonpickle.encode(timbre),
                    },
                    phase="singing_performance",
                )
                for interval in [-1.3, -2.6, 0.0, 1.3, 2.6]
                for register in ["low", "high"]
            ],
        )


# gender split
gender_choices = ["male", "female", "other / prefer not to say"]
gender_markup = [Markup(f"<DIV style=text-align:left;>{choice}</DIV>") for choice in gender_choices]

GenderSplit = Module(
    "gender_split",
    ModularPage(
        "Gender",
        Prompt(
            Markup(f"<h4>Please select the gender with which you identify.</h4>")
        ),
        PushButtonControl(
            labels=gender_markup,
            choices=gender_choices,
            arrange_vertically=True),
        time_estimate=5,
    ),
    conditional(
        "set_participant_gender",
        lambda experiment, participant: (Markup(participant.answer).striptags() == "male") or (
                Markup(participant.answer).striptags() == "female"),
        CodeBlock(
            lambda experiment, participant: participant.var.set("gender", Markup(participant.answer).striptags())),
        join(
            CodeBlock(lambda experiment, participant: logger.info(
                "********** ANSWER: {} **********".format(Markup(participant.answer).striptags()))),
            ModularPage(
                "voice_range",
                Prompt(
                    Markup(f"<h4> Select which voice range best matches you.</h4>")
                ),
                PushButtonControl(
                    labels=["Low Voice Range", "High Voice Range"],
                    choices=["Low Voice Range", "High Voice Range"],
                    arrange_vertically=True),
                time_estimate=5,
            ),
            conditional(
                label="assign-to-high-or-low",
                condition=lambda experiment, participant: participant.answer == "Low Voice Range",
                logic_if_true=CodeBlock(lambda experiment, participant: participant.var.set("gender", "male")),
                logic_if_false=CodeBlock(lambda experiment, participant: participant.var.set("gender", "female")),
                fix_time_credit=False
            )
        ),
        fix_time_credit=True
    ),
    # PageMaker(
    #     lambda participant: InfoPage(f"Your gender is: {participant.var.gender}"),
    #     time_estimate=5,
    # ),
)


# self-reported register
SelfReportedRegister = Module(
    "selfreported_register",
    ModularPage(
        "selfreported_register",
        Prompt(
            Markup(
                f"""
                <h4>What is your vocal range?</h4>
                Everyone is different in the range of notes they can sing.
                This is known as <b><b>vocal range</b></b>.
                <br> For example, male voice is generally lower than female voice.
                <br><br>
                <b><b>Please select which vocal range best matches you:</b></b>
                """)
        ),
        PushButtonControl(
            labels=["Low Vocal Range", "High Vocal Range"],
            choices=["Low Vocal Range", "High Vocal Range"],
            arrange_vertically=True),
            time_estimate=5,
            ),
    conditional(
        label="assign-to-high-or-low",
        condition=lambda experiment, participant: participant.answer == "Low Vocal Range",
        logic_if_true=CodeBlock(lambda experiment, participant: participant.var.set("selfreported_register", "low")),
        logic_if_false=CodeBlock(lambda experiment, participant: participant.var.set("selfreported_register", "high")),
        fix_time_credit=False
        )
)


