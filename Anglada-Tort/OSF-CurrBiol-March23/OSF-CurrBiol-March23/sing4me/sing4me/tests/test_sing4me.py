import os

from sing4me.singing_extract import (
	read_file_and_detect,
	analyze,
	PlotOptions,
	generate_sine_tones,
	compute_stats
)
from sing4me.utils import (
	freq2midi,
	ms_to_samples
)
from sing_experiments import params
singing_config = params.singing_2intervals

here = os.path.abspath(os.path.dirname(__file__))


def generate_pilot_analysis_suite(audio_dir, list_file= []):
	import json
	for fname in os.listdir(audio_dir):
		if list_file:
			if not (fname in list_file):
				continue
		if fname.endswith(".wav"):
			fp = os.path.join(audio_dir, fname)
			plot_fp = "/".join([audio_dir, fname[:-4]]) + ".png"
			analysis = analyze(
				filename=fp,
				extract_config=singing_config,
				plot_options=PlotOptions(save=True, path=plot_fp)
			)

			mjson = json.dumps(analysis)
			f = open(plot_fp.replace(".png",".txt"),"w")
			f.write(mjson)
			f.close()


generate_pilot_analysis_suite(
	#audio_dir=os.path.dirname(here) + "/tests/old_pilot_examples/"
	#audio_dir=os.path.dirname(here) + "/tests/bad_2int/"
	audio_dir=os.path.dirname(here) + "/tests/good_2int/"
	#audio_dir=os.path.dirname(here) + "/tests/bad_2int/"
	#audio_dir=os.path.dirname(here) + "/tests/bad_4int/"
	#audio_dir=os.path.dirname(here) + "/tests/bad_4int/" #
	#audio_dir=os.path.dirname(here) + "/tests/bad_india/"
	#audio_dir=os.path.dirname(here) + "/tests/good_4int/"
	#audio_dir=os.path.dirname(here) + "/tests/problem_cases/"
	
)

# then do automatically everything!
todos=['problem_cases','bad_4int','bad_2int','good_4int','good_2int','bad_india','old_pilot_examples'] 
for todo in todos:
	dname=os.path.dirname(here) + "/tests/" + todo + "/"
	generate_pilot_analysis_suite(audio_dir=dname)



# generate_pilot_analysis_suite(
# 	audio_dir="/Users/jacoby/cap/sing4me/sing4me/tests/old_pilot_examples/"
# )


# old

# def test_example1():
# 	ex2 = os.path.dirname(here) + "/tests/example_audio/example1.wav"
# 	analysis = analyze(ex2, extract_config=singing_config)
# 	assert len(analysis) == 2

# def test_example2():
# 	ex2 = os.path.dirname(here) + "/tests/example_audio/example2.wav"
# 	analysis = analyze(ex2, extract_config=singing_config)
# 	assert len(analysis) == 1

# 	target = [49]
# 	extracted_median_f0 = [round(x["median_f0"]) for x in analysis]
# 	assert target == extracted_median_f0

# def test_example4():
# 	ex4 = os.path.dirname(here) + "/tests/example_audio/example4.wav"
# 	analysis = analyze(ex4, extract_config=singing_config)
# 	targets = [55, 51] # hehe both are flat.
# 	extracted_median_f0s = [round(x["median_f0"]) for x in analysis]
# 	assert targets == extracted_median_f0s

# def test_example6():
# 	ex6 = os.path.dirname(here) + "/tests/example_audio/example6.wav"
# 	analysis = analyze(ex6, extract_config=singing_config)

# 	assert len(analysis) == 3
# 	targets = [52, 57, 57]
# 	extracted_median_f0s = [round(x["median_f0"]) for x in analysis]
# 	assert targets == extracted_median_f0s

# def test_example8():
# 	ex8 = os.path.dirname(here) + "/tests/example_audio/example8.wav"
# 	analysis = analyze(ex8, extract_config=singing_config)
# 	assert len(analysis) == 3

# def test_example9():
# 	ex8 = os.path.dirname(here) + "/tests/example_audio/example9.wav"
# 	analysis = analyze(ex8, extract_config=singing_config)
# 	assert len(analysis) == 3

# def test_empty_recording():
# 	empty = os.path.dirname(here) + "/tests/example_audio/empty.wav"
# 	analysis = analyze(empty, target_pitches=[60, 61], extract_config=singing_config)
# 	assert len(analysis) == 0

