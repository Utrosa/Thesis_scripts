# Singing Experiments

![Logo](logo.png)


[Computational Auditory Perception Group](https://www.aesthetics.mpg.de/en/research/research-group-computational-auditory-perception.html), \
Max Planck Institute for Empirical Aesthetics.


_sing4me_ is a python package for singing experiments in laboratory and online settings. In currently installs two modules:
1. sing4me: the code supporting singing extract with tests.
2. sing_experiments: resources to run singing experiments with sing4me, including methods to generate and transform melodies, pre-screening tests, paramaters, and questionaires. 

## Documentation

Sing4Me documentation: to build

Source code: https://gitlab.com/computational-audition/sing4me



## Installation (macOS)
_sing4me_ needs Python 3.x (tested using Python 3.7 and 3.9). It has only been tested in macOS.


### Virtual environment

1. Set up a virtual environment. For example:

```
pip3 install virtualenv

pip3 install virtualenvwrapper

export WORKON_HOME=$HOME/.virtualenvs

mkdir -p $WORKON_HOME

export VIRTUALENVWRAPPER_PYTHON=$(which python3)

source $(which virtualenvwrapper.sh)

mkvirtualenv sing4me --python $(which python3)
```

This will automatically activate the virtual environment, so you should see it between brackets in your terminal.

Note that for this example, we call the virtual environment __sing4me__, but you can give it any other name.
<br><br>
The following commands are optional, but they will be useful to easily activate your virtual environment from the terminal.

```
echo "export VIRTUALENVWRAPPER_PYTHON=$(which python3)" >> ~/.zshrc

echo "source $(which virtualenvwrapper.sh)" >> ~/.zshrc
```

In the future, you can activate this virtual environment by typing in the terminal:
```
workon sing4me
```


### Install Sing4Me

2. Download or clone `sing4me` from the repository. For example:
```
git clone git@gitlab.com:computational-audition-lab/sing4me.git
```

3. Then make sure you are working from that folder in the terminal:
```
cd sing4me
```
4. Install the requirements:
```
pip3 install -r requirements.txt
```
5. Install Sing4Me:
```
pip3 install -e .
```
The -e flag makes the Sing4Me code editable.

### Verify successful installation:
```
sing4me --version
```
You are now done with the installation and you can begin using Sing4Me.



## Tests

To run tests cd into sing4me and type `pytest tests`. The code will look for all filenames in the tests directory that begin with the word “test”. You can add tests by adding functions like this to the tests/test_sing4me.py (note the assertion):
```
def test_example8():
	ex8 = os.path.dirname(here) + "/tests/6c355b9f-6f21-43b4-b451-6641831eb529.wav"
	analysis = analyze(ex8, extract_config=singing_config)
	assert len(analysis) == 3
```

There is also a simple function called `generate_pilot_analysis_suite` that will generate plots and analysis for a set of recordings in a folder. Just uncomment this, add an appropriate directory, and run the file (or just import and run it in python, of course).
```
generate_pilot_analysis_suite(
	audio_dir=os.path.dirname(here) + "/tests/good_2int/"
)
```



## License

 MIT License

