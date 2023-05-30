####################################################################################################
# File:     setup.py
# Purpose:  Setup of the package.
####################################################################################################
import os

from setuptools import find_packages

with open(os.path.join("sing4me", "VERSION")) as version_file:
    version = version_file.read().strip()

from setuptools import setup, Command


class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')


__MODULES__ = [
    "sing4me.detect_peaks",
    "sing4me.singing_extract",
    "sing4me.utils",
    "sing_experiments.params",
    "sing_experiments.primes",
    "sing_experiments.questionnaire",
    "sing_experiments.resources",
    "sing_experiments.melodies",

]

setup(
    name="sing4me",
    version=version,
    py_modules=__MODULES__,
    author="Manuel Anglada-Tort, Peter Harrison, Nori Jacoby",
    author_email="manel.anglada.tort@gmail.com",
    description="Python package for singing experiments",
    long_description="Python package for singing experiments",
    long_description_content_type="text/markdown",
    url="https://gitlab.com/computational-audition-lab/sing4me",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    include_package_data=True,
    package_data={"sing4me": ["VERSION"]},
    cmdclass={'clean': CleanCommand},
    install_requires=[
        "click>=7.1.2",
        "flask>=1.0.2",
        "docopt>=0.6.2",
        "matplotlib>=3.3.1",
        "numpy>=1.18.4",
        "scipy>=1.5.3",
        "praat-parselmouth>=0.4.0",
        "textgrid",
        # requirements for iterated_singing_demo
        # "sounddevice>=0.4.1",
        # "soundfile>=0.10.3.post1",
        # "jupyter>=1.0.0"
    ],
    entry_points={
        "console_scripts": [
            "sing4me = sing4me.cli:sing4me"
        ]
    }
)

# python3.7 setup.py sdist bdist_wheel
