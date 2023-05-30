from setuptools import setup

setup(
    name = 'pdf_tonality',
    version = '0.1',
    author = 'Michael W. Weiss',
    email = 'michael.weiss@umontreal.ca',
    packages = ['pdf_tonality'],
    python_requires='>3.6',
    install_requires = [
        'pandas>=1.3.4',
        'numpy>=1.21.3',
        'matplotlib>=3.4.3',
        'pytest>=6.2.5',
        'scipy>=1.7.1',
        'IPython>=7.28.0',
        'rich==10.12.0',
    ],
)