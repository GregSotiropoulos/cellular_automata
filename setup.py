# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from pathlib import Path

here = Path(__file__).parent.resolve()

with open(here/'README.md') as fh:
    long_description = fh.read()

setup(
    name='cellular_automata',
    version='1.0.0',
    author='Greg Sotiropoulos',
    author_email='greg.sotiropoulos@gmail.com',
    description='Cellular Automata (Game of Life) simulator and GUI',
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords=[
        'cellular',
        'automaton',
        'game of life',
        'conway',
        'algorithms',
        'games'
    ],
    url='https://github.com/gregsotiropoulos/cellular_automata',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Games/Entertainment :: Life Games',
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy != 1.19.4',
        'scipy',
        'matplotlib',
        'keyboard'
    ],
    package_data={
        '': ['LICENSE', 'README.md'],
        # package name : [filename1, filename2, ...]
        'cellular_automata': ['glider.png']
    },
    project_urls={
        'Source': 'https://github.com/gregsotiropoulos/cellular_automata'
    }
)
