#!/usr/bin/env python3
# encoding: utf-8

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.md') as file:
    readme = file.read()

# Setup
setup(
    name='BindingSitesFromFragments',
    version='0.1',
    author='James Lucas',
    author_email='james.lucas@berkeley.edu',
    description='',
    long_description=readme,
    url='https://github.com/jaaamessszzz/BindingSitesFromFragments',
    keywords=[
        'Binding Sites',
        'Fragments'
    ],
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    packages=[
        'BindingSitesFromFragments',
    ],
    install_requires=[
        'docopt',
        'matplotlib',
        'numpy',
        'decorator',
        'pubchempy',
        'pypdb',
        'pandas',
        'biopython',
        'seaborn',
        'multiprocess',
        'networkx',
        'pathos',
        'prody',
        'scipy',
        'xmltodict',
        'appdirs'
    ],
    entry_points={
        'console_scripts': [
            'bsff = BindingSitesFromFragments.commands.bsff:main',
            'bsff_clean = BindingSitesFromFragments.commands.bsff:main'
        ],
    },
    include_package_data=True,
    zip_safe=False,
)