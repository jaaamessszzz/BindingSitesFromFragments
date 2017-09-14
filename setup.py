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
    name='Binding_Sites_From_Fragments',
    version='0.1',
    author='James Lucas',
    author_email='james.lucas@berkeley.edu',
    description='',
    long_description=readme,
    url='https://github.com/jaaamessszzz/Binding-Sites-From-Fragments',
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
        'Binding_Sites_From_Fragments',
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
        'multiprocess',
        'mysqlclient',
        'networkx',
        'pathos',
        'prody',
        'scipy',
        'xmltodict',
        'appdirs'
    ],
    entry_points={
        'console_scripts': [
            'bsff = Binding_Sites_From_Fragments.usage:main',
        ],
    },
    include_package_data=True,
    zip_safe=False,
)