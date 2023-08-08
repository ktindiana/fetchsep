#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'astropy',
    'python-dateutil',
    'matplotlib',
    'numpy',
    'scipy',
    'pandas',
    'wget',
    'zulu',
    'asciitable',
    'lmfit',
    'cycler',
    'scikit-learn',
    'seaborn',
    'netCDF4'
]


test_requirements = [ ]

setup(
    author="Kathryn Whitman",
    author_email='kathryn.whitman@nasa.gov',
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    description="Identify SEP elevations above background in a time series (idsep) and analyze events individually (opsep)",
    entry_points={
        'console_scripts': [
            'fetchsep=fetchsep.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='fetchsep',
    name='fetchsep',
    packages=find_packages(include=['fetchsep', 'fetchsep.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/rickyegeland/fetchsep',
    version='0.1.0',
    zip_safe=False,
)
