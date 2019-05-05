"""Setup script for ``neutcurve``.

Written by Jesse Bloom.
"""

import re
import sys
try:
    from setuptools import setup
except ImportError:
    raise ImportError("You must install `setuptools`")

if not (sys.version_info[0] == 3 and sys.version_info[1] >= 6):
    raise RuntimeError(
                'neutcurve requires Python 3.6 or higher.\n'
                'You are using Python {0}.{1}'.format(
                    sys.version_info[0], sys.version_info[1])
                )

# get metadata from package `__init__.py` file as here:
# https://packaging.python.org/guides/single-sourcing-package-version/
metadata = {}
init_file = 'neutcurve/__init__.py'
with open(init_file) as f:
    init_text = f.read()
for dataname in ['version', 'author', 'email', 'url']:
    matches = re.findall(
            '__' + dataname + r'__\s+=\s+[\'"]([^\'"]+)[\'"]',
            init_text)
    if len(matches) != 1:
        raise ValueError(f"found {len(matches)} matches for {dataname} "
                         f"in {init_file}")
    else:
        metadata[dataname] = matches[0]

with open('README.rst') as f:
    readme = f.read()

# main setup command
setup(
    name='neutcurve',
    version=metadata['version'],
    author=metadata['author'],
    author_email=metadata['email'],
    url=metadata['url'],
    download_url='https://github.com/jbloomlab/neutcurve/tarball/' +
                 metadata['version'],  # tagged version on GitHub
    description='Plot and fit neutralization curves.',
    long_description=readme,
    license='GPLv3',
    install_requires=[
        'dmslogo>=0.2.1',
        'matplotlib>=3.0.0',
        'pandas>=0.24',
        'pyyaml>=3.13',
        'scipy>=1.1.0',
        'xlrd>=1.2',
        ],
    platforms='Linux and Mac OS X.',
    packages=['neutcurve'],
    package_dir={'neutcurve': 'neutcurve'},
)
