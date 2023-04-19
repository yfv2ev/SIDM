"""Simple setup module to allow for pip installation"""

import setuptools

setuptools.setup(
    name='sidm',
    version='0.1.0',
    description='Self-Interacting Dark Matter analysis',
    url='https://github.com/btcardwell/SIDM',
    author='Bryan Cardwell',
    author_email='bryan.cardwell@cern.ch',
    license='BSD 3-clause',
    packages=['sidm', 'sidm.tools', 'sidm.definitions'],
    include_package_data=True,
)
