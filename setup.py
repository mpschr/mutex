from distutils.core import setup
from setuptools import find_packages

setup(name='mutex',
      version='0.1.0',
      description='A simple algorithm to estimate the significance of a mutual exclusive pattern or co-occurrence',
      author='Michael P Schroeder',
      author_email='michael.p.schroeder@gmail.com',
      url='https://github.com/mpschr/mutex',
      packages=find_packages(),
      install_requires=['pandas', 'numpy']
     )
