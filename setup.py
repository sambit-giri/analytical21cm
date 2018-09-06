'''
Created on 06 September 2018
@author: Raghunath Ghara, Sambit Giri
Setup script
'''

#from setuptools import setup, find_packages
from distutils.core import setup
import os, sys

os.chdir('src/')
os.system('bash makelibfile.sh')
os.chdir('../')

setup(name='analytical21cm',
      version='0.1',
      author='Raghunath Ghara, Sambit Giri',
      author_email='raghunath.ghara@astro.su.se, sambit.giri@astro.su.se',
      package_dir = {'analytical21cm' : 'src'},
      packages=['analytical21cm'],
      #package_data={'src' : ['ANRAGpylib.so'],},
      #include_package_data=True,
)

dst = sys.argv[2].split('=')[-1]
if '~' in dst: dst = os.path.expanduser("~")+'/'
from shutil import copyfile
copyfile('src/ANRAGpylib.so', dst+'/lib/python/analytical21cm/ANRAGpylib.so')

