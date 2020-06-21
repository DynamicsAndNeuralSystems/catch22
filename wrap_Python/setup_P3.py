from distutils.core import setup, Extension
import os

sourceDir = "../C/";

sourceFileList = [sourceDir+file for file in os.listdir(sourceDir) if file.endswith(".c") and not 'main' in file]; #  and not (file == "sampen.c" or file == "run_features.c")];

# the c++ extension module
extension_mod = Extension("catch22_C", 
	sources=["catch22_wrap_P3.c"] + sourceFileList, 
	include_dirs=[sourceDir])

# setup(name = "catch22", ext_modules=[extension_mod], packages=['catch22'])

setup(
name="catch22",
version="0.1.0",
author="Carl H Lubba",
url="https://github.com/chlubba/catch22",
description="CAnonical Time-series Features, see description and license on GitHub.",
ext_modules=[extension_mod],
packages=['catch22'],
classifiers=[
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Operating System :: OS Independent",
    ],
)
