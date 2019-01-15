from distutils.core import setup, Extension
import os

sourceDir = "../C/";

sourceFileList = [sourceDir+file for file in os.listdir(sourceDir) if file.endswith(".c") and not 'main' in file]; #  and not (file == "sampen.c" or file == "run_features.c")];

# the c++ extension module
extension_mod = Extension("catch22_C", 
	sources=["catch22_wrap.c"] + sourceFileList, 
	include_dirs=[sourceDir])

setup(name = "catch22", ext_modules=[extension_mod])