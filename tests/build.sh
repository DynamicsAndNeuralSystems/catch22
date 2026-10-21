#!/bin/sh
# Compile sources into one executable with $CC and $CFLAGS: tests/build.sh <out> <source>...
# Takes GCC-style compilers and MSVC's cl, whose output and library flags differ; CI's
# build matrix calls this so each entry is just a compiler and a flag string.
out=$1; shift
case "$CC" in
  cl) cl -nologo $CFLAGS -Fe:"$out.exe" "$@" ;;
  *)  $CC $CFLAGS "$@" -o "$out" -lm ;;
esac
