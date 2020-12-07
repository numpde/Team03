#!/usr/bin/env bash
env=$(pwd)/environment.yml
echo Writing to $env
conda env export -n Team03 | grep -v "^prefix:" > $env
