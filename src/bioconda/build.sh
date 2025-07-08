#!/bin/bash

# Compile the project
make all

# Install the binaries
mkdir -p "$PREFIX/bin"
cp bin/* "$PREFIX/bin"
