#!/bin/bash

# Activate virtual environment
source .venv/bin/activate

# Set any necessary environment variables for WSL
export DISPLAY=:0
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Install or update pyinstaller if needed
pip install --upgrade pyinstaller

# Build the application using the WSL-specific spec file
pyinstaller build-windows.spec --distpath build/build-wsl