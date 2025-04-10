#!/bin/bash

# Exit on error
set -e

# Setup Frontend
echo "Setting up frontend dependencies..."
if [ -d "frontend" ]; then
    pushd frontend
    if [ -f "package.json" ]; then
        npm install
    else
        echo "Error: package.json not found in frontend directory"
        popd
        exit 1
    fi
    popd
else
    echo "Error: frontend directory not found"
    exit 1
fi

# Setup Backend
echo "Setting up backend dependencies..."
if [ -d "backend" ]; then
    pushd backend
    if [ -d ".venv" ]; then
        echo "Virtual environment already exists. Skipping creation."
    else
        echo "Creating virtual environment..."
        python3 -m venv .venv
        if [ $? -ne 0 ]; then
            echo "Error: Failed to create virtual environment"
            popd
            exit 1
        fi
    fi

    echo "Activating virtual environment..."
    source .venv/bin/activate
    if [ $? -ne 0 ]; then
        echo "Error: Failed to activate virtual environment"
        popd
        exit 1
    fi

    if [ -f "requirements.txt" ]; then
        echo "Installing requirements..."
        pip install -r requirements.txt
        if [ $? -ne 0 ]; then
            echo "Error: Failed to install requirements"
            deactivate
            popd
            exit 1
        fi
    else
        echo "Error: requirements.txt not found"
        deactivate
        popd
        exit 1
    fi

    if [ -f "requirements_dev.txt" ]; then
        echo "Installing development requirements..."
        pip install -r requirements_dev.txt
        if [ $? -ne 0 ]; then
            echo "Error: Failed to install development requirements"
            deactivate
            popd
            exit 1
        fi
    else
        echo "Error: requirements_dev.txt not found"
        deactivate
        popd
        exit 1
    fi

    deactivate
    popd
else
    echo "Error: backend directory not found"
    exit 1
fi

echo "Setup complete!"
