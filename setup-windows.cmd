@echo off

REM Setup Frontend
echo Setting up frontend dependencies...
if exist package.json (
    call npm install
) else (
    echo Error: package.json not found in root directory
    exit /b 1
)

if exist frontend (
    pushd frontend
    if exist package.json (
        call npm install
    ) else (
        echo Error: package.json not found in frontend directory
        popd
        exit /b 1
    )
    popd
) else (
    echo Error: frontend directory not found
    exit /b 1
)

REM Setup Backend
echo Setting up backend dependencies...
if exist backend (
    pushd backend
    if exist .venv (
        echo Virtual environment already exists. Skipping creation.
    ) else (
        echo Creating virtual environment...
        py -m venv .venv
        if %ERRORLEVEL% NEQ 0 (
            echo Error: Failed to create virtual environment
            popd
            exit /b 1
        )
    )
    echo Activating virtual environment...
    call .venv\Scripts\activate.bat
    if %ERRORLEVEL% NEQ 0 (
        echo Error: Failed to activate virtual environment
        popd
        exit /b 1
    )
    if exist requirements.txt (
        echo Installing requirements...
        pip install -r requirements.txt
        if %ERRORLEVEL% NEQ 0 (
            echo Error: Failed to install requirements
            call deactivate
            popd
            exit /b 1
        )
    ) else (
        echo Error: requirements.txt not found
        call deactivate
        popd
        exit /b 1
    )
    if exist requirements_dev.txt (
        echo Installing development requirements...
        pip install -r requirements_dev.txt
        if %ERRORLEVEL% NEQ 0 (
            echo Error: Failed to install development requirements
            call deactivate
            popd
            exit /b 1
        )
    ) else (
        echo Error: requirements_dev.txt not found
        call deactivate
        popd
        exit /b 1
    )
    call deactivate
    popd
) else (
    echo Error: backend directory not found
    exit /b 1
)

echo Setup complete!
