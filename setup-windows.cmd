@REM @echo off
@REM setlocal enabledelayedexpansion

@REM echo ===== Starting Project Setup =====

@REM REM Check if Python is installed
@REM py --version > nul 2>&1
@REM if %ERRORLEVEL% NEQ 0 (
@REM     echo Error: Python is not installed or not in PATH
@REM     echo Please install Python from https://www.python.org/downloads/
@REM     exit /b 1
@REM )

@REM REM Setup Frontend (uncomment if needed)
@REM REM echo Setting up frontend dependencies...
@REM REM if exist package.json (
@REM REM     npm install
@REM REM ) else (
@REM REM     echo Warning: package.json not found in root directory
@REM REM )
@REM REM 
@REM REM if exist frontend (
@REM REM     cd frontend
@REM REM     if exist package.json (
@REM REM         npm install
@REM REM     ) else (
@REM REM         echo Warning: package.json not found in frontend directory
@REM REM     )
@REM REM     cd ..
@REM REM )

@REM REM Setup Backend
@REM echo Setting up backend environment...
@REM cd backend || (
@REM     echo Error: backend directory not found
@REM     exit /b 1
@REM )

@REM if exist .venv (
@REM     echo Virtual environment already exists. Skipping creation.
@REM ) else (
@REM     echo Creating virtual environment...
@REM     py -m venv .venv
@REM     if %ERRORLEVEL% NEQ 0 (
@REM         echo Error: Failed to create virtual environment
@REM         exit /b 1
@REM     )
@REM )

@REM echo Activating virtual environment and installing requirements...
@REM call .venv\Scripts\activate.bat
@REM if %ERRORLEVEL% NEQ 0 (
@REM     echo Error: Failed to activate virtual environment
@REM     exit /b 1
@REM )

@REM if exist requirements.txt (
@REM     echo Installing main requirements...
@REM     pip install -r requirements.txt
@REM     if %ERRORLEVEL% NEQ 0 (
@REM         echo Error: Failed to install requirements
@REM         exit /b 1
@REM     )
@REM ) else (
@REM     echo Warning: requirements.txt not found
@REM )

@REM if exist requirements_dev.txt (
@REM     echo Installing development requirements...
@REM     pip install -r requirements_dev.txt
@REM     if %ERRORLEVEL% NEQ 0 (
@REM         echo Error: Failed to install development requirements
@REM         exit /b 1
@REM     )
@REM ) else (
@REM     echo Warning: requirements_dev.txt not found
@REM )

@REM call deactivate
@REM cd ..

@REM echo ===== Setup complete! =====
@REM echo To activate the virtual environment, run: 
@REM echo     cd backend ^& .venv\Scripts\activate.bat

@REM exit /b 0

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
