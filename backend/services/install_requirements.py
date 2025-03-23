import sys
import subprocess
import argparse
from pathlib import Path
import platform
import os

# Base directory for the project
ROOT_DIR = Path(__file__).resolve().parent.parent.parent
BASE_DIR = Path(__file__).resolve().parent

# Configuration for all services (only the "path" key is used in this script)
SERVICES = {
    "api-gateway": {
        "path": "api_gateway",
    },
    "download": {
        "path": "download_service",
    },
    "aligning": {
        "path": "aligning_service",
    },
    "analysis": {
        "path": "analysis_service", 
    },
    "conversion": {
        "path": "conversion_service",
    },
    "database": {
        "path": "database_service",
    }
}


def run_command(cmd, cwd=None, check=True):
    """Run a shell command and print output"""
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(
        cmd,
        cwd=cwd,
        check=False,
        text=True,
        capture_output=True
    )
    
    if result.stdout:
        print(result.stdout)
    
    if result.stderr:
        print(f"Error: {result.stderr}", file=sys.stderr)
        
    if check and result.returncode != 0:
        print(f"Command failed with exit code {result.returncode}")
        return False
        
    return result.returncode == 0


def install_python_requirements():
    """Install Python requirements for all services"""
    print("\n=== Installing Python Requirements ===")
    
    success = True
    for service_name, service_info in SERVICES.items():
        service_path = BASE_DIR / service_info["path"]
        print(service_path)
        req_file = service_path / "requirements.txt"
        
        if req_file.exists():
            print(f"\nInstalling requirements for {service_name}...")
            if not run_command([sys.executable, "-m", "pip", "install", "-r", str(req_file)]):
                success = False
                print(f"Failed to install requirements for {service_name}")
        else:
            print(f"No requirements.txt found for {service_name}")
    
    return success


def install_blast_plus():
    """Install BLAST+ package based on the operating system"""
    print("\n=== Installing BLAST+ ===")
    
    system = platform.system().lower()
    
    if system == "linux":
        # For Ubuntu/Debian
        if os.path.exists("/usr/bin/apt"):
            print("Installing BLAST+ using apt...")
            return run_command(["sudo", "apt", "update"]) and \
                   run_command(["sudo", "apt", "install", "-y", "ncbi-blast+"])
        # For CentOS/RHEL/Fedora
        elif os.path.exists("/usr/bin/yum"):
            print("Installing BLAST+ using yum...")
            return run_command(["sudo", "yum", "install", "-y", "ncbi-blast+"])
        # For Arch Linux
        elif os.path.exists("/usr/bin/pacman"):
            print("Installing BLAST+ using pacman...")
            return run_command(["sudo", "pacman", "-Sy", "blast"])
        else:
            print("Unsupported Linux distribution. Please install BLAST+ manually.")
            return False
            
    elif system == "darwin":  # macOS
        # Check if homebrew is installed
        if os.path.exists("/usr/local/bin/brew") or os.path.exists("/opt/homebrew/bin/brew"):
            print("Installing BLAST+ using Homebrew...")
            return run_command(["brew", "install", "blast"])
        else:
            print("Homebrew not found. Please install Homebrew or install BLAST+ manually.")
            return False
            
    elif system == "windows":
        print("For Windows, we recommend installing BLAST+ manually from NCBI:")
        print("https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")
        print("\nAlternatively, if you're using WSL, run this script within the WSL environment.")
        return False
        
    else:
        print(f"Unsupported operating system: {system}")
        return False


def check_blast_installation():
    """Check if BLAST+ is properly installed"""
    print("\n=== Checking BLAST+ Installation ===")
    
    try:
        result = subprocess.run(
            ["blastn", "-version"],
            check=False,
            text=True,
            capture_output=True
        )
        
        if result.returncode == 0:
            print(f"BLAST+ is installed: {result.stdout.strip()}")
            return True
        else:
            print("BLAST+ is not properly installed or not in the PATH")
            return False
            
    except FileNotFoundError:
        print("BLAST+ is not installed or not in the PATH")
        return False


def create_venv():
    """Create a virtual environment if it doesn't exist"""
    print("\n=== Setting up Python Virtual Environment ===")
    
    venv_path = ROOT_DIR / ".venv"
    if venv_path.exists():
        print("Virtual environment already exists")
        return True
        
    print(f"Creating virtual environment at {venv_path}...")
    try:
        import venv
        venv.create(venv_path, with_pip=True)
        
        # Print activation instructions
        print("\nTo activate the virtual environment:")
        if platform.system().lower() == "windows":
            print(f"    {venv_path}\\Scripts\\activate")
        else:
            print(f"    source {venv_path}/bin/activate")
            
        return True
    except Exception as e:
        print(f"Failed to create virtual environment: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Setup script for KATH-PVP project")
    parser.add_argument("--no-venv", action="store_true", help="Skip virtual environment creation")
    parser.add_argument("--no-blast", action="store_true", help="Skip BLAST+ installation")
    parser.add_argument("--no-requirements", action="store_true", help="Skip Python requirements installation")
    args = parser.parse_args()
    
    print("=== KATH-PVP Setup Script ===")
    
    # Create virtual environment
    if not args.no_venv:
        if not create_venv():
            print("Failed to create virtual environment")
    
    # Install Python requirements
    if not args.no_requirements:
        if not install_python_requirements():
            print("Warning: Some Python requirements failed to install")
    
    # Install BLAST+
    if not args.no_blast:
        if not install_blast_plus():
            print("Warning: Failed to install BLAST+")
    
    # Check BLAST+ installation
    blast_installed = check_blast_installation()
    
    print("\n=== Setup Summary ===")
    print(f"- Directories: Created")
    print(f"- BLAST+: {'Installed' if blast_installed else 'Not installed properly'}")
    print(f"- Python Requirements: {'Installed' if not args.no_requirements else 'Skipped'}")
    
    print("\nSetup complete!")
    
    if not blast_installed:
        print("\nWARNING: BLAST+ is not properly installed or not in the PATH.")
        print("The aligning service may not work without BLAST+.")
        print("WSL may be required for Windows users. (python does not support windows installation)")
        print("Please install BLAST+ manually from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())