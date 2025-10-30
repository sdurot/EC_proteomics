#!/usr/bin/env python3
"""
Virtual Environment Setup Script for EC Proteomics Analysis

This script provides multiple options for setting up your analysis environment:
1. Standard Python virtual environment
2. UV-based setup (if UV is available)
3. Conda environment setup

Usage:
    python setup_environment.py [--method venv|uv|conda] [--jupyter]
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path

def run_command(cmd, description, check=True):
    """Run a command with error handling."""
    print(f" {description}")
    print(f"   Running: {cmd}")
    
    try:
        if isinstance(cmd, str):
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        else:
            result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0 and check:
            print(f"❌ Error: {result.stderr}")
            return False
        else:
            if result.stdout:
                print(f"   Output: {result.stdout.strip()}")
            print("✅ Success!")
            return True
    except Exception as e:
        print(f"❌ Exception: {e}")
        return False

def setup_venv():
    """Set up standard Python virtual environment."""
    print("\n Setting up Python Virtual Environment")
    print("=" * 50)
    
    # Create virtual environment
    if not run_command([sys.executable, "-m", "venv", ".venv"], 
                      "Creating virtual environment"):
        return False
    
    # Determine activation script path
    if os.name == 'nt':  # Windows
        activate_script = ".venv\\Scripts\\activate"
        pip_path = ".venv\\Scripts\\pip.exe"
        python_path = ".venv\\Scripts\\python.exe"
    else:  # Unix/Linux/macOS
        activate_script = "source .venv/bin/activate"
        pip_path = ".venv/bin/pip"
        python_path = ".venv/bin/python"
    
    print(f"\n Installing packages with {pip_path}")
    
    # Upgrade pip first
    if not run_command([python_path, "-m", "pip", "install", "--upgrade", "pip"], 
                      "Upgrading pip"):
        print("⚠️ Warning: Could not upgrade pip, continuing anyway...")
    
    # Install requirements
    if os.path.exists("requirements.txt"):
        if run_command([pip_path, "install", "-r", "requirements.txt"], 
                      "Installing requirements"):
            print("\n Virtual environment setup complete!")
            print(f"\n To activate the environment:")
            if os.name == 'nt':
                print(f"   .venv\\Scripts\\activate")
            else:
                print(f"   source .venv/bin/activate")
            return True
    else:
        print("❌ requirements.txt not found!")
        return False

def setup_uv():
    """Set up UV-based environment."""
    print("\n Setting up UV Environment")
    print("=" * 50)
    
    # Check if UV is available
    if not run_command(["uv", "--version"], "Checking UV installation", check=False):
        print("❌ UV not found. Installing UV...")
        if not run_command([sys.executable, "-m", "pip", "install", "uv"], 
                          "Installing UV"):
            return False
    
    # Initialize UV project
    if not run_command(["uv", "init", "--no-readme"], "Initializing UV project", check=False):
        print("⚠️ UV init failed, trying alternative approach...")
    
    # Create UV environment and install packages
    if os.path.exists("requirements.txt"):
        # Try to install from requirements.txt
        if run_command(["uv", "pip", "install", "-r", "requirements.txt"], 
                      "Installing packages with UV"):
            print("\n✅ UV environment setup complete!")
            print(f"\n To use UV environment:")
            print(f"   uv run python your_script.py")
            print(f"   uv run jupyter notebook")
            return True
        else:
            print("❌ UV package installation failed")
            print(" This might be due to network drive issues")
            print("   Trying fallback to standard virtual environment...")
            return setup_venv()
    else:
        print("❌ requirements.txt not found!")
        return False

def setup_conda():
    """Set up Conda environment."""
    print("\n Setting up Conda Environment")
    print("=" * 50)
    
    env_name = "ec-proteomics"
    
    # Check if conda is available
    if not run_command(["conda", "--version"], "Checking Conda installation", check=False):
        print("❌ Conda not found. Please install Anaconda or Miniconda first.")
        return False
    
    # Create conda environment
    if not run_command(["conda", "create", "-n", env_name, "python=3.11", "-y"], 
                      f"Creating conda environment '{env_name}'"):
        return False
    
    # Install packages
    if os.path.exists("requirements.txt"):
        # Use pip within conda environment
        if run_command([
            "conda", "run", "-n", env_name, "pip", "install", "-r", "requirements.txt"
        ], "Installing packages in conda environment"):
            print(f"\n✅ Conda environment '{env_name}' setup complete!")
            print(f"\n To activate the environment:")
            print(f"   conda activate {env_name}")
            return True
    else:
        print("❌ requirements.txt not found!")
        return False

def launch_jupyter(method):
    """Launch Jupyter notebook."""
    print("\n Launching Jupyter Notebook")
    print("=" * 50)
    
    if method == "venv":
        if os.name == 'nt':
            cmd = [".venv\\Scripts\\jupyter.exe", "notebook"]
        else:
            cmd = [".venv/bin/jupyter", "notebook"]
    elif method == "uv":
        cmd = ["uv", "run", "jupyter", "notebook"]
    elif method == "conda":
        cmd = ["conda", "run", "-n", "ec-proteomics", "jupyter", "notebook"]
    else:
        cmd = ["jupyter", "notebook"]
    
    print(f" Starting Jupyter with command: {' '.join(cmd)}")
    try:
        subprocess.run(cmd)
    except KeyboardInterrupt:
        print("\n Jupyter notebook stopped.")
    except Exception as e:
        print(f"❌ Error launching Jupyter: {e}")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Setup EC Proteomics Analysis Environment")
    parser.add_argument("--method", choices=["venv", "uv", "conda"], default="venv",
                       help="Environment setup method")
    parser.add_argument("--jupyter", action="store_true",
                       help="Launch Jupyter notebook after setup")
    
    args = parser.parse_args()
    
    print(" EC Proteomics Analysis - Environment Setup")
    print("=" * 60)
    print(f" Method: {args.method}")
    print(f" Launch Jupyter: {args.jupyter}")
    
    # Check if we're in the right directory
    if not os.path.exists("proteomics_analysis.py"):
        print("❌ Error: proteomics_analysis.py not found!")
        print("   Please run this script from the EC_proteomics directory.")
        return 1
    
    # Setup environment based on method
    success = False
    if args.method == "venv":
        success = setup_venv()
    elif args.method == "uv":
        success = setup_uv()
    elif args.method == "conda":
        success = setup_conda()
    
    if not success:
        print(f"\n❌ Environment setup failed with method '{args.method}'")
        if args.method != "venv":
            print(" Try fallback: python setup_environment.py --method venv")
        return 1
    
    print(f"\n Setup complete! Environment ready for proteomics analysis.")
    
    # Launch Jupyter if requested
    if args.jupyter:
        launch_jupyter(args.method)
    
    
    return 0

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)