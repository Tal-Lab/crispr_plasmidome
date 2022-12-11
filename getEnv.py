"""
Created on 11/12/2022

Author: Lucy Androsiuk
"""
import platform
import os
from dotenv import load_dotenv
from pathlib import Path

platform = platform.platform()
print(platform)

if "macOS" in platform:
    dotenv_path = Path('.env-mac')

if "Windows" in platform:
    dotenv_path = Path('.env-win')

load_dotenv(dotenv_path=dotenv_path)

print(os.getenv("PROJECT_PATH"))