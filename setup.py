"""
Setup configuration for FaCC (Fibonacci Cosmology)
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
readme_file = Path(__file__).parent / "Readme.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

# Read requirements
requirements_file = Path(__file__).parent / "requirements.txt"
if requirements_file.exists():
    with open(requirements_file, "r") as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]
else:
    requirements = []

setup(
    name="facc",
    version="0.1.0",
    author="Bryan David Persaud",
    author_email="bryan@imediacorp.com",
    description="Fibonacci Cosmology: Testing φ-modulation in cosmic structure",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/imediacorp/FaCC",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.9",
    install_requires=requirements,
    extras_require={
        "dev": [
            "flake8>=6.0.0",
            "pytest>=7.0.0",
        ],
    },
)

