from setuptools import setup, find_packages
from pathlib import Path
import versioneer

DESCRIPTION = 'A Python-based converter for sonar logs used in PINGMapper'
LONG_DESCRIPTION = Path('README.md').read_text()

setup(
    name="pingverter",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Cameron Bodine",
    author_email="bodine.cs@gmail.email",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        ],
    python_requires=">=3.6",
    install_requires=["versioneer"],
)