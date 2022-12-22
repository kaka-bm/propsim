from setuptools import setup
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='propsim',
    version='0.0.5',
    description='A simple and intuitive tool for simulating different types of aircraft engines.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=['ambiance','pyyaml'],
    py_modules=["propsim"],
    package_dir={'':'src'}
)