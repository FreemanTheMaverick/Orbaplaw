from setuptools import setup, find_packages

__version__ = "1.0.2"

setup(
		name = "Orbaplaw",
		version = __version__,
		author = "FreemanTheMaverick",
		description = "Orbital alignment analysis for plane wave basis sets",
		packages = find_packages(),
		url = "https://github.com/FreemanTheMaverick/Orbaplaw.git",
		install_requires = ["numpy", "scipy", "pyscf", "Maniverse"],
		classifiers = ["Programming Language :: Python :: 3"]
)
