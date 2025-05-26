from setuptools import setup, find_packages

__version__ = "1.2.0"

setup(
		name = "Orbaplaw",
		version = __version__,
		author = "FreemanTheMaverick",
		description = "Orbital alignment analysis for plane wave basis sets",
		long_description = open("README.md").read(),
		long_description_content_type = "text/markdown",
		packages = find_packages(),
		entry_points = { "console_scripts": [
			"orbaplaw = Orbaplaw.main:main"
		]},
		url = "https://github.com/FreemanTheMaverick/Orbaplaw.git",
		install_requires = ["numpy", "scipy", "pyscf", "Maniverse"],
		classifiers = ["Programming Language :: Python :: 3"]
)
