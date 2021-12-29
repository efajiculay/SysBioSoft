import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
	name="BioSANS2020",
	version="0.2.7",
    author="Erickson Erigio Fajiculay",
    author_email="efajiculay@yahoo.com",
    description="Symbolic and Numeric Software for Systems Biology",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/efajiculay/SysBioSoft",
	package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
	license="GPLv3",
    classifiers=[
		"Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 "
		"(GPLv3)",
        "Operating System :: OS Independent",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
	install_requires=[
		"DateTime",
		"pathlib",
		"Pillow",
		"numpy",
		"thread6",
		"pandas",
		"matplotlib==3.3.3",
		"python-libsbml",
		"scipy",
		"func_timeout",
		"sympy",
		"applescript"
	],
    entry_points={
		'gui_scripts' : ['BioSANS=BioSANS2020.RunBioSANS:BioSANS'], 
        'console_scripts': ['BioSSL=BioSANS2020.RunBioSANS:BioSSL']
    },
    python_requires='>=3.7',
)

#https://python-packaging.readthedocs.io/en/latest/metadata.html
#https://pypi.org/pypi?%3Aaction=list_classifiers




