import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
	#name="BioSANS2020",
    name="BioSANS2020-efajiculay",
    #version="0.1.3",
	#version="0.0.2",
    author="Erickson Erigio Fajiculay",
    author_email="efajiculay@yahoo.com",
    description="Symbolic and Numeric Software for Systems Biology",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/efajiculay/SysBioSoft",
    packages=setuptools.find_packages(),
	license="GPLv3",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 "
		"(GPLv3)",
        "Operating System :: OS Independent",
    ],
	install_requires=[
		"DateTime",
		"Pillow",
		"numpy",
		"pandas",
		"matplotlib==3.3.3",
		"python-libsbml",
		"scipy",
		"sdeint",
		"sympy",
		"applescript"
	],
    python_requires='>=3.7',
)

