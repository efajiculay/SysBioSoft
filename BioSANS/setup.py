import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="BioSANS2020-efajiculay", # Replace with your own username
    version="0.0.1",
    author="Erickson Erigio Fajiculay",
    author_email="efajiculay@yahoo.com",
    description="Symbolic and Numeric Software for Systems Biology",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)