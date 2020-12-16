import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyrotein",
    version="0.1.2",
    author="Cong Wang",
    author_email="wangimagine@gmail.com",
    description="A tiny package for structure analysis of macromolecules.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/carbonscott/pyrotein",
    keywords = ['PDB', 'structure biology', 'protein', 'analysis'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
