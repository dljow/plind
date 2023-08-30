import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="plind",
    version="0.1.0",
    author="Dylan Jow, Fang Xi Lin, Emily Tyhurst",
    description="Picard-Lefschetz integrator for highly oscillatory functions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    py_modules=["plind"],
    package_dir={'':'plind'},
    install_requires=['numpy','scipy','matplotlib']
)