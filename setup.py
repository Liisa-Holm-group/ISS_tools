from setuptools import setup, find_packages

setup(
    name="ISS_tools",
    version="0.1.0",
    packages=find_packages(include=["ISS_tools", "ISS_tools.*"]),
    install_requires=[
        "pandas",
        # add any other dependencies here
    ],
    entry_points={
        "console_scripts": [
            "pfamtool=scripts.pfamtool:main",  # if pfamtool.py defines a main()
        ],
    },
    author="Liisa Holm",
    description="Tools for parsing and annotating ISS/DALI output using Pfam",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3.7",
)
