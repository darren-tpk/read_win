from setuptools import setup, find_packages

setup(
    name="read_win",
    version="0.1.0",
    description="Read WIN seismic/acoustic data into ObsPy Streams",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Darren Tan, Gilles Seropian",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "obspy",
    ],
    python_requires=">=3.10",
)