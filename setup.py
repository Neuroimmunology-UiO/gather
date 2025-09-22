from setuptools import setup, find_packages

about = {}
with open("gather/__init__.py") as f:
    exec(f.read(), about)

setup(
    name="gather",
    version=about["__version__"],
    description=about["__description__"],
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author=about["__author__"],
    author_email="seyedmos@uio.no",
    url=about["__url__"],
    packages=find_packages(exclude=["tests*", "docs*"]),
    install_requires=[
        "numpy>=1.19",
        "pandas",
        "networkx",
        "joblib",
        "tqdm",
        "biopython",
        "argcomplete",
        "glob2",
        "rich",
    ],
    entry_points={
    "console_scripts": [
        "10x-asm=gather.tenx_asm:main",
        "sc-asm=gather.sc_asm:main",
        "postproc=gather.postproc:main",
    ],
    },

    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

