from setuptools import setup, find_packages
from pathlib import Path

about = {}
here = Path(__file__).parent

with open(here / "gather" / "__init__.py") as f:
    exec(f.read(), about)

setup(
    name="gather",
    version=about["__version__"],
    description=about["__description__"],
    long_description=(here / "README.md").read_text(),
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
            "tenx-asm=gather.tenx_asm:main",
            "sc-asm=gather.sc_asm:main",
            "postproc=gather.postproc:main",
        ],
    },
    include_package_data=True,
    package_data={
        "gather": [
            "database/**/*",
        ]
    },
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)
