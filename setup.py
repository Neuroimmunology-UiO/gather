from setuptools import setup, find_packages

setup(
    name="gather",
    version="0.1.0",
    description="A short description of the gbat package",
    author="Neuroimmunology UiO",
    author_email="seyedmos@uio.no",
    url="https://github.com/Neuroimmunology-UiO/GATHER.git",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.19",
        "pandas",
        "networkx",
        "joblib",
        "tqdm",
        "biopython",
    ],
    entry_points={
        'console_scripts': [
            'G_assembler=gather.G_assembler:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
