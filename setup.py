from setuptools import setup, find_packages

setup(
    name="gather",
    version="0.1.0",
    description="Assembly of heavy and light chain BCR using De Bruijn graph",
    author="Neuroimmunology UiO",
    author_email="seyedmos@uio.no",
    url="https://github.com/Neuroimmunology-UiO/gather.git",
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
            'gassembler=gather.gassembler:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
