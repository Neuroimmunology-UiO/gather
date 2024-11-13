# Gather

This repository provides the environment setup for `gather`, including all necessary dependencies. Follow the instructions below to create and activate the environment using `conda`.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Prerequisites

Ensure you have Conda installed. If you don't have Conda, download and install it from the [official Anaconda website](https://www.anaconda.com/products/individual) or [Miniconda site](https://docs.conda.io/en/latest/miniconda.html).

## Installation

### Using the Environment File

To create and activate the `gather-env` environment using the provided `environment.yml` file, follow these steps:

1. **Clone the Repository**:

    ```bash
    git clone https://github.com/Neuroimmunology-UiO/gather.git
    cd gather
    ```

2. **Create the Environment**:

    Use the `environment.yml` file to create the environment. This will install all necessary dependencies.

    ```bash
    conda env create -f environment.yml
    ```

3. **Activate the Environment**:

    Once the environment is created, activate it using the following command:

    ```bash
    conda activate gather-env
    ```

### Installing Dependencies Separately

If you prefer to install dependencies manually or using `requirements.txt`, you can follow these instructions:

1. **Create and activate a conda environment** (optional but recommended):

    ```bash
    conda create --name gather-env python=3.8
    conda activate gather-env
    ```

2. **Install the dependencies using pip**:

    ```bash
    pip install -r requirements.txt
    ```

3. **Install `bcalm` (if not available on PyPI)**:

    ```bash
    conda install -c bioconda bcalm
    ```

## Usage

Once you have the environment set up and activated, you can use the `gather` environment for your project. Refer to the specific usage instructions or documentation for your project.

## Contributing

If you wish to contribute to this project, please follow these steps:

1. **Fork the repository**.
2. **Create a new branch** for your feature or bug fix.
3. **Make your changes** and commit them with clear, descriptive messages.
4. **Push your changes** to your forked repository.
5. **Create a pull request** to have your changes reviewed and merged into the main repository.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
