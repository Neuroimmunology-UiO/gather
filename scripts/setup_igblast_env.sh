#!/bin/bash

# Set the IgBLAST version
VERSION="1.22.0"

# Determine the operating system
OS_TYPE=$(uname -s)
ARCH_TYPE=$(uname -m)

# Set download URL based on OS
case "$OS_TYPE" in
    Linux)
        if [ "$ARCH_TYPE" = "x86_64" ]; then
            FILE="ncbi-igblast-${VERSION}-x64-linux.tar.gz"
        else
            echo "Unsupported architecture: $ARCH_TYPE"
            exit 1
        fi
        ;;
    Darwin)
        if [ "$ARCH_TYPE" = "x86_64" ]; then
            FILE="ncbi-igblast-${VERSION}-x64-macosx.tar.gz"
        else
            echo "Unsupported architecture: $ARCH_TYPE"
            exit 1
        fi
        ;;
    *)
        echo "Unsupported OS: $OS_TYPE"
        exit 1
        ;;
esac

# Download and extract IgBLAST
wget "https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/${VERSION}/${FILE}"
tar -zxf "${FILE}"

# Create necessary directories
mkdir -p ~/bin
mkdir -p ~/share/igblast
mkdir -p ~/share/germlines/imgt

# Copy executables to ~/bin
cp ncbi-igblast-${VERSION}/bin/* ~/bin

# Download reference databases and setup IGDATA directory
fetch_igblastdb.sh -o ~/share/igblast
cp -r ncbi-igblast-${VERSION}/internal_data ~/share/igblast
cp -r ncbi-igblast-${VERSION}/optional_file ~/share/igblast

# Build IgBLAST database from IMGT reference sequences
fetch_imgtdb.sh -o ~/share/germlines/imgt
imgt2igblast.sh -i ~/share/germlines/imgt -o ~/share/igblast

# Set IGDATA environment variable
echo 'export IGDATA=~/share/igblast' >> ~/.bashrc
export IGDATA=~/share/igblast

echo "IgBLAST setup completed successfully."

