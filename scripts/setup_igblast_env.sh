#!/bin/bash
set -euo pipefail

# IgBLAST installer (Linux & macOS)
VERSION="1.22.0"
OS_TYPE=$(uname -s)
ARCH_TYPE=$(uname -m)

# Pick correct tarball
case "$OS_TYPE" in
    Linux)
        [ "$ARCH_TYPE" = "x86_64" ] || { echo "Unsupported architecture: $ARCH_TYPE"; exit 1; }
        FILE="ncbi-igblast-${VERSION}-x64-linux.tar.gz"
        ;;
    Darwin)
        [ "$ARCH_TYPE" = "x86_64" ] || { echo "Unsupported architecture: $ARCH_TYPE"; exit 1; }
        FILE="ncbi-igblast-${VERSION}-x64-macosx.tar.gz"
        ;;
    *)
        echo "Unsupported OS: $OS_TYPE"
        exit 1
        ;;
esac

echo "==> Installing IgBLAST ${VERSION} ..."

# Download & extract
wget -q "https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/${VERSION}/${FILE}"
tar -zxf "${FILE}"
rm -f "${FILE}"

# Create data dirs
mkdir -p ~/share/igblast ~/share/germlines/imgt

# Copy executables
cp ncbi-igblast-${VERSION}/bin/* ~/bin

# Copy support files
cp -r ncbi-igblast-${VERSION}/internal_data ~/share/igblast
cp -r ncbi-igblast-${VERSION}/optional_file ~/share/igblast

# Fetch reference databases (requires Immcantation helper scripts in PATH)
fetch_igblastdb.sh -o ~/share/igblast
fetch_imgtdb.sh -o ~/share/germlines/imgt
imgt2igblast.sh -i ~/share/germlines/imgt -o ~/share/igblast

# Export IGDATA
if ! grep -q 'export IGDATA=~/share/igblast' ~/.bashrc; then
    echo 'export IGDATA=~/share/igblast' >> ~/.bashrc
fi
export IGDATA=~/share/igblast

echo "==> IgBLAST setup completed."
echo "   Test with: igblastn -version"
