#!/usr/bin/env bash
set -euo pipefail

### ================================
### Gather external deps installer
### - BLAST+, SPAdes, Immcantation scripts, IgPhyML
### - Linux and macOS (Darwin)
### - No sudo required; installs into $HOME
### - Excludes IgBLAST/IMGT and all R packages
### ================================

# ---- Versions (edit here if needed) ----
BLAST_VER="2.16.0"
SPADES_VER="4.2.0"
IGPHYML_REPO="https://github.com/immcantation/igphyml.git"
IMMCANTATION_REPO="https://github.com/immcantation/immcantation.git"

# ---- Directories ----
INSTALL_DIR="${HOME}/opt"
BIN_DIR="${HOME}/bin"
mkdir -p "${INSTALL_DIR}" "${BIN_DIR}"

OS_TYPE="$(uname -s)"
ARCH_TYPE="$(uname -m)"

echo "==> Detected system: ${OS_TYPE} ${ARCH_TYPE}"
echo "==> Install dir: ${INSTALL_DIR}"
echo "==> User bin dir: ${BIN_DIR}"

need_cmd () { command -v "$1" >/dev/null 2>&1 || return 1; }
fetch () {
  # fetch URL to file ($1=url, $2=output)
  if need_cmd wget; then
    wget -q -O "$2" "$1"
  elif need_cmd curl; then
    curl -LfsS -o "$2" "$1"
  else
    echo "ERROR: Need either wget or curl." >&2
    exit 1
  fi
}

append_path_once () {
  local line='export PATH="${HOME}/bin:$PATH"'
  # Update both bash and zsh profiles as appropriate
  for rc in "${HOME}/.bashrc" "${HOME}/.zshrc"; do
    [ -f "$rc" ] || continue
    grep -Fq "$line" "$rc" || echo "$line" >> "$rc"
  done
  export PATH="${HOME}/bin:$PATH"
}

### --------------------------
### 1) BLAST+ (Linux/macOS)
### --------------------------
install_blast () {
  echo "==> Installing NCBI BLAST+ ${BLAST_VER} ..."
  local tar=""
  case "$OS_TYPE" in
    Linux)  tar="ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz" ;;
    Darwin)
      # NCBI provides macOS x64 builds; on Apple Silicon we try anyway
      tar="ncbi-blast-${BLAST_VER}+-x64-macosx.tar.gz"
      if [ "${ARCH_TYPE}" = "arm64" ]; then
        echo "   [Note] macOS arm64 detected. The BLAST+ macOS tarball is x86_64."
        echo "          If it fails to run, use a conda fallback:  conda install -c bioconda blast"
      fi
      ;;
    *) echo "Unsupported OS for BLAST+: $OS_TYPE"; return 1 ;;
  esac

  local url="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VER}/${tar}"
  local tmp="${PWD}/${tar}"
  fetch "$url" "$tmp"
  tar -xzf "$tmp" -C "${INSTALL_DIR}"
  rm -f "$tmp"
  mv "${INSTALL_DIR}/ncbi-blast-${BLAST_VER}+" "${INSTALL_DIR}/ncbi-blast"
  ln -sf "${INSTALL_DIR}/ncbi-blast/bin/"* "${BIN_DIR}/"
  echo "   BLAST+ installed → ${INSTALL_DIR}/ncbi-blast"
}

### --------------------------
### 2) SPAdes 4.2.0 (Linux/macOS)
### --------------------------
install_spades () {
  echo "==> Installing SPAdes ${SPADES_VER} ..."
  local tar=""
  local url=""
  if [ "$OS_TYPE" = "Linux" ]; then
    tar="SPAdes-${SPADES_VER}-Linux.tar.gz"
    url="https://github.com/ablab/spades/releases/download/v${SPADES_VER}/${tar}"
  elif [ "$OS_TYPE" = "Darwin" ]; then
    tar="SPAdes-${SPADES_VER}-Darwin.tar.gz"
    url="https://github.com/ablab/spades/releases/download/v${SPADES_VER}/${tar}"
    if [ "${ARCH_TYPE}" = "arm64" ]; then
      echo "   [Note] macOS arm64 detected. SPAdes Darwin tarball may be x86_64."
      echo "          If binaries fail, build from source or run via Rosetta."
    fi
  else
    echo "Unsupported OS for SPAdes: $OS_TYPE"; return 1
  fi

  local tmp="${PWD}/${tar}"
  fetch "$url" "$tmp"
  tar -xzf "$tmp" -C "${INSTALL_DIR}"
  rm -f "$tmp"
  # Move/rename to stable path
  local src_dir_linux="${INSTALL_DIR}/SPAdes-${SPADES_VER}-Linux"
  local src_dir_mac="${INSTALL_DIR}/SPAdes-${SPADES_VER}-Darwin"
  if [ -d "$src_dir_linux" ]; then
    mv "$src_dir_linux" "${INSTALL_DIR}/SPAdes"
  elif [ -d "$src_dir_mac" ]; then
    mv "$src_dir_mac" "${INSTALL_DIR}/SPAdes"
  fi
  ln -sf "${INSTALL_DIR}/SPAdes/bin/"* "${BIN_DIR}/"
  echo "   SPAdes installed → ${INSTALL_DIR}/SPAdes"
  echo "   Try: spades.py --test"
}

### ------------------------------------
### 3) Immcantation helper scripts
### ------------------------------------
install_immcantation_scripts () {
  echo "==> Installing Immcantation helper scripts ..."
  local tmp_dir="${INSTALL_DIR}/immcantation_tmp"
  rm -rf "$tmp_dir"
  git clone --depth 1 "${IMMCANTATION_REPO}" "$tmp_dir" >/dev/null 2>&1
  mkdir -p "${INSTALL_DIR}/immcantation_scripts"
  cp -r "$tmp_dir/scripts/"* "${INSTALL_DIR}/immcantation_scripts/"
  rm -rf "$tmp_dir"
  chmod +x "${INSTALL_DIR}/immcantation_scripts/"* || true
  ln -sf "${INSTALL_DIR}/immcantation_scripts/"* "${BIN_DIR}/"
  echo "   Immcantation scripts → ${INSTALL_DIR}/immcantation_scripts"
}

### --------------------------
### 4) IgPhyML
###   - Linux: compile (OpenMP build)
###   - macOS: best-effort via Homebrew LLVM
###   - If macOS compile is painful, print Docker hint
### --------------------------
install_igphyml_linux () {
  echo "==> Installing IgPhyML (Linux, compile) ..."
  local dir="${INSTALL_DIR}/igphyml"
  rm -rf "$dir"
  git clone --depth 1 "${IGPHYML_REPO}" "$dir" >/dev/null 2>&1
  ( cd "$dir" && chmod +x make_phyml_omp && ./make_phyml_omp )
  # Primary binary is src/igphyml; link both igphyml and phyml if present
  [ -x "${dir}/src/igphyml" ] && ln -sf "${dir}/src/igphyml" "${BIN_DIR}/igphyml"
  [ -x "${dir}/src/phyml" ]   && ln -sf "${dir}/src/phyml"   "${BIN_DIR}/phyml"
  echo "   IgPhyML built → ${dir}/src"
}

install_igphyml_macos () {
  echo "==> Installing IgPhyML (macOS, best-effort compile) ..."
  if ! need_cmd brew; then
    echo "   Homebrew not found. Install from https://brew.sh and re-run."
    echo "   Alternatively, use Docker (recommended on macOS):"
    echo "     docker pull immcantation/suite:4.5.0"
    return 1
  fi
  # Prereqs for build (non-sudo, user-space)
  brew list autoconf >/dev/null 2>&1 || brew install autoconf
  brew list automake >/dev/null 2>&1 || brew install automake
  brew list llvm     >/dev/null 2>&1 || brew install llvm

  local LLVM_PREFIX
  LLVM_PREFIX="$(brew --prefix llvm)"
  local CC_BIN="${LLVM_PREFIX}/bin/clang"
  # Find omp.h inside LLVM tree
  local OMP_H
  OMP_H="$(/usr/bin/find "${LLVM_PREFIX}/lib/clang" -type f -name omp.h 2>/dev/null | head -n1 || true)"
  local LLVM_LIB="${LLVM_PREFIX}/lib"

  if [ ! -x "$CC_BIN" ] || [ -z "$OMP_H" ]; then
    echo "   Could not locate LLVM clang or omp.h. Check your Homebrew LLVM install."
    echo "   LLVM prefix: ${LLVM_PREFIX}"
    echo "   Try: brew reinstall llvm"
    return 1
  fi

  local dir="${INSTALL_DIR}/igphyml"
  rm -rf "$dir"
  git clone --depth 1 "${IGPHYML_REPO}" "$dir" >/dev/null 2>&1

  # Prepend settings to Makefile.am files if not already present
  prepend_if_missing () {
    local file="$1"
    local text="$2"
    if ! grep -Fq "$text" "$file"; then
      # Prepend (portable)
      tmpfile="$(mktemp)"
      printf "%s\n%s" "$text" "$(cat "$file")" > "$tmpfile"
      mv "$tmpfile" "$file"
    fi
  }

  local header_main="CC=${CC_BIN}"
  local header_src=$(cat <<EOF
CC=${CC_BIN}
MACOMP=${OMP_H}
MACLLVM=${LLVM_LIB}
EOF
)

  prepend_if_missing "${dir}/Makefile.am" "${header_main}"
  prepend_if_missing "${dir}/src/Makefile.am" "${header_src}"

  ( cd "$dir" && chmod +x make_phyml_omp && ./make_phyml_omp || true )

  if [ -x "${dir}/src/igphyml" ]; then
    ln -sf "${dir}/src/igphyml" "${BIN_DIR}/igphyml"
  fi
  if [ -x "${dir}/src/phyml" ]; then
    ln -sf "${dir}/src/phyml" "${BIN_DIR}/phyml"
  fi

  if ! command -v igphyml >/dev/null 2>&1; then
    echo "   IgPhyML binary was not produced. On macOS, Docker is recommended:"
    echo "     docker run -it --workdir /data -v \"\$(pwd)\":/data:z immcantation/suite:4.5.0 bash"
  else
    echo "   IgPhyML built → ${dir}/src"
  fi
}

install_igphyml () {
  case "$OS_TYPE" in
    Linux)  install_igphyml_linux ;;
    Darwin) install_igphyml_macos ;;
    *) echo "Unsupported OS for IgPhyML: $OS_TYPE" ;;
  esac
}

### --------------------------
### Run installers
### --------------------------
install_blast
install_spades
install_immcantation_scripts
install_igphyml
append_path_once

echo
echo "=============================================="
echo " All done. Binaries linked into: ${BIN_DIR}"
echo " New shells will pick up PATH automatically."
echo
echo " Quick checks:"
echo "   - blastn -version"
echo "   - spades.py --test"
echo "   - igphyml --version   (or see Docker note on macOS)"
echo "=============================================="
