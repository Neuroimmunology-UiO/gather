#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Gather dependency installer (Linux & macOS)
# - Installs: BLAST+, SPAdes, Immcantation helper scripts, IgPhyML
# - Installs R deps: BiocManager, treeio, ggtree, dowser, dplyr, ggrepel
# - No sudo; installs under $HOME by default
# - Excludes: IgBLAST/IMGT and any Python packages (handled by conda)
# ============================================================

### ---------- Versions (edit if needed) ----------
BLAST_VER="${BLAST_VER:-2.16.0}"
SPADES_VER="${SPADES_VER:-4.2.0}"
IGPHYML_REPO="${IGPHYML_REPO:-https://github.com/immcantation/igphyml.git}"
IMMCANTATION_REPO="${IMMCANTATION_REPO:-https://github.com/immcantation/immcantation.git}"

### ---------- Install locations ----------
PREFIX="${PREFIX:-$HOME/opt}"
BIN_DIR="${PREFIX_BIN:-$HOME/bin}"
mkdir -p "$PREFIX" "$BIN_DIR"

OS_TYPE="$(uname -s)"
ARCH_TYPE="$(uname -m)"

echo "==> Gather deps installer"
echo "    OS: ${OS_TYPE} ${ARCH_TYPE}"
echo "    PREFIX: ${PREFIX}"
echo "    BIN: ${BIN_DIR}"

need_cmd() { command -v "$1" >/dev/null 2>&1; }
fetch() {
  # $1=url $2=out
  if need_cmd wget; then
    wget -q -O "$2" "$1"
  elif need_cmd curl; then
    curl -LfsS -o "$2" "$1"
  else
    echo "ERROR: Need 'wget' or 'curl' to download files." >&2
    exit 1
  fi
}
append_path_once () {
  local line='export PATH="${HOME}/bin:$PATH"'
  for rc in "${HOME}/.bashrc" "${HOME}/.zshrc"; do
    [ -f "$rc" ] || continue
    grep -Fq "$line" "$rc" || echo "$line" >> "$rc"
  done
  export PATH="${HOME}/bin:$PATH"
}

### ------------------------------------------------
### 1) BLAST+ (Linux/macOS)
### ------------------------------------------------
install_blast () {
  echo "==> Installing NCBI BLAST+ ${BLAST_VER} ..."
  local tar url
  case "$OS_TYPE" in
    Linux)  tar="ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz" ;;
    Darwin)
      tar="ncbi-blast-${BLAST_VER}+-x64-macosx.tar.gz"
      if [ "${ARCH_TYPE}" = "arm64" ]; then
        echo "   [Note] macOS arm64 detected. BLAST+ tarball is x86_64."
        echo "          If binaries fail under Rosetta, consider: conda install -c bioconda blast"
      fi
      ;;
    *) echo "Unsupported OS for BLAST+: $OS_TYPE"; return 1 ;;
  esac
  url="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VER}/${tar}"
  tmp="${PWD}/${tar}"
  fetch "$url" "$tmp"
  tar -xzf "$tmp" -C "$PREFIX"
  rm -f "$tmp"
  mv "${PREFIX}/ncbi-blast-${BLAST_VER}+" "${PREFIX}/ncbi-blast"
  ln -sf "${PREFIX}/ncbi-blast/bin/"* "$BIN_DIR/"
  echo "   BLAST+ → ${PREFIX}/ncbi-blast"
}

### ------------------------------------------------
### 2) SPAdes (Linux/macOS)
### ------------------------------------------------
install_spades () {
  echo "==> Installing SPAdes ${SPADES_VER} ..."
  local tar url
  if [ "$OS_TYPE" = "Linux" ]; then
    tar="SPAdes-${SPADES_VER}-Linux.tar.gz"
  elif [ "$OS_TYPE" = "Darwin" ]; then
    tar="SPAdes-${SPADES_VER}-Darwin.tar.gz"
    if [ "${ARCH_TYPE}" = "arm64" ]; then
      echo "   [Note] macOS arm64 detected. Darwin tarball may be x86_64."
      echo "          If binaries fail, build from source or run via Rosetta."
    fi
  else
    echo "Unsupported OS for SPAdes: $OS_TYPE"; return 1
  fi
  url="https://github.com/ablab/spades/releases/download/v${SPADES_VER}/${tar}"
  tmp="${PWD}/${tar}"
  fetch "$url" "$tmp"
  tar -xzf "$tmp" -C "$PREFIX"
  rm -f "$tmp"
  if [ -d "${PREFIX}/SPAdes-${SPADES_VER}-Linux" ]; then
    mv "${PREFIX}/SPAdes-${SPADES_VER}-Linux" "${PREFIX}/SPAdes"
  elif [ -d "${PREFIX}/SPAdes-${SPADES_VER}-Darwin" ]; then
    mv "${PREFIX}/SPAdes-${SPADES_VER}-Darwin" "${PREFIX}/SPAdes"
  fi
  ln -sf "${PREFIX}/SPAdes/bin/"* "$BIN_DIR/"
  echo "   SPAdes → ${PREFIX}/SPAdes (try: spades.py --test)"
}

### ------------------------------------------------
### 3) Immcantation helper scripts
### ------------------------------------------------
install_immcantation_scripts () {
  echo "==> Installing Immcantation helper scripts ..."
  local dst="${PREFIX}/immcantation_scripts"
  rm -rf "$dst"
  mkdir -p "$dst"
  if [ -n "${IMMCANTATION_SCRIPTS_LOCAL:-}" ]; then
    echo "   Copying from local: $IMMCANTATION_SCRIPTS_LOCAL"
    cp -r "${IMMCANTATION_SCRIPTS_LOCAL}/"* "$dst/"
  else
    echo "   Cloning from ${IMMCANTATION_REPO}"
    local tmp="${PREFIX}/immcantation_tmp"
    rm -rf "$tmp"
    git clone --depth 1 "$IMMCANTATION_REPO" "$tmp" >/dev/null 2>&1
    cp -r "$tmp/scripts/"* "$dst/"
    rm -rf "$tmp"
  fi
  chmod +x "$dst/"* || true
  ln -sf "$dst/"* "$BIN_DIR/"
  echo "   Immcantation scripts → ${dst}"
}

### ------------------------------------------------
### 4) IgPhyML
###     - Linux: compile (OpenMP build if available)
###     - macOS: best-effort via Homebrew LLVM
### ------------------------------------------------
install_igphyml_linux () {
  echo "==> Installing IgPhyML (Linux compile) ..."
  local dir="${PREFIX}/igphyml"
  rm -rf "$dir"
  git clone --depth 1 "$IGPHYML_REPO" "$dir" >/dev/null 2>&1
  ( cd "$dir" && chmod +x make_phyml_omp && ./make_phyml_omp || ./make_phyml || true )
  [ -x "${dir}/src/igphyml" ] && ln -sf "${dir}/src/igphyml" "${BIN_DIR}/igphyml"
  [ -x "${dir}/src/phyml" ]   && ln -sf "${dir}/src/phyml"   "${BIN_DIR}/phyml"
  if ! command -v igphyml >/dev/null 2>&1; then
    echo "   IgPhyML build failed. Install toolchain (autoconf, automake, OpenMP) or use Docker:"
    echo "     docker run -it --workdir /data -v \"\$(pwd)\":/data:z immcantation/suite:4.5.0 bash"
  else
    echo "   IgPhyML built → ${dir}/src"
  fi
}

install_igphyml_macos () {
  echo "==> Installing IgPhyML (macOS best-effort) ..."
  if ! need_cmd brew; then
    echo "   Homebrew not found. Install from https://brew.sh or use Docker:"
    echo "     docker run -it --workdir /data -v \"\$(pwd)\":/data:z immcantation/suite:4.5.0 bash"
    return 1
  fi
  brew list autoconf >/dev/null 2>&1 || brew install autoconf
  brew list automake >/dev/null 2>&1 || brew install automake
  brew list llvm     >/dev/null 2>&1 || brew install llvm

  local LLVM_PREFIX CC_BIN OMP_H LLVM_LIB
  LLVM_PREFIX="$(brew --prefix llvm)"
  CC_BIN="${LLVM_PREFIX}/bin/clang"
  OMP_H="$(/usr/bin/find "${LLVM_PREFIX}/lib/clang" -type f -name omp.h 2>/dev/null | head -n1 || true)"
  LLVM_LIB="${LLVM_PREFIX}/lib"

  if [ ! -x "$CC_BIN" ] || [ -z "$OMP_H" ]; then
    echo "   Could not locate LLVM clang or omp.h. Try: brew reinstall llvm"
    return 1
  fi

  local dir="${PREFIX}/igphyml"
  rm -rf "$dir"
  git clone --depth 1 "$IGPHYML_REPO" "$dir" >/dev/null 2>&1

  prepend_if_missing () {
    local file="$1" text="$2"
    if ! grep -Fq "$text" "$file"; then
      local tmpf; tmpf="$(mktemp)"
      printf "%s\n%s" "$text" "$(cat "$file")" > "$tmpf"
      mv "$tmpf" "$file"
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

  ( cd "$dir" && chmod +x make_phyml_omp && ./make_phyml_omp || ./make_phyml || true )
  [ -x "${dir}/src/igphyml" ] && ln -sf "${dir}/src/igphyml" "${BIN_DIR}/igphyml"
  [ -x "${dir}/src/phyml" ]   && ln -sf "${dir}/src/phyml"   "${BIN_DIR}/phyml"

  if ! command -v igphyml >/dev/null 2>&1; then
    echo "   IgPhyML not built. Recommended on macOS:"
    echo "     docker pull immcantation/suite:4.5.0"
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

### ------------------------------------------------
### 5) R packages (Bioc + CRAN)
###     - Uses current Rscript or a Conda env (GATHER_CONDA_ENV)
### ------------------------------------------------
install_r_deps () {
  local CONDA_ENV="${GATHER_CONDA_ENV:-}"
  local RSCRIPT_CMD=()

  echo "==> Installing R dependencies (BiocManager, treeio, ggtree, dowser, dplyr, ggrepel)"

  if [ -n "$CONDA_ENV" ] && need_cmd conda; then
    if ! conda run -n "$CONDA_ENV" Rscript -e "q()" >/dev/null 2>&1; then
      echo "   Conda env '$CONDA_ENV' not ready; creating with r-base + r-essentials..."
      if need_cmd mamba; then
        mamba create -y -n "$CONDA_ENV" -c conda-forge r-base r-essentials
      else
        conda create -y -n "$CONDA_ENV" -c conda-forge r-base r-essentials
      fi
    fi
    RSCRIPT_CMD=(conda run -n "$CONDA_ENV" Rscript)
  elif need_cmd Rscript; then
    RSCRIPT_CMD=(Rscript)
  elif need_cmd conda; then
    CONDA_ENV="${CONDA_ENV:-gather_r}"
    echo "   No system Rscript; creating conda env '$CONDA_ENV' with r-base + r-essentials..."
    if need_cmd mamba; then
      mamba create -y -n "$CONDA_ENV" -c conda-forge r-base r-essentials
    else
      conda create -y -n "$CONDA_ENV" -c conda-forge r-base r-essentials
    fi
    RSCRIPT_CMD=(conda run -n "$CONDA_ENV" Rscript)
  else
    echo "   ERROR: No Rscript and no conda found. Install R or Conda and re-run." >&2
    return 1
  fi

  echo "   Using: ${RSCRIPT_CMD[*]}"

  "${RSCRIPT_CMD[@]}" - <<'RS_EOF'
options(warn=1)
repos <- getOption("repos")
if (is.null(repos) || is.na(repos["CRAN"]) || repos["CRAN"] == "@CRAN@") {
  repos["CRAN"] <- "https://cloud.r-project.org"; options(repos=repos)
}
lib_user <- Sys.getenv("R_LIBS_USER")
if (!nzchar(lib_user)) {
  lib_user <- file.path(path.expand("~"), "R", paste0("library-", paste(R.version$major, R.version$minor, sep=".")))
  Sys.setenv(R_LIBS_USER = lib_user)
}
if (!dir.exists(lib_user)) dir.create(lib_user, recursive=TRUE, showWarnings = FALSE)
.libPaths(c(lib_user, .libPaths()))
message("LibPaths: ", paste(.libPaths(), collapse=" | "))

ncpus <- 1L
try(ncpus <- parallel::detectCores(), silent = TRUE)

cran_pkgs <- c("dplyr", "ggrepel")
bioc_pkgs <- c("treeio", "ggtree", "dowser")

install_cran <- function(pkgs) {
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(need)) {
    message("Installing CRAN: ", paste(need, collapse=", "))
    install.packages(need, Ncpus = ncpus)
  } else {
    message("CRAN packages already installed: ", paste(pkgs, collapse=", "))
  }
}

install_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) {
    message("Installing BiocManager...")
    install.packages("BiocManager", Ncpus = ncpus)
  }
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(need)) {
    message("Installing Bioconductor: ", paste(need, collapse=", "))
    BiocManager::install(need, ask=FALSE, update=TRUE)
  } else {
    message("Bioconductor packages already installed: ", paste(pkgs, collapse=", "))
  }
}

install_cran(cran_pkgs)
install_bioc(bioc_pkgs)

cat("\nInstalled versions:\n")
for (p in c(cran_pkgs, bioc_pkgs, "BiocManager")) {
  if (requireNamespace(p, quietly=TRUE)) {
    cat(sprintf("  %-12s %s\n", p, as.character(utils::packageVersion(p))))
  } else {
    cat(sprintf("  %-12s (NOT INSTALLED)\n", p))
  }
}
RS_EOF
}

### ------------------------------------------------
### Run selected steps
### ------------------------------------------------
append_path_once

: "${SKIP_BLAST:=0}"
: "${SKIP_SPADES:=0}"
: "${SKIP_IMMCANTATION:=0}"
: "${SKIP_IGPHYML:=0}"
: "${SKIP_R:=0}"

[ "$SKIP_BLAST" -eq 1 ]        || install_blast
[ "$SKIP_SPADES" -eq 1 ]       || install_spades
[ "$SKIP_IMMCANTATION" -eq 1 ] || install_immcantation_scripts
[ "$SKIP_IGPHYML" -eq 1 ]      || install_igphyml
[ "$SKIP_R" -eq 1 ]            || install_r_deps

echo
echo "=============================================="
echo " Done. Binaries are symlinked to: ${BIN_DIR}"
echo " New shells will pick up PATH automatically."
echo
echo " Quick checks:"
echo "   - blastn -version"
echo "   - spades.py --test"
echo "   - igphyml --version   (or see Docker note on macOS)"
echo "   - R -q -e \"library(ggtree); library(dowser)\""
echo "=============================================="
