#!/usr/bin/env bash
#
# One-time setup. Run it from the repository root:
#
#     bash install.sh              # create the three conda environments, then install gtftools
#     bash install.sh gtftools     # just gtftools
#     bash install.sh envs         # just the environments
#
# Re-running is safe: an environment that already exists is left alone, and gtftools is
# reinstalled over itself.
#
# gtftools is the one dependency that is not on conda-forge or bioconda, so
# envs/processing.yml cannot cover it and `run.sh annotation` needs it. It is a single
# MIT-licensed Python script; this installs it as `gtftools` on the environment's PATH.

set -euo pipefail

[[ -d processing_scripts && -f README.md ]] || {
    echo "ERROR: run this from the repository root (see README)." >&2; exit 2; }

GTFTOOLS_URL="${GTFTOOLS_URL:-https://www.genemine.org/codes/GTFtools_0.9.0.zip}"

install_envs() {
    command -v conda >/dev/null || {
        echo "ERROR: conda not on PATH. Install miniforge/miniconda first." >&2; exit 1; }
    local yml name
    for yml in envs/notebooks.yml envs/processing.yml envs/r.yml; do
        name="$(awk '/^name:/{print $2; exit}' "$yml")"
        if conda env list | awk '{print $1}' | grep -qx "$name"; then
            echo "== $name already exists, skipping" >&2
        else
            echo "== creating $name from $yml" >&2
            conda env create -f "$yml"
        fi
    done
}

install_gtftools() {
    # Default to the processing environment's bin, since that is where the annotation
    # stage runs. Override with GTFTOOLS_BIN=<dir> to put it anywhere else on your PATH.
    local bindir
    if [[ -n "${GTFTOOLS_BIN:-}" ]]; then
        bindir="$GTFTOOLS_BIN"
    elif command -v conda >/dev/null \
         && [[ -d "$(conda info --base)/envs/h3k27me2-processing/bin" ]]; then
        bindir="$(conda info --base)/envs/h3k27me2-processing/bin"
    else
        bindir="$HOME/.local/bin"
        echo "NOTE: h3k27me2-processing not found; installing to $bindir instead." >&2
        echo "      Make sure that directory is on your PATH." >&2
    fi
    mkdir -p "$bindir"

    command -v unzip >/dev/null || {
        echo "ERROR: unzip not on PATH; it is needed to unpack the gtftools archive." >&2; exit 1; }

    local tmp; tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' RETURN
    echo "== downloading $GTFTOOLS_URL" >&2
    if command -v curl >/dev/null; then
        curl -fsSL --retry 3 -o "$tmp/gtftools.zip" "$GTFTOOLS_URL"
    else
        wget -q -O "$tmp/gtftools.zip" "$GTFTOOLS_URL"
    fi
    unzip -q -o "$tmp/gtftools.zip" -d "$tmp/unpacked"

    local src; src="$(find "$tmp/unpacked" -name gtftools.py -print -quit)"
    [[ -n "$src" ]] || { echo "ERROR: gtftools.py not found in the archive." >&2; exit 1; }

    # The shipped shebang is #!/usr/bin/python, which is python2 or absent on most systems.
    # gtftools needs python3 + numpy, both of which envs/processing.yml provides.
    { echo '#!/usr/bin/env python3'; tail -n +2 "$src"; } > "$bindir/gtftools"
    chmod +x "$bindir/gtftools"
    echo "== installed $bindir/gtftools" >&2
    "$bindir/gtftools" --help >/dev/null \
      && echo "== gtftools runs" >&2 \
      || { echo "ERROR: installed gtftools does not run (needs python3 + numpy)." >&2; exit 1; }
    # Also copy the upstream licence next to it, since gtftools is third-party MIT code.
    local lic; lic="$(find "$tmp/unpacked" -name 'license*' -print -quit)"
    [[ -n "$lic" ]] && cp "$lic" "$bindir/gtftools.LICENSE"
    return 0
}

case "${1:-all}" in
  all)       install_envs; install_gtftools ;;
  envs)      install_envs ;;
  gtftools)  install_gtftools ;;
  -h|--help|help)
    sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//' ;;
  *) echo "usage: bash install.sh [all|envs|gtftools]" >&2; exit 2 ;;
esac
