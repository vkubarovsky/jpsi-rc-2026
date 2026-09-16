#!/usr/bin/env bash
# Copy the repo's deliverables into OneDrive, where the paper is compiled
# and where the older analysis scripts live.
#
#   rc_vector_mesons.py  ->  Work/Python_Lib/RC_vector_mesons.py   (the library)
#   paper/Jpsi_RC.tex    ->  Tex/2024_rad_corr_paper/
#   paper/Figures/*.pdf  ->  Tex/2024_rad_corr_paper/Figures/
#
# The repo is the source of truth; OneDrive is a publication target.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OD="$HOME/OneDrive"
LIB="$OD/Work/Python_Lib"
TEX="$OD/Tex/2024_rad_corr_paper"

[ -d "$LIB" ] || { echo "missing $LIB" >&2; exit 1; }
[ -d "$TEX" ] || { echo "missing $TEX" >&2; exit 1; }

echo "library -> $LIB"
cp "$REPO/rc_vector_mesons.py" "$LIB/RC_vector_mesons.py"
cp "$REPO/validate_rc.py"      "$LIB/validate_rc.py"

echo "paper   -> $TEX"
cp "$REPO/paper/Jpsi_RC.tex" "$TEX/"
mkdir -p "$TEX/Figures"
cp "$REPO"/paper/Figures/*.pdf "$TEX/Figures/"

echo "done"
