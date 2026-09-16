#!/usr/bin/env bash
# Copy the repo's deliverables into OneDrive, where the paper is compiled
# and where the older analysis scripts live.
#
#   the whole repo       ->  Work/2026_Jpsi_radiative_corrections/  (full mirror,
#                            so the paper and code are reachable from any machine)
#   rc_vector_mesons.py  ->  Work/Python_Lib/RC_vector_mesons.py   (the library)
#   paper/Jpsi_RC.tex    ->  Tex/2024_rad_corr_paper/
#   paper/Jpsi_RC.pdf    ->  Tex/2024_rad_corr_paper/   (the built paper)
#   paper/Figures/*.pdf  ->  Tex/2024_rad_corr_paper/Figures/
#
# The mirror deliberately omits .git: OneDrive syncing a git directory can
# corrupt it.  To get the repo itself on another machine, clone it from
# github.com/vkubarovsky/jpsi-rc-2026.
#
# The repo is the source of truth; OneDrive is a publication target.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OD="$HOME/OneDrive"
LIB="$OD/Work/Python_Lib"
TEX="$OD/Tex/2024_rad_corr_paper"

[ -d "$LIB" ] || { echo "missing $LIB" >&2; exit 1; }
[ -d "$TEX" ] || { echo "missing $TEX" >&2; exit 1; }

MIRROR="$OD/Work/2026_Jpsi_radiative_corrections"
echo "mirror  -> $MIRROR"
mkdir -p "$MIRROR"
rsync -a --delete --exclude '.git' --exclude '__pycache__' --exclude '.DS_Store' \
      "$REPO"/ "$MIRROR"/

echo "library -> $LIB"
cp "$REPO/rc_vector_mesons.py" "$LIB/RC_vector_mesons.py"
cp "$REPO/validate_rc.py"      "$LIB/validate_rc.py"

echo "paper   -> $TEX"
cp "$REPO/paper/Jpsi_RC.tex" "$TEX/"
[ -f "$REPO/paper/Jpsi_RC.pdf" ] && cp "$REPO/paper/Jpsi_RC.pdf" "$TEX/"
mkdir -p "$TEX/Figures"
cp "$REPO"/paper/Figures/*.pdf "$TEX/Figures/"

echo "done"
