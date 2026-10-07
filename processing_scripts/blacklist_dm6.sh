# Builds data/Kc_merged_blacklist_hg38.bed, the Drosophila spike-in blacklist. That output
# ALREADY SHIPS in this repository, so you do not need to run this — which is just as well,
# because its four inputs are lab-internal and are not redistributed here. Pass them as
# arguments.
#
# KNOWN DISCREPANCY: which Drosophila cell line the spike-in came from.
# The manuscript Methods say "hg38 positions mapped to by reads from S2 profiles were
# blacklisted", but the four inputs below are named for Kc and BT_Dm, not S2. Both Kc167
# and S2 are Drosophila lines and the blacklisting method is identical either way, so this
# does not affect the shipped blacklist or anything downstream of it — but the cell line
# named in the manuscript may not be the one these profiles came from. The shipped
# blacklist is the authoritative artefact; the naming is unreconciled.

# Usage (from the repository root):
#     bash processing_scripts/blacklist_dm6.sh <dm1.bed> <dm2.bed> <dm3.bed> <dm4.bed>
# Each input is a BED of hg38 positions hit by reads from a Drosophila-only CUT&Tag
# profile (one Kc167 profile and three BT_Dm profiles in the original build).
[[ $# -eq 4 ]] || {
    echo "usage: bash processing_scripts/blacklist_dm6.sh <dm1.bed> <dm2.bed> <dm3.bed> <dm4.bed>" >&2
    exit 2; }
Kc1=$1
Kc2=$2
Kc3=$3
Kc4=$4
OUT=data
mkdir -p "$OUT"
command -v bedtools >/dev/null || module load BEDTools/2.31.0-GCC-12.3.0

# Clean, combine, sort, and merge
cat "$Kc1" "$Kc2" "$Kc3" "$Kc4" | \
  awk 'NF>=3 {print $1, $2, $3}' OFS='\t' | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i stdin | \
  awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $3 - $2}' > "${OUT}/Kc_merged_blacklist_hg38.bed"