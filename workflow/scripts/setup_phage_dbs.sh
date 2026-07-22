#!/usr/bin/env bash
#
# Download the six databases the phage workflow needs (~20 GB total).
#
# These are deliberately NOT wired into the Snakemake DAG: they are large,
# slow, and shared between projects, so they belong in one place on the cluster
# rather than being re-fetched per run. Run this once, then point the
# `databases:` block of your phage config at the resulting directories.
#
# Each download runs in its tool's own conda environment, so create those first:
#
#   for e in genomad cenotetaker checkv pharokka phold phynteny; do
#       conda env create -n "$e" -f workflow/envs/"$e".yaml
#   done
#
# Usage:
#   bash workflow/scripts/setup_phage_dbs.sh /path/to/phage_dbs [tool ...]
#
# With no tool names, all six are downloaded. Otherwise only the named ones,
# which is useful for resuming after a failure.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    sed -n '2,22p' "$0"
    exit 1
fi

DB_ROOT="$1"
shift
TOOLS=("$@")
if [[ ${#TOOLS[@]} -eq 0 ]]; then
    TOOLS=(genomad cenotetaker checkv pharokka phold phynteny)
fi

mkdir -p "$DB_ROOT"
echo "Installing phage databases under: $DB_ROOT"

run_in_env() {
    # Run a command inside a named conda env without requiring `conda activate`
    # to work in a non-interactive shell.
    local env_name="$1"; shift
    conda run --no-capture-output -n "$env_name" "$@"
}

for tool in "${TOOLS[@]}"; do
    echo ""
    echo "=== $tool"
    case "$tool" in
        genomad)
            run_in_env genomad genomad download-database "$DB_ROOT"
            echo "    -> $DB_ROOT/genomad_db"
            ;;
        cenotetaker)
            run_in_env cenotetaker get_ct3_dbs -o "$DB_ROOT/ct3_dbs" \
                --hmm T --hallmark_tax T --refseq_tax T --mmseqs_cdd T --domain_list T
            echo "    -> $DB_ROOT/ct3_dbs"
            ;;
        checkv)
            run_in_env checkv checkv download_database "$DB_ROOT"
            echo "    -> $DB_ROOT/checkv-db-v* (note the version suffix)"
            ;;
        pharokka)
            run_in_env pharokka install_databases.py -o "$DB_ROOT/pharokka_db"
            echo "    -> $DB_ROOT/pharokka_db"
            ;;
        phold)
            run_in_env phold phold install -d "$DB_ROOT/phold_db"
            echo "    -> $DB_ROOT/phold_db"
            ;;
        phynteny)
            run_in_env phynteny install_models -o "$DB_ROOT/phynteny_models"
            echo "    -> $DB_ROOT/phynteny_models"
            ;;
        *)
            echo "Unknown tool: $tool" >&2
            exit 1
            ;;
    esac
done

cat <<EOF

Done. Add the resulting paths to your phage config, for example:

databases:
  genomad: $DB_ROOT/genomad_db
  cenotetaker: $DB_ROOT/ct3_dbs
  checkv: $DB_ROOT/checkv-db-v1.5
  pharokka: $DB_ROOT/pharokka_db
  phold: $DB_ROOT/phold_db
  phynteny: $DB_ROOT/phynteny_models

CheckV's directory carries a version suffix that changes between releases, so
check the actual name before copying the line above.
EOF
