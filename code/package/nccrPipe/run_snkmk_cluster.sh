#!/bin/bash
if [ $# -eq 0 ]
  then
      echo 'No arguments supplied'
      exit 1
fi

snakemake -s Snakefile_test --use-conda -k --cluster "DIR=\$(dirname {params.qoutfile}); mkdir -p \"\${{DIR}}\"; qsub -S /bin/bash -V -cwd -o {params.qoutfile} -e {params.qerrfile} -pe smp {threads} -l h_vmem={params.mem}M" -p -j 4 --max-jobs-per-second 1 $1

