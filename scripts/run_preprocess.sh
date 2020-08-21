snakemake -s rules/preprocess.smk --configfile configs/test_config.yaml --use-conda -k --cluster "DIR=\$(dirname {params.qoutfile}); mkdir -p \"\${{DIR}}\"; qsub -S /bin/bash -V -cwd -o {params.qoutfile} -e {params.qerrfile} -pe smp {threads} -l h_vmem={params.mem}M" -p -j 4 --max-jobs-per-second 1