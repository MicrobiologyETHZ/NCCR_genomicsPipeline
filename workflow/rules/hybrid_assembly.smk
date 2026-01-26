

rule qc_short_reads:
    input:
        fq1 = f"{DATADIR}/{{sample}}{SHORT_READ_SUFFIX_R1}",
        fq2 = f"{DATADIR}/{{sample}}{SHORT_READ_SUFFIX_R2}",
        adapters = Path(config['adapters']),
        phix = Path(config['phix'])
    output:
        fq1_clean = OUTDIR /'clean_reads/{sample}/{sample}.1.fq.gz',
        fq2_clean = OUTDIR /'clean_reads/{sample}/{sample}.2.fq.gz',
        adapter_matched = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.adapter.matched.fq.gz',
        adapter_singletons = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.adapter.singletons.fq.gz',
        adapter_stats = OUTDIR /'clean_reads/{sample}/{sample}.adapter.stats',
        phix_matched = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.phix.matched.fq.gz',
        phix_singletons = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.phix.singletons.fq.gz',
        phix_stats = OUTDIR /'clean_reads/{sample}/{sample}.phix.stats',
        qc_failed = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.qc.failed.fq.gz',
        qc_singletons = OUTDIR /'clean_reads/{sample}/{sample}.s.fq.gz',
        qc_stats = OUTDIR /'clean_reads/{sample}/{sample}.qc.stats',
        marker = touch(OUTDIR /'clean_reads/{sample}/{sample}.qc.done')
    params:
        trimq = config['trimq'],
        maq = config['mapq'],
        minlen = config['minlen'],
        qoutfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qout',
        qerrfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qerr',
        scratch = 500,
        mem = 8000,
        time = 235
    conda:
        "preprocessing"

    log:
        log = OUTDIR/'logs/qc/{sample}.qc.log'
    threads:
        8
    shell:
        "bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t "
        "in={input.fq1} in2={input.fq2} "
        "out=stdout.fq outm={output.adapter_matched} "
        "outs={output.adapter_singletons} "
        "refstats={output.adapter_stats} statscolumns=5 "
        "overwrite=t ref={input.adapters} "
        "ktrim=r k=23 mink=11 hdist=1  2>> {log.log} | "
        "bbduk.sh -Xmx1G usejni=t pigz=t bgzip=f "
        "interleaved=true overwrite=t "
        "in=stdin.fq out=stdout.fq "
        "outm={output.phix_matched} outs={output.phix_singletons} "
        "ref={input.phix} k=31 hdist=1 "
        "refstats={output.phix_stats} statscolumns=5 2>> {log.log}| "
        "bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t  "
        "overwrite=t interleaved=true "
        "in=stdin.fq fastawrap=10000 "
        "out1={output.fq1_clean} out2={output.fq2_clean} "
        "outm={output.qc_failed} outs={output.qc_singletons} "
        "minlength={params.minlen} qtrim=rl maq={params.maq} maxns=1  "
        "stats={output.qc_stats} statscolumns=5 "
        "trimq={params.trimq}  2>> {log.log};"


rule qc_hifi_reads:
    input:
        long = f"{LONG_DATA_DIR}/{{sample}}{LONG_SUFFIX}"
    output:
        long_clean = f"{OUTDIR}/{{sample}}/{{sample}}_hifi_clean.fastq.gz",
        stats = f"{OUTDIR}/{{sample}}/hifi_stats.txt"
    params:
        qoutfile = lambda wildcards: OUTDIR /f'logs/{wildcards.sample}.qc_hifi.qout',
        qerrfile = lambda wildcards: OUTDIR /f'logs/{wildcards.sample}.qc_hifi.qerr',
        scratch = 2000,
        mem = 8000,
        time = 1400
    conda:
        "metaflye"
    threads: 4
    shell:
        """
        # Minimal filtering for HiFi - just remove very short reads
        seqtk seq -L 5000 {input.long} | gzip > {output.long_clean}
        
        # Generate comprehensive statistics
        echo "HiFi Read Statistics:" > {output.stats}
        seqkit stats -T {output.long_clean} >> {output.stats}
        
        # Calculate N50
        #echo "Read length N50:" >> {output.stats}
        #seqkit fx2tab -l {output.long_clean} | \
        #awk '{{print $2}}' | sort -rn | \
        #awk '{{a[i++]=$1; s+=$1}} END {{print a[int(i/2)]}}' >> {output.stats}
        
        # Quality distribution sample
        #echo "Quality sample (first 100 reads):" >> {output.stats}
        #seqkit fx2tab -q {output.long_clean} | head -100 | \
        #awk '{{print $3}}' | awk '{{sum+=$1; count++}} END {{print "Mean quality:", sum/count}}' >> {output.stats}
        """


# MetaFlye optimized for HiFi
rule metaflye_hifi_assembly:
    input:
        #long = f"{OUTDIR}/{{sample}}/{{sample}}_hifi_clean.fastq.gz"
        long = f"{LONG_DATA_DIR}/{{sample}}{LONG_SUFFIX}"
    output:
        assembly = f"{OUTDIR}/assemblies/metaflye/{{sample}}/assembly.fasta",
        assembly_info = f"{OUTDIR}/assemblies/metaflye/{{sample}}/assembly_info.txt",
        assembly_graph = f"{OUTDIR}/assemblies/metaflye/{{sample}}/assembly_graph.gfa"
    params:
        outdir = f"{OUTDIR}/assemblies/metaflye/{{sample}}/flye_tmp",
        genome_size = config.get('estimated_genome_size', '40m'),
        qoutfile = lambda wildcards: OUTDIR /f'logs/{wildcards.sample}.metaflye.qout',
        qerrfile = lambda wildcards: OUTDIR /f'logs/{wildcards.sample}.metaflye.qerr',
        scratch = 8000,
        time = 10440,
        mem = 4000
    conda:
        "metaflye"
    threads: 128
    shell:
        """
        # HiFi-optimized metaFlye
        flye --pacbio-hifi {input.long} \
             --scaffold \
             --out-dir {params.outdir} \
             --genome-size {params.genome_size} \
             # --min-overlap 3000 \
             --meta \
             --iterations 1 \
            --resume \
             --threads 128
        
        # Copy outputs
        cp {params.outdir}/assembly.fasta {output.assembly}
        cp {params.outdir}/assembly_info.txt {output.assembly_info}
        cp {params.outdir}/assembly_graph.gfa {output.assembly_graph}
        
        # Check for circular contigs (common with HiFi)
        echo "Checking for circular contigs..." >> {output.assembly_info}
        grep -c "circular=true" {params.outdir}/assembly_info.txt >> {output.assembly_info} || echo "0" >> {output.assembly_info}

        #rm -rf {params.outdir}
        """



# Map short reads to assembly and generate BAM file
rule map_reads_for_polishing:
    input:
        assembly = f"{OUTDIR}/assemblies/{{method}}/{{sample}}/assembly.fasta",
        r1 = f"{OUTDIR}/clean_reads/{{sample}}/{{sample}}.1.fq.gz",
        r2 = f"{OUTDIR}/clean_reads/{{sample}}/{{sample}}.2.fq.gz"
    output:
        bam = f"{OUTDIR}/assemblies/{{method}}/mapping/{{sample}}/mapped_reads.bam",
        bai = f"{OUTDIR}/assemblies/{{method}}/mapping/{{sample}}/mapped_reads.bam.bai",
        stats = f"{OUTDIR}/assemblies/{{method}}/mapping/{{sample}}/mapping_stats.txt"
    params:
        outdir = f"{OUTDIR}/assemblies/{{method}}/mapping/{{sample}}",
        qoutfile = lambda wildcards: OUTDIR /f'logs/{wildcards.method}_{wildcards.sample}.bwa.qout',
        qerrfile = lambda wildcards: OUTDIR /f'logs/{wildcards.method}_{wildcards.sample}.bwa.qerr',
        scratch = 8000,
        time = 10440,
        mem = 4000
    conda:
        "metaflye"
    threads: THREADS
    resources:
        mem_mb=16000,
        runtime=360  # 6 hours
    shell:
        """
        mkdir -p {params.outdir}
        
        echo "Mapping short reads to assembly for {wildcards.sample}" > {output.stats}
        echo "Started: $(date)" >> {output.stats}
        
        # Index the assembly
        bwa index {input.assembly}
        
        # Map reads
        bwa mem -t 32 {input.assembly} {input.r1} {input.r2} | \
        samtools view -bS - | \
        samtools sort -@ 32 -o {output.bam}
        
        # Index BAM file
        samtools index {output.bam}
        
        # Generate mapping statistics
        echo "Mapping statistics:" >> {output.stats}
        samtools flagstat {output.bam} >> {output.stats}
        
        echo "Coverage statistics:" >> {output.stats}
        samtools depth {output.bam} | \
        awk '{{sum+=$3; count++}} END {{
            if (count > 0) {{
                print "Mean coverage:", sum/count
                print "Total bases covered:", count
            }} else {{
                print "No coverage data"
            }}
        }}' >> {output.stats}
        
        echo "Mapping completed: $(date)" >> {output.stats}
        
        # Clean up BWA index files
        rm -f {input.assembly}.*
        """


rule pilon_polish:
    input:
        assembly = f"{OUTDIR}/assemblies/metaflye/{{sample}}/assembly.fasta",
        bam = f"{OUTDIR}/assemblies/metaflye/mapping/{{sample}}/mapped_reads.bam",
        bai = f"{OUTDIR}/assemblies/metaflye/mapping/{{sample}}/mapped_reads.bam.bai"
    output:
        polished = f"{OUTDIR}/assemblies/metaflye_polished/{{sample}}/pilon.fasta",
        #changes = f"{OUTDIR}/assemblies/metaflye_polished/{{sample}}/pilon.changes",
        
    params:
        outdir = f"{OUTDIR}/assemblies/metaflye_polished/{{sample}}",
        pilon_memory = config.get('pilon_memory', '16g'),
        min_coverage = config.get('min_polish_coverage', 5),
        qoutfile = lambda wildcards: OUTDIR /f'logs/{wildcards.sample}.pilon.qout',
        qerrfile = lambda wildcards: OUTDIR /f'logs/{wildcards.sample}.pilon.qerr',
        conda_env = "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/Projects_NCCR/conda_envs/metaflye",
        scratch = 8000,
        time = 10440,
        mem = 48000
    conda:
        "metaflye"
    threads: THREADS
    log:
        log = OUTDIR/'logs/{sample}.pilon.log'
    shell:
        """
        mkdir -p {params.outdir}
        
        echo "Starting Pilon polishing for {wildcards.sample}" > {log.log}
        echo "Started: $(date)" >> {log.log}
        echo "Using BAM file: {input.bam}" >> {log.log}

        # Run Pilon with the BAM file
        echo "Running Pilon..." >> {log.log}
        java -Xmx200G -jar {params.conda_env}/share/pilon-1.24-0/pilon.jar \
              --genome {input.assembly} \
              --frags {input.bam} \
              --output pilon \
              --outdir {params.outdir} \
              --changes \
              --fix all \
              --mindepth {params.min_coverage} \
              --minmq 20 \
              2>&1 | tee -a {log.log}
        
       
        """


rule enhanced_binning:
    input:
        assembly = f"{OUTDIR}/assemblies/{{method}}/{{sample}}/assembly.fasta",
        bam = f"{OUTDIR}/assemblies/{{method}}/mapping/{{sample}}/mapped_reads.bam",
        bai = f"{OUTDIR}/assemblies/{{method}}/mapping/{{sample}}/mapped_reads.bam.bai"

    output:
        bin_dir = directory(f"{OUTDIR}/binning/{{method}}/{{sample}}/metabat2_bins"),
        depth_file = f"{OUTDIR}/binning/{{method}}/{{sample}}/depth.txt",
        bin_summary = f"{OUTDIR}/binning/{{method}}/{{sample}}/bin_summary.txt"
    conda:
        "metaflye"
    params:
        qoutfile = lambda wildcards: OUTDIR /f'logs/{wildcards.method}_{wildcards.sample}.binning.qout',
        qerrfile = lambda wildcards: OUTDIR /f'logs/{wildcards.method}_{wildcards.sample}.binning.qerr',
    threads: THREADS
    resources:
        mem_mb=24000
    shell:
        """
        mkdir -p {output.bin_dir}
        
        # Calculate depth
        jgi_summarize_bam_contig_depths --outputDepth {output.depth_file} {input.bam}
        
        # MetaBAT2 with settings optimized for 8 strains
        metabat2 -i {input.assembly} \
                 -a {output.depth_file} \
                 -o {output.bin_dir}/bin \
                 -t {threads} \
                 --minContig 1000 \
                 --maxP 95 \
                 --minS 60 \
                 --maxEdges 200 \
                 --ptMaxEdges 50 \
                 --seed 42
        
        # Analyze bins
        bin_count=$(ls {output.bin_dir}/bin.*.fa 2>/dev/null | wc -l)
        echo "Binning Summary for {wildcards.method} {wildcards.sample}" > {output.bin_summary}
        echo "Total bins generated: $bin_count" >> {output.bin_summary}
        echo "Expected bins: 8" >> {output.bin_summary}
        
        if [[ $bin_count -eq 8 ]]; then
            echo "✓ Achieved expected 8 bins" >> {output.bin_summary}
        elif [[ $bin_count -lt 8 ]]; then
            echo "⚠ Fewer bins than expected (may indicate merged strains)" >> {output.bin_summary}
        else
            echo "⚠ More bins than expected (may indicate fragmented strains)" >> {output.bin_summary}
        fi
        
        # Calculate bin sizes
        echo "Bin sizes:" >> {output.bin_summary}
        for bin in {output.bin_dir}/bin.*.fa; do
            if [[ -f $bin ]]; then
                size=$(seqkit stats -T $bin | tail -1 | cut -f5)
                echo "$(basename $bin): $size bp" >> {output.bin_summary}
            fi
        done
        """