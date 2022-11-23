# Adapted from https://github.com/nf-core/rnaseq/blob/master/bin/fastq_dir_to_samplesheet.py

# !/usr/bin/env python

import os
import sys
import glob


def fastq_dir_to_samplesheet(
        fastq_dir,
        samplesheet_file,
        read1_extension="_R1_001.fastq.gz",
        read2_extension="_R2_001.fastq.gz",
        sanitise_name=False,
        sanitise_name_delimiter="_",
        sanitise_name_index=1,
):
    def sanitize_sample(path, extension):
        """Retrieve sample id from filename"""
        sample = os.path.basename(path).replace(extension, "")
        unit = ""
        if sanitise_name:
            sample = sanitise_name_delimiter.join(
                os.path.basename(path).split(sanitise_name_delimiter)[
                :sanitise_name_index
                ]
            )
            unit = sanitise_name_delimiter.join(
                os.path.basename(path).split(sanitise_name_delimiter)[
                sanitise_name_index:
                ]).replace(extension, "")
        return sample, unit

    def get_fastqs(extension):
        """
        Needs to be sorted to ensure R1 and R2 are in the same order
        when merging technical replicates. Glob is not guaranteed to produce
        sorted results.
        See also https://stackoverflow.com/questions/6773584/how-is-pythons-glob-glob-ordered
        """
        return sorted(
            glob.glob(os.path.join(fastq_dir, f"**/*{extension}"), recursive=True)
        )

    read_dict = {}
    ## Get read 1 files
    for read1_file in get_fastqs(read1_extension):
        sample, unit = sanitize_sample(read1_file, read1_extension)
        if sample not in read_dict:
            read_dict[sample] = {"unit": [], "R1": [], "R2": []}
        read_dict[sample]["R1"].append(read1_file)
        read_dict[sample]['unit'].append(unit)
    ## Get read 2 files
    for read2_file in get_fastqs(read2_extension):
        sample, _ = sanitize_sample(read2_file, read2_extension)
        read_dict[sample]["R2"].append(read2_file)
    ## Write to file
    if len(read_dict) > 0:
        out_dir = os.path.dirname(samplesheet_file)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

        with open(samplesheet_file, "w") as fout:
            header = ["sample", "unit", "fastq_1", "fastq_2"]
            fout.write(",".join(header) + "\n")
            for sample, reads in sorted(read_dict.items()):
                for idx, read_1 in enumerate(reads["R1"]):
                    read_2 = ""
                    unit = reads['unit'][idx]
                    if idx < len(reads["R2"]):
                        read_2 = reads["R2"][idx]
                    sample_info = ",".join([sample, unit, read_1, read_2])
                    fout.write(f"{sample_info}\n")
    else:
        error_str = (
            "\nWARNING: No FastQ files found so samplesheet has not been created!\n\n"
        )
        error_str += "Please check the values provided for the:\n"
        error_str += "  - Path to the directory containing the FastQ files\n"
        error_str += "  - '--read1_extension' parameter\n"
        error_str += "  - '--read2_extension' parameter\n"
        print(error_str)
        sys.exit(1)
