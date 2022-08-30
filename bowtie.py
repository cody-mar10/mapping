"""bowtie2 mapping"""
import logging
import os
import shlex
import subprocess
import time
from pathlib import Path
from typing import List

from samtools import index_bam, samtools_view_and_sort
from utils import get_read_pairs


def build_index(reference: Path, threads: int, __log: Path) -> None:
    """Build bowtie2 index

    Args:
        reference (Path): file path to reference genome
        threads (int): number of parallel worker threads
        __log (Path): path to store stderr log output of commands ran
    """
    command = shlex.split(f"bowtie2-build --threads {threads} {reference} {reference}")
    with __log.open("ab") as logfile:
        subprocess.run(command, stdout=logfile, stderr=subprocess.STDOUT)


def bowtie2_map(
    reference: Path, r1: str, r2: str, threads: int, output: Path, __log: Path
) -> None:
    """Map read pairs to the reference index and convert the resulting
    SAM alignment to the sorted BAM output. Only the final BAM file is kept.

    Args:
        reference (Path): file path to reference genome
        r1 (str): path to forward reads
        r2 (str): path to reverse reads
        threads (int): number of parallel worker threads
        output (Path): name of output BAM file
        __log (Path): path to store stderr log output of commands ran
    """
    # TODO: if multiple sets of reads, use -mm ?
    # TODO: accept other types of reads formats? -> set that logic off to get_read_pairs
    mapping = shlex.split(
        f"bowtie2 -p {threads} --no-unal -x {reference} -1 {r1} -2 {r2}"
    )

    with __log.open("ab") as logfile:
        mapping_process = subprocess.Popen(
            mapping, stdout=subprocess.PIPE, stderr=logfile
        )
        # sleep so there is some output for samtools view
        time.sleep(0.1)

        samtools_view_and_sort(mapping_process, threads, output, logfile)


def main(reference: Path, reads: List[str], threads: int, outdir: Path,) -> None:
    """Run bowtie mapping and samtools bam file sorting"""
    commands_log = outdir.joinpath("bowtie2.log")
    ext = reference.suffix
    if not reference.with_suffix(f"{ext}.1.bt2").exists():
        logging.info(f"Building index for {reference}")
        build_index(reference, threads, commands_log)
    else:
        logging.info(
            f"SKIPPING building bowtie2 a index for {reference} since it already exists."
        )

    for r1, r2 in get_read_pairs(*reads):
        output_basename = f'{os.path.basename(r1).rsplit("_1", 1)[0]}.bam'
        output = outdir.joinpath(output_basename)
        logging.info(f"Mapping {r1} and {r2} against the reference.")
        bowtie2_map(reference, r1, r2, threads, output, commands_log)
        index_bam(output, commands_log)
        logging.info(f"Mapping completed - saved to {output}")
