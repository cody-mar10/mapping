"""Common samtools workflows"""
from __future__ import annotations
import subprocess
import shlex
from pathlib import Path
from typing import BinaryIO


def samtools_view_and_sort(
    mapping: subprocess.Popen[bytes], threads: int, output: Path, logfile: BinaryIO
) -> None:
    """Given an open subprocess that is piping mapped reads to stdout, use
    samtools to convert the raw alignment to binary and sort the alignment.

    Args:
        mapping (subprocess.Popen[bytes]): opened subprocess piping mapped reads to stdout
        threads (int): number of parallel worker threads
        output (Path): name of output BAM file
        logfile (BinaryIO): opened log file to store stderr log
    """
    # TODO: supposedly -u is better for piping to another samtools command?
    sam2bam = shlex.split(f"samtools view -uS -F4 -")
    sort_bam = shlex.split(f"samtools sort -@{threads-1} - -o {output}")

    bam_process = subprocess.Popen(
        sam2bam, stdin=mapping.stdout, stdout=subprocess.PIPE, stderr=logfile,
    )
    sorting_process = subprocess.Popen(
        sort_bam, stdin=bam_process.stdout, stderr=logfile
    )
    sorting_process.wait()


def index_bam(bamfile: Path, __log: Path) -> None:
    """Index a sorted BAM file

    Args:
        bamfile (Path): file path to sorted BAM file
        __log (Path): path to store stderr log output of commands ran
    """
    command = shlex.split(f"samtools index -b {bamfile}")
    with __log.open("ab") as logfile:
        subprocess.run(command, stderr=logfile)
