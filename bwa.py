"""bwa mem mapping"""
import argparse
import logging
import os
import shlex
import subprocess
import time
from pathlib import Path
from typing import List

from samtools import index_bam, samtools_view_and_sort
from utils import get_read_pairs


def get_bwa_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="bwa mem mapper argument parser -- only common args implemented"
    )

    alg_opts = parser.add_argument_group("ALGORITHM OPTIONS")
    scoring_opts = parser.add_argument_group("SCORING OPTIONS")
    io_opts = parser.add_argument_group("INPUT/OUTPUT OPTIONS")

    #########----------ALGORITHM OPTIONS----------#########
    alg_opts.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="number of parallel worker threads (default: %(default)s)",
    )
    alg_opts.add_argument(
        "-k",
        "--min-seed-len",
        metavar="K",
        type=int,
        default=19,
        help="minimum seed length (default: %(default)s)",
    )
    alg_opts.add_argument(
        "-w",
        "--width",
        type=int,
        default=100,
        help="band width for banded alignment (default: %(default)s)",
    )
    alg_opts.add_argument(
        "-d",
        "--diagonal",
        type=int,
        default=100,
        help="off-diagonal X-dropoff (default: %(default)s)",
    )
    alg_opts.add_argument(
        "-r",
        type=float,
        default=1.5,
        help="look for internal seeds inside a seed longer than {-k} * [%(default)s]",
    )
    alg_opts.add_argument(
        "-y",
        type=int,
        default=20,
        help="seed occurrence for the 3rd round seeding (default: %(default)s)",
    )
    alg_opts.add_argument(
        "-c",
        "--skip-seed",
        metavar="C",
        type=int,
        default=500,
        help="skip seeds with more than this many occurrences (default: %(default)s)",
    )
    alg_opts.add_argument(
        "-D",
        "--drop",
        type=float,
        default=0.5,
        help="drop chains shorter than this fraction of the longest overlapping chain (default: %(default)s)",
    )
    alg_opts.add_argument(
        "-W",
        type=int,
        default=0,
        help="discard a chain if seeded bases shorter than this (default: %(default)s)",
    )
    alg_opts.add_argument(
        "-m",
        "--mate-rescues",
        metavar="m",
        type=int,
        default=50,
        help="perform as most this many rounds of mate rescues for each read (default: %(default)s)",
    )
    alg_opts.add_argument(
        "-S",
        "--skip-mate-rescue",
        default=False,
        action="store_true",
        help="skip mate rescue",
    )
    alg_opts.add_argument(
        "-P",
        "--skip-pairing",
        default=False,
        action="store_true",
        help="skip pairing; mate rescue performed unless -S also in use",
    )
    #########----------SCORING OPTIONS----------#########
    scoring_opts.add_argument(
        "-A",
        "--scale",
        metavar="A",
        type=int,
        default=1,
        help="score for a sequence match, which scales options -TdBOELU unless overridden (default: %(default)s)",
    )
    scoring_opts.add_argument(
        "-B",
        "--mismatch-penalty",
        metavar="B",
        type=int,
        default=4,
        help="penalty for a mismatch (default: %(default)s)",
    )
    scoring_opts.add_argument(
        "-O",
        "--gap-open",
        metavar="O",
        nargs=2,
        type=int,
        default=(6, 6),
        help="gap open penalties for indels (default: %(default)s)",
    )
    scoring_opts.add_argument(
        "-E",
        "--gap-ext",
        metavar="E",
        nargs=2,
        type=int,
        default=(1, 1),
        help="gap extension penalties; a gap of size k costs {-O} + {-E}*k (default: %(default)s)",
    )
    scoring_opts.add_argument(
        "-L",
        "--clipping-penalty",
        metavar="L",
        nargs=2,
        type=int,
        default=(5, 5),
        help="penalty for 5'- and 3'-end clipping (default: %(default)s)",
    )
    scoring_opts.add_argument(
        "-U",
        "--unpaired-penality",
        metavar="U",
        type=int,
        default=17,
        help="penalty for unpaired read pair (default: %(default)s)",
    )
    scoring_opts.add_argument(
        "-x",
        "--read-type",
        choices={"pacbio", "ont2d", "intractg", "null"},
        default="null",
        help="read type. Setting -x changes multiple parameters (default: %(default)s)",
    )
    #########----------IO OPTIONS----------#########
    io_opts.add_argument(
        "-p",
        "--interleaved",
        default=False,
        action="store_true",
        help="use if reads are interleaved into a single file",
    )
    io_opts.add_argument(
        "-R",
        "--read-group-header",
        metavar="R",
        default="null",
        help="read group header line such as @RG\tID:foo\tSM:bar (default: %(default)s)",
    )
    io_opts.add_argument(
        "-H",
        "--header",
        metavar="H",
        default="null",
        help="insert this into header if it starts with @; or insert lines into output file (default: %(default)s)",
    )
    io_opts.add_argument(
        "-o",
        "--output",
        default="stdout",
        help="sam file to output results to (default: %(default)s)",
    )
    io_opts.add_argument(
        "-j",
        default=False,
        action="store_true",
        help="treat ALT contigs as part of the primary assembly (ie ignore <idxbase>.alt file",
    )
    # TODO: -5, -q, -K, -v, -T, -h, -a, -C, -V, -Y, -M, -I
    return parser


def build_index(reference: Path, __log: Path) -> None:
    """Build a bwa index

    Args:
        reference (Path): file path to reference genome
        __log (Path): path to store stderr log output of commands ran
    """
    command = shlex.split(f"bwa index {reference}")
    with __log.open("ab") as logfile:
        subprocess.run(command, stderr=logfile)


def bwa_mem(
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
    # TODO: samtools view and sort also are using 1.5 CPUs, should I decrease input here?
    mapping = shlex.split(f"bwa mem -t {threads} {reference} {r1} {r2}")

    with __log.open("ab") as logfile:
        mapping_process = subprocess.Popen(
            mapping, stdout=subprocess.PIPE, stderr=logfile
        )
        # sleep so there is some output for samtools view
        time.sleep(0.1)

        samtools_view_and_sort(mapping_process, threads, output, logfile)


# TODO: this code is almost identical to bowtie....need a common base func
def main(reference: Path, reads: List[str], threads: int, outdir: Path,) -> None:
    """Run bwa mem mapping and samtools bam file sorting"""
    commands_log = outdir.joinpath("bwa-mem.log")
    ext = reference.suffix
    ref_basename = reference.stem
    if not reference.with_suffix(f"{ext}.amb").exists():
        logging.info(f"Building index for {reference}")
        build_index(reference, commands_log)
    else:
        logging.info(
            f"SKIPPING building bwa index for {reference} since it already exists."
        )

    for r1, r2 in get_read_pairs(*reads):
        output_basename = (
            f'{ref_basename}_{os.path.basename(r1).rsplit("_1", 1)[0]}.sorted.bam'
        )
        output = outdir.joinpath(output_basename)
        logging.info(f"Mapping {r1} and {r2} against the reference.")
        bwa_mem(reference, r1, r2, threads, output, commands_log)
        index_bam(output, commands_log)
        logging.info(f"Mapping completed - saved to {output}")
