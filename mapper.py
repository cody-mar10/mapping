#!/usr/bin/env python3
import argparse
import logging
import sys
from pathlib import Path
from typing import List
import resource

import bowtie
import bwa
from utils import check_installs

DISPATCH_TABLE = {"bowtie2": bowtie.main, "bwa": bwa.main}

def parse_args() -> argparse.Namespace:
    # TODO: probably should have subparsers for each mapper...
    parser = argparse.ArgumentParser(
        description=(
            "Pipeline for mapping reads to a reference, converting sam to bam,"
            " and sorting the bam file. Example usage:\n\n./mapping.py -i ref.fna"
            " -r reads/*.fastq.gz -t 25 -o outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # ubparsers = parser.add_subparsers(dest="mapper")
    # bwa_parser = bwa.get_bwa_argparser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="input reference genome file in fasta format",
    )
    parser.add_argument(
        "-r",
        "--reads",
        required=True,
        nargs="+",
        help="path or glob pattern to the reads. Ex: /path/to/reads/*.fastq.gz",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default=Path.cwd().resolve(),
        help="output directory for bam files (default: %(default)s)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=15,
        help="number of worker threads to use (default: %(default)s)",
    )
    parser.add_argument(
        "-m",
        "--mapper",
        choices={"bwa", "bowtie2"},
        default="bwa",
        help="mapper to use (default: %(default)s)",
    )
    parser.add_argument(
        "-l",
        "--log",
        default="mapping.log",
        help="name of log file in output directory (default: %(default)s)",
    )

    return parser.parse_args()


def main(
    reference: Path,
    reads: List[str],
    mapper: str,
    threads: int,
    outdir: Path,
    log: Path,
) -> None:
    """_summary_

    Args:
        reference (Path): _description_
        reads (List[str]): _description_
        mapper (str): _description_
        threads (int): _description_
        outdir (Path): _description_
        log (Path): _description_

    Raises:
        NotImplementedError: _description_
    """
    outdir.mkdir(exist_ok=True, parents=True)
    logging.basicConfig(
        filename=log,
        level=logging.DEBUG,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.info(f"Command: {' '.join(sys.argv)}")
    if mapper == "bwa":
        logging.info(f"Running bwa mem for mapping")
    elif mapper == "bowtie2":
        logging.info(f"Running bowtie2 for mapping")

    mapping_fn = DISPATCH_TABLE[mapper]
    mapping_fn(reference, reads, threads, outdir)

    maxmem_GB = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / (1024 ** 2)
    logging.info(f"Max memory usage: {maxmem_GB:.3f} GB")


if __name__ == "__main__":
    args = parse_args()
    check_installs(args.mapper, "samtools")

    outdir = Path(args.outdir)
    main(
        reference=Path(args.input),
        reads=args.reads,
        mapper=args.mapper,
        threads=args.threads,
        outdir=outdir,
        log=outdir.joinpath(args.log),
    )
