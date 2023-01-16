#!/usr/bin/env python3
import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional
import resource

import bowtie
import bwa
from utils import check_installs, read_batch_file

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
    io_args = parser.add_argument_group("IO ARGS")
    reads_args = io_args.add_mutually_exclusive_group(required=True)
    batch_args = parser.add_argument_group("BATCH MODE ARGS")

    io_args.add_argument(
        "-i",
        "--input",
        required=True,
        help="input reference genome file in fasta format",
    )
    reads_args.add_argument(
        "-r",
        "--reads",
        nargs="+",
        help="path or glob pattern to the reads. Ex: /path/to/reads/*.fastq.gz -- Mutually exclusive with --batch",
    )
    reads_args.add_argument(
        "--batch",
        default=False,
        action="store_true",
        help="use if -i file is a tab-delimited file that maps assemblies (fasta files) to all reads that need to be mapped to the reference -- Mutually exclusive with -r/--reads",
    )
    io_args.add_argument(
        "-o",
        "--outdir",
        default=Path.cwd().resolve(),
        help="output directory for bam files (default: %(default)s)",
    )

    batch_args.add_argument(
        "-fd",
        "--ref-dir",
        help="directory that contains all reference fasta files (assemblies). If not provided, assumes that the batch file has the correct paths",
    )
    batch_args.add_argument(
        "-rd",
        "--reads-dir",
        help="directory that contains all reads. If not provided, assumes that the batch file has the correct paths",
    )
    batch_args.add_argument(
        "-no-header",
        help="use if your batch file provided in -i with --batch has no header",
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


def setup(outdir: Path, log: Path, mapper: str):
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


def _main(mapper: str, reference: Path, reads: List[str], threads: int, outdir: Path):
    mapping_fn = DISPATCH_TABLE[mapper]
    mapping_fn(reference, reads, threads, outdir)

    maxmem_GB = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / (1024**2)
    logging.info(f"Max memory usage: {maxmem_GB:.3f} GB")


def main(
    reference: Path,
    reads: List[str],
    mapper: str,
    threads: int,
    outdir: Path,
    log: Path,
) -> None:
    """single reference mode

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
    setup(outdir, log, mapper)
    _main(mapper, reference, reads, threads, outdir)


def main_batch(
    batchfile: Path,
    mapper: str,
    threads: int,
    outdir: Path,
    log: Path,
    header: bool = True,
    refdir: Optional[Path] = None,
    readsdir: Optional[Path] = None,
) -> None:
    """batch reference mode

    Args:
        batchfile (Path): _description_
        mapper (str): _description_
        threads (int): _description_
        outdir (Path): _description_
        log (Path): _description_

    Raises:
        NotImplementedError: _description_
    """
    setup(outdir, log, mapper)

    for reference, reads in read_batch_file(batchfile, header):
        # samtools_view_and_sort calls proc.wait() so won't open all at once
        # actually need to make the mapping step wait
        if refdir is not None:
            reference = refdir.joinpath(reference)

        if readsdir is not None:
            reads = [readsdir.joinpath(read).as_posix() for read in reads]
        _main(mapper, reference, reads, threads, outdir)


if __name__ == "__main__":
    args = parse_args()
    check_installs(args.mapper, "samtools")

    outdir = Path(args.outdir)

    if args.batch:

        main_batch(
            batchfile=Path(args.input),
            mapper=args.mapper,
            threads=args.threads,
            outdir=outdir,
            log=outdir.joinpath(args.log),
            header=(not args.no_header),
            refdir=None if args.ref_dir is None else Path(args.ref_dir),
            readsdir=None if args.reads_dir is None else Path(args.reads_dir),
        )
    else:
        main(
            reference=Path(args.input),
            reads=args.reads,
            mapper=args.mapper,
            threads=args.threads,
            outdir=outdir,
            log=outdir.joinpath(args.log),
        )
