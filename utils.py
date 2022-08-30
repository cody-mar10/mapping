import re
import shutil
from collections import defaultdict
from typing import DefaultDict, List, Tuple


def check_install(name: str) -> bool:
    """Check if a program is installed and in $PATH

    Args:
        name (str): name of command line program or executable

    Returns:
        bool: True if program is installed and in $PATH
    """
    return shutil.which(name) is not None


def check_installs(*programs: str) -> None:
    """Check multiple programs to see if they are installed

    Raises:
        RuntimeError: if any program dependency is not installed
    """
    if any(not check_install(program) for program in programs):
        raise RuntimeError(
            f"Check that the following dependencies are installed: {programs}"
        )


def get_read_pairs(*reads: str) -> List[Tuple[str, str]]:
    """Given all reads file paths, identify the paired reads that share
    the same sample basename.

    >>> reads = ["ABC_1.fastq.gz", "ABC_2.fastq.gz", \
        "DEF_1.fastq.gz", "DEF_2.fastq.gz"]
    >>> get_read_pairs(*reads)
    >>> [("ABC_1.fastq.gz", "ABC_2.fastq.gz"), ("DEF_1.fastq.gz", "DEF_2.fastq.gz")]

    Returns:
        List[Tuple[str, str]]: list of read pairs
    """
    _read_pairs: DefaultDict[str, List[str]] = defaultdict(list)

    # works for .fastq.gz .fastq .fq.gz .fq
    pattern = re.compile("_[12].f(?:ast)?q(?:.gz)?")
    for read in reads:
        sample = pattern.sub("", read)
        _read_pairs[sample].append(read)

    return [tuple(pair) for pair in _read_pairs.values()]
