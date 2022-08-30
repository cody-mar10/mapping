import shutil


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

