"""Read and compare ZSoil version identifiers.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`ZSoil_version` helper parses version strings and supports comparisons used when
output formats depend on the ZSoil release.
"""

import os


class ZSoil_version ():
    def __init__(self, project):
        """Initialize the zsoil version object for the supplied project data.

        :param project: Path or filename prefix of the ZSoil project.
        :type project: str
        """
        self.my_dat_file = project + '.dat'
        self.dat_version = self.instanciate()

    def get_dat_version():
        """Return input-data version.

        :return: The input-data version.
        :rtype: Optional[int]
        """
        return self.dat_version

    def instanciate(self):
        """Load this object's data from the associated ZSoil project input.

        :return: Version number read from the DAT file, or None when unavailable.
        :rtype: Optional[int]
        """
        filename = self.my_dat_file

        if os.path.exists(filename):
            dat = open(filename)
            done = False
            while not done:
                line = dat.readline()
                if "DATFILE_VERSION " in line:
                    texts = line.split("DATFILE_VERSION")
                    while texts[0] == "" and len(texts) > 0:
                        del texts[0]
                    return int(texts[0].strip())
                    done = True
        else:
            return None

        f.close()


def main():
    """Run this module's example or command-line workflow."""
    #prj = "D:\\v16\\inp\\bugs\\18-12-dyn-strange\\test-pushover"
    #prj = "d:/v17/bench_TMP/benchmarks/boxd5"
    #prj = "d:/v17/bench_TMP/benchmarks/filldrawdown2d"
    prj = "d:/v17/bench_ref_1414/benchmarks/mlcut"
    v = ZSoil_version(prj)
    print(v.get_dat_version())


if __name__ == '__main__':
    main()
