import os

#=====================================================
class ZSoil_version ():
#=====================================================

    #=====================================================
    def __init__ (self,project):
    #=====================================================
        self.my_dat_file = project + '.dat'
        self.dat_version = self.instanciate ()

    #=====================================================
    def get_dat_version (self):
    #=====================================================
        return self.dat_version

    #=====================================================
    def instanciate (self):
    #=====================================================
        filename = self.my_dat_file

        if os.path.exists (filename):
            dat = open (filename)
            done = False
            while not done:
                line = dat.readline()
                if "DATFILE_VERSION " in line:
                    texts = line.split("DATFILE_VERSION")
                    while texts[0] =="" and len(texts) > 0:
                        del texts[0]
                    return int(texts[0].strip())
                    done = True
        else:
            return None

        f.close ()


#=====================================================
def main():
#=====================================================
    #prj = "D:\\v16\\inp\\bugs\\18-12-dyn-strange\\test-pushover"
    #prj = "d:/v17/bench_TMP/benchmarks/boxd5"
    #prj = "d:/v17/bench_TMP/benchmarks/filldrawdown2d"
    prj = "d:/v17/bench_ref_1414/benchmarks/mlcut"
    v = ZSoil_version (prj)
    print(v.get_dat_version())

if __name__ == '__main__':
    main()
