#=====================================================
class RCF_layers_info ():
#=====================================================

    #=====================================================
    def __init__ (self,project):
    #=====================================================
        self.my_project = project
        f = open (project+".lay","rt")
        lines = f.readlines ()

        n_ele_groups = int (lines[0])
        count = 1

        self.ele_rcf = []

        for i in range (n_ele_groups):

            a = lines [count]
            a = a.split ()

            group_id = a [0]
            nitems   = int ( a [1] )

            aux1 = []

            for j in range (nitems):
                count = count + 1
                a = lines [count]
                a = a.split ()

                rsl_item  = a [0]

                #correcting misprints in older versions
                if rsl_item == "STRESESS":
                    rsl_item = "STRESSES"

                n_comp    = int ( a [1] )
                rsl_items = []

                #print 'item = ',j,rsl_item,n_comp

                if n_comp > 1:
                    for i in range (n_comp):
                        rsl_items.append (a [i+2])
                else:
                    if len(a) > 2:
                        rsl_items.append (a [2])
                    else:
                        rsl_items.append ('')

                aux2 = [rsl_item,n_comp,rsl_items]
                aux1.append (aux2)
            count = count + 1

            aux3 = [group_id,aux1]
            self.ele_rcf.append (aux3)

        f.close ()


    #=====================================================
    def give_col_for_ele_rsl (self,group_id_,rsl_item_,comp_id_):
    #=====================================================
        for i in range (len(self.ele_rcf)):
            group_rcf = self.ele_rcf [i]
            group_id  = (group_rcf [0]).strip()
            if group_id == group_id_.strip():
                list_items = group_rcf [1]
                col = 0
                for j in range (len(list_items)):
                    item  = list_items [j]
##                    print item [0]
##                    print item [1]
##                    print item [2]
                    if item [0].find (rsl_item_) != -1:
                        n_comps = item [1]
                        comps   = item [2]
                        for k in range (n_comps):
                            if comps [k].find (comp_id_) != -1:
                                return col,n_comps
                            else:
                                col = col + 1
                    else:
                        n_comps = item [1]
                        col = col + n_comps
        return -1,0

    #========================================
    def give_one_lay_size (self,group_id_):
    #========================================
        for i in range (len(self.ele_rcf)):
            group_rcf = self.ele_rcf [i]
            group_id  = (group_rcf [0]).strip()
            if group_id == group_id_.strip():
                list_items = group_rcf [1]
                count = 0
                for j in range (len(list_items)):
                    item  = list_items [j]
                    n_comps = item [1]
                    count = count + n_comps
                return count
        return -1

    #=====================================================
    def is_equivalent_to (self,rcf_):
    #=====================================================
        if self.ele_rcf != rcf_.ele_rcf:
            return False
        else:
            return True

#=====================================================
def main():
#=====================================================
    rcf = RCF_layers_info ('D:\\v17\\bench_TMP\\benchmarks\\eplbeam')
    rcf1= RCF_layers_info ('D:\\v17\\bench_TMP\\benchmarks\\eplbeam')
    print(rcf.is_equivalent_to(rcf1))


if __name__ == '__main__':
    main()