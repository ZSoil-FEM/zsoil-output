


#=====================================================
class ExecutionLOG ():
#=====================================================

    #=====================================================
    def __init__ (self,project):
    #=====================================================
        if ".log" in project:
            self.my_logfile = project
        else:
            self.my_logfile = project + ".log"

        self.INIT_STATE = 0
        self.STABILITY  = 1
        self.TIME_DEPNT = 2
        self.ARC_LENGTH = 3
        self.DYNAMICS   = 4
        self.PUSHOVER   = 5
        self.EIGENMODES = 6

        self.drivers_dict = {self.INIT_STATE:"Initial State" ,self.STABILITY:"Stability",\
                             self.TIME_DEPNT:"Time Dependent",self.ARC_LENGTH:"Displ.Driven ( with arc length )",\
                             self.DYNAMICS  :"Dynamics"      ,self.PUSHOVER  :"Static pushover",\
                             self.EIGENMODES:"Eigenvalue analysis"}

        self.drivers_dict_inv = {v: k for k, v in self.drivers_dict.items()}
        self.echo = []
        self.instanciate ()
        self.verbose = False


    #=====================================================
    def instanciate (self):
    #=====================================================
        filename = self.my_logfile
        try:
            f = open (filename,"rt")
        except IOError as e:
            print("log file:",filename)
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            return []
        except:
            print("log file:",filename)
            print("I/O error({0}): {1}")
            raise
            return []

        log_content = f.readlines ()
        done = False
        count = 0

        while not done:
            # seek for the next driver
            while not done and not "Driver type" in log_content [count]:
                count = count + 1
                if count >= len(log_content):
                    done = True

            if not done:
                # driver is found
                # aa = log_content [count]
                driver_id = (log_content [count].split(":")) [-1]
                driver_id = driver_id.strip()

                if self.drivers_dict_inv [driver_id] == self.INIT_STATE:
                    count = self.analyze_init_state_driver (log_content,count,self.echo)
                elif self.drivers_dict_inv [driver_id] == self.STABILITY:
                    count = self.analyze_stability_driver  (log_content,count,self.echo)
                elif self.drivers_dict_inv [driver_id] == self.TIME_DEPNT:
                    count = self.analyze_timedependent_driver (log_content,count,self.echo)
                elif self.drivers_dict_inv [driver_id] == self.ARC_LENGTH:
                    count = self.analyze_arclength_driver(log_content, count, self.echo)
                elif self.drivers_dict_inv [driver_id] == self.DYNAMICS:
                    count = self.analyze_dynamics_driver (log_content,count,self.echo)
                elif self.drivers_dict_inv [driver_id] == self.PUSHOVER:
                  count = self.analyze_pushover_driver(log_content, count, self.echo)
                elif self.drivers_dict_inv[driver_id] == self.EIGENMODES:
                  count = self.analyze_eigenmodes_driver(log_content, count, self.echo)

            if count >= len(log_content):
                done = True

        f.close ()

    #=====================================================
    def analyze_timedependent_driver (self,log_content,count,ret):
    # =====================================================
        done = False
        count = count + 1

        while not done:

            if count >= len(log_content)-2:
                return count

            ok = False
            while not ok:
                if "Time" in log_content [count]:
                    ok = True
                else:
                    count = count +1
                    if count < len(log_content):
                        if "Driver type" in log_content [count]:
                            ok = True
                            done = True
                    else:
                        ok   = True
                        done = True

            if ok and not done:
                aux_dict = {"REACTIONS": []}
                aux_dict.update ({"DRIVER":self.TIME_DEPNT})
                step =  float((log_content [count].split (":"))[1])
                aux_dict.update ({"STEP":step})

                aux_dict.update({"EXTRA_DATA": []})


                done_aux = False

                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "----> 1 Stage : Solution of Steady State Fluid Flow" in log_content [count]:
                        aux_dict.update ({"SOLVING_STAGE":"FLOW_SOLVING"})
                        done_aux = True
                    elif "----> 2 Stage : Solution of the Deformation State" in log_content [count]:
                        aux_dict.update ({"SOLVING_STAGE":"DEF_SOLVING"})
                        done_aux = True
                    elif "Nonlinear solver" in log_content [count]:
                        aux_dict.update ({"SOLVING_STAGE":"STD_SOLVING"})
                        done_aux = True
                        count = count - 1

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Nonlinear solver" in log_content [count]:
                        texts = (log_content [count]).split(":")
                        aux_dict.update ({"NL_SOLVER":texts[-1].strip()})
                        done_aux = True

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Linear solver" in log_content[count]:
                        texts = (log_content[count]).split(":")
                        aux_dict.update({"LIN_SOLVER": texts[-1].strip()})
                        done_aux = True

                #add iterations

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "ITER" in log_content [count]:
                        done_aux = True
                    if count >= len(log_content):
                        done_aux = True

                echo_iter = []
                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    line  = log_content [count].strip()
                    if line == "":
                        texts = [""]
                    else:
                        texts = line.split()
                    try:
                        iter = int(texts[0])
                        size = len (texts)
                        aux = []
                        aux.append (iter)
                        for i in range (1,size-1):
                            aux.append (float(texts [i]))
                        aux.append (texts [-1])
                        echo_iter.append (aux)
                    except ValueError:
                        done_aux = True
                aux_dict.update ({"ITERATIONS":echo_iter})
                ret.append (aux_dict)

                if aux_dict ["SOLVING_STAGE"] != "FLOW_SOLVING":

                    done_aux = False
                    while not done_aux:
                        count = count + 1
                        if count >= len(log_content):
                            return count
                        if "SUM OF REACTIONS" in log_content [count]:
                            done_aux = True
                        if count >= len(log_content):
                            done_aux = True

                    if "SUM OF REACTIONS" in log_content[count]:
                        reactions = aux_dict ["REACTIONS"]
                        done_aux = False
                        while not done_aux:
                            count = count + 1
                            if count >= len(log_content):
                                return count
                            line  = log_content [count]
                            line  = line.replace ("|","")
                            line  = line.strip()
                            texts = line.split ()
                            try:
                                size = len(texts)
                                aux = []
                                for i in range(size):
                                    aux.append(float(texts[i]))
                                for v in aux:
                                    reactions.append (v)
                                done_aux = True
                            except ValueError:
                                done_aux = False

        return count


    #=====================================================
    def analyze_dynamics_driver (self,log_content,count,ret):
    # =====================================================
        done = False
        count = count + 1

        while not done:

            if count >= len(log_content)-2:
                return count

            ok = False
            while not ok:
                if "Dynamic time integration step" in log_content [count]:
                    ok = True
                else:
                    count = count +1
                    if count < len(log_content):
                        if "Driver type" in log_content [count]:
                            ok = True
                            done = True
                    else:
                        ok   = True
                        done = True

            if ok and not done:
                aux_dict = {"REACTIONS": []}
                aux_dict.update ({"DRIVER":self.DYNAMICS})
                step =  float((log_content [count].split (":"))[1])
                aux_dict.update ({"STEP":step})
                aux_dict.update ({"SOLVING_STAGE":"STD_SOLVING"})

                aux_dict.update({"EXTRA_DATA": []})

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Nonlinear solver" in log_content [count]:
                        texts = (log_content [count]).split(":")
                        aux_dict.update ({"NL_SOLVER":texts[-1].strip()})
                        done_aux = True

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Linear solver" in log_content[count]:
                        texts = (log_content[count]).split(":")
                        aux_dict.update({"LIN_SOLVER": texts[-1].strip()})
                        done_aux = True

                #add iterations

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "ITER" in log_content [count]:
                        done_aux = True
                    if count >= len(log_content):
                        done_aux = True

                echo_iter = []
                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    line  = log_content [count].strip()
                    if line == "":
                        texts = [""]
                    else:
                        texts = line.split()
                    try:
                        iter = int(texts[0])
                        size = len (texts)
                        aux = []
                        aux.append (iter)
                        for i in range (1,size-1):
                            aux.append (float(texts [i]))
                        aux.append (texts [-1])
                        echo_iter.append (aux)
                    except ValueError:
                        done_aux = True
                aux_dict.update ({"ITERATIONS":echo_iter})
                ret.append (aux_dict)

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "SUM OF REACTIONS" in log_content [count]:
                        done_aux = True
                    if count >= len(log_content):
                        done_aux = True

                if "SUM OF REACTIONS" in log_content[count]:
                    reactions = aux_dict ["REACTIONS"]
                    done_aux = False
                    while not done_aux:
                        count = count + 1
                        if count >= len(log_content):
                            return count
                        line  = log_content [count]
                        line  = line.replace ("|","")
                        line  = line.strip()
                        texts = line.split ()
                        try:
                            size = len(texts)
                            aux = []
                            for i in range(size):
                                aux.append(float(texts[i]))
                            for v in aux:
                                reactions.append (v)
                            done_aux = True
                        except ValueError:
                            done_aux = False

        return count

    # =====================================================
    def analyze_init_state_driver(self, log_content, count, ret):
        # =====================================================

        # for i in range (120):
        #     print i+1,"---",log_content [i].strip()


        done = False
        count = count + 1

        def_flow_flag = 0

        while not done:

            if count >= len(log_content)-2:
                return count

            aux_dict = {"REACTIONS": []}
            ok = False
            while not ok:

                if "----> 1 Stage : Solution of Steady State Fluid Flow" in log_content[count]:
                    aux_dict.update ({"SOLVING_STAGE":"FLOW_SOLVING"})
                    def_flow_flag = 1
                    ok = True
                elif "Actual Load Factor" in log_content[count]:
                    aux_dict.update ({"SOLVING_STAGE":"DEF_SOLVING"})
                    ok = True
                else:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "SUM OF REACTIONS" in log_content[count] or "Nonlinear solver" in log_content[count]:
                        ok = True
                        done = True
            if ok and not done:

                if not "SOLVING_STAGE" in list(aux_dict.keys()):
                    aux_dict.update({"SOLVING_STAGE": "STD_SOLVING"})

                aux_dict.update({"DRIVER": self.INIT_STATE})
                if aux_dict ["SOLVING_STAGE"] == "DEF_SOLVING":
                    step = float((log_content[count].split(":"))[1])
                else:
                    step = 1.0
                aux_dict.update({"STEP": step})

                aux_dict.update({"EXTRA_DATA": []})

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Nonlinear solver" in log_content[count]:
                        texts = (log_content[count]).split(":")
                        aux_dict.update({"NL_SOLVER": texts[-1].strip()})
                        done_aux = True

                # print "55kddddd",log_content[55]

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Linear solver" in log_content[count]:
                        texts = (log_content[count]).split(":")
                        aux_dict.update({"LIN_SOLVER": texts[-1].strip()})
                        done_aux = True

                # add iterations

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "ITER" in log_content[count]:
                        done_aux = True
                    if count >= len(log_content):
                        done_aux = True

                echo_iter = []
                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    line = log_content[count].strip()
                    # print count+1,line
                    if line == "":
                        texts = [""]
                    else:
                        texts = line.split()
                    try:
                        iter = int(texts[0])
                        size = len(texts)
                        aux = []
                        aux.append(iter)
                        for i in range(1, size - 1):
                            aux.append(float(texts[i]))
                        aux.append(texts[-1])
                        echo_iter.append(aux)
                    except ValueError:
                        done_aux = True
                aux_dict.update({"ITERATIONS": echo_iter})
                ret.append(aux_dict)

        if "SUM OF REACTIONS" in log_content[count]:
            aux_dict = ret[-1]  # take last record
            reactions = aux_dict["REACTIONS"]
            done_aux = False
            while not done_aux:
                count = count + 1
                if count >= len(log_content):
                    return count
                line = log_content[count]
                line = line.replace("|", "")
                line = line.strip()
                texts = line.split()
                try:
                    size = len(texts)
                    aux = []
                    for i in range(size):
                        aux.append(float(texts[i]))
                    for v in aux:
                        reactions.append (v)
                    done_aux = True
                    count = count + 1
                except ValueError:
                    done_aux = False

        return count

    #=====================================================
    def analyze_stability_driver (self,log_content,count,ret):
    # =====================================================
        done = False
        count = count + 1

        while not done:

            if count >= len(log_content)-2:
                return count

            ok = False
            while not ok:
                if "Safety Factor" in log_content [count]:
                    ok = True
                else:
                    count = count +1
                    if count < len(log_content):
                        if "Driver type" in log_content [count]:
                            ok = True
                            done = True
                    else:
                        ok   = True
                        done = True

            if ok and not done:
                aux_dict = {"REACTIONS": []}
                aux_dict.update ({"DRIVER":self.STABILITY})
                step =  float((log_content [count].split (":"))[1])
                aux_dict.update ({"STEP":step})
                aux_dict.update({"SOLVING_STAGE": "DEF_SOLVING"})
                aux_dict.update({"EXTRA_DATA": []})

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Nonlinear solver" in log_content [count]:
                        texts = (log_content [count]).split(":")
                        aux_dict.update ({"NL_SOLVER":texts[-1].strip()})
                        done_aux = True

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Linear solver" in log_content[count]:
                        texts = (log_content[count]).split(":")
                        aux_dict.update({"LIN_SOLVER": texts[-1].strip()})
                        done_aux = True

                #add iterations

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "ITER" in log_content [count]:
                        done_aux = True
                    if count >= len(log_content):
                        done_aux = True

                echo_iter = []
                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    line  = log_content [count].strip()
                    if line == "":
                        texts = [""]
                    else:
                        texts = line.split()
                    try:
                        iter = int(texts[0])
                        size = len (texts)
                        aux = []
                        aux.append (iter)
                        for i in range (1,size-1):
                            aux.append (float(texts [i]))
                        aux.append (texts [-1])
                        echo_iter.append (aux)
                    except ValueError:
                        done_aux = True
                aux_dict.update ({"ITERATIONS":echo_iter})
                ret.append (aux_dict)

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "SUM OF REACTIONS" in log_content [count]:
                        done_aux = True
                    if count >= len(log_content):
                        done_aux = True

                if "SUM OF REACTIONS" in log_content[count]:
                    reactions = aux_dict ["REACTIONS"]
                    done_aux = False
                    while not done_aux:
                        count = count + 1
                        if count >= len(log_content):
                            return count
                        line  = log_content [count]
                        line  = line.replace ("|","")
                        line  = line.strip()
                        texts = line.split ()
                        try:
                            size = len(texts)
                            aux = []
                            for i in range(size):
                                aux.append(float(texts[i]))
                            for v in aux:
                                reactions.append (v)
                            done_aux = True
                        except ValueError:
                            done_aux = False


        return count


    #=====================================================
    def analyze_eigenmodes_driver (self,log_content,count,ret):
    # =====================================================
        done = False
        count = count + 1

        text_w = "w"
        eigen_index = 1
        text_wi = text_w + str(eigen_index)

        ok = False
        while not ok:
            if text_wi in log_content[count]:
                ok = True
            else:
                count = count + 1
                if count >= len(log_content):
                    return count


        aux_dict = {"REACTIONS": []}
        aux_dict.update ({"DRIVER":self.EIGENMODES})
        step =  0.0
        aux_dict.update ({"STEP":step})
        aux_dict.update ({"SOLVING_STAGE": "DEF_SOLVING"})
        aux_dict.update ({"NL_SOLVER" : "-"})
        aux_dict.update ({"LIN_SOLVER": "-"})
        aux_dict.update ({"ITERATIONS": []})

        wi_list = []
        eigen_index = 0
        done = False

        while not done:

            if count >= len(log_content)-2:
                return count

            eigen_index = eigen_index + 1
            text_wi = text_w + str(eigen_index)
            ok = False
            while not ok:
                if text_wi in log_content [count]:
                    ok = True
                else:
                    line = log_content [count].strip()
                    if line == "":
                        done = True
                        aux_dict.update({"EXTRA_DATA": wi_list})
                        ret.append (aux_dict)
                        return count

                    count = count + 1
                    if count >= len(log_content):
                        ret.append (aux_dict)
                        return count


            line  = log_content [count].strip()
            if line == "":
                texts = [""]
            else:
                texts = line.split("=")
                texts = texts[-1].split("Hz")
            try:
                wi_Hz = float(texts[0])
                wi_list.append (wi_Hz)
            except ValueError:
                a = 1 #nothing

            count = count + 1

        return count


    #=====================================================
    def analyze_arclength_driver (self,log_content,count,ret):
    # =====================================================
        done = False
        count = count + 1

        while not done:

            if count >= len(log_content)-2:
                return count

            ok = False
            while not ok:
                if "Arc Length Step Nr" in log_content [count]:
                    ok = True
                else:
                    count = count +1
                    if count < len(log_content):
                        if "Driver type" in log_content [count]:
                            ok = True
                            done = True
                    else:
                        ok   = True
                        done = True

            if ok and not done:
                aux_dict = {"REACTIONS": []}
                aux_dict.update ({"DRIVER":self.ARC_LENGTH})
                step =  int((log_content [count].split (":"))[1])
                aux_dict.update ({"STEP":step})
                aux_dict.update({"SOLVING_STAGE": "DEF_SOLVING"})
                aux_dict.update({"EXTRA_DATA": []})

                done_aux = False

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Nonlinear solver" in log_content [count]:
                        texts = (log_content [count]).split(":")
                        aux_dict.update ({"NL_SOLVER":texts[-1].strip()})
                        done_aux = True

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Linear solver" in log_content[count]:
                        texts = (log_content[count]).split(":")
                        aux_dict.update({"LIN_SOLVER": texts[-1].strip()})
                        done_aux = True

                #add iterations
                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "ITER" in log_content [count]:
                        done_aux = True
                    if count >= len(log_content):
                        done_aux = True

                echo_iter = []
                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        texts = [""]
                        done_aux = True
                    else:
                        line = log_content[count].strip()
                    if line == "":
                        texts = [""]
                    else:
                        texts = line.split()
                    try:
                        iter = int(texts[0])
                        size = len (texts)
                        aux = []
                        aux.append (iter)
                        for i in range (1,size-1):
                            aux.append (float(texts [i]))
                        aux.append (texts [-1])
                        echo_iter.append (aux)
                    except ValueError:
                        done_aux = True
                aux_dict.update ({"ITERATIONS":echo_iter})
                ret.append (aux_dict)

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "SUM OF REACTIONS" in log_content [count]:
                        done_aux = True
                    if count >= len(log_content):
                        done_aux = True

                if "SUM OF REACTIONS" in log_content[count]:
                    reactions = aux_dict ["REACTIONS"]
                    done_aux = False
                    while not done_aux:
                        count = count + 1
                        if count >= len(log_content):
                            return count
                        line  = log_content [count]
                        line  = line.replace ("|","")
                        line  = line.strip()
                        texts = line.split ()
                        try:
                            size = len(texts)
                            aux = []
                            for i in range(size):
                                aux.append(float(texts[i]))
                            for v in aux:
                                reactions.append (v)
                            done_aux = True
                        except ValueError:
                            done_aux = False

                done_aux = False
                while not done_aux:
                    if "Variable load factor =" in log_content[count]:
                        done_aux = True
                    else:
                        count = count + 1
                        if count >= len(log_content):
                            done_aux = True

                if count < len(log_content):
                    line  = log_content[count]
                    texts = line.split ()
                    try:
                        size = len(texts)
                        aux = []
                        aux.append (float(texts[3]))
                        aux.append (float(texts[8]))
                        for v in aux:
                            aux_dict["EXTRA_DATA"].append(v)
                    except ValueError:
                        a = 1

        return count

    #=====================================================
    def analyze_pushover_driver (self,log_content,count,ret):
    # =====================================================
        done = False
        count = count + 1

        while not done:

            if count >= len(log_content)-2:
                return count

            ok = False
            while not ok:
                if "Pushover   Step Nr" in log_content [count]:
                    ok = True
                else:
                    count = count +1
                    if count < len(log_content):
                        if "Driver type" in log_content [count]:
                            ok = True
                            done = True
                    else:
                        ok   = True
                        done = True

            if ok and not done:
                aux_dict = {"REACTIONS": []}
                aux_dict.update ({"DRIVER":self.ARC_LENGTH})
                step =  int((log_content [count].split (":"))[1])
                aux_dict.update ({"STEP":step})
                aux_dict.update({"SOLVING_STAGE": "DEF_SOLVING"})
                aux_dict.update({"EXTRA_DATA": []})

                done_aux = False

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Nonlinear solver" in log_content [count]:
                        texts = (log_content [count]).split(":")
                        aux_dict.update ({"NL_SOLVER":texts[-1].strip()})
                        done_aux = True

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "Linear solver" in log_content[count]:
                        texts = (log_content[count]).split(":")
                        aux_dict.update({"LIN_SOLVER": texts[-1].strip()})
                        done_aux = True

                #add iterations
                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "ITER" in log_content [count]:
                        done_aux = True
                    if count >= len(log_content):
                        done_aux = True

                echo_iter = []
                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    line  = log_content [count].strip()
                    if line == "":
                        texts = [""]
                    else:
                        texts = line.split()
                    try:
                        iter = int(texts[0])
                        size = len (texts)
                        aux = []
                        aux.append (iter)
                        for i in range (1,size-1):
                            aux.append (float(texts [i]))
                        aux.append (texts [-1])
                        echo_iter.append (aux)
                    except ValueError:
                        done_aux = True
                aux_dict.update ({"ITERATIONS":echo_iter})
                ret.append (aux_dict)

                done_aux = False
                while not done_aux:
                    count = count + 1
                    if count >= len(log_content):
                        return count
                    if "SUM OF REACTIONS" in log_content [count]:
                        done_aux = True
                    if count >= len(log_content):
                        done_aux = True

                if "SUM OF REACTIONS" in log_content[count]:
                    reactions = aux_dict ["REACTIONS"]
                    done_aux = False
                    while not done_aux:
                        count = count + 1
                        if count >= len(log_content):
                            return count
                        line  = log_content [count]
                        line  = line.replace ("|","")
                        line  = line.strip()
                        texts = line.split ()
                        try:
                            size = len(texts)
                            aux = []
                            for i in range(size):
                                aux.append(float(texts[i]))
                            for v in aux:
                                reactions.append (v)
                            done_aux = True
                        except ValueError:
                            done_aux = False

                done_aux = False
                while not done_aux:
                    if "Load multiplier =" in log_content[count]:
                        done_aux = True
                    else:
                        count = count + 1
                        if count >= len(log_content):
                            done_aux = True

                if count < len(log_content):
                    line  = log_content[count]
                    texts = line.split ()
                    try:
                        size = len(texts)
                        aux = []
                        aux.append (float(texts[3]))
                        aux.append (float(texts[7]))
                        for v in aux:
                            aux_dict["EXTRA_DATA"].append(v)
                    except ValueError:
                        a = 1

        return count

#=====================================================
    def compare_one_record (self,rec1,rec2):
#=====================================================
        err_keys = []
        if rec1 ["DRIVER"] != rec2 ["DRIVER"]:
            err_keys.append ("DRIVER")
            return err_keys

        rel_val = max(abs(rec1 ["STEP"]),abs(rec2 ["STEP"]))
        rel_val = max(rel_val,1.0e-8)
        if abs(rec1 ["STEP"] - rec2 ["STEP"])/rel_val > 1.0e-6:
            err_keys.append ("STEP")

        if rec1 ["SOLVING_STAGE"] != rec2 ["SOLVING_STAGE"]:
            err_keys.append ("SOLVING_STAGE")

        if rec1 ["NL_SOLVER"] != rec2 ["NL_SOLVER"]:
            err_keys.append ("NL_SOLVER")

        if rec1 ["LIN_SOLVER"] != rec2 ["LIN_SOLVER"]:
            err_keys.append ("LIN_SOLVER")

        if len(rec1 ["EXTRA_DATA"]) == len(rec2 ["EXTRA_DATA"]):

            size = len(rec1 ["EXTRA_DATA"])
            aux1 = rec1 ["EXTRA_DATA"]
            aux2 = rec2 ["EXTRA_DATA"]

            skip = False
            for i in range (size):
                ref_value = max(abs(aux1[i]),abs(aux2[i]))
                ref_value = max(1.0e-7,ref_value)
                dv = aux1[i]-aux2[i]
                if abs(dv)/ref_value > 1.0e-5:
                    if not skip:
                        err_keys.append("EXTRA_DATA")
                        skip = True
        else:
            err_keys.append("EXTRA_DATA")


        if len(rec1["REACTIONS"]) == len(rec2["REACTIONS"]):

            size = len(rec1["REACTIONS"])
            aux1 = rec1["REACTIONS"]
            aux2 = rec2["REACTIONS"]

            skip = 0
            for i in range(size):
                ref_value = max(abs(aux1[i]), abs(aux2[i]))
                ref_value = max(1.0, ref_value)
                dR = abs(aux1[i] - aux2[i])
                if dR / ref_value > 1.0e-5:
                    if not skip:
                        err_keys.append("REACTIONS")
                        skip = 1
        else:
            err_keys.append("REACTIONS")

        if len(rec1["ITERATIONS"]) == len(rec2["ITERATIONS"]):

            size = len(rec1["ITERATIONS"])
            aux1 = rec1["ITERATIONS"]
            aux2 = rec2["ITERATIONS"]

            skip = False
            for i in range(size):
                if i == 0:
                    row1    = aux1 [i]
                    row_ref = row1 [:]
                if len(aux1 [i]) != len(aux2 [i]) and not skip:
                    err_keys.append("ITERATIONS")
                    skip = True
                if not skip:
                    row1 = aux1 [i]
                    row2 = aux2 [i]
                    for j in range (len(row1)-1):
                        ref_value = row_ref [j] # max(abs(row1[j]), abs(row2[j]))
                        ref_value = max(1.0e-4, ref_value)
                        dF = abs(row1 [j] - row2 [j])
                        if dF / ref_value > 2.0e-4:
                            if not skip:
                                err_keys.append("ITERATIONS")
                                skip = True
                    if row1[-1] != row2[-1]:
                        if not skip:
                            err_keys.append("ITERATIONS")
                            skip = True

        else:
            err_keys.append("ITERATIONS")

        return err_keys

#=====================================================
    def is_equivalent_to (self,log_,log_aux=None):
# =====================================================
        err = 0
        if log_aux != None and len(self.echo) != len(log_.echo):
            log_aux.write ("size(" + self.my_logfile + ")=" + str(len(self.echo)) +"\n")
            if self.verbose:
                print("size(" + self.my_logfile + ")=" + str(len(self.echo)))
            log_aux.write ("size(" + log_.my_logfile + ")=" + str(len(log_.echo)) +"\n")
            if self.verbose:
                print("size(" + log_.my_logfile + ")=" + str(len(log_.echo)))

        size = min (len(self.echo),len(log_.echo))
        for i in range (size):
            rec1 = self.echo [i]
            rec2 = log_.echo [i]
            ret = self.compare_one_record (rec1,rec2)
            if len(ret) != 0:
                if log_aux != None:
                    log_aux.write ("difference in log record "  + str(i+1)  + "\n")
                    if self.verbose:
                        print("difference in log record "  + str(i+1))
                    log_aux.write ("conflicts in keywords : ")
                    if self.verbose:
                        print("conflicts in keywords : ")
                    for item in ret:
                        log_aux.write ( item )
                        log_aux.write ("\n")
                        if self.verbose:
                            print(item)
                return False

        if len(self.echo) != len(log_.echo):
            return False

        return True

#=====================================================
def main():
#=====================================================
    #prj = "D:\\v16\\inp\\bugs\\18-12-dyn-strange\\test-pushover"
    #prj = "d:/v17/bench_TMP/benchmarks/boxd5"
    #prj = "d:/v17/bench_TMP/benchmarks/filldrawdown2d"
    prj = 'D:\\V17\BENCH_REF\\concrete\\CPDM-KUPFER--1-0'
    log1 = ExecutionLOG (prj)
    prj = 'D:\\V17\BENCH_TMP\\concrete\\CPDM-KUPFER--1-0'
    log2 = ExecutionLOG (prj)
    #for item in log1.echo:
    #    print item
    aux_log = open ("d:\\jaja.txt","wt")
    log1.verbose = True
    status = log1.is_equivalent_to (log2,aux_log)
    aux_log.close ()
    print(status)

if __name__ == '__main__':
    main()
