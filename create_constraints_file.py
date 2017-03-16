def create_constraints_file(galaxy,single_or_double,target,filter_objs,mask):
    file = open("./"+galaxy+"/constraints"+"_"+mask+".txt", "w")

    if single_or_double == 1:
        file.write("# Component parameter constraint (Optional and written with #)\n")
        file.write("# Component number: single fit\n")
        file.write("2 n 0.1 to 10\n")
        file.write("2 mag -2. +2.\n")
        file.write("2 re 0.01 to 100\n")
        file.write("2 x -2 2\n")
        file.write("2 y -2 2\n")
        file.write("2 q 0 to 1\n")
        file.write("2 pa -90.01 to 90.01\n")
        file.write("\n")
    else:
        file.write("# Component parameter constraint (Optional and written with #)\n")
        file.write("# Component number: first component fit\n")
        file.write("2 n 0.1 to 10\n")
        file.write("2 mag -2. +2.\n")
        file.write("2 re 0.01 to 100\n")
        file.write("2 x -2 2\n")
        file.write("2 y -2 2\n")
        file.write("2 q 0 to 1\n")
        file.write("2 pa -90.01 to 90.01\n")
        file.write("\n")
        file.write("# Component number: second component fit\n")
        file.write("3 n 0.1 to 10\n")
        file.write("3 mag -2. +2.\n")
        file.write("3 re 0.01 to 100\n")
        file.write("3 x -2 2\n")
        file.write("3 y -2 2\n")
        file.write("3 q 0 to 1\n")
        file.write("3 pa -90.01 to 90.01\n")
        file.write("\n")

    order = single_or_double+1 #the first +1 is for the sky
    for ii in range(len(filter_objs)):
        if (ii != target) & (filter_objs[ii]==True) & (mask.find("normal")!=-1):
            order = order+1
            file.write("# Component number: "+str(order)+" - fit\n")
            file.write(str(order)+" n 0.1 to 10\n")
            file.write(str(order)+" mag -2. +2.\n")
            file.write(str(order)+" re 0.01 to 100\n")
            file.write(str(order)+" x -2 2\n")
            file.write(str(order)+" y -2 2\n")
            file.write(str(order)+" q 0 to 1\n")
            file.write(str(order)+" pa -90.01 to 90.01\n")
            file.write("\n")

    file.close()
