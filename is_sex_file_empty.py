def is_sex_file_empty(filename):

    ff = open(filename, 'r')
    empty = True
    for line in ff:
        if line.startswith('#'):
            continue
        else:
            chuncks = line.split()
            if chuncks[0] == '1':
                empty = False
                break
            else:
                break
    return(empty)
