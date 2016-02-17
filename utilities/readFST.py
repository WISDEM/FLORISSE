def readFST(caseLocation):

    from os import listdir, path


    files = [file for file in listdir(caseLocation) if file.endswith(".fst")]
    files.sort()

    properties = {'gearbox ratio': list()}

    for file in files:

        f = open(path.join(caseLocation, file), "r")
        lines = f.readlines()
        f.close()

        lineI = 0

        def getfromline(line):
            x = line[0]
            x = float(x)
            return x

        while lineI < len(lines):
            line = lines[lineI].strip().split()
            if len(line) > 1:
                if line[1] == "GBRatio":
                    properties['gearbox ratio'].append(getfromline(line))
            lineI += 1
        
    return properties
