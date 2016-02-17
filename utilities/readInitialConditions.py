def readInitialConditions(file = 'initialConditions'):

    f = open(file, "r")
    lines = f.readlines()
    f.close()
    
    turbineIndices = list()
    properties = dict()
    lineI = 0

    def getfromline(line):
        x = line.split(';')
        x = x[0].split()
        x = float(x[1])
        return x

    while lineI < len(lines):
        line = lines[lineI].strip()
        if "Ug" in line:
            properties['wind speed'] = getfromline(line)
        if "Udir" in line:
            properties['wind direction'] = getfromline(line)
        lineI += 1

    return properties
