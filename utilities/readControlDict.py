def readControlDict(file = 'controlDict'):

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
        if "startTime" in line and not "startFrom" in line:
            properties['start time'] = getfromline(line)
        lineI += 1

    return properties
