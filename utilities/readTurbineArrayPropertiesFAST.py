def readTurbineArrayPropertiesFAST(file = 'TurbineArrayPropertiesFAST'):

    f = open(file, "r")
    lines = f.readlines()
    f.close()

    turbineIndices = list()
    properties = {'x': list(), 'y': list(), 'z': list(), 'hubZ': list(), 'rotorDiameter': list(), 'yaw': list()}
    lineI = 0

    def getfromline(line):
        x = line.split(';')
        x = x[0].split()
        x = float(x[1])
        return x

    while lineI < len(lines):
        line = lines[lineI].strip()
        if "turbine" in line and line[7].isdigit():
            turbineI = int(filter(str.isdigit, line))
            turbineIndices.append(turbineI)
        if "refx" in line:
            properties['x'].append(getfromline(line))
        if "refy" in line:
            properties['y'].append(getfromline(line))
        if "refz" in line:
            properties['z'].append(getfromline(line))
        if "hubz" in line:
            properties['hubZ'].append(getfromline(line))
        if "rotorDiameter" in line:
            properties['rotorDiameter'].append(getfromline(line))
        if "yawAngle" in line:
            properties['yaw'].append(getfromline(line))
        lineI += 1

    nTurbines = len(turbineIndices)

    # resort if turbines not in order, extend to each turbine if in general
    for key in properties.keys():
        if nTurbines > 1 and len(properties[key]) == 1:
            properties[key] = properties[key] * nTurbines
        properties[key] = [properties[key][i] for i in turbineIndices]

    return properties
