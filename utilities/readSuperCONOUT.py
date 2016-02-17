import numpy as np


class SuperCONOUT:
    def __init__(self, time, data):
        self.time = time
        self.data = data
        self.nTurbines = len(self.data)
        self.nOutputs = len(self.data[0].keys())
        self.outputNames = self.data[0].keys()


def readSuperCONOUT(filename='superCONOUT.csv', createDictionary = True):
    """imports standard SOWFA superCONOUT.csv files

    input: file = location of .csv-file
    output:
    SCO is object with following attributes:
    SCO.nTurbines = number of turbines
    SCO.nOutputs = number of output signals at each turbine
    SCO.outputNames = names of output signals
    SCO.time = time vector (s)
    SCO.data = list of dictionaries, one dictionary with all signals per turbine, example:
                SCO.data[0]['Power W'] = power of first turbine (numpy vector)
               see SCO.outputNames for list of all signal names

    see example below

    Pieter Gebraad, 2015"""

    file = open(filename, 'r')
    lines = file.readlines()
    file.close()

    labels = filter(None,lines[0].strip().split(',')) # read column labels
    labels = labels[1:]  # remove "Time (s) from labels"
    nTurbines = int(labels[-1].split(':')[0][1:]) # read number of turbines from last label
    labels = labels[::nTurbines] # one label per output
    nOutputs = len(labels)

    # remove turbine names from labels
    labels = [label.split(':')[1:] for label in labels]
    outputNames = [''.join(label) for label in labels]

    # convert read lines to data matrix and time vector
    lines = np.array([line.split(',') for line in lines[1:]]) # split columns to make matrix
    lines = np.delete(lines, -1, 1)                           # remove last column ('\n')
    #lines.astype(np.float)                                    # convert to floats PF edit, this line doesn't do anything...

    time = lines[:, 0].astype(float)
    data = lines[:, 1:].astype(float)

    # split data into several matrices, each corresponding to a turbine
    indicesList = [turbineI + np.array(range(0, nOutputs))*nTurbines for turbineI in range(0, nTurbines)]
    data = np.array([data.take(indices, 1) for indices in indicesList])

    if createDictionary:
        SCO = SuperCONOUT(time, [dict(zip(outputNames, [data[turbineI][:, i] for i in range(0, nOutputs)])) for turbineI in range(0, nTurbines)])
        return SCO
    else:
        return time, data, outputNames, nTurbines, nOutputs

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    signalsToPlot = ['Power (W)', 'Blade 1 Pitch (Rad)', 'Yaw Angle (Degrees)']

    SCO = readSuperCONOUT()
    fig, axes = plt.subplots(ncols=SCO.nTurbines, nrows=len(signalsToPlot))

    formatter = ticker.ScalarFormatter()
    formatter.set_powerlimits((-3, 3))

    for signalI in range(0, len(signalsToPlot)):
        signal = signalsToPlot[signalI]
        for turbineI in range(0,SCO.nTurbines):
            axes[signalI, turbineI].plot(SCO.time, SCO.data[turbineI][signal])
            axes[signalI, turbineI].set_xlabel('time (s)')
            axes[signalI, turbineI].set_ylabel(signal)

    for turbineI in range(0, SCO.nTurbines):
        axes[0, turbineI].yaxis.set_major_formatter(formatter)

    plt.show()




