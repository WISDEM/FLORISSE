from matplotlib import pyplot as plt
import numpy as np
import re
from readSuperCONOUT import readSuperCONOUT

filenames=['SScon_T1.txt', 'CONOUT/SScon_T1.txt']

ssmat = {'A':list(), 'B':list(), 'C':list(), 'D':list()}
file = open(filenames[0], 'r')
lines = file.readlines()
file.close()

mat = ''
for line in lines:
    if line[:21] ==   'A (nStates * nStates)':
        mat = 'A'
    elif line[:21] == 'B (nStates * nInputs)':
        mat = 'B'
    elif line[:22] == 'C (nOutputs * nStates)':
        mat = 'C'
    elif line[:22] == 'D (nOutputs * nInputs)':
        mat = 'D'
    elif mat != '':
        line = np.array([np.float(v) for v in re.split(' |\t',line)])
        ssmat[mat].append(line)

for mat in ssmat.keys():
    ssmat[mat] = np.array(ssmat[mat])
    print mat
    print ssmat[mat]

file = open(filenames[1], 'r')
lines = file.readlines()
file.close()

time = list()
opI = list()    # timeseries of operating point
opO = list()    # timeseries of operating point
s = list()    # timeseries of state
i = list()    # timeseries of input
o = list()    # timeseries of output

count = -1

for line in lines:
    if line[0] == 't':
        time.append(np.float(line.split('t ')[1].split('\r\n')[0]))
        opI.append(None);
        opO.append(None);
        s.append(None);
        i.append(None);
        o.append(None);
        count += 1
    elif line[:4] == 'op i':
        opI[count] = line.split('op i ')[1].split('\r\n')[0].split(',')[:-1]
        opI[count] = np.array([np.float(v) for v in opI[count]])
    elif line[:4] == 'op o':
        opO[count] = line.split('op o ')[1].split('\r\n')[0].split(',')[:-1]
        opO[count] = np.array([np.float(v) for v in opO[count]])
    elif line[0] == 's':
        s[count] = line.split('s ')[1].split('\r\n')[0].split(',')[:-1]
        s[count] = np.array([np.float(v) for v in s[count]])
    elif line[0] == 'i':
        i[count] = line.split('i ')[1].split('\r\n')[0].split(',')[:-1]
        i[count] = np.array([np.float(v) for v in i[count]])
    elif line[0] == 'o':
        o[count] = line.split('o ')[1].split('\r\n')[0].split(',')[:-1]
        o[count] = np.array([np.float(v) for v in o[count]])

time = np.array(time)
opI = np.array(opI)
opO = np.array(opO)

timeSS = np.array([time[j] for j in range(len(time)) if s[j]!=None])
opIss = np.array([opI[j] for j in range(len(time)) if s[j]!=None])
opOss = np.array([opO[j] for j in range(len(time)) if s[j]!=None])
o = np.array([o[j] for j in range(len(time)) if s[j]!=None])
i = np.array([i[j] for j in range(len(time)) if s[j]!=None])
s = np.array([s[j] for j in range(len(time)) if s[j]!=None])

# mimic state-space system
(nk,ns) = s.shape
(nk,no) = o.shape
s_m = np.zeros((nk+1,ns))
o_m = np.zeros((nk,no))
for k in range(nk):
    o_m[k] = np.dot(ssmat['C'],s_m[k])+np.dot(ssmat['D'],i[k])
    s_m[k+1] = np.dot(ssmat['A'],s_m[k])+np.dot(ssmat['B'],i[k])
s_m = s_m[1:,:]

plt.figure()
plt.plot(timeSS,s)
plt.plot(timeSS,s_m,'--')
plt.title('s')
plt.figure()
plt.plot(timeSS,i+opIss)
plt.plot(timeSS,i+opIss)
plt.title('i')

plt.figure()
plt.plot(time,opI)
plt.title('op i')
plt.figure()
plt.plot(time,opO)
plt.title('op o')

plt.figure()
plt.plot(timeSS,o[:,0]+opOss[:,0])
plt.plot(timeSS,o_m[:,0]+opOss[:,0],'k--')
SCO = readSuperCONOUT('CONOUT/superCONOUT.csv')
plt.plot(SCO.time, SCO.data[0]['Blade 1 Pitch (Rad)'],'r--')


plt.figure()
plt.plot(timeSS,o[:,1]+opOss[:,1])
plt.plot(timeSS,o_m[:,1]+opOss[:,1],'k--')
SCO = readSuperCONOUT('CONOUT/superCONOUT.csv')
plt.plot(SCO.time, SCO.data[0]['TorqueScaling (-)'],'r--')

plt.title('o')



plt.show()