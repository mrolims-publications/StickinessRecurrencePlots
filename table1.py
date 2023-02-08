import numpy as np
import pandas as pd
import os
from functions import corr_coef

cc = np.zeros(4)

path = '/home/matheus/Doutorado/StickinessRecurrencePlots/Data/'

datafile = path + 'fig3.dat'
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig3.py?\nStopping execution...' % datafile)
    sys.exit()
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
k = np.array(df[0])
lyap = np.array(df[1])
rte = np.array(df[2])

print('Evaluating the correlation coefficient of Figure 3...')
cc[0] = corr_coef(lyap, rte)
print('Correlation coefficient: %.2f\n' % cc[0])

datafile = path + 'fig4a.dat'
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig4ab.py a?\nStopping execution...' % datafile)
    sys.exit()
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
lyap = np.array(df[2])
datafile = path + 'fig4d.dat'
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig4de.py d?\nStopping execution...' % datafile)
    sys.exit()
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
rte = np.array(df[2])

print('Evaluating the correlation coefficient of Figure 4(a) and 4(d)...')
cc[1] = corr_coef(lyap, rte)
print('Correlation coefficient: %.2f\n' % cc[1])

datafile = path + 'fig4b.dat'
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig4ab.py b?\nStopping execution...' % datafile)
    sys.exit()
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
lyap = np.array(df[2])
datafile = path + 'fig4e.dat'
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig4de.py e?\nStopping execution...' % datafile)
    sys.exit()
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
rte = np.array(df[2])

print('Evaluating the correlation coefficient of Figure 4(b) and 4(e)...')
cc[2] = corr_coef(lyap, rte)
print('Correlation coefficient: %.2f\n' % cc[2])

datafile = path + 'fig4c.dat'
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig4c.py?\nStopping execution...' % datafile)
    sys.exit()
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
lyap = np.array(df[2])
datafile = path + 'fig4f.dat'
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig4f.py?\nStopping execution...' % datafile)
    sys.exit()
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
rte = np.array(df[2])

print('Evaluating the correlation coefficient of Figure 4(c) and 4(f)...')
cc[2] = corr_coef(lyap, rte)
print('Correlation coefficient: %.2f\n' % cc[2])