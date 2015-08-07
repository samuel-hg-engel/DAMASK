#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

def normalize(vec):
    return vec/np.sqrt(np.inner(vec,vec))

def E_hkl(stiffness,vec):   # stiffness = (c11,c12,c44)
    v = normalize(vec)
    S11 = (stiffness[0]+stiffness[1])/(stiffness[0]*stiffness[0]+stiffness[0]*stiffness[1]-2.0*stiffness[1]*stiffness[1])
    S12 = (            -stiffness[1])/(stiffness[0]*stiffness[0]+stiffness[0]*stiffness[1]-2.0*stiffness[1]*stiffness[1])
    S44 = 1.0/stiffness[2]

    invE = S11-(S11-S12-0.5*S44)* (1.0 - \
                                         (v[0]**4+v[1]**4+v[2]**4) \
                                    /#------------------------------------
                                         np.inner(v,v)**2 \
                                  )

    return 1.0/invE

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing directional stiffness based on given cubic stiffness values C11, C12, and C44 in consecutive columns.

""", version = scriptID)

parser.add_option('-c','--stiffness',
                  dest = 'stiffness',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'heading of column containing C11 (followed by C12, C44) field values')
parser.add_option('-d','--direction','--hkl',
                  dest = 'hkl',
                  type = 'int', nargs = 3, metavar = 'int int int',
                  help = 'direction of elastic modulus [%default]')
parser.set_defaults(hkl = (1,1,1),
                    )

(options,filenames) = parser.parse_args()

if options.stiffness == None:
  parser.error('no data column specified...')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = ['STDIN']

for name in filenames:
  if not (name == 'STDIN' or os.path.exists(name)): continue
  table = damask.ASCIItable(name = name, outname = name+'_tmp',
                            buffered = False)
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name != 'STDIN' else ''))

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  remarks = []
  columns = []
  
  for i,column in enumerate(table.label_index(options.stiffness)):
    if   column <  0: remarks.append('column {} not found.'.format(options.stiffness[i]))
    else:
      columns.append(column)
      table.labels_append(['E{}{}{}({})'.format(*options.hkl,options.stiffness[i]))                 # extend ASCII header with new labels

  if remarks != []: table.croak(remarks)

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for column in columns:
      table.data_append(E_hkl(map(float,table.data[column:column+3]),options.hkl))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
  if name != 'STDIN': os.rename(name+'_tmp',name)                                                   # overwrite old one with tmp new
