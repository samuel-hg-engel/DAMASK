#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add data in column(s) of second ASCIItable selected from row that is given by the value in a mapping column.

""", version = scriptID)

parser.add_option('-c','--map',
                  dest = 'map',
                  type = 'string', metavar = 'string',
                  help = 'heading of column containing row mapping')
parser.add_option('-o','--offset',
                  dest = 'offset',
                  type = 'int', metavar = 'int',
                  help = 'offset between mapping column value and actual row in mapped table [%default]')
parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help='heading of column(s) to be mapped')
parser.add_option('-a','--asciitable',
                  dest = 'asciitable',
                  type = 'string', metavar = 'string',
                  help = 'mapped ASCIItable')

parser.set_defaults(offset = 0,
                   )

(options,filenames) = parser.parse_args()

if options.label == None:
  parser.error('no data columns specified.')
if options.map == None:
  parser.error('no mapping column given.')

# ------------------------------------------ process mapping ASCIItable ---------------------------

if options.asciitable != None and os.path.isfile(options.asciitable):

  mappedTable = damask.ASCIItable(name = options.asciitable,buffered = False, readonly = True) 
  mappedTable.head_read()                                                                           # read ASCII header info of mapped table
  missing_labels = mappedTable.data_readArray(options.label)

  if len(missing_labels) > 0:
    mappedTable.croak('column{} {} not found...'.format('s' if len(missing_labels) > 1 else '',', '.join(missing_labels)))

else:
  parser.error('no mapped ASCIItable given.')

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

  errors = []

  mappedColumn = table.label_index(options.map)  
  if mappedColumn <  0: errors.append('mapping column {} not found.'.format(options.map))

  if errors != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(mappedTable.labels)                                                           # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    table.data_append(mappedTable.data[int(table.data[mappedColumn])+options.offset-1])             # add all mapped data types
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
  if name != 'STDIN': os.rename(name+'_tmp',name)                                                   # overwrite old one with tmp new

mappedTable.close()                                                                                 # close mapped input ASCII table
