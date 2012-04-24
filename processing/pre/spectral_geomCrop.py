#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP


# -----------------------------
class extendedOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
    
    ACTIONS = Option.ACTIONS + ("extend",)
    STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
    TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
    ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

    def take_action(self, action, dest, opt, value, values, parser):
        if action == "extend":
            lvalue = value.split(",")
            values.ensure_value(dest, []).extend(lvalue)
        else:
            Option.take_action(self, action, dest, opt, value, values, parser)


# ----------------------- MAIN -------------------------------

identifiers = {
        'resolution': ['a','b','c'],
        'dimension':  ['x','y','z'],
          }
mappings = {
        'resolution': lambda x: int(x),
        'dimension':  lambda x: float(x),
          }


parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Crops a spectral geometry description at given offset and by given points.
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-r', '--resolution', dest='resolution', type='int', nargs = 3, \
                  help='a,b,c resolution of cropped box')
parser.add_option('-o', '--offset', dest='offset', type='int', nargs = 3, \
                  help='a,b,c offset of cropped box')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement')

parser.set_defaults(resolution = [0,0,0])
parser.set_defaults(offset = [0,0,0])
parser.set_defaults(twoD = False)

(options, filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w')})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': print file['name']

  #  get labels by either read the first row, or - if keyword header is present - the last line of the header

  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    headers  = [firstline]+[file['input'].readline() for i in range(headerlines)]
  else:
    headerlines = 1
    headers = firstline

  content = file['input'].readlines()
  file['input'].close()
  if 0 in options.resolution:  
    for header in headers:
      headitems = header.split()
      if headitems[0] == 'resolution':         # located resolution entry
        for i in xrange(3):
          options.resolution[i] = mappings['resolution'](headitems[headitems.index(identifiers['resolution'][i])+1])

  resolution = [0,0,0]
  dimension = [0.0,0.0,0.0]
  new_header = []
  for header in headers:
    headitems = header.split()
    if headitems[0] == 'resolution':         # located resolution entry
      for i in xrange(3):
        resolution[i] = mappings['resolution'](headitems[headitems.index(identifiers['resolution'][i])+1])
      header = "resolution\ta %i\tb %i\tc %i\n"%(options.resolution[0],options.resolution[1],options.resolution[2],)
    if headitems[0] == 'dimension':          # located dimension entry
      for i in xrange(3):
        dimension[i] = mappings['dimension'](headitems[headitems.index(identifiers['dimension'][i])+1])
      header = "dimension\tx %f\ty %f\tz %f\n"%(dimension[0]/resolution[0]*options.resolution[0],
                                              dimension[1]/resolution[1]*options.resolution[1],
                                              dimension[2]/resolution[2]*options.resolution[2],)
    
    new_header.append(header)
    
  if resolution == [0,0,0]:
    print 'no resolution info found.'
    sys.exit(1)
  if dimension == [0.0,0.0,0.0]:
    print 'no dimension info found.'
    sys.exit(1)

  if file['name'] != 'STDIN':
    print 'resolution: %s'%(' x '.join(map(str,resolution)))
    print 'dimension:  %s'%(' x '.join(map(str,dimension)))
    
  microstructure = numpy.zeros(resolution,'i')
  i = 0
  for line in content:  
    for item in map(int,line.split()):
      microstructure[i%resolution[0],(i/resolution[0])%resolution[1],i/resolution[0]/resolution[1]] = item
      i += 1
  
  formatwidth = int(math.floor(math.log10(microstructure.max())))
  

            
# ------------------------------------------ assemble header ---------------------------------------  

  output = ''.join(new_header)

# ------------------------------------- regenerate texture information ----------------------------------  

  for z in xrange(options.offset[2],options.offset[2]+options.resolution[2]):
    for y in xrange(options.offset[1],options.offset[1]+options.resolution[1]):
      output += {True:' ',False:'\n'}[options.twoD].join(map(lambda x: ('%%%ii'%formatwidth)%x, microstructure[options.offset[0]:options.offset[0]+options.resolution[0],y,z])) + '\n'
    
    #output += '\n'
    
# ------------------------------------------ output result ---------------------------------------  

  file['output'].write(output)

  if file['name'] != 'STDIN':
    file['output'].close
    os.rename(file['name']+'_tmp',file['name'])
    