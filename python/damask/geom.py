import math
from io import StringIO

import numpy as np

from . import util


class Geom():
  """Geometry definition for grid solvers"""

  def __init__(self,microstructure,size,origin=[0.0,0.0,0.0],homogenization=1,comments=[]):
    """New geometry definition from array of microstructures and size"""
    self.set_size(size)
    self.set_origin(origin)
    self.set_microstructure(microstructure)
    self.set_homogenization(homogenization)
    self.set_comments(comments)
    

  def __repr__(self):
    """Basic information on geometry definition"""
    return util.srepr([
           'grid     a b c:      {}'.format(' x '.join(map(str,self.get_grid  ()))),
           'size     x y z:      {}'.format(' x '.join(map(str,self.get_size  ()))),
           'origin   x y z:      {}'.format(' x '.join(map(str,self.get_origin()))),
           'homogenization:      {}'.format(self.get_homogenization()),
           '# microstructures:   {}'.format(len(np.unique(self.microstructure))),
           'max microstructure:  {}'.format(np.nanmax(self.microstructure)),
          ])

  def update(self,microstructure=None,size=None,origin=None,rescale=False):
    """Updates microstructure and size"""
    grid_old    = self.get_grid()
    size_old    = self.get_size()
    origin_old  = self.get_origin()
    unique_old  = len(np.unique(self.microstructure))
    max_old     = np.nanmax(self.microstructure)
    
    if size is not None and rescale:
      raise ValueError('Either set size explicitly or rescale automatically')

    self.set_microstructure(microstructure)
    self.set_size(self.get_grid()/grid_old*self.size if rescale else size)
    self.set_origin(origin)
   
    message = ['grid     a b c:      {}'.format(' x '.join(map(str,grid_old)))]
    if np.any(grid_old != self.get_grid()):
      message[-1] = util.delete(message[-1])
      message.append('grid     a b c:      {}'.format(' x '.join(map(str,self.get_grid()))))

    message.append('size     x y z:      {}'.format(' x '.join(map(str,size_old))))
    if np.any(size_old != self.get_size()):
      message[-1] = util.delete(message[-1])
      message.append('size     x y z:      {}'.format(' x '.join(map(str,self.get_size()))))

    message.append('origin   x y z:      {}'.format(' x '.join(map(str,origin_old))))
    if np.any(origin_old != self.get_origin()):
      message[-1] = util.delete(message[-1])
      message.append('origin   x y z:      {}'.format(' x '.join(map(str,self.get_origin()))))

    message.append('homogenization:      {}'.format(self.get_homogenization()))

    message.append('# microstructures:   {}'.format(unique_old))
    if unique_old != len(np.unique(self.microstructure)):
      message[-1] = util.delete(message[-1])
      message.append('# microstructures:   {}'.format(len(np.unique(self.microstructure))))

    message.append('max microstructure:  {}'.format(max_old))
    if max_old != np.nanmax(self.microstructure):
      message[-1] = util.delete(message[-1])
      message.append('max microstructure:  {}'.format(np.nanmax(self.microstructure)))
    
    return util.return_message(message)

  def set_comments(self,comments):
    self.comments = []
    self.add_comments(comments)
    
  def add_comments(self,comments):
    self.comments += [str(c) for c in comments] if isinstance(comments,list) else [str(comments)]

  def set_microstructure(self,microstructure):
    if microstructure is not None:
      if len(microstructure.shape) != 3:
        raise ValueError('Invalid microstructure shape {}'.format(*microstructure.shape))
      elif microstructure.dtype not in [int,float]:
        raise TypeError('Invalid data type {} for microstructure'.format(microstructure.dtype))
      else:
        self.microstructure = np.copy(microstructure)

  def set_size(self,size):
    if size is not None:
      if len(size) != 3 or any(np.array(size)<=0):
        raise ValueError('Invalid size {}'.format(*size))
      else:
        self.size = np.array(size)

  def set_origin(self,origin):
    if origin is not None:
      if len(origin) != 3:
        raise ValueError('Invalid origin {}'.format(*origin))
      else:
        self.origin = np.array(origin)

  def set_homogenization(self,homogenization):
    if homogenization is not None:
      if not isinstance(homogenization,int) or homogenization < 1:
        raise TypeError('Invalid homogenization {}'.format(homogenization))
      else:
        self.homogenization = homogenization


  def get_microstructure(self):
    return np.copy(self.microstructure)

  def get_size(self):
    return np.copy(self.size)

  def get_origin(self):
    return np.copy(self.origin)

  def get_grid(self):
    return np.array(self.microstructure.shape)

  def get_homogenization(self):
    return self.homogenization

  def get_comments(self):
    return self.comments[:]

  def get_header(self):
    header =  ['{} header'.format(len(self.comments)+4)] + self.comments
    header.append('grid   a {} b {} c {}'.format(*self.get_grid()))
    header.append('size   x {} y {} z {}'.format(*self.get_size()))
    header.append('origin x {} y {} z {}'.format(*self.get_origin()))
    header.append('homogenization {}'.format(self.get_homogenization()))
    return header
      
  @classmethod
  def from_file(cls,fname):
    """Reads a geom file"""
    with (open(fname) if isinstance(fname,str) else fname) as f:
      f.seek(0)
      header_length,keyword = f.readline().split()[:2]
      header_length = int(header_length)
      content = f.readlines()

    if not keyword.startswith('head') or header_length < 3:
      raise TypeError('Header length information missing or invalid')

    comments = [] 
    for i,line in enumerate(content[:header_length]):
      items = line.lower().strip().split()
      key = items[0] if len(items) > 0 else ''
      if   key == 'grid':
        grid   = np.array([  int(dict(zip(items[1::2],items[2::2]))[i]) for i in ['a','b','c']])
      elif key == 'size':
        size   = np.array([float(dict(zip(items[1::2],items[2::2]))[i]) for i in ['x','y','z']])
      elif key == 'origin':
        origin = np.array([float(dict(zip(items[1::2],items[2::2]))[i]) for i in ['x','y','z']])
      elif key == 'homogenization':
        homogenization = int(items[1])
      else:
        comments.append(line.strip())

    microstructure = np.empty(grid.prod())                                                          # initialize as flat array
    i = 0
    for line in content[header_length:]:
      items = line.split()
      if len(items) == 3:
        if   items[1].lower() == 'of':
          items = np.ones(int(items[0]))*float(items[2])
        elif items[1].lower() == 'to':
          items = np.linspace(int(items[0]),int(items[2]),
                              abs(int(items[2])-int(items[0]))+1,dtype=float)
        else:                          items = list(map(float,items))
      else:                            items = list(map(float,items))
      
      microstructure[i:i+len(items)] = items
      i += len(items)
    
    if i != grid.prod():
      raise TypeError('Invalid file: expected {} entries,found {}'.format(grid.prod(),i))
    
    microstructure = microstructure.reshape(grid,order='F')
    if not np.any(np.mod(microstructure.flatten(),1) != 0.0):                                     # no float present
      microstructure = microstructure.astype('int')
    
    return cls(microstructure.reshape(grid),size,origin,homogenization,comments)

  def to_file(self,fname):
    """Writes to file"""
    header = self.get_header()
    grid   = self.get_grid()
    format_string = '%{}i'.format(int(math.floor(math.log10(self.microstructure.max())+1))) if self.microstructure.dtype == int \
               else '%g'
    np.savetxt(fname,
               self.microstructure.reshape([grid[0],np.prod(grid[1:])],order='F').T,
               header='\n'.join(header), fmt=format_string, comments='')
               
  def show(self):
    """Show raw content (as in file)"""
    f=StringIO()
    self.to_file(f)
    f.seek(0)
    return ''.join(f.readlines())
