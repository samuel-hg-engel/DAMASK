import numpy as np
import h5py
from typing import Union, Sequence, Dict, Any, Collection

from ._typehints import FileHandle
from . import Config
from . import Rotation
from . import Orientation
from . import util
from . import Table


class ConfigMaterial(Config):
    """
    Material configuration.

    Manipulate material configurations for storage in YAML format.
    A complete material configuration file has the entries 'material',
    'phase', and 'homogenization'. For use in DAMASK, it needs to be
    stored as 'material.yaml'.

    """

    def __init__(self,
                 d: Dict[str, Any] = None,
                 **kwargs):
        """
        New material configuration.

        Parameters
        ----------
        d : dictionary or YAML string, optional
            Initial content. Defaults to None, in which case empty entries for
            any missing material, homogenization, and phase entry are created.
        kwargs : key=value pairs, optional
            Initial content specified as pairs of key=value.

        """
        default: Collection
        if d is None:
            for section, default in {'material':[],'homogenization':{},'phase':{}}.items():
                if section not in kwargs: kwargs.update({section:default})

        super().__init__(d,**kwargs)


    def save(self,
             fname: FileHandle = 'material.yaml',
             **kwargs):
        """
        Save to yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            Filename or file for writing. Defaults to 'material.yaml'.
        **kwargs
            Keyword arguments parsed to yaml.dump.

        """
        super().save(fname,**kwargs)


    @classmethod
    def load(cls,
             fname: FileHandle = 'material.yaml') -> 'ConfigMaterial':
        """
        Load from yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            Filename or file to read from. Defaults to 'material.yaml'.

        Returns
        -------
        loaded : damask.ConfigMaterial
            Material configuration from file.

        """
        return super(ConfigMaterial,cls).load(fname)


    @staticmethod
    def load_DREAM3D(fname: str,
                     grain_data: str = None,
                     cell_data: str = None,
                     cell_ensemble_data: str = 'CellEnsembleData',
                     phases: str = 'Phases',
                     Euler_angles: str = 'EulerAngles',
                     phase_names: str = 'PhaseName',
                     base_group: str = None) -> 'ConfigMaterial':
        """
        Load DREAM.3D (HDF5) file.

        Data in DREAM.3D files can be stored per cell ('CellData')
        and/or per grain ('Grain Data'). Per default, cell-wise data
        is assumed.

        damask.Grid.load_DREAM3D allows to get the corresponding geometry
        for the grid solver.

        Parameters
        ----------
        fname : str
            Filename of the DREAM.3D (HDF5) file.
        grain_data : str
            Name of the group (folder) containing grain-wise data. Defaults
            to None, in which case cell-wise data is used.
        cell_data : str
            Name of the group (folder) containing cell-wise data. Defaults to
            None in wich case it is automatically detected.
        cell_ensemble_data : str
            Name of the group (folder) containing data of cell ensembles. This
            group is used to inquire the name of the phases. Phases will get
            numeric IDs if this group is not found. Defaults to 'CellEnsembleData'.
        phases : str
            Name of the dataset containing the phase ID (cell-wise or grain-wise).
            Defaults to 'Phases'.
        Euler_angles : str
            Name of the dataset containing the crystallographic orientation as
            Euler angles in radians (cell-wise or grain-wise). Defaults to 'EulerAngles'.
        phase_names : str
            Name of the dataset containing the phase names. Phases will get
            numeric IDs if this dataset is not found. Defaults to 'PhaseName'.
        base_group : str
            Path to the group (folder) that contains geometry (_SIMPL_GEOMETRY),
            and grain- or cell-wise data. Defaults to None, in which case
            it is set as the path that contains _SIMPL_GEOMETRY/SPACING.

        Notes
        -----
        Homogenization and phase entries are emtpy and need to be defined separately.

        Returns
        -------
        loaded : damask.ConfigMaterial
            Material configuration from file.

        """
        b = util.DREAM3D_base_group(fname) if base_group is None else base_group
        c = util.DREAM3D_cell_data_group(fname) if cell_data is None else cell_data
        f = h5py.File(fname,'r')

        if grain_data is None:
            phase = f['/'.join([b,c,phases])][()].flatten()
            O = Rotation.from_Euler_angles(f['/'.join([b,c,Euler_angles])]).as_quaternion().reshape(-1,4) # noqa
            _,idx = np.unique(np.hstack([O,phase.reshape(-1,1)]),return_index=True,axis=0)
            idx = np.sort(idx)
        else:
            phase = f['/'.join([b,grain_data,phases])][()]
            O = Rotation.from_Euler_angles(f['/'.join([b,grain_data,Euler_angles])]).as_quaternion() # noqa
            idx = np.arange(phase.size)

        if cell_ensemble_data is not None and phase_names is not None:
            try:
                names = np.array([s.decode() for s in f['/'.join([b,cell_ensemble_data,phase_names])]])
                phase = names[phase]
            except KeyError:
                pass


        base_config = ConfigMaterial({'phase':{k if isinstance(k,int) else str(k): None for k in np.unique(phase)},
                                      'homogenization':{'direct':{'N_constituents':1}}})
        constituent = {k:np.atleast_1d(v[idx].squeeze()) for k,v in zip(['O','phase'],[O,phase])}

        return base_config.material_add(**constituent,homogenization='direct')


    @staticmethod
    def from_table(table: Table,
                   **kwargs) -> 'ConfigMaterial':
        """
        Generate from an ASCII table.

        Parameters
        ----------
        table : damask.Table
            Table that contains material information.
        **kwargs
            Keyword arguments where the key is the property name and
            the value specifies either the label of the data column in the table
            or a constant value.

        Returns
        -------
        new : damask.ConfigMaterial
            Material configuration from values in table.

        Examples
        --------
        >>> import damask
        >>> import damask.ConfigMaterial as cm
        >>> t = damask.Table.load('small.txt')
        >>> t
            3:pos  pos  pos   4:qu   qu    qu    qu   phase    homog
        0     0    0    0     0.19  0.8   0.24 -0.51  Aluminum SX
        1     1    0    0     0.8   0.19  0.24 -0.51  Steel    SX
        2     1    1    0     0.8   0.19  0.24 -0.51  Steel    SX
        >>> cm.from_table(t,O='qu',phase='phase',homogenization='homog')
        material:
          - constituents:
              - O: [0.19, 0.8, 0.24, -0.51]
                v: 1.0
                phase: Aluminum
            homogenization: SX
          - constituents:
              - O: [0.8, 0.19, 0.24, -0.51]
                v: 1.0
                phase: Steel
            homogenization: SX
        homogenization: {SX: null}
        phase: {Aluminum: null, Steel: null}

        >>> cm.from_table(t,O='qu',phase='phase',homogenization='single_crystal')
        material:
          - constituents:
              - O: [0.19, 0.8, 0.24, -0.51]
                v: 1.0
                phase: Aluminum
            homogenization: single_crystal
          - constituents:
              - O: [0.8, 0.19, 0.24, -0.51]
                v: 1.0
                phase: Steel
            homogenization: single_crystal
        homogenization: {single_crystal: null}
        phase: {Aluminum: null, Steel: null}

        """
        kwargs_ = {k:table.get(v) if v in table.labels else np.atleast_2d([v]*len(table)).T for k,v in kwargs.items()}

        _,idx = np.unique(np.hstack(list(kwargs_.values())),return_index=True,axis=0)
        idx = np.sort(idx)
        kwargs_ = {k:np.atleast_1d(v[idx].squeeze()) for k,v in kwargs_.items()}
        for what in ['phase','homogenization']:
            if what not in kwargs_: kwargs_[what] = what+'_label'

        return ConfigMaterial().material_add(**kwargs_)


    @property
    def is_complete(self) -> bool:
        """
        Check for completeness.

        Only the general file layout is considered.
        This check does not consider whether specific parameters for
        a particular phase/homogenization model are missing.

        Returns
        -------
        complete : bool
            Whether the material.yaml definition is complete.

        """
        def LabeledList(label,items):
            return f'{label.capitalize()}{"s" if len(items)>1 else ""} {util.srepr(items,",",quote=True)}'

        ok = True
        msg = []
        all = set(['homogenization','phase','material'])
        miss = set([item for item in all if item not in self])
        empty = set([item for item in all-miss if self[item] is None])

        if miss:
            msg.append(f'{LabeledList("top-level",miss)} missing')
            ok = False
        if empty:
            msg.append(f'{LabeledList("top-level",empty)} empty')

        if ok:
            ok &= len(self['material']) > 0
            if len(self['material']) < 1: msg.append('No materials defined')

            homogenization = set()
            phase          = set()
            for i,v in enumerate(self['material']):
                if 'homogenization' in v:
                    homogenization.add(v['homogenization'])
                else:
                    msg.append(f'No homogenization specified for material {i}')
                    ok = False

                if 'constituents' in v:
                    for ii,vv in enumerate(v['constituents']):
                        if 'O' not in vv:
                            msg.append(f'No orientation specified for constituent {ii} of material {i}')
                            ok = False
                        if 'phase' in vv:
                            phase.add(vv['phase'])
                        else:
                            msg.append(f'No phase specified for constituent {ii} of material {i}')
                            ok = False

            for v,other in {'phase':phase,
                            'homogenization':homogenization}.items():
                me = set([] if v in empty else self[v])
                if _miss := other - me:
                    msg.append(f'{LabeledList(v,_miss)} missing')
                    ok = False
                if len(_empty := [item for item in me if self[v][item] is None]) > 0:
                    msg.append(f'{LabeledList(v,_empty)} undefined')
                    ok = False

        print(util.srepr(msg))
        return ok


    @property
    def is_valid(self) -> bool:
        """
        Check for valid content.

        Only the generic file content is considered.
        This check does not consider whether parameters for a
        particular phase/homogenization mode are out of bounds.

        Returns
        -------
        valid : bool
            Whether the material.yaml definition is valid.

        """
        ok = True

        if 'phase' in self:
            for k,v in self['phase'].items():
                if v is not None and 'lattice' in v:
                    try:
                        Orientation(lattice=v['lattice'])
                    except KeyError:
                        print(f"Invalid lattice '{v['lattice']}' in phase '{k}'")
                        ok = False

        if 'material' in self:
            for i,m in enumerate(self['material']):
                if 'constituents' in m:
                    v = 0.0
                    for c in m['constituents']:
                        v += float(c['v'])
                        if 'O' in c:
                            try:
                                Rotation.from_quaternion(c['O'])
                            except ValueError:
                                print(f"Invalid orientation '{c['O']}' in material '{i}'")
                                ok = False
                    if not np.isclose(v,1.0):
                        print(f"Total fraction v = {v} ≠ 1 in material '{i}'")
                        ok = False

        return ok


    def material_rename_phase(self,
                              mapping: Dict[str, str],
                              ID: Sequence[int] = None,
                              constituent: Sequence[int] = None) -> 'ConfigMaterial':
        """
        Change phase name in material.

        Parameters
        ----------
        mapping: dictionary
            Mapping from old name to new name
        ID: list of ints, optional
            Limit renaming to selected material IDs.
        constituent: list of ints, optional
            Limit renaming to selected constituents.

        Returns
        -------
        updated : damask.ConfigMaterial
            Updated material configuration.

        """
        dup = self.copy()
        for i,m in enumerate(dup['material']):
            if ID is not None and i not in ID: continue
            for c in m['constituents']:
                if constituent is not None and c not in constituent: continue
                try:
                    c['phase'] = mapping[c['phase']]
                except KeyError:
                    continue
        return dup


    def material_rename_homogenization(self,
                                       mapping: Dict[str, str],
                                       ID: Sequence[int] = None) -> 'ConfigMaterial':
        """
        Change homogenization name in material.

        Parameters
        ----------
        mapping: dictionary
            Mapping from old name to new name
        ID: list of ints, optional
            Limit renaming to selected homogenization IDs.

        Returns
        -------
        updated : damask.ConfigMaterial
            Updated material configuration.

        """
        dup = self.copy()
        for i,m in enumerate(dup['material']):
            if ID is not None and i not in ID: continue
            try:
                m['homogenization'] = mapping[m['homogenization']]
            except KeyError:
                continue
        return dup


    def material_add(self,
                     **kwargs: Any) -> 'ConfigMaterial':
        """
        Add material entries.

        Parameters
        ----------
        **kwargs
            Key-value pairs.
            First index of array-like values runs over materials,
            whereas second index runs over constituents.

        Returns
        -------
        updated : damask.ConfigMaterial
            Updated material configuration.

        Examples
        --------
        Create two grains of ferrite and one grain of martensite, each with random orientation:

        >>> import damask
        >>> m = damask.ConfigMaterial()
        >>> m = m.material_add(phase = ['Ferrite','Martensite','Ferrite'],
        ...                    O = damask.Rotation.from_random(3),
        ...                    homogenization = 'SX')
        >>> m
        material:
          - constituents:
              - O: [0.577764, -0.146299, -0.617669, 0.513010]
                v: 1.0
                phase: Ferrite
            homogenization: SX
          - constituents:
              - O: [0.184176, 0.340305, 0.737247, 0.553840]
                v: 1.0
                phase: Martensite
            homogenization: SX
          - constituents:
              - O: [0.47925185, -0.04294454, 0.78760173, -0.3849116 ]
                v: 1.0
                phase: Ferrite
            homogenization: SX
        homogenization: {SX: null}
        phase: {Ferrite: null, Martensite: null}

        Create hundred materials that each approximate a duplex stainless steel microstructure
        with three austenite and one relatively bigger ferrite grain of random orientation each:

        >>> import damask
        >>> m = damask.ConfigMaterial()
        >>> m = m.material_add(phase = np.array(['Austenite']*3+['Ferrite']),
        ...                    O = damask.Rotation.from_random((100,4)),
        ...                    v = np.array([0.2]*3+[0.4]),
        ...                    homogenization = 'Taylor')
        >>> m
        material:
          - constituents:
              - v: 0.2
                phase: Austenite
                O: [0.46183665006602664, 0.2215160420973196, -0.5594313187331139, 0.6516702781083836]
              - v: 0.2
                phase: Austenite
                O: [0.11321658382410027, 0.6354079414360444, 0.00562701344273936, 0.7638108992590535]
              - v: 0.2
                phase: Austenite
                O: [0.050991978809077604, 0.8069522034362003, -0.11352928955610851, -0.5773552285027659]
              - v: 0.4
                phase: Ferrite
                O: [0.9460076150721788, 0.15880754622367604, -0.0069841062241482385, -0.28249066842661014]
            homogenization: Taylor
          .
          .
          .
          - constituents:
              - v: 0.2
                phase: Austenite
                O: [0.12531400788494199, -0.18637769037997565, 0.31737548053338394, -0.9213210951197429]
              - v: 0.2
                phase: Austenite
                O: [0.37453930577161404, -0.33529507696450805, -0.3266564259130028, -0.800370601162502]
              - v: 0.2
                phase: Austenite
                O: [0.035776891752713764, -0.720706371010592, -0.4540438656728926, -0.5226342017569017]
              - v: 0.4
                phase: Ferrite
                O: [0.6782596727966124, -0.20800082041703685, -0.138636083554039, 0.6909989227925536]
            homogenization: Taylor

        homogenization: {Taylor: null}

        phase: {Austenite: null, Ferrite: null}

        """
        _constituent_properties = ['phase','O','v','V_e']
        _dim = {'O':(4,),'V_e':(3,3,)}
        _ex = dict((k, -len(v)) for k, v in _dim.items())

        N,n = 1,1
        shaped : Dict[str, Union[None,np.ndarray]] = \
                 {'v': None,
                  'phase': None,
                  'homogenization': None,
                  }

        for k,v in kwargs.items():
            shaped[k] = np.array(v)
            s = shaped[k].shape[:_ex.get(k,None)]                               # type: ignore
            N = max(N,s[0]) if len(s)>0 else N
            n = max(n,s[1]) if len(s)>1 else n

        shaped['v'] = np.array(1./n) if shaped['v'] is None else shaped['v']

        mat: Sequence[dict] = [{'constituents':[{} for _ in range(n)]} for _ in range(N)]

        for k,v in shaped.items():
            target = (N,n) + _dim.get(k,())
            obj = np.broadcast_to(np.array(v).reshape(util.shapeshifter(() if v is None else v.shape,
                                                                        target,
                                                                        mode = 'right')),
                                  target)
            for i in range(N):
                if k in _constituent_properties:
                    for j in range(n):
                        mat[i]['constituents'][j][k] = obj[i,j].item() if isinstance(obj[i,j],np.generic) else obj[i,j]
                else:
                    mat[i][k] = obj[i,0].item() if isinstance(obj[i,0],np.generic) else obj[i,0]

        dup = self.copy()
        dup['material'] = dup['material'] + mat if 'material' in dup else mat

        for what in [item for item in ['phase','homogenization'] if shaped[item] is not None]:
            for k in np.unique(shaped[what]):                                   # type: ignore
                if k not in dup[what]: dup[what][str(k)] = None

        return dup
