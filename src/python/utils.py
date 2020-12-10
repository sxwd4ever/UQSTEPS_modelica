'''
    Utility class and method for Tests
'''
from numpy.core.fromnumeric import transpose
from numpy.core.function_base import linspace
from numpy.core.getlimits import _title_fmt
import pandas as pd
import xlwings as xw
import numpy as np
from xlwings.main import Range
from typing import Dict, Tuple

class ExcelHelper(object):
    """
    Helper class for reading and writing excels
    """

    def __init__(self, sht:xw.Sheet):
        super().__init__()
        self._sht: xw.Sheet = sht # the Sheet this helper working on

    def reset(self, sht=None):
        if sht != None:
            self._sht = sht

    def write_dict(self, dict_:dict, ul=(1, 1)) -> Tuple[int, int]:
        (u, l) = ul
        df:pd.DataFrame = pd.DataFrame.from_dict(dict_, orient='index') # (dict_, index=[0])
        rng:xw.Range = self._sht.range(u, l)

        rng.options(pd.DataFrame, header=False).value = df
        b = u + len(dict_)

        width = 0

        for val in dict_.values():
            if isinstance(val, (list, tuple, np.ndarray)) and not isinstance(val, (str)):
                num_elem = len(val)
            else:
                num_elem = 1

            if num_elem > width :
                width = num_elem

        r = l + width

        return (b, r)

    def read_dict(self, ul, inlcudehead=True) -> Dict:
        (u, l) = ul
        if inlcudehead:
            u += 1 # jump over the head line

        rng:xw.Range = self._to_range((u, l))

        dict_ = {}
        b = rng.end('down').row
        r = rng.end('right').column

        # value range for config, omit the title
        rng = self._sht.range((u, l), (b, r))

        # for VBA-based xlwings, index starts from 1

        for i in range(u, b + 1):
            key = str(self._sht.cells(i, l).value)
            values = []
            for j in range(l + 1, r + 1):
                cur_val = self._sht.cells(i,j).value
                values.append(cur_val)
            if len(values) == 1:
                dict_[key] = values[0]
            else:
                dict_[key] = values

        return dict_

    def write_table(self, table:dict, up=1, left=1, title:dict=None, linespacing=False) -> Tuple[int, int]:
        """
        """
        from xlwings import constants 
        ul = (up, left)
        (b, r) = (up, left + 1)

        if title != None:            
            (b, r) = self.write_dict(dict_=title, ul=ul)
            rng = self._to_range(ul, (b, r))
            self._draw_border(rng)

        (b, r) = self.write_dict(dict_=table, ul=(b,left))
        rng = self._to_range((up, left),(b - 1, r))
        self._draw_border(rng)

        # adjust the key column
        rng.columns[0].autofit()

        # merge for the table title
        rng = self._to_range((up, left + 1), (up, r))
        rng.merge()
        rng.api.HorizontalAlignment = constants.HAlign.xlHAlignCenter        
        rng.color = (200, 200, 200)

        if linespacing:
            b += 1

        return (b, r)

    def _draw_border(self, rng:xw.Range):
        # bid is the enumeration value for XLBordersIndex
        # xlEdgeBottom, xlEdgeLeft, xlEdgeRight, xlEdgeTop, xlInsideHorizontal, xlInsideVertical
        border_set = {9, 7, 10, 8, 12, 11}
        
        for bid in border_set: 
            rng.api.Borders(bid).weight = 2

    def _to_range(self, ul, br=None):
        if br == None:
            return self._sht.range(ul)

        return self._sht.range(ul, br)

def from_degC(T_c) -> float :
    return T_c + 273.15

def from_bar(p) -> float:
    return p * 1e5

def mkdir_filepath(file_name:str):
    '''
    make containing dir for file named file_name
    '''

    from pathlib import Path
    import os.path    

    dirs = file_name.split(os.path.pathsep)
    cur_dir = ""
    for dir_name in dirs[0:-1]:
        cur_dir = os.path.join(cur_dir, dir_name)
        dir = Path(cur_dir)

        if not dir.exists():
            dir.mkdir()

    

