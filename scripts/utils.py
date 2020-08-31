'''
    Utility class and method for Tests
'''
from numpy.core.fromnumeric import transpose
from numpy.core.function_base import linspace
from numpy.core.getlimits import _title_fmt
import pandas as pd
import xlwings as xw

class ExcelHelper(object):
    """
    Helper class for reading and writing excels
    """

    def __init__(self, sht:xw.Sheet, col_start = 1):
        super().__init__()
        self._row_idx: int = 1
        self._col_start = col_start
        self._sht: xw.Sheet = sht # the Sheet this helper working on

    def reset(self, sht=None):
        self._row_idx = 1
        if sht != None:
            self._sht = sht

    def _do_write_dict(self, dict_:dict, linespacing=False):

        df:pd.DataFrame = pd.DataFrame.from_dict(dict_, orient='index') # (dict_, index=[0])
        rng:xw.Range = self._sht.range(self._row_idx, self._col_start)

        rng.options(pd.DataFrame, header=False).value = df
        self._row_idx: int = rng.expand().last_cell.row + 1

        if linespacing:
            self._row_idx += 1

    def write_dict(self, dict_:dict, title_dict:dict=None, linespacing=False):
        """
        """
        if title_dict != None:            
            self._do_write_dict(dict_=title_dict)

        self._do_write_dict(dict_=dict_, linespacing=True)

class TestDataIO(object):
    """

    """

    pass