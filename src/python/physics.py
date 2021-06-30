'''
definition of basic Quantity, relations in physics
'''

import json
import jsonpickle

class Unit(object):
    '''
    class for unit conversion
    # IMPORTANT: SIunits will be used in this script : [T] = K, [p] = pa, [length] = m, [angle] = rad, [mass] = kg, [time] = s   
    '''
    from pint import UnitRegistry
    ureg = UnitRegistry()

    @classmethod
    def convert(cls, val, src_unit, dest_unit):
        ureg = Unit.ureg
        ureg.Unit
        src = val * ureg[src_unit]
        return src.to(ureg[dest_unit]).magnitude

class Quantity(object):
    '''
        pyhsics quantity : magnitude [unit], such as 1 m, 10 s or 1e6 pa 
    '''
    from pint import UnitRegistry
    ureg = UnitRegistry()

    # Non-Dimensional Unit for Integer, Text or constant
    ND_UNIT = '1'

    def __init__(self, mag = 0, def_unit = ND_UNIT, unit = None, disp_prec = 2):
        super().__init__()
        
        self._def_unit = def_unit        
        self._unit = unit

        if unit == None:
            self._unit = def_unit

        if not self.unit == self.def_unit:
            self._mag = self.__convert(mag, self.unit, self.def_unit) 
        else:
            self._mag = mag 
        
        self._disp_prec = disp_prec

    def __convert(self, val, src, dst):
        return Quantity.ureg.convert(val, src, dst)

    @property
    def unit(self):
        """The unit property."""
        return self._def_unit
    @unit.setter
    def def_unit(self, value):
        self._def_unit = value

    @property
    def mag(self):
        """The magnitude property."""
        return self._mag
        
    @mag.setter
    def mag(self, value):
        self._mag = value

    @property
    def unit(self):
        """
        The disp_unit property. for disp purpose
        like pa could be displayed in kPa or mPa for short
        """
        return self._unit

    @unit.setter
    def unit(self, value):
        self._unit = value

    @property
    def disp_prec(self):
        """The disp precison of magnitude."""
        return self._disp_prec

    @disp_prec.setter
    def disp_prec(self, value):
        self._disp_prec = value

    def __repr__(self):
        val = self.mag
        if not self.unit == self.def_unit:
            val = Quantity.ureg.convert(self.mag, self.def_unit, self.unit)
        fmt = "<0:1.{0}f> <1>".format(self.disp_prec).replace('<','{').replace('>', '}')

        if self.unit == Quantity.ND_UNIT:
            return fmt.format(val)

        return fmt.format(val, self.unit)

class Length(Quantity):
    '''
        Quantity length in meter
    '''
    def __init__(self, arg, unit=None):
        super().__init__(arg, 'm', unit)

    @classmethod
    def mm(cls, val):
        return Length(val, 'mm')
        
    @classmethod
    def km(cls, val):
        return Length(val, 'km')

class Mass(Quantity):
    '''
        Quantity length in kg
    '''
    def __init__(self, arg, unit=None):
        super().__init__(arg, 'kg', unit)

class Angle(Quantity):
    '''
    Quantity - Density in rad
    '''

    def __init__(self, arg, unit=None):
        super().__init__(arg,'rad', unit=unit)

    @classmethod    
    def deg(cls, val):
        return Angle(val, '°')
class Pressure(Quantity):

    '''
    Quantity Pressure
    '''
    def __init__(self, val, unit=None):
        super().__init__(val, 'Pa', unit)

    @classmethod
    def bar(cls, val):
        return Pressure(val, 'bar')

    @classmethod
    def MPa(cls, val):
        return Pressure(val, 'MPa')

class Temperature(Quantity):
    '''
    Quantity - Temperature in K
    '''

    def __init__(self, val, unit=None):
        super().__init__(val, 'K', unit)
    
    @classmethod
    def degC(cls, val):
        return Temperature(val, 'degC')

class MDot(Quantity):
    '''
    Quantity - mass flow rate in kg/s
    '''

    def __init__(self, val, unit=None):
        super().__init__(val, 'kg/s', unit)

class Velocity(Quantity):
    '''
    Quantity - Velocity in m/s
    '''

    def __init__(self, val, unit=None):
        super().__init__(val, 'm/s', unit=unit)

class Density(Quantity):
    '''
    Quantity - Density in kg/m^3
    '''

    def __init__(self, val, unit=None):
        super().__init__(val, 'kg/m**3', unit=unit)

def main():

    p = Pressure.bar(90)

    print(p)

    p = Pressure(80)

    print(p)

    p = Pressure(90, 'bar')

    print(p)

    p = Pressure(90, 'MPa')

    p.unit = 'Pa'

    print(p)

    T = Temperature(500)

    print(T)

    T = Temperature.degC(30)

    print(T)

    T = Temperature(50, 'degC')
    T.unit = 'K'
    
    print(T)

    mdot = MDot(10)

    print(mdot)

    rho = Density(12.3)

    print(rho)

    rho.unit = 'kg/cm^3'
    rho.disp_prec = 10

    print(rho)

    a = Angle(45)

    print(a)

    a = Angle.deg(30)

    print(a)

    a.unit = 'rad'

    print(a)

    val = Length(2, 'mm')

    print(val)   

    val.unit = 'm'
    val.disp_prec = 5
    print(val)

    empJSON = jsonpickle.encode(val, unpicklable=False)

    data = json.dumps(empJSON, indent=4)

    print(data)

if __name__ == "__main__":
    main()