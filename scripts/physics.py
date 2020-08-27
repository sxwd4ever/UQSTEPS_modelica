'''
definition of basic Quantity, relations in physics
'''
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

    @classmethod
    def from_bar(cls, val):
        return Unit.convert(val, 'bar', 'Pa')

    @classmethod
    def from_degC(cls, val):
        return Unit.convert(val, '°C','K')

    @classmethod
    def from_deg(cls, val): 
        '''
        deg to rad for angle
        '''
        return Unit.convert(val, '°', 'rad')

class Quantity(object):
    '''
        pyhsics quantity : magnitude + unit, such 1 m, 10 s or 1e6 pa 
    '''
    from pint import UnitRegistry
    ureg = UnitRegistry()

    def __init__(self, mag = 0, unit = 1, disp_unit = [], disp_prec = 2):
        super().__init__()
        
        self._unit = unit        
        self._disp_unit = disp_unit

        if disp_unit == []:
            self._disp_unit = unit

        if not self.disp_unit == self.unit:
            self._magnitude = self.__convert(mag, self.disp_unit, self.unit) 
        else:
            self._magnitude = mag 
        
        self._disp_prec = disp_prec

    def __convert(self, val, src, dst):
        return Quantity.ureg.convert(val, src, dst)

    @property
    def unit(self):
        """The unit property."""
        return self._unit
    @unit.setter
    def unit(self, value):
        self._unit = value

    @property
    def magnitude(self):
        """The magnitude property."""
        return self._magnitude
    @magnitude.setter
    def magnitude(self, value):
        self._magnitude = value

    @property
    def disp_unit(self):
        """
        The disp_unit property. for disp purpose
        like pa could be displayed in kPa or mPa for short
        """
        return self._disp_unit

    @disp_unit.setter
    def disp_unit(self, value):
        self._disp_unit = value

    @property
    def disp_prec(self):
        """The disp precison of magnitude."""
        return self._disp_prec

    @disp_prec.setter
    def disp_prec(self, value):
        self._disp_prec = value

    def __repr__(self):
        val = self.magnitude
        if not self.unit == self.disp_unit:
            val = Quantity.ureg.convert(self._magnitude, self.unit, self.disp_unit)
        fmt = "<0:1.{0}f> <1>".format(self.disp_prec).replace('<','{').replace('>', '}')

        return fmt.format(val, self.disp_unit)

class Pressure(Quantity):
    '''
    Quantity Pressure
    '''
    def __init__(self, val, disp_unit=[]):
        super().__init__(val, 'Pa', disp_unit)

    @classmethod
    def from_bar(cls, val):
        return Pressure(val, 'bar')

class Temperature(Quantity):
    '''
    Quantity - Temperature in K
    '''

    def __init__(self, val, disp_unit=[]):
        super().__init__(val, 'K', disp_unit)
    
    @classmethod
    def from_degC(cls, val):
        return Temperature(val, 'degC')

class MDot(Quantity):
    '''
    Quantity - mass flow rate in kg/s
    '''

    def __init__(self, val, disp_unit=[]):
        super().__init__(val, 'kg/s', disp_unit)

class Velocity(Quantity):
    '''
    Quantity - Velocity in m/s
    '''

    def __init__(self, val, disp_unit=[]):
        super().__init__(val, 'm/s', disp_unit=disp_unit)

class Density(Quantity):
    """docstring for Density in kg/s."""

    def __init__(self, val, disp_unit=[]):
        super().__init__(val, 'kg/m**3', disp_unit=disp_unit)

class Angle(Quantity):
    """docstring for Degree."""
    def __init__(self, arg, disp_unit=[]):
        super().__init__(arg,'rad', disp_unit=disp_unit)
        self.arg = arg

    @classmethod    
    def from_deg(cls, val):
        return Angle(val, '°')

def main():

    p = Pressure.from_bar(90)

    print(p)

    p = Pressure(80)

    print(p)

    p = Pressure(90, 'bar')

    print(p)

    p = Pressure(90, 'MPa')

    p.disp_unit = 'Pa'

    print(p)

    T = Temperature(500)

    print(T)

    T = Temperature.from_degC(30)

    print(T)

    T = Temperature(50, 'degC')
    T.disp_unit = 'K'
    
    print(T)

    mdot = MDot(10)

    print(mdot)

    rho = Density(12.3)

    print(rho)

    rho.disp_unit = 'kg/cm**3'
    rho.disp_prec = 10

    print(rho)

    a = Angle(45)

    print(a)

    a = Angle.from_deg(30)

    print(a)

    a.disp_unit = 'rad'

    print(a)

if __name__ == "__main__":
    main()