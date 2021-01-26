from sympy import symbols
from sympy.physics.optics import TWave
from sympy.physics.units import speed_of_light, meter, second 

class HeNe_laser:

    def __init__(self,transition=None,lifetime=None,power=None,color=None,Spectral_width=None,amplitude=None,frequency=None,phase=S.Zero,time_period=None):
        self.transition = transition #fixed
        self.lifetime = lifetime #fixed
        self.power = power #fixed
        self.color = color #fixed
        self.Spectral_width = Spectral_width
        self._frequency = frequency
        self._amplitude = amplitude
        self._phase = phase
        self._time_period = time_period
        # if time_period is not None:
        #     self._frequency = 1/self._time_period
        # if frequency is not None:
        #     self._time_period = 1/self._frequency
        #     if time_period is not None:
        #         if frequency != 1/time_period:
        #             raise ValueError("frequency and time_period should be consistent.")
        # if frequency is None and time_period is None:
        #     raise ValueError("Either frequency or time period is needed.")

    class red(HeNe_laser):
        def __init__(self,transition,lifetime,power,GHz,α)
        pass

    class infrared_one(HeNe_laser):
        pass

    class infrared_three(HeNe_laser):
        pass

    @property
    def frequency(self):

        return self._frequency

    @property
    def time_period(self):
        return self._time_period

    @property
    def wavelength(self):
        return c/(self._frequency*self._n)

    @property
    def amplitude(self):
        return self._amplitude

    @property
    def phase(self):
        return 2*pi*self._frequency

    @property
    def wavenumber(self):
        return 2*pi/self.wavelength

    def __str__(self):
        """String representation of a TWave."""
        from sympy.printing import sstr
        return type(self).__name__ + sstr(self.args)

    __repr__ = __str__

    def __add__(self, other):
        """
        Addition of two waves will result in their superposition.
        The type of interference will depend on their phase angles.
        """
        if isinstance(other, TWave):
            if self._frequency == other._frequency and self.wavelength == other.wavelength:
                return TWave(sqrt(self._amplitude**2 + other._amplitude**2 + 2 *
                                  self.amplitude*other.amplitude*cos(
                                      self._phase - other.phase)),
                             self.frequency,
                             atan2(self._amplitude*cos(self._phase)
                             +other._amplitude*cos(other._phase),
                             self._amplitude*sin(self._phase)
                             +other._amplitude*sin(other._phase))
                             )
            else:
                raise NotImplementedError("Interference of waves with different frequencies"
                    " has not been implemented.")
        else:
            raise TypeError(type(other).__name__ + " and TWave objects can't be added.")

    def _eval_rewrite_as_sin(self, *args, **kwargs):
        return self._amplitude*sin(self.wavenumber*Symbol('x')
            - self.angular_velocity*Symbol('t') + self._phase + pi/2, evaluate=False)

    def _eval_rewrite_as_cos(self, *args, **kwargs):
        return self._amplitude*cos(self.wavenumber*Symbol('x')
            - self.angular_velocity*Symbol('t') + self._phase)

    def _eval_rewrite_as_pde(self, *args, **kwargs):
        from sympy import Function
        mu, epsilon, x, t = symbols('mu, epsilon, x, t')
        E = Function('E')
        return Derivative(E(x, t), x, 2) + mu*epsilon*Derivative(E(x, t), t, 2)

    def _eval_rewrite_as_exp(self, *args, **kwargs):
        from sympy import exp, I
        return self._amplitude*exp(I*(self.wavenumber*Symbol('x')
            - self.angular_velocity*Symbol('t') + self._phase))
