import numpy as np
import os
from scipy.optimize import fsolve, curve_fit
from scipy.interpolate import interp1d
from .utils import n_air, n_glass, n_water, parabola, Index
from scipy import constants as sc
from uncertainties import umath


class Angles:
    '''
    Angles class to obtain internal angles from the external ones 
    and viceversa.

    :param alpha: prism's angle [rad]
    :param waveleght: light's wavelenght [um]
    :param theta_ext: external angle (optional) [rad]
    :param theta_int: internal angle (optional) [rad]
    
    '''

    def __init__(self, alpha: float, wavelength: float | list, theta_ext: float| list | None = None, 
                 theta_int: float | list | None = None) -> None:
        self.alpha = alpha 
        self.wavelength = wavelength
        self.theta_ext = theta_ext
        self.theta_int = theta_int

        
        self.n_air = n_air(self.wavelength)
        self.n_glass = n_glass(self.wavelength)
        self.n_water = n_water(self.wavelength)

        return None 

    def int_angle(self, uncertainty: bool = False):
        '''
        method to compute the internal angle from the external one
        (the variable theta_ext needs to be defined)

        it allows you to compute internal angle with its uncertainty given an external 
        angle in ufloat format
        '''
        try:
            if uncertainty:
                theta_1 = umath.asin((self.n_air *umath.sin(self.theta_ext))/self.n_glass)
                theta_2 = self.alpha + theta_1
                return theta_2

            else:
                theta_1 = np.arcsin((self.n_air *np.sin(self.theta_ext))/self.n_glass)
                theta_2 = self.alpha + theta_1
                return theta_2
        
        except Exception:
            print('The variable theta_ext is not defined.')
            return None
    
    def int_critical_angle(self, water : bool = False):
        '''
        method to compute the internal critical angle
        '''

        if water: 
            theta_c = np.arcsin(self.n_water/self.n_glass)
        else:
            theta_c = np.arcsin(self.n_air/self.n_glass)

        return theta_c

    def ext_angle(self):
        '''
        method to compute the external angle from the internal one
        (the variable theta_int needs to be defined)
        '''

        try:
            theta_1 = self.theta_int - self.alpha
            theta_ext = np.arcsin((self.n_glass * np.sin(theta_1))/self.n_air)
        except Exception: 
            print('The variable theta_int is not defined.')

        return theta_ext
    
    def critical_lambda(self):
        '''
        if correct method to obtain the critical wavelength
        '''

        def f(wavelength):

            theta_int =  Angles(self.alpha, wavelength, theta_ext = self.theta_ext).int_angle()
            obj = Angles(self.alpha, wavelength, self.theta_ext, theta_int)
            
            theta_2 = obj.int_angle()
            theta_c = obj.int_critical_angle()

            if theta_2 is None: 
                return float('inf')
            else:
                return theta_2 - theta_c

        lambda_c = fsolve(f,x0 = 0.3)

        return lambda_c[0]


class Transmittance: 
    '''
    class to compute the theoretical transmittance and reflection.

    :param distance: distance between prisms [m]
    :param alpha: prism's angle [rad]
    :param theta_ext: external angle [rad]
    '''

    def __init__(self, distance, alpha, theta_ext, water : bool = False) -> None:

        self.distance = distance
        self.alpha = alpha
        self.theta_ext = theta_ext

        self.wavelength = 0.6328

        self.theta_int = Angles(self.alpha, self.wavelength, self.theta_ext).int_angle()
        self.theta_c = Angles(self.alpha, self.wavelength, self.theta_ext).int_critical_angle(water = water)

        
        self.n_glass = n_glass(self.wavelength)

        if water:
            self.n_air = n_water(self.wavelength) #if the medium in the gap is water we use self.n_air as the label to n_water
        else: 
            self.n_air = n_air(self.wavelength)


        return None
    
    def _fabry_perot(self, _distance, _theta_int, _n_air, _n_glass):
        '''
        internal method to compute transmission outside FTIR range
        '''
        k0 = (2*np.pi)/(self.wavelength *10**(-6))
        beta = _n_glass * k0 * np.sin(_theta_int)

        kp = np.sqrt((_n_glass**2*k0**2)-beta**2)
        kc = np.sqrt((_n_air**2*k0**2)-beta**2)


        T = (1+(1/4)*(((kp**2-kc**2)/(kp*kc))**2)*np.sin(kc*_distance)**2)**(-1)
        R = 1 - T

        return T, R
    
    def _ftir(self, _distance, _theta_int, _n_air, _n_glass):
        '''
        internal method to compute transmission in FTIR range
        '''
        k0 = (2*np.pi)/(self.wavelength *10**(-6))
        beta = _n_glass * k0 * np.sin(_theta_int)

        kp = np.sqrt((_n_glass**2*k0**2)-beta**2)
        kappa = np.sqrt(beta**2-(_n_air**2*k0**2))


        T = (1+((kappa**2+kp**2)**2/(4*(kp*kappa)**2))*np.sinh(kappa*_distance)**2)**(-1)
        R = 1 - T

        return T, R

    def transmittance(self):
        '''
        method to compute the total transmittance outside and inside FTIR range
        '''

        T = []
        for i in range(len(self.theta_int)):
            if self.theta_int[i] < self.theta_c:
                T.append(self._fabry_perot(self.distance, self.theta_int[i],
                                           self.n_air,self.n_glass)[0])
            else:
                T.append(self._ftir(self.distance, self.theta_int[i],
                                    self.n_air,self.n_glass)[0])
                
        return T

    def reflectance(self):
        '''
        method to compute the total reflectance outside and inside FTIR range
        '''

        R = []
        for i in range(len(self.theta_int)):
            if self.theta_int[i] < self.theta_c:
                R.append(self._fabry_perot(self.distance, self.theta_int[i],
                                           self.n_air,self.n_glass)[1])
            else:
                R.append(self._ftir(self.distance, self.theta_int[i],
                                    self.n_air,self.n_glass)[1])
                
        return R


class Maximums: 

    '''
    Method to compute the maximas of T and R and the corresponding angles/wavelength 

    :param T: transmittance
    :param R: reflectance
    '''
    def __init__(self, T : None | list = None, R: None | list = None ) -> None:
        
        self.T = T
        self.R = R

        self.angle : float | None = None  # value of theta associated to T/R maximum
        self.angle_R : float | None = None

        self.T_max : float | None = None #value of T / R associated to the maximum depending on whether T or R is provided
        self.R_max : float | None = None

        self.T_max_index : int | None = None
        self.R_max_index : int | None = None

        self.wavelength = 0.6328

  
    def interpolation(self, angles : np.array, value_1 : float, value_2 : float) -> [float, float]:
        '''
        Method to obtain maximum value of T,R and the angle associated to that 
        maximum by doing a quadratic interpolation.

        value_1 and value_2 are the angle values that delimiter the maximum
        we are studying.  

        theta is a continuosly increasing function so we can use get_index()
        '''

        i_1 = Index(angles).get_index(value_1)
        i_2 = Index(angles).get_index(value_2)

        sorted_indices = np.sort([i_1, i_2])
        i_1, i_2  = sorted_indices


        vect=[]
        if self.T is not None:
            vect = [self.T]
        elif self.R is not None:
            vect.append(self.R)
        else:
            'Neither T nor R are defined.'

        for vector in vect:

            vector_1 = np.array(vector[i_1:i_2])

            if 'R' in vector:
                index_1 = np.where(vector == vector_1.min())[0][0]
                self.R_max_index = index_1
            else:
                index_1 = np.where(vector == vector_1.max())[0][0]
                self.T_max_index = index_1


            v_1 = vector[index_1-10:index_1+10] ; angles_1 = angles[index_1-10:index_1+10]

            interp_parabolic_1 = interp1d(angles_1, v_1, kind='quadratic')

            # maximum: 
            x_int_1 = np.linspace(angles[index_1-10],angles[index_1+9], 10000)
            v_int_1 = interp_parabolic_1(x_int_1)


            #store data:

            if 'R' in vector:
                self.R_max = v_int_1.min()
                self.angle_R = x_int_1[np.where(v_int_1 ==  v_int_1.min())[0][0]] 

                return self.R_max, self.angle_R
            
            else:
                self.T_max = v_int_1.max()
                self.angle_T = x_int_1[np.where(v_int_1 ==  v_int_1.max())[0][0]] 

                return self.T_max, self.angle_T


    def parabolic(self, angles : np.array , value_1 : float , value_2 : float) -> [float, float]:
        
        '''
        Method to obtain maximum value of T,R and the angle associated to that 
        maximum by doing a parabolic adjustment.

        value_1 and value_2 are the angle values that delimiter the maximum
        we are studying.  

        theta is a continuosly increasing function so we can use get_index()
        '''

        i_1 = Index(angles).get_index(value_1)
        i_2 = Index(angles).get_index(value_2)
        

        sorted_indices = np.sort([i_1, i_2])

        # Assign the sorted indices
        i_1, i_2  = sorted_indices

        vect=[]
        if self.T is not None:
            vect = [self.T]
        elif self.R is not None:
            vect.append(self.R)
        else:
            'Neither T nor R are defined.'


        for vector in vect:
            vector_1 = np.array(vector[i_1:i_2])

            if 'R' in vector:
                index_1 = np.where(vector == vector_1.min())[0][0]
            else:
                index_1 = np.where(vector == vector_1.max())[0][0]

            v_1 = vector[index_1-10:index_1+10] ; angles_1 = angles[index_1-10:index_1+10]

            params1, _ = curve_fit(parabola, angles_1, v_1)

            # maximum: 
            x_int_1 = np.linspace(angles[index_1-10],angles[index_1+9], 10000)
            v_int_1 = parabola(x_int_1, *params1)


            #store data:
            if 'R' in vector:
                self.R_max = v_int_1.min()
                self.angle_R = x_int_1[np.where(v_int_1 ==  v_int_1.min())[0][0]] 

                return self.R_max, self.angle_R
            
            else:
                self.T_max = v_int_1.max()
                self.angle_T = x_int_1[np.where(v_int_1 ==  v_int_1.max())[0][0]] 

                return self.T_max, self.angle_T
        
    def distance(self, theta_1: float , theta_2: float,  water: bool = False ) -> float:
        
        '''
        method to compute the distance. 
        
        :water paraeter: if True the gap is filled with water
        :theta_1 parameter: angle corresponding to the 1st maximum
        :theta_2 parameter: angle corresponding to the 2nd maximum
        '''


        if water:
            n_air_ = n_water(self.wavelength)
        else:
            n_air_ = n_air(self.wavelength)
        
        n_glass_ = n_glass(self.wavelength)

        d = (1/2) * ( (( n_air_**2 - (n_glass_*np.sin(theta_2))**2 )**(1/2))/self.wavelength - (( n_air_**2 - (n_glass_*np.sin(theta_1))**2 )**(1/2))/self.wavelength  )**(-1) 
            
        return abs(d) 

    
