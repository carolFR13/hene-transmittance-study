import numpy as np
import matplotlib.pyplot as plt
from optics import Monochromatic, Angles, Maximums, Transmittance


'''
code to find the distance between prisms with He-Ne data. 

we use experimental graph to find the relative angular distance between 
the two maximums. the distance would be the one that provides the same
relative angular distance in the theorical graph.

we start assuming a particular value of the prism's angle.
'''

# read experimental data

obj = Monochromatic('sources/2023-10-10/set3')
obj._get_transmittance()

T_exp, theta_e = obj.data['b_r_TE']['T']['1']

theta_ext = theta_e * np.pi/180

alpha = 45 * np.pi/180
wavelength = 0.6328 #HeNe wavelength

theta_int = Angles(alpha = alpha, wavelength = wavelength, theta_ext = theta_ext).int_angle() # [rad]
theta_c = Angles(alpha = alpha, wavelength = wavelength, theta_ext = theta_ext).int_critical_angle(water = True)

print(theta_c * 180/np.pi)


# the best values to delimiter the maximums for this alpha are: 

values = [53.8, 56.2, 58.1, 59.6, 60.6]


maximums = []

for i in range(len(values)):
    try:
        _ , angle_1  = Maximums( T = T_exp).interpolation(angles = (theta_int*180/np.pi), 
                                                      value_1 = values[i],value_2 = values[i+1])

        if round(angle_1, 3) not in [round(val, 3) for val in maximums]: 
            maximums.append(angle_1)
    except Exception:
        break

print(maximums)

relat_dist_exp = np.abs(np.diff(maximums)) #[°]

plt.figure()
plt.plot((theta_int*180/np.pi), T_exp, '.')
plt.grid(alpha = 0.7)
plt.title(r'Experimental transmittance for $\alpha$ = %2.2f°' % (alpha*180/np.pi) )
plt.xlabel(r'$\theta_{int} (°)$')
plt.ylabel('T')
plt.show()
# print(relat_dist_exp)

'''
we compute the distance by choosing the value that minimizes the 
difference between the relative angular distance for the experimental data 
and the one we obtain from theoretical data.

the relative distance is computed with the values 1,2,3 found for the experimental
graph since we want both graphs to look alike (?)
'''

d = np.linspace(1,3.5,100) * 10**(-6)  #we expect the distance to be near 3 um

min1 = 1000 ; min2 = 1000
possible_d_1 = [] ; possible_d_2 = []
min_vector_1 = [] ; min_vector_2 = []

for dist in d:

    T_teo = Transmittance(dist, alpha, theta_ext, water = True).transmittance()

    # plt.figure()
    # plt.plot((theta_int*180/np.pi), T_teo, '.')
    # plt.grid(alpha = 0.7)
    # plt.title(r'Experimental transmittance for $\alpha$ = %2.2f°' % (alpha*180/np.pi) )
    # plt.xlabel(r'$\theta_{int} (°)$')
    # plt.ylabel('T')
    # plt.show()

    maximums = []
    for i in range(len(values)):
        try:
            _ , angle  = Maximums( T = T_teo).parabolic(angles = (theta_int*180/np.pi),
                                                        value_1 = values[i], value_2 = values[i+1])

            if round(angle, 3) not in [round(val, 3) for val in maximums]: 
                maximums.append(angle)
            
        except Exception:
            break
    
    # print('len max', len(maximums))

    relat_dist_teo = np.abs(np.diff(maximums)) #[°]

    min_vector_1.append(abs(relat_dist_exp[0] - relat_dist_teo[0]))
    min_vector_2.append(abs(relat_dist_exp[1] - relat_dist_teo[1]))

    if abs(relat_dist_exp[0] - relat_dist_teo[0]) < min1:
        min1 = abs(relat_dist_exp[0] - relat_dist_teo[0])
        final_d_1 = dist
    if abs(relat_dist_exp[1] - relat_dist_teo[1]) < min2:
        min2 = abs(relat_dist_exp[1] - relat_dist_teo[1])
        final_d_2 = dist

    else:
        continue

print('Final distances:', final_d_1 , final_d_2)

plt.figure()
plt.plot(d, min_vector_1,'.', label='1st maximum')
plt.plot(d, min_vector_2, '.', label = '2nd maximum')
plt.grid(alpha = 0.7)
plt.title(r'Obtention of distance between prisms' )
plt.xlabel(r' d ($\mu m$)')
plt.ylabel(r'|$\Delta d_{exp}$-$\Delta d_{teo}$|')
plt.legend(loc = 'best')
plt.show()



T_teo = Transmittance(final_d_1, alpha, theta_ext, water = True).transmittance()

plt.figure()
plt.plot((theta_int*180/np.pi), T_teo, '.', label='theoretical')
plt.plot((theta_int*180/np.pi), T_exp, '.', label = 'experimental')
plt.grid(alpha = 0.7)
plt.title(r'Transmittance for $\alpha$ = %2.2f° and d = %2.2f $\mu$m' % (alpha*180/np.pi,final_d_1*1e06) )
plt.xlabel(r'$\theta_{int} (°)$')
plt.legend(loc='best')
plt.ylabel('T')
plt.show()



# AIR MEAUSUREMENTS: sources/2023-09-29/set2

# for alpha = 45 : d = 1.5555555555555556e-06 /3
# for alpha = 44.833333333333336 : d = 1.5050505050505051e-06 / 3
# for alpha = 44.59090909090909 : d = 1.5050505050505051e-06 / 1.5555555555555556e-06

# WATER MEASUREMENTS: sources/2023-10-10/set3

# for alpha = 44.59090909090909 : d = 3.3484848484848483e-06 / 3.3484848484848483e-06
# for alpha = 45 : d = 3.3484848484848483e-06 3.3484848484848483e-06