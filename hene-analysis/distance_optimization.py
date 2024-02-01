import numpy as np
import matplotlib.pyplot as plt
from optics import Index, Monochromatic, Angles, Maximums, Transmittance


'''
code to find the distance between prisms with He-Ne data. 

we use  the experimental graph to find the relative angular distance between 
the two maximums. the distance would be the one that provides the same
relative angular distance in the theorical graph.

we start assuming a particular value of the prism's angle.
'''

# read experimental data

# air 
obj = Monochromatic('user_codes/sources/2023-09-29/set1')
obj2 = Monochromatic('user_codes/sources/2023-09-29/set2')

# water

obj3 = Monochromatic('user_codes/sources/2023-10-10/set3')
obj4 = Monochromatic('user_codes/sources/2023-10-10/set4')

water = True

obj4._get_transmittance()

T_exp, theta_e = obj4.data['b_r_TE']['T']['1']


alpha = 44.88333333333333 * np.pi/180
wavelength = 0.6328 #HeNe wavelength


plt.figure()
plt.plot((theta_e), T_exp, '.')
plt.grid(alpha = 0.7)
plt.title(r'Experimental transmittance for $\alpha$ = %2.2f°' % (alpha*180/np.pi) )
plt.xlabel(r'$\theta_{int} (°)$')
plt.ylabel('T')
plt.show()

start_value = 380.8

start = Index(theta_e).get_index(start_value)

T_exp = T_exp[:start]
theta_e = theta_e[:start]

theta_ext = theta_e * np.pi/180

theta_int = Angles(alpha = alpha, wavelength = wavelength, theta_ext = theta_ext).int_angle() # [rad]
theta_c = Angles(alpha = alpha, wavelength = wavelength, theta_ext = theta_ext).int_critical_angle(water=water)

print(theta_c * 180/np.pi)

# we shorten the vectors to get the limited range of interest


plt.figure()
plt.plot((theta_int*180/np.pi), T_exp, '.')
plt.grid(alpha = 0.7)
plt.title(r'Experimental transmittance for $\alpha$ = %2.2f°' % (alpha*180/np.pi) )
plt.xlabel(r'$\theta_{int} (°)$')
plt.ylabel('T')
plt.show()


# the best values to delimiter the maximums for this alpha are: 



values = [39.5,40.75,41.4]
values = [34.9,39.3,41.6]
# values = [53.8, 56.2, 58.1, 59.6, 60.6]
values = [60.1,61,61.6]


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

print(relat_dist_exp)

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
distances_1 = [] ; distances_2 = []

for dist in d:

    T_teo = Transmittance(dist, alpha, theta_ext, water = water).transmittance()

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
            _ , angle  = Maximums( T = T_teo).interpolation(angles = (theta_int*180/np.pi),
                                                        value_1 = values[i], value_2 = values[i+1])

            if round(angle, 3) not in [round(val, 3) for val in maximums]: 
                maximums.append(angle)
            
        except Exception:
            break
    
    # print('len max', len(maximums))

    try:

        relat_dist_teo = np.abs(np.diff(maximums)) #[°]

        min_vector_1.append(abs(relat_dist_exp[0] - relat_dist_teo[0]))
        distances_1.append(dist)

        if abs(relat_dist_exp[0] - relat_dist_teo[0]) < min1:
            min1 = abs(relat_dist_exp[0] - relat_dist_teo[0])
            final_d_1 = dist

        # min_vector_2.append(abs(relat_dist_exp[1] - relat_dist_teo[1]))
        # distances_2.append(dist)
            
        # if abs(relat_dist_exp[1] - relat_dist_teo[1]) < min2:
        #     min2 = abs(relat_dist_exp[1] - relat_dist_teo[1])
        #     final_d_2 = dist

    except Exception:
        continue

print('Final distance 1:', final_d_1 )
# print('Final distance 2:', final_d_2 )


plt.figure()
plt.plot(distances_1, min_vector_1,'.')
# plt.plot(distances_2, min_vector_2,'.')
plt.grid(alpha = 0.7)
plt.title(r'Obtention of distance between prisms' )
plt.xlabel(r' d ($\mu m$)')
plt.ylabel(r'|$\Delta d_{exp}$-$\Delta d_{teo}$|')
plt.legend(loc = 'best')
plt.show()



T_teo = Transmittance(final_d_1, alpha, theta_ext, water = water).transmittance()

plt.figure()
plt.plot((theta_int*180/np.pi), T_teo, label='theoretical')
plt.plot((theta_int*180/np.pi), T_exp, '.', label = 'experimental')
plt.grid(alpha = 0.7)
plt.title(r'Transmittance for $\alpha$ = %2.2f° and d = %2.2f $\mu$m' % (alpha*180/np.pi,final_d_1*1e06) )
plt.xlabel(r'$\theta_{int} (°)$')
plt.legend(loc='best')
plt.ylabel('T')
plt.show()



# AIR MEASUREMENTS: sources/2023-09-29/set1

# for alpha = 45 ; d = 2.8939393939393934e-06
# for alpha = 44.94393939393939 ; d = 2.8939393939393934e-06
 


# AIR MEAUSUREMENTS: sources/2023-09-29/set2

# considering 3 peaks:

# for alpha = 45 : d = 1.5555555555555556e-06 / 3
# for alpha = 44.833333333333336 : d = 1.5050505050505051e-06 / 3
# for alpha = 44.59090909090909 : d = 1.5050505050505051e-06 / 1.5555555555555556e-06

# considering 2 peaks:

# for alpha = 45 : d = 1.5555555555555556e-06
# for alpha = 44.87323232323232 : d = 1.5555555555555556e-06



# WATER MEASUREMENTS: sources/2023-10-10/set3

# using 3 peaks:

# for alpha = 44.59090909090909 : d = 3.3484848484848483e-06 / 3.3484848484848483e-06
# for alpha = 45 : d = 3.3484848484848483e-06 3.3484848484848483e-06

# using the 2 peaks at the right side:

# for alpha = 45 : d = 3.424242424242424e-06
# for alpha = 44.92878787878787 ; d = 3.424242424242424e-06

# WATER MEASUREMENTS: sources/2023-10-10/set4

# using the 2 peaks at the right side:

# for alpha = 45 : d = 3.424242424242424e-06
# for alpha = 44.92878787878787 ; d = 3.424242424242424e-06