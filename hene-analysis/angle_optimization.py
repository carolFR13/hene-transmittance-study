import numpy as np
import matplotlib.pyplot as plt
from optics import Monochromatic, Angles, Maximums, Transmittance, Index


'''
code to find the best value for the prim's angle from He-Ne data.

we will proceed as follows: first we find the value of alpha that
adjusts the experimental angle associated to the first maximum of T to suit the 
theoretical angle for the same maximum. Next we will do the same for the second one. 
The optimal value of alpha we will choose it to be the mean of both.
'''

# read experimental data
# air 
obj = Monochromatic('user_codes/sources/2023-09-29/set1') #a_a
obj2 = Monochromatic('user_codes/sources/2023-09-29/set2') #a_a

# water

obj3 = Monochromatic('user_codes/sources/2023-10-10/set3') #b_r
obj4 = Monochromatic('user_codes/sources/2023-10-10/set4') #b_r

obj4._get_transmittance()

T_exp, theta_e = obj4.data['b_r_TE']['T']['1']

plt.figure()
plt.plot((theta_e), T_exp, '.')
plt.grid(alpha = 0.7)
plt.title(r'Experimental transmittance')
plt.xlabel(r'$\theta_{exp} (°)$')
plt.ylabel('T')
plt.show()

start_value = 380.7

start = Index(theta_e).get_index(start_value)

T_exp = T_exp[:start]
theta_e = theta_e[:start]

theta_ext = theta_e * np.pi/180

plt.figure()
plt.plot((theta_e), T_exp, '.')
plt.grid(alpha = 0.7)
plt.title(r'Experimental transmittance')
plt.xlabel(r'$\theta_{exp} (°)$')
plt.ylabel('T')
plt.show()


'''
the values to delimiter tha maximums will vary significantly with the prism's angle

the way to obtain the values of the internal angle for the maximums will be to 
search for the maximum values of T in each graph. Since those values won't vary
with the prism's angle, we will find the angle associated to those values by 
finding the index of the transmittance array.

interpolation seems to work better when you need maximums values of T. parabolic 
adjustement doesn't return the right value.
'''

#finding theoretical and experimental values of the transmittance maximums

wavelength = 0.6328
distance =  3.424242424242424e-06
alpha = 45 * np.pi/180 #arbitrary prism's angle to find maximums
water = True


theta_int = Angles(alpha = alpha, wavelength = wavelength, theta_ext = theta_ext).int_angle() # [rad]


# experimental maximums 

plt.figure()
plt.plot((theta_int * 180/np.pi) ,T_exp,'.')
plt.grid(alpha=0.7)
plt.title(r'Experimental transmittance for $\alpha$ = %2.2f°' % (alpha*180/np.pi))
plt.xlabel(r'$\theta_{int} (°)$')
plt.ylabel('T')
plt.show()

values = [60.1,61,61.6]

maximums = []
max_indexes = []
print(len(values))
for i in range(len(values)):
    try:
        obj_max_exp = Maximums( T = T_exp)
        T_max_exp , _  = obj_max_exp.interpolation(angles = (theta_int*180/np.pi), 
                                                    value_1 = values[i], value_2 = values[i+1])
        index = obj_max_exp.T_max_index
        if T_max_exp not in [val for val in maximums]: 
            maximums.append(T_max_exp)
        if index not in [val for val in max_indexes]:
            max_indexes.append(index)
    
    except Exception:
        break

print('Maximum values of T for experimental graph:', maximums)
print('Maximum indeces:', max_indexes) #index asociated to that maximums in order to find the maximum located near that area

# theoretical maximums

T_teo = Transmittance(distance, alpha, theta_ext, water = water).transmittance()

plt.figure()
plt.plot((theta_int*180/np.pi), T_teo, '.')
plt.grid(alpha = 0.7)
plt.title(r'Theoretical transmittance for $\alpha$ = %2.2f°' % (alpha*180/np.pi) )
plt.xlabel(r'$\theta_{int} (°)$')
plt.ylabel('T')
plt.show()


maximums_teo = []
max_indexes_teo = []
print(len(values))
for i in range(len(values)):
    try:
        obj_max_teo = Maximums( T = T_teo)
        T_max_teo , _  = obj_max_teo.interpolation(angles = (theta_int*180/np.pi), 
                                                    value_1 = values[i], value_2 = values[i+1])
        index_teo = obj_max_teo.T_max_index
        if T_max_teo not in [val for val in maximums_teo]: 
            maximums_teo.append(T_max_teo)
        if index_teo not in [val for val in max_indexes_teo]:
            max_indexes_teo.append(index_teo)
    
    except Exception:
        break


print('Maximum values of T for theoretical graph:', maximums_teo)
print('Maximum indeces:', max_indexes_teo) 


'''
the following loop finds the values of alpha that minimizes the distance between 
the theoretical angle where the maximum is and the experimental angle for the 
same maximum
'''

a = np.linspace(44.3, 45.8, 100) * np.pi/180

min1 = 1000
min2 = 1000
# min3 = 1000
# min4 = 1000

for alpha in a:

    print('alpha:', alpha*180/np.pi)
    theta_int = Angles(alpha = alpha, wavelength = wavelength, 
                   theta_ext = theta_ext).int_angle() # [rad]
    i_1, i_2  = max_indexes
    i_1_exp = Index(T_exp).get_delimited_index(maximums[0], i_1) #the inferior maximum is bivaluated
    i_2_exp =  Index(T_exp).get_delimited_index(maximums[1], i_2)

    angle_1_exp = (theta_int[i_1_exp]* 180/np.pi) #highest value of theta_int
    angle_2_exp = (theta_int[i_2_exp]* 180/np.pi)

    print('experimental angles',angle_1_exp, angle_2_exp)


    T_teo = Transmittance(distance, alpha, theta_ext, water = water).transmittance()

    ind_1, ind_2 = max_indexes_teo

    i_2_teo = Index(T_teo).get_delimited_index(maximums_teo[0],ind_1)
    i_1_teo = Index(T_teo).get_delimited_index(maximums_teo[1],ind_2)
    
    angle_1_teo = (theta_int[i_2_teo]* 180/np.pi)
    angle_2_teo = (theta_int[i_1_teo]* 180/np.pi)

    print('theoretical angles:',angle_1_teo, angle_2_teo)

    # plt.figure()
    # plt.plot((theta_int * 180/np.pi) ,T_exp,'.',label = 'exp')
    # plt.plot((theta_int*180/np.pi), T_teo, '.')
    # plt.grid(alpha = 0.7)
    # plt.title(r'Theoretical transmittance for $\alpha$ = %2.2f°' % (alpha*180/np.pi) )
    # plt.xlabel(r'$\theta_{int} (°)$')
    # plt.ylabel('T')
    # plt.legend(loc = 'best')
    # plt.show()

    if abs(angle_1_exp-angle_1_teo) < min1: 
        final_a_1 = alpha
        min1 = abs(angle_1_exp-angle_1_teo)
    else:
        pass

    if abs(angle_2_exp-angle_2_teo) < min2:
        final_a_2 = alpha
        min2 = abs(angle_2_exp-angle_2_teo)
    else:
        pass
    
    # if abs(angle_3_exp-angle_3_teo) < min3:
    #     final_a_3 = alpha
    #     min3 = abs(angle_3_exp-angle_3_teo)
    # else:
    #     pass

print(final_a_1*180/np.pi)
print(final_a_2*180/np.pi)


print('Mean value: ', (final_a_1*180/np.pi+final_a_2*180/np.pi)/2)

mean = (final_a_1+final_a_2)/2

T_teo = Transmittance(distance, mean, theta_ext, water = water).transmittance()
theta_int = Angles(alpha = mean, wavelength = wavelength, theta_ext = theta_ext).int_angle() # [rad]

plt.figure()
plt.plot((theta_int*180/np.pi),T_teo,'.', label = 'teo')
plt.plot((theta_int*180/np.pi),T_exp,'.', label = 'exp')
plt.grid(alpha=0.7)
plt.title(r'Theoretical transmittance for $\alpha = %1.3f°$' % (final_a_2*180/np.pi))
plt.xlabel(r'$\theta_{int} (°)$')
plt.ylabel('T')
plt.show()



# AIR MEASUREMENTS: sources/2023-09-29/set1

# for d = 2.8939393939393934e-06 : alpha = 44.93636363636363 / alpha = 44.951515151515146
 


# AIR MEAUSUREMENTS: sources/2023-09-29/set2

# considering 3 peaks:

# for d = 1.5555555555555556e-06 : alpha = 45.333333333333336 / 44.833333333333336 / 44.833333333333336
# for d = 1.5050505050505051e-06 : alpha = 44.59090909090909 / 44.57575757575758 / 44.81818181818182

# considering only 2 peaks:

# for d =  1.5555555555555556e-06 : alpha = 44.865656565656565 / 44.88080808080808



# WATER MEASUREMENTS: sources/2023-10-10/set3

# considering only 2 peaks:

# for d = 3.424242424242424e-06 ; alpha = 44.92878787878787


# WATER MEASUREMENTS: sources/2023-10-10/set4

# considering only 2 peaks:

# for d = 3.424242424242424e-06 ; alpha = 44.88333333333333