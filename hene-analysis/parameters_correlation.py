import numpy as np
import matplotlib.pyplot as plt
from optics import Monochromatic, Angles, Transmittance


# read experimental data

obj = Monochromatic('sources/2023-10-10/set3')
obj._get_transmittance()

T_exp, theta_e = obj.data['b_r_TE']['T']['1']

theta_ext = theta_e * np.pi/180 #[rad]
alpha = 44.87878787878788 * np.pi/180
wavelength = 0.6328 #HeNe wavelength


theta_int = Angles(alpha = alpha, wavelength = wavelength, theta_ext = theta_ext).int_angle() # [rad]

distance = 3.3484848484848483e-06

T_teo = Transmittance(distance, alpha, theta_ext, water = True).transmittance()

plt.figure()
plt.plot((theta_int*180/np.pi), T_teo, '.', label='theoretical')
plt.plot((theta_int*180/np.pi), T_exp, '.', label = 'experimental')
plt.grid(alpha = 0.7)
plt.title(r'Transmittance for $\alpha$ = %2.2f° and d = %2.2f $\mu$m' % (alpha*180/np.pi,distance*1e06) )
plt.xlabel(r'$\theta_{int} (°)$')
plt.legend(loc='best')
plt.ylabel('T')
plt.show()

num_points = 100
distance_values = np.linspace(2e-6, 6e-6, num_points)  # Adjust the range based on your problem
alpha_values = np.linspace(44, 46, num_points) * np.pi/180  # Assuming alpha is between 0 and pi/2

# Initialize variables to store the best parameters and correlation
best_distance = None
best_alpha = None
best_correlation = -1  # Initialize with a low value

# Loop through the parameter space
for dist in distance_values:
    for a in alpha_values:
        # Calculate theoretical transmittance
        T_theoretical = Transmittance(dist, a, theta_ext, water = True).transmittance()

        # Calculate correlation coefficient
        correlation = np.corrcoef(T_theoretical, T_exp)[0, 1]

        # Update best parameters if the correlation is improved
        if correlation > best_correlation:
            best_correlation = correlation
            best_distance = distance
            best_alpha = alpha

print(f"Best Distance: {best_distance}")
print(f"Best Alpha: {best_alpha* 180/np.pi}")

T_teo = Transmittance(best_distance, best_alpha, theta_ext, water = True).transmittance()
theta_int = Angles(alpha = best_alpha, wavelength = wavelength, theta_ext = theta_ext).int_angle() # [rad]


plt.figure()
plt.plot((theta_int*180/np.pi), T_teo, '.', label='theoretical')
plt.plot((theta_int*180/np.pi), T_exp, '.', label = 'experimental')
plt.grid(alpha = 0.7)
plt.title(r'Transmittance for $\alpha$ = %2.2f° and d = %2.2f $\mu$m' % (best_alpha*180/np.pi,best_distance*1e06) )
plt.xlabel(r'$\theta_{int} (°)$')
plt.legend(loc='best')
plt.ylabel('T')
plt.show()