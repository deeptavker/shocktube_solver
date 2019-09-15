from matplotlib import pyplot as plt
import json 
import numpy as np

exact_v_x = [-0.5,-0.2,0,0.4,0.4,0.5]
exact_v_y = [0,0,1,1,0,0]

exact_e_x = [-0.5,-0.2,0,0.2, 0.2,0.4,0.4,0.5]
exact_e_y = [2.5,2.5,1.75,1.75,2.75,2.75,2,2]

exact_rho_x = [-0.5,-0.2,0,0.2,0.2,0.4,0.4,0.5]
exact_rho_y = [1,1,0.4,0.4,0.22,0.22, 0.175, 0.175]

exact_p_x = [-0.5,-0.2,0,0.4,0.4,0.5]
exact_p_y = [1, 1, 0.3, 0.3, 0.1, 0.1]

with open('shocktube_ce_output.json', 'r') as shock2:
    data = json.load(shock2)
    
for key in data.keys():
    data[key] = np.array(data[key])
    
fig, ax = plt.subplots(2,2, figsize=(10,10))

mask = np.where(abs(np.array(data['x'])) <= 0.5)
mask = mask[0]

ax[0,0].plot(exact_e_x, exact_e_y)
ax[0,0].scatter(data['x'][mask], data['e'][mask], s=5)
ax[0,0].set_title('Parameter : e')

ax[0,1].plot(exact_p_x, exact_p_y)
ax[0,1].scatter(data['x'][mask], data['p'][mask], s=5)
ax[0,1].set_title('Parameter : p')

ax[1,0].plot(exact_rho_x, exact_rho_y)
ax[1,0].scatter(data['x'][mask], data['rho'][mask], s=5)
ax[1,0].set_title('Parameter : rho')

ax[1,1].plot(exact_v_x, exact_v_y)
ax[1,1].scatter(data['x'][mask], data['v'][mask], s=5)
ax[1,1].set_title('Parameter : v')

plt.savefig('figure_ce_output.png')


with open('shocktube_sd_output.json', 'r') as shock2:
    data = json.load(shock2)
    
for key in data.keys():
    data[key] = np.array(data[key])
    
fig, ax = plt.subplots(2,2, figsize=(10,10))

mask = np.where(abs(np.array(data['x'])) <= 0.5)
mask = mask[0]

ax[0,0].plot(exact_e_x, exact_e_y)
ax[0,0].scatter(data['x'][mask], data['e'][mask], s=5)
ax[0,0].set_title('Parameter : e')

ax[0,1].plot(exact_p_x, exact_p_y)
ax[0,1].scatter(data['x'][mask], data['p'][mask], s=5)
ax[0,1].set_title('Parameter : p')

ax[1,0].plot(exact_rho_x, exact_rho_y)
ax[1,0].scatter(data['x'][mask], data['rho'][mask], s=5)
ax[1,0].set_title('Parameter : rho')

ax[1,1].plot(exact_v_x, exact_v_y)
ax[1,1].scatter(data['x'][mask], data['v'][mask], s=5)
ax[1,1].set_title('Parameter : v')

plt.savefig('figure_sd_output.png')