from shocktube import *
import json

left = np.linspace(-0.5,0,320)
dxl = left[1] - left[0]
right = np.linspace(0,0.5,40)[1:]
dxr = right[1] - right[0]

left_boundary = np.linspace(-0.5 - (35 * dxl), -0.5 - dxl, 35)
right_boundary = np.linspace(0.5 + dxr, 0.5 + (35 * dxr), 35)

h = 2*(right[1] - right[0])

left = np.append(left_boundary, left)
right = np.append(right, right_boundary)

x = np.append(left, right)

rho = np.append(np.ones_like(left), np.ones_like(right)*0.125)
p = np.append(np.ones_like(left), np.ones_like(right)*0.1)
v = np.zeros_like(x)
gamma = 1.4
epsilon = 0.5
eta = 1e-04
m = 0.0015625000000000

st = Shocktube(x = x, p = p, rho = rho, \
               v = v, m = m, h = h, gamma = gamma, \
               epsilon = epsilon, eta = eta, kernel = cubic_spline)


for i in range(2000):
    st.update_euler(dt=1e-04)
    print(i)

data = dict()
data['rho'] = st.rho
data['p'] = st.p
data['v'] = st.v
data['x'] = st.x
data['e'] = st.e

for key in data.keys():
    data[key] = list(data[key])

with open('shocktube_ce_output.json', 'w') as shock:
    json.dump(data, shock)


### Using summation density 

left = np.linspace(-0.5,0,320)
dxl = left[1] - left[0]
right = np.linspace(0,0.5,40)[1:]
dxr = right[1] - right[0]

left_boundary = np.linspace(-0.5 - (35 * dxl), -0.5 - dxl, 35)
right_boundary = np.linspace(0.5 + dxr, 0.5 + (35 * dxr), 35)

h = 2*(right[1] - right[0])

left = np.append(left_boundary, left)
right = np.append(right, right_boundary)

x = np.append(left, right)

rho = np.append(np.ones_like(left), np.ones_like(right)*0.125)
p = np.append(np.ones_like(left), np.ones_like(right)*0.1)
v = np.zeros_like(x)
gamma = 1.4
epsilon = 0.5
eta = 1e-04
m = 0.0015625000000000

st = Shocktube(x = x, p = p, rho = rho, \
               v = v, m = m, h = h, gamma = gamma, \
               epsilon = epsilon, eta = eta, kernel = cubic_spline)

    
for i in range(2000):
    st.update_euler_SD(dt=1e-04)   
    print(i)


data = dict()
data['rho'] = st.rho
data['p'] = st.p
data['v'] = st.v
data['x'] = st.x
data['e'] = st.e

for key in data.keys():
    data[key] = list(data[key])

with open('shocktube_sd_output.json', 'w') as shock:
    json.dump(data, shock)
