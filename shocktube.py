import numpy as np

@np.vectorize 
def cubic_spline(r, h, order):
    
    q = float(abs(r))/h
    sigma = float(2.0 / (3.0 * h))
    if order == 0:
        if q <=1.0 and q >=0.0: 
            return sigma * (1.0 - (1.5 * q * q) + (0.75 * q * q * q))
        elif q > 1.0 and q <= 2.0:
            return sigma * 0.25 * ((2.0 - q) ** 3.0)
        else:
            return 0.0
    else:
        diff_multiplier = float(np.sign(r) / h)
        if q <=1.0 and q >=0.0: 
            return float(sigma * ((-3.0 * q) + (2.25 * q * q)) * diff_multiplier)
        elif q > 1.0 and q <= 2.0:
            return float(sigma * -0.75 * ((2 - q) ** 2) * diff_multiplier)
        else:
            return 0.0

class Shocktube(object):
    
    ### Need to add solid boundaries ? ###
    
    def __init__(self, x, p, rho, v, m, h, gamma, epsilon, eta, kernel):
        
        self.x = x
        self.p = p
        self.rho = rho
        self.e = p / ((gamma - 1) * rho)
        self.m = m
        self.h = h
        self.gamma = gamma
        self.epsilon = epsilon
        self.eta = eta
        self.kernel = kernel
        self.v = v
        self.num = len(self.x)
        
    #### One Step Euler Integrator #### 

    def update_euler(self, dt):
        
        ## Define temp arrays for storing increments ##
        
        x_inc = np.zeros_like(self.x, dtype='float32')
        rho_inc = np.zeros_like(self.x, dtype='float32')
        v_inc = np.zeros_like(self.x, dtype='float32')
        e_inc = np.zeros_like(self.x, dtype='float32')
        
        ## Iterate over all particles and store the accelerations in temp arrays ##
        ## 10 particles to the left and 10 particles to the right as boundary ###
        for i in range(35, self.num - 35):

            ### Evaluating variables ###

            vij = self.v[i] - self.v
            dwij = self.kernel(self.x[i] - self.x, self.h, order = 1)
            wij = self.kernel(self.x[i] - self.x, self.h, order = 0)
            xij = self.x[i] - self.x

            p_rho_i = self.p[i] / (self.rho[i] * self.rho[i])
            p_rho_j = self.p / (self.rho * self.rho)
            
            p_rho_ij =  p_rho_i + p_rho_j

            ### Artificial viscosity ###

            numerator = self.h * vij * xij
            denominator = (xij ** 2.0) + (self.eta ** 2.0)
            mu_ij = numerator / denominator
            mu_ij[mu_ij > 0] = 0.0 #only activated for approaching particles

            ci = (self.gamma * self.p[i] / self.rho[i]) ** 0.5
            cj = (self.gamma * self.p / self.rho) ** 0.5
            cij = 0.5 * (ci + cj)
            rhoij = 0.5 * (self.rho[i] + self.rho)
            numerator = (-1 * cij * mu_ij) + (mu_ij ** 2)
            denominator = rhoij
            pi_ij = numerator / denominator

            ### Evaluating gradients ###

            grad_rho = self.rho[i] * np.sum(self.m * vij * dwij / self.rho)
            grad_v = -1 * np.sum(self.m * (p_rho_ij + pi_ij) * dwij)            
            grad_e = 0.5 * np.sum(self.m * (p_rho_ij + pi_ij) * vij * dwij)

            rho_inc[i] = dt * grad_rho
            v_inc[i] = dt * grad_v
            e_inc[i] = dt * grad_e
            
            ### Get XSPH Velocity and calculate increment ###

            rho_avg = 0.5 * (self.rho[i] + self.rho)
            correction = self.epsilon * np.sum(self.m * -1 * vij * wij / rho_avg)
            xsph_velocity =  self.v[i] + correction

            x_inc[i] = dt * xsph_velocity

            
        ## Update the original arrays using the increment arrays ##
            
        self.rho += rho_inc
        self.v += v_inc
        self.e += e_inc

        ### Update pressure ###
        
        self.p = (self.gamma - 1) * self.e * self.rho

        ### Update positions ###

        self.x += x_inc

    def update_euler_SD(self, dt):
        
        ## Define temp arrays for storing increments ##
        
        x_inc = np.zeros_like(self.x, dtype='float32')
        rho_inc = np.copy(self.rho) # not increment, direct corection
        v_inc = np.zeros_like(self.x, dtype='float32')
        e_inc = np.zeros_like(self.x, dtype='float32')
        
        ## Iterate over all particles and store the accelerations in temp arrays ##
        ## 10 particles to the left and 10 particles to the right as boundary ###
        for i in range(35, self.num - 35):

            ### Evaluating variables ###

            vij = self.v[i] - self.v
            dwij = self.kernel(self.x[i] - self.x, self.h, order = 1)
            wij = self.kernel(self.x[i] - self.x, self.h, order = 0)
            xij = self.x[i] - self.x

            p_rho_i = self.p[i] / (self.rho[i] * self.rho[i])
            p_rho_j = self.p / (self.rho * self.rho)
            
            p_rho_ij =  p_rho_i + p_rho_j

            ### Artificial viscosity ###

            numerator = self.h * vij * xij
            denominator = (xij ** 2.0) + (self.eta ** 2.0)
            mu_ij = numerator / denominator
            mu_ij[mu_ij > 0] = 0.0 #only activated for approaching particles

            ci = (self.gamma * self.p[i] / self.rho[i]) ** 0.5
            cj = (self.gamma * self.p / self.rho) ** 0.5
            cij = 0.5 * (ci + cj)
            rhoij = 0.5 * (self.rho[i] + self.rho)
            numerator = (-1 * cij * mu_ij) + (mu_ij ** 2)
            denominator = rhoij
            pi_ij = numerator / denominator

            ### Evaluating gradients ###

            grad_v = -1 * np.sum(self.m * (p_rho_ij + pi_ij) * dwij)            
            grad_e = 0.5 * np.sum(self.m * (p_rho_ij + pi_ij) * vij * dwij)

            rho_inc[i] = np.sum(self.m * wij)
            v_inc[i] = dt * grad_v
            e_inc[i] = dt * grad_e
            
            ### Get XSPH Velocity and calculate increment ###

            rho_avg = 0.5 * (self.rho[i] + self.rho)
            correction = self.epsilon * np.sum(self.m * -1 * vij * wij / rho_avg)
            xsph_velocity =  self.v[i] + correction

            x_inc[i] = dt * xsph_velocity

            
        ## Update the original arrays using the increment arrays ##
            
        self.rho = rho_inc #not incremented, corrected directly
        self.v += v_inc
        self.e += e_inc

        ### Update pressure ###
        
        self.p = (self.gamma - 1) * self.e * self.rho

        ### Update positions ###

        self.x += x_inc
        
if __name__ == '__main__':

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

