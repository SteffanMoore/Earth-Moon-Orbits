# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 16:10:08 2020

@author: steff
"""
import numpy as np
import math
import matplotlib.pyplot as plt

M_E = 5.9722E24
G = 6.67408E-11
d_E_M = 4e8 #average distance from centre of the Earth to the centre of the Moon
R_E = 6.3714e6
R_M = 1.7381e6
M_M = 7.347E22


class massless_object:                                 #creates objects of negligible mass
    def __init__(self, init_x, init_y, init_vel_x, init_vel_y):
        self.x, self.y = init_x, init_y                             #sets up the initial position
        self.temp_x, self.temp_y = self.x, self.y                   #tempory position for use with k vectors
        self.vel_x, self.vel_y = init_vel_x, init_vel_y             #sets up the initial velocity
        self.history = [[init_x, init_y]]
        
    def record_position(self):                              #adds coordinates to the history list
        self.history.append([self.x, self.y])
        
    def position_history(self):                     #returns coordinate history
        return self.history
    
    def obj_dist(self, x, y):                 #returns x and y distance from object to a point
        distance = [self.temp_x - x, self.temp_y - y]
        return distance
    
    def forces(self, x_E, y_E, x_M, y_M):       #returns acting force components
        #sets up distances to celestial bodies
        obj_E = self.obj_dist(x_E, y_E)
        obj_M = self.obj_dist(x_M, y_M)
        
        #forces acting on the object
        F_x = (-G*M_E*obj_E[0])/(((obj_E[0])**2 + (obj_E[1])**2)**(3/2)) + (-G*M_M*obj_M[0])/(((obj_M[0])**2 + (obj_M[1])**2)**(3/2))
        F_y = (-G*M_E*obj_E[1])/(((obj_E[0])**2 + (obj_E[1])**2)**(3/2)) + (-G*M_M*obj_M[1])/(((obj_M[0])**2 + (obj_M[1])**2)**(3/2))
        
        force_components = [F_x, F_y]
        return force_components
    
    def next_iteration(self, h, x_E, y_E, x_M, y_M):            #calulates the next position of object
        #runge kutta K_vectors
        self.temp_x, self.temp_y = self.x, self.y
        K_matrix = np.empty((4,4))
        
        for n in range(4):                  #runs through rows of the k matrix and assigns them values
            if n < 1:
                acting_force = self.forces(x_E, y_E, x_M, y_M)
                K_matrix[n] = np.array([self.vel_x, self.vel_y, acting_force[0], acting_force[1]])
            else:                                                   #after setting up the first row, the rest use this
                self.temp_x = self.x + (h/2)*K_matrix[n-1][0]
                self.temp_y = self.y + (h/2)*K_matrix[n-1][1]       #temporary postition values so real position not preemptively changed
                acting_force = self.forces(x_E, y_E, x_M, y_M)
                K_matrix[n] = np.array([self.vel_x + (h/2)*K_matrix[n-1][2], self.vel_y + (h/2)*K_matrix[n-1][3], acting_force[0], acting_force[1]])

        iteration = []            #sets up a list of changes for next iteration
        for n in range(4):
            iteration.append((h/6)*(K_matrix[0][n] + K_matrix[1][n] + K_matrix[2][n] + K_matrix[3][n]))
        
        
        #iteration is in the order(x, y, vel_x, vel_y)
        return iteration
    
    def keep_iteration(self, changes):              #updates and records positions and velocities
        self.x += changes[0]; self.y += changes[1]; self.vel_x += changes[2]; self.vel_y += changes[3]
        self.record_position()
        

    
class massive_object(massless_object):
    def __init__(self, init_x, init_y, init_vel_x, init_vel_y, radius, mass):
        super().__init__(init_x, init_y, init_vel_x, init_vel_y)
        self.mass = mass
        self.radius = radius
        

def polar_to_cartesian(angle, radius):                  #quick polar to cartesian coordinate converter
    x = radius*math.sin(angle*math.pi/180)
    y = radius*math.cos(angle*math.pi/180)
    return x, y


def find_hypotoneuse(x, y):
    hyp = math.sqrt((x)**2 + (y)**2)
    return hyp

    
class simulation:
    def __init__(self):
        self.massive_object_list = []        #used for collision detection
        self.massless_object_list = []
        
        self.Earth = massive_object(0, 0, 0, 0, R_E, M_E)
        self.massive_object_list.append(self.Earth)
        
        Moon_input = "n"
        while Moon_input != "y":
            Moon_input = input("Would you like the moon to be simulated? (y/n):")
            if Moon_input == "y":
                self.Moon = massive_object(d_E_M, 0, 0, 0, R_M, M_M)
                self.massive_object_list.append(self.Moon)
            elif Moon_input == "n":
                self.Moon = massive_object(d_E_M, 0, 0, 0, R_M, 0)
                Moon_input = "y"
            else:
                print("Please enter a correct answer")
        
        menu_loop = "repeat"
        while menu_loop != "end repeat":            #menu for starting conditions
            menu_loop = "end repeat"
            orbit_height_input = input(f"Input starting orbit height in metres (Earth's radius is {R_E}m):")
            angular_position_input = input("Input the angle at which the spacecraft is in orbit (0 is North with angle increasing clockwise to 359):")
            speed_input = input("Input the initial speed of the spacecraft (m/s):")
            direction_input = input("Input the initial direction of the spacecraft using the same bearings as before:")
            
            #if a mistake is made the user can use this to redo
            redo_input = input("If you aren't happy with your starting conditions, enter 'y' to restart:")
            if redo_input == "y":
                menu_loop = "repeat"
                
            #error prevention    
            print("-"*60)
            if float(orbit_height_input) < R_E:             #prevents immediate collision
                print("Initial orbit must be greater than the radius of the Earth!")
                menu_loop = "repeat"
            elif float(speed_input) >= 3e8:                 #prevents faster than light travel (shouln't even be APPROACHING these speeds)
                print("Easy tiger, this isn't a relativistic model; no need to go breaking the laws of Physics.")
                menu_loop = "repeat"    
            print("-"*60)
            
        #converts the inputs into usable coordinates
        init_x, init_y = polar_to_cartesian(float(angular_position_input), float(orbit_height_input))
        init_vel_x, init_vel_y = polar_to_cartesian(float(direction_input), float(speed_input))
        
        #initialisation of the spacecraft object
        spacecraft = massless_object(init_x, init_y, init_vel_x, init_vel_y)
        self.massless_object_list.append(spacecraft)
        
    def path(self):
        menu_loop = "repeat"
        while menu_loop != "end repeat":            #menu for timestep
            menu_loop = "end repeat"
            time_input = input("Input the time to run the simulation:")
            timestep_input = input("Input the timestep of the simulation:")
            
            #converts inputs to give a number of iterations and an actual simulation time
            t_given, h = float(time_input), float(timestep_input)
            iteration_no = int(t_given/h)
            t_max = iteration_no * h
            
            redo_input = input(f"The simulation will be run for {t_max}s over {iteration_no} iterations with a timestep of {h}. If you wish to change this, press 'n':")
            if redo_input == "n":
                menu_loop = "repeat"
            
        for t in range(iteration_no):
            collision = "no"
            
            for i in range(len(self.massless_object_list)):
                iteration_changes = self.massless_object_list[i].next_iteration(h, self.Earth.x, self.Earth.y, self.Moon.x, self.Moon.y)
               
                self.massless_object_list[i].temp_x = self.massless_object_list[i].x + iteration_changes[0]
                self.massless_object_list[i].temp_y = self.massless_object_list[i].y + iteration_changes[1]
                
                Earth_distances = self.massless_object_list[i].obj_dist(self.Earth.x, self.Earth.y)
                Moon_distances = self.massless_object_list[i].obj_dist(self.Moon.x, self.Moon.y)
                
                obj_to_Earth = find_hypotoneuse(Earth_distances[0], Earth_distances[1])
                obj_to_Moon = find_hypotoneuse(Moon_distances[0], Moon_distances[1])
                
                if obj_to_Earth > self.Earth.radius and obj_to_Moon > self.Moon.radius:
                    self.massless_object_list[i].keep_iteration(iteration_changes)
                else:
                    print(f"Spacecraft has crashed at {(t+1)*h}s")
                    collision = "yes"
                    break
                
            if collision == "yes":
                break
            
    def plot(self):
        for n in range(len(self.massless_object_list)):         #cycles through all massless objects
            path = self.massless_object_list[n].position_history()
            x = np.zeros(len(path))
            y = np.zeros(len(path))
            for i in range(len(path)):
                x[i] = path[i][0]
                y[i] = path[i][1]
                
        E_circle = plt.Circle((0, 0), self.Earth.radius, color = "xkcd:marine blue")
        M_circle = plt.Circle((4e8, 0), self.Moon.radius, color = "xkcd:steel grey")
        fig, ax = plt.subplots()
        ax.add_artist(E_circle)
        ax.add_artist(M_circle)
        
        plt.plot(x,y, color = "xkcd:rust red")
        plt.axis("equal")
        plt.xlabel("x position (m)")
        plt.ylabel("y position (m)")
        plt.show()



model = simulation()
model.path()
model.plot()