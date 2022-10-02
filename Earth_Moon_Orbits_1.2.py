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



def find_info(planet_input, term):               #function to retrieve planetary characteristics stored in dictionaries
    Earth_dictionary = {"radius": 6.3714E6, "mass": 5.9722E24, "start_x": 0, "start_y": 0}
    Moon_dictionary = {"radius": 1.7381E6, "mass": 7.347E22, "start_x": 4E8, "start_y": 0}
    Bodies = {"Earth": Earth_dictionary, "Moon": Moon_dictionary}
    
    try:                        #checks whether theres a dictionary for the planet
        Bodies[planet_input]
    except KeyError:
        return f"The body {planet_input} is unlisted"
    try:                        #checks whether the information wanted is stored
        search = Bodies[planet_input][term]
        return search
    except KeyError:
        return f"The body {planet_input} has no {term} characteristic"


def polar_to_cartesian(angle, radius):                  #quick polar to cartesian coordinate converter
    x = radius*math.sin(angle*math.pi/180)
    y = radius*math.cos(angle*math.pi/180)
    return x, y


def find_hypotoneuse(x, y):
    hyp = math.sqrt((x)**2 + (y)**2)
    return hyp


class massless_object:                                 #creates objects of negligible mass
    def __init__(self, init_x, init_y, init_vel_x, init_vel_y, name):
        self.name = name
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
    
    def forces(self, mass_list):       #returns acting force components
        F_x, F_y = 0, 0
        for i in range(len(mass_list)):         #cycles through masses and adds their contribution to force components
            obj_to_mass = self.obj_dist(mass_list[i].x, mass_list[i].y)     #finds distance from mass
            
            #adds each masses contribution to total force
            F_x += (-G*mass_list[i].mass*obj_to_mass[0])/(((obj_to_mass[0])**2 + (obj_to_mass[1])**2)**(3/2))
            F_y += (-G*mass_list[i].mass*obj_to_mass[1])/(((obj_to_mass[0])**2 + (obj_to_mass[1])**2)**(3/2))
        
        force_components = [F_x, F_y]
        return force_components
    
    def next_iteration(self, h, mass_list):            #calulates the next position of object
        self.temp_x, self.temp_y = self.x, self.y
        K_matrix = np.empty((4,4))              #sets up the empty k matrix for runge-kutta
        
        #runs through rows of the k matrix and assigns them values
        for n in range(4):
            if n < 1:
                acting_force = self.forces(mass_list)
                K_matrix[n] = np.array([self.vel_x, self.vel_y, acting_force[0], acting_force[1]])
            else:                                                   #after setting up the first row, the rest use this
                self.temp_x = self.x + (h/2)*K_matrix[n-1][0]
                self.temp_y = self.y + (h/2)*K_matrix[n-1][1]       #temporary postition values so real position not preemptively changed
                acting_force = self.forces(mass_list)
                K_matrix[n] = np.array([self.vel_x + (h/2)*K_matrix[n-1][2], self.vel_y + (h/2)*K_matrix[n-1][3], acting_force[0], acting_force[1]])

        iteration = []            #sets up a list of changes for next iteration
        for n in range(4):
            iteration.append((h/6)*(K_matrix[0][n] + K_matrix[1][n] + K_matrix[2][n] + K_matrix[3][n]))
        
        
        #iteration is in the order(x, y, vel_x, vel_y)
        return iteration
    
    def keep_iteration(self, changes):              #updates and records positions and velocities
        self.x += changes[0]; self.y += changes[1]; self.vel_x += changes[2]; self.vel_y += changes[3]
        self.record_position()
        

    
class massive_object(massless_object):              #child class of massless object with addition of mass and radius
    def __init__(self, init_x, init_y, init_vel_x, init_vel_y, name, radius, mass):
        super().__init__(init_x, init_y, init_vel_x, init_vel_y, name)
        self.mass = mass
        self.radius = radius
        

    
class simulation:
    def __init__(self):
        possible_massive_objects = ["Earth", "Moon"]
        
        self.massive_object_list = []        #used for collision detection
        self.massless_object_list = []
        
        
        #startup text explaining the simulation
        print("The Earth is automatically generated however other celestial bodies are not.")
        print("The bodies which can currently be simulated are as follows:")
        for i in range(len(possible_massive_objects) - 1):
            print(possible_massive_objects[i + 1])
        
        
        #generates the Earth
        self.massive_object_list.append(massive_object(find_info("Earth", "start_x"), find_info("Earth", "start_y"), 0, 0, "Earth", find_info("Earth", "radius"), find_info("Earth", "mass")))
        bodies_generated = ["Earth"]
        
        
        #body generation section
        body_generation = "y"
        while body_generation != "n":                       #generates masses until told otherwise
            body_input = input("Which, if any, massive objects do you wish to simulate? (press enter if none):")
            
            if body_input == "":
                body_generation = "n"
            elif bodies_generated.count(body_input) >= 1:       #checks whether body already made
                print("Body already generated")
            else:
                try:                        #checks if the input planet already has stored characteristics
                    possible_massive_objects.index(body_input)
                    self.massive_object_list.append(massive_object(find_info(body_input, "start_x"), find_info(body_input, "start_y"), 0, 0, body_input, find_info(body_input, "radius"), find_info(body_input, "mass")))
                    bodies_generated.append(body_input)
                except ValueError:          #if planet isn't a predefined one, checks with user that there's not a syntax error
                    double_check = input(f"Did you make a mistake? If not a massive body of name {body_input} will be generated. (y/n)")
                    if double_check == "y":
                        pass
                    elif double_check == "n":       #once user comfirms, they can generate a new mass
                        new_x = input("Input new body x position:")
                        new_y = input("Input new body y position:")
                        new_r = input("Input new body radius:")
                        new_m = input("Input new body mass:")
                        self.massive_object_list.append(massive_object(new_x, new_y, 0, 0, body_input, new_r, new_m))
                    
        


        #menu for setting up spacecraft starting conditions        
        menu_loop = "repeat"
        while menu_loop != "end repeat":
            menu_loop = "end repeat"
            orbit_height_input = input(f"Input starting orbit height in metres (Earth's radius is {R_E}m):")
            angular_position_input = input("Input the angle at which the spacecraft is in orbit (0 is North with angle increasing clockwise to 359):")
            speed_input = input("Input the initial speed of the spacecraft (m/s):")
            direction_input = input("Input the initial direction of the spacecraft using the same bearings as before:")
            
            #if a mistake is made the user can use this to redo
            redo_input = input("If you aren't happy with your starting conditions, enter 'y' to restart:")
            if redo_input == "y":
                menu_loop = "repeat"
             
                
            #converts the inputs into usable coordinates
            init_x, init_y = polar_to_cartesian(float(angular_position_input), float(orbit_height_input))
            init_vel_x, init_vel_y = polar_to_cartesian(float(direction_input), float(speed_input))
            
            
            #error prevention    
            print("-"*60)
            for i in range(len(self.massive_object_list) - 1):
                distance = find_hypotoneuse(init_x - self.massive_object_list[i + 1].x, init_y - self.massive_object_list[i + 1].y) 
                
                if distance <= self.massive_object_list[i + 1].radius:      #checks for collisions with generated masses
                    print("Initial position collides with an object")
                    menu_loop = "repeat"
                
            if float(orbit_height_input) < find_info("Earth", "radius"):             #prevents immediate collision with Earth
                print("Initial orbit must be greater than the radius of the Earth!")
                menu_loop = "repeat"
            elif float(speed_input) >= 3e8:                 #prevents faster than light travel (shouln't even be APPROACHING these speeds)
                print("Easy tiger, this isn't a relativistic model; no need to go breaking the laws of Physics.")
                menu_loop = "repeat"    
            print("-"*60)
            
        
        #initialisation of the spacecraft object
        spacecraft = massless_object(init_x, init_y, init_vel_x, init_vel_y, "spacecraft")
        self.massless_object_list.append(spacecraft)
        
        
        
        
    def path(self):             #sets up the timescale of the simulation and then carries it out
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
                
                
        #carries out the selected number of runge-kutta iterations
        for t in range(iteration_no):
            collision = "no"                #as soon as a collision is detected, the simulation stops
            
            for i in range(len(self.massless_object_list)):
                iteration_changes = self.massless_object_list[i].next_iteration(h, self.massive_object_list)
                
                #stores the next iteration temporarily to allow testing before accepting the values
                self.massless_object_list[i].temp_x = self.massless_object_list[i].x + iteration_changes[0]
                self.massless_object_list[i].temp_y = self.massless_object_list[i].y + iteration_changes[1]
                
                #runs through stored masses to check for collisions
                for n in range(len(self.massive_object_list)):
                    distances = self.massless_object_list[i].obj_dist(self.massive_object_list[n].x, self.massive_object_list[n].y)
                    distance_to_mass = find_hypotoneuse(distances[0], distances[1])
                    
                    if distance_to_mass > self.massive_object_list[n].radius:
                        pass
                    else:
                        collision_speed = find_hypotoneuse(self.massless_object_list[i].vel_x, self.massless_object_list[i].vel_y)
                        print(f"\nSpacecraft has crashed into {self.massive_object_list[n].name} at {(t+1)*h}s going at a speed of {collision_speed}m/s")
                        collision = "yes"
                
                if collision == "no":       #if the next iteration passes without a collision, it's kept
                    self.massless_object_list[i].keep_iteration(iteration_changes)
                else:
                    break
                
            if collision == "yes":
                break
            
            
            
            
    def plot(self):         #graphs the simulation once complete
        for n in range(len(self.massless_object_list)):         #cycles through all massless objects setting up their paths over time
            path = self.massless_object_list[n].position_history()
            x = np.zeros(len(path))
            y = np.zeros(len(path))
            for i in range(len(path)):
                x[i] = path[i][0]
                y[i] = path[i][1]
        
        #stores colours for a number of predefined planets
        colours = {"Mercury": "sand brown", "Venus": "gold", "Earth": "marine blue", "Moon": "steel grey", "Mars": "rust red", "Jupiter": "apricot", "Saturn": "buff", "Uranus": "eggshell blue", "Neptune": "cerulean blue", "Pluto": "light grey"}
        fig, ax = plt.subplots()
        
        for n in range(len(self.massive_object_list)):      #selects the colour used to display planets and adds them to the plot
            try:
                colour = "xkcd:" + colours[self.massive_object_list[n].name]
            except KeyError:
                colour = "xkcd:violet"
                
            mass_circle = plt.Circle((self.massive_object_list[n].x, self.massive_object_list[n].y), self.massive_object_list[n].radius, color = colour)
            ax.add_artist(mass_circle)
        
        plt.plot(x,y, color = "xkcd:rust red")
        plt.axis("equal")
        plt.xlabel("x position (m)")
        plt.ylabel("y position (m)")
        plt.show()



model = simulation()
model.path()
model.plot()