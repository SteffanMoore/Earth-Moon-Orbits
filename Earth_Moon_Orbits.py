import numpy as np
import math
import matplotlib.pyplot as plt

M_E = 5.9722E24
G = 6.67408E-11
d_E_M = 4e8 #average distance from centre of the Earth to the centre of the Moon
R_E = 6.3714e6
R_M = 1.7381e6
M_M = 7.347E22


def main():
    model = simulation()
    model.path()
    model.plot()


def find_info(planet_input, term):
    """
    Retrieve planetary characteristics stored in dictionaries
    """

    Earth_dictionary = {"radius": 6.3714E6, "mass": 5.9722E24, "start_x": 0, "start_y": 0}
    Moon_dictionary = {"radius": 1.7381E6, "mass": 7.347E22, "start_x": 4E8, "start_y": 0}
    Bodies = {"Earth": Earth_dictionary, "Moon": Moon_dictionary}
    
    # Checks whether there is a dictionary and term for the input planet and term
    if planet_input in Bodies:
        try:
            search = Bodies[planet_input][term]
            return search
        except KeyError:
            return f"The body {planet_input} has no {term} characteristic"
    else:
        raise KeyError(f"The body {planet_input} is unlisted")


def polar_to_cartesian(angle, radius):
    """
    Input polar coordinates to return cartesian (x, y) coordinates.
    """

    x = radius*math.sin(angle * math.pi / 180)
    y = radius*math.cos(angle * math.pi / 180)

    return x, y


def find_hypotoneuse(x, y):
    """
    Returns the hypoteneuse of a triangle with sides of length x and y.
    """

    hyp = math.sqrt((x)**2 + (y)**2)

    return hyp


class massless_object:
    """
    The massless object class creates a object of negligible mass on astronomical scale (eg. a spaceship). This object is
    capable of calculating the current forces acting on it and finding its position after the next time iteration.
    """

    def __init__(self, init_x, init_y, init_vel_x, init_vel_y, name):
        self.name = name
        self.x, self.y = init_x, init_y                             #sets up the initial position
        self.temp_x, self.temp_y = self.x, self.y                   #tempory position for use with k vectors
        self.vel_x, self.vel_y = init_vel_x, init_vel_y             #sets up the initial velocity
        self.history = [[init_x, init_y]]
        
    def record_position(self):
        """
        Adds current coordinates to the history list.
        """

        self.history.append([self.x, self.y])
        
    def position_history(self):
        """
        Returns all previous coordinates at which the object 
        """
        
        return self.history
    
    def obj_dist(self, x, y):
        """
        Returns the x and y distance between the object and a given point
        """

        distance = [self.temp_x - x, self.temp_y - y]

        return distance
    
    def forces(self, mass_list):
        """
        Returns the x and y force components acting on the object
        """

        x_force, y_force = 0, 0

        # Cycles through the masses and adds their contributions to the force components
        for mass in mass_list:
            mass_distance = self.obj_dist(mass.x, mass.y)
            force_eq_denominator = ((mass_distance[0])**2 + (mass_distance[1])**2)**(3/2)
            
            x_force += (-G * mass.mass * mass_distance[0]) / force_eq_denominator
            y_force += (-G * mass.mass * mass_distance[1]) / force_eq_denominator
        
        force_components = [x_force, y_force]

        return force_components
    
    def next_iteration(self, h, mass_list):
        """
        Calculates the position of the object after the next iteration
        """

        self.temp_x, self.temp_y = self.x, self.y
        K_matrix = np.empty((4, 4))
        
        # Assigns each row of the K matrix values to be used in Runge-Kutta
        for n in range(4):
            if n < 1:
                acting_force = self.forces(mass_list)
                K_matrix[n] = np.array([self.vel_x, self.vel_y, acting_force[0], acting_force[1]])

            # After the first row, temporary position values are used to find forces without altering the actual position of the object
            else:
                self.temp_x = self.x + (h / 2) * K_matrix[n - 1][0]
                self.temp_y = self.y + (h / 2) * K_matrix[n - 1][1]
                acting_force = self.forces(mass_list)
                K_matrix[n] = np.array([self.vel_x + (h / 2) * K_matrix[n - 1][2], self.vel_y + (h / 2) * K_matrix[n - 1][3],
                                        acting_force[0], acting_force[1]])

        # Next iteration position and velocity is calculated using the K matrix
        iteration = []
        for n in range(4):
            iteration.append((h / 6) * (K_matrix[0][n] + K_matrix[1][n] + K_matrix[2][n] + K_matrix[3][n]))
        
        # Iteration returned as [x, y, vel_x, vel_y]
        return iteration
    
    def keep_iteration(self, changes):
        """
        Updates and records positions and velocities
        """

        self.x += changes[0]; self.y += changes[1]; self.vel_x += changes[2]; self.vel_y += changes[3]
        self.record_position()
        

    
class massive_object(massless_object):
    """
    Adds mass and radius to simulate a astronomically significant object
    """

    def __init__(self, init_x, init_y, init_vel_x, init_vel_y, name, radius, mass):
        super().__init__(init_x, init_y, init_vel_x, init_vel_y, name)
        self.mass = mass
        self.radius = radius
        

class simulation:
    """
    Contains all objects to be simulated and 
    """

    def __init__(self):
        possible_massive_objects = ["Earth", "Moon"]
        
        self.massive_object_list = []        # Used for collision detection
        self.massless_object_list = []
        
        
        # Startup text explaining the simulation
        print("The Earth is automatically generated however other celestial bodies are not.")
        print("The bodies which can currently be simulated are as follows:")
        for i in range(len(possible_massive_objects) - 1):
            print(possible_massive_objects[i + 1])
        
        
        # Generates the Earth
        self.massive_object_list.append(massive_object(find_info("Earth", "start_x"), find_info("Earth", "start_y"), 0, 0, "Earth", find_info("Earth", "radius"), find_info("Earth", "mass")))
        bodies_generated = ["Earth"]
        
        
        # Body generation section
        body_generation = "y"
        while body_generation != "n":
            body_input = input("Which, if any, massive objects do you wish to simulate? (press enter if none):")
            
            if body_input == "":
                body_generation = "n"
            elif bodies_generated.count(body_input) >= 1:
                print("Body already generated")
            else:
                # Checks if the input planet already has stored characteristics throwing error if not
                try:
                    possible_massive_objects.index(body_input)
                    self.massive_object_list.append(massive_object(find_info(body_input, "start_x"), find_info(body_input, "start_y"), 0, 0, body_input, find_info(body_input, "radius"), find_info(body_input, "mass")))
                    bodies_generated.append(body_input)
                except ValueError:
                    name_mispelled = input(f"Did you make a mistake? If not a massive body of name {body_input} will be generated. (y/n)")
                    if name_mispelled == "y":
                        pass
                    elif name_mispelled == "n":
                        new_x = input("Input new body x position:")
                        new_y = input("Input new body y position:")
                        new_r = input("Input new body radius:")
                        new_m = input("Input new body mass:")
                        self.massive_object_list.append(massive_object(new_x, new_y, 0, 0, body_input, new_r, new_m))
                    
        


        # Menu for setting up spacecraft starting conditions        
        menu_loop = "repeat"
        while menu_loop != "end repeat":
            menu_loop = "end repeat"
            orbit_height_input = input(f"Input starting orbit height in metres (Earth's radius is {R_E}m):")
            angular_position_input = input("Input the angle at which the spacecraft is in orbit (0 is North with angle increasing clockwise to 359):")
            speed_input = input("Input the initial speed of the spacecraft (m/s):")
            direction_input = input("Input the initial direction of the spacecraft using the same bearings as before:")
            
            # If a mistake is made the user can use this to redo
            redo_input = input("If you aren't happy with your starting conditions, enter 'y' to restart:")
            if redo_input == "y":
                menu_loop = "repeat"
             
                
            # Converts the inputs into usable coordinates
            init_x, init_y = polar_to_cartesian(float(angular_position_input), float(orbit_height_input))
            init_vel_x, init_vel_y = polar_to_cartesian(float(direction_input), float(speed_input))
            
            
            # Error prevention    
            print("-"*60)
            for i in range(len(self.massive_object_list) - 1):
                distance = find_hypotoneuse(init_x - self.massive_object_list[i + 1].x, init_y - self.massive_object_list[i + 1].y) 
                
                # Checks for collisions with generated masses
                if distance <= self.massive_object_list[i + 1].radius:
                    print("Initial position collides with an object")
                    menu_loop = "repeat"
                
            # Prevents immediate collision with Earth
            if float(orbit_height_input) < find_info("Earth", "radius"):
                print("Initial orbit must be greater than the radius of the Earth!")
                menu_loop = "repeat"
            # Prevents faster than light travel (shouldn't even be APPROACHING these speeds)
            elif float(speed_input) >= 3e8:
                print("Easy tiger, this isn't a relativistic model; no need to go breaking the laws of Physics.")
                menu_loop = "repeat"    
            print("-"*60)
            
        
        # Initialisation of the spacecraft object
        spacecraft = massless_object(init_x, init_y, init_vel_x, init_vel_y, "spacecraft")
        self.massless_object_list.append(spacecraft)
        
    def path(self):
        """
        Sets up the timescale of the simulation and then carries it out
        """

        menu_loop = "repeat"
        while menu_loop != "end repeat":            # Menu for timestep
            menu_loop = "end repeat"
            time_input = input("Input the time to run the simulation:")
            timestep_input = input("Input the timestep of the simulation:")
            
            # Converts inputs to give a number of iterations and an actual simulation time
            t_given, h = float(time_input), float(timestep_input)
            iteration_no = int(t_given/h)
            t_max = iteration_no * h
            
            redo_input = input(f"The simulation will be run for {t_max}s over {iteration_no} iterations with a timestep of {h}. If you wish to change this, press 'n':")
            if redo_input == "n":
                menu_loop = "repeat"
                
                
        # Carries out the selected number of runge-kutta iterations
        for t in range(iteration_no):
            collision = "no"
            
            for massless_object in self.massless_object_list:
                iteration_changes = massless_object.next_iteration(h, self.massive_object_list)
                
                # Stores the next iteration temporarily to allow testing before accepting the values
                massless_object.temp_x = massless_object.x + iteration_changes[0]
                massless_object.temp_y = massless_object.y + iteration_changes[1]
                
                # Runs through stored masses to check for collisions
                for massive_object in self.massive_object_list:
                    distances = massless_object.obj_dist(massive_object.x, massive_object.y)
                    distance_to_mass = find_hypotoneuse(distances[0], distances[1])
                    
                    if distance_to_mass <= massive_object.radius:
                        collision_speed = find_hypotoneuse(massless_object.vel_x, massless_object.vel_y)
                        print(f"\nSpacecraft has crashed into {massive_object.name} at {(t + 1) * h}s going at a speed of {collision_speed}m/s")
                        collision = "yes"
                
                # If there's no collision, the iteration is kept, otherwise the simulation ends
                if collision == "no":
                    massless_object.keep_iteration(iteration_changes)
                else:
                    break
                
            if collision == "yes":
                break
             
    def plot(self):
        """
        Plot the simulation once completed
        """

        for massless_object in self.massless_object_list:         #cycles through all massless objects setting up their paths over time
            path = massless_object.position_history()
            x = np.zeros(len(path))
            y = np.zeros(len(path))
            for i, position in enumerate(path):
                x[i] = position[0]
                y[i] = position[1]
        
        # Stores colours for a number of predefined planets
        colours = {"Mercury": "sand brown", "Venus": "gold", "Earth": "marine blue", "Moon": "steel grey",
                   "Mars": "rust red", "Jupiter": "apricot", "Saturn": "buff", "Uranus": "eggshell blue",
                   "Neptune": "cerulean blue", "Pluto": "light grey"}
        fig, ax = plt.subplots()
        
        # Selects the colour used to display planets and adds them to the plot
        for n in range(len(self.massive_object_list)):
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

if __name__ == "__main__":
    main()