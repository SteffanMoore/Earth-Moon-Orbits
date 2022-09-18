# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 16:10:08 2020

@author: steff
"""
import numpy as np
import math

M_E = 5.9722E24
G = 6.67408E-11
d_E_M = 4e8 #average distance from centre of the Earth to the centre of the Moon
R_E = 6.3714e6
R_M = 1.7381e6


def x_M(x_E):               #finds the x distance to the moon
    return x_E - d_E_M

def y_M(y_E):               #finds the y distance to the moon
    return y_E
    
def F_x(x_E, y_E):
    return (-G*M_E*x_E)/(((x_E)**2+(y_E)**2)**(3/2)) + (-G*M_M*x_M(x_E))/(((x_M(x_E))**2+(y_M(y_E))**2)**(3/2))

def F_y(x_E, y_E):
    return (-G*M_E*y_E)/(((x_E)**2+(y_E)**2)**(3/2)) + (-G*M_M*y_M(y_E))/(((x_M(x_E))**2+(y_M(y_E))**2)**(3/2))
     
def K_vectors(t_max, h, v_x, v_y, x_E, y_E): #this finds the K values and uses them to return the position of the rocket
    for t in range(1, t_max + 1, h): 
        K_1 = np.array([[v_x], [v_y], [F_x(x_E, y_E)], [F_y(x_E, y_E)]])
        K_2 = np.array([[v_x + (h/2)*K_1[2]], [v_y + (h/2)*K_1[3]], [F_x(x_E + (h/2)*K_1[0], y_E + (h/2)*K_1[1])], [F_y(x_E + (h/2)*K_1[0], y_E + (h/2)*K_1[1])]])
        K_3 = np.array([[v_x + (h/2)*K_2[2]], [v_y + (h/2)*K_2[3]], [F_x(x_E + (h/2)*K_2[0], y_E + (h/2)*K_2[1])], [F_y(x_E + (h/2)*K_2[0], y_E + (h/2)*K_2[1])]])
        K_4 = np.array([[v_x + (h/2)*K_3[2]], [v_y + (h/2)*K_3[3]], [F_x(x_E + (h/2)*K_3[0], y_E + (h/2)*K_3[1])], [F_y(x_E + (h/2)*K_3[0], y_E + (h/2)*K_3[1])]])
        
        x_E += (h/6)*(float(K_1[0])+float(K_2[0])+float(K_3[0])+float(K_4[0]))
        y_E += (h/6)*(float(K_1[1])+float(K_2[1])+float(K_3[1])+float(K_4[1]))
        v_x += (h/6)*(float(K_1[2])+float(K_2[2])+float(K_3[2])+float(K_4[2]))
        v_y += (h/6)*(float(K_1[3])+float(K_2[3])+float(K_3[3])+float(K_4[3]))
        
        if math.sqrt((x_E)**2 + (y_E)**2) < R_E:
            break
        elif math.sqrt((x_M(x_E))**2 + (y_M(x_E))**2) < R_M:
            break

        initial_x_y_position.append(x_E)
        initial_x_y_position.append(y_E) 
    return initial_x_y_position

def x_values():
    for n in range(len(x_y_positions)):
        if n % 2 == 0:
            x_list.append(x_y_positions[n])
    return x_list

def y_values():
    for n in range(len(x_y_positions)):
        if n % 2 != 0:
            y_list.append(x_y_positions[n])
    return y_list
 
MyInput = '0'
while MyInput != 'q':
    MyInput = input('Enter a choice: "a" for orbit around the Earth; "b" for orbit around the Earth and the Moon or "q" to quit: ')
    print('You entered the choice: ',MyInput)

    if MyInput == 'a':
        print('You have chosen part (a)')
        
        Input_v_x = input('Enter a value for initial x velocity (m/s): ')
        v_x = float(Input_v_x)
        Input_v_y = input('Enter a value for initial y velocity (m/s): ')
        v_y = float(Input_v_y)
        Input_x_E = input('Enter a value for initial orbit height (m): ')
        x_E = float(Input_x_E)
        y_E = 0
        Input_t_max = input('Enter maximum time (s): ')
        t_max = int(Input_t_max)
        Input_h = input('Enter an integer step size (s): ')
        h = int(Input_h)
        M_M = 0
        
        initial_x_y_position = [x_E, y_E]
        x_list = []
        y_list = []
        
        x_y_positions = K_vectors(t_max, h, v_x, v_y, x_E, y_E)
        x_positions = x_values()
        y_positions = y_values()          

        import matplotlib.pyplot as plt

        x = np.zeros(len(x_positions))
        y = np.zeros(len(y_positions))

        for i in range(len(x_positions)):
            x[i] = x_positions[i]
            y[i] = y_positions[i]

        E_circle = plt.Circle((0, 0), R_E, color = "xkcd:marine blue")
        fig, ax = plt.subplots()
        ax.add_artist(E_circle)
        
        plt.plot(x, y, color = "xkcd:rust red")
        plt.axis("equal")
        plt.xlabel("x position (m)")
        plt.ylabel("y position (m)")
        plt.show()


    elif MyInput == 'b':
        print('You have chosen part (b)')
        
        Input_v_x = input('Enter a value for initial x velocity (m/s): ')
        v_x = float(Input_v_x)
        Input_v_y = input('Enter a value for initial y velocity (m/s): ')
        v_y = float(Input_v_y)
        Input_x_E = input('Enter a value for initial x position (m): ')
        x_E = float(Input_x_E)
        Input_y_E = input('Enter a value for initial y position (m): ')
        y_E = float(Input_y_E)
        Input_t_max = input('Enter maximum time (s): ')
        t_max = int(Input_t_max)
        Input_h = input('Enter a step size (s): ')
        h = int(Input_h)
        M_M = 7.347E22
        
        initial_x_y_position = [x_E, y_E]
        x_list = []
        y_list = []
        
        x_y_positions = K_vectors(t_max, h, v_x, v_y, x_E, y_E)
        x_positions = x_values()
        y_positions = y_values()          

        import matplotlib.pyplot as plt

        x = np.zeros(len(x_positions))
        y = np.zeros(len(y_positions))

        for i in range(len(x_positions)):
            x[i] = x_positions[i]
            y[i] = y_positions[i]
        
        
        E_circle = plt.Circle((0, 0), R_E, color = "xkcd:marine blue")
        M_circle = plt.Circle((4e8, 0), R_M, color = "xkcd:steel grey")
        fig, ax = plt.subplots()
        ax.add_artist(E_circle)
        ax.add_artist(M_circle)
        
        plt.plot(x,y, color = "xkcd:rust red")
        plt.axis("equal")
        plt.xlabel("x position (m)")
        plt.ylabel("y position (m)")
        plt.show()
       

    elif MyInput != 'q':
        print('This is not a valid choice')

print('You have chosen to finish - goodbye.')