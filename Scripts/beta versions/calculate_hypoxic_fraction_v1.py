#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun May  9 22:09:58 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Calculate Hypoxic Fraction (Dichotomous Network)

Tested in Python 3.7.4.

"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
from collections import deque
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import scipy
import time
import vtk
from vtk.util.numpy_support import vtk_to_numpy

# Starts stopwatch to clock execution time
start_time = time.time()

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read a .vti file and return the field distribution
def get_vti_data(input_filename):
    
    # Import the file
    extension = input_filename.split('.').pop()
    reader = None
    if extension == 'vtk':
        reader = vtk.vtkDataSetReader() 
    elif extension == 'vti':
        reader = vtk.vtkXMLImageDataReader() 
    else:
        raise RuntimeError('Unknown File Type: %s ' % input_filename)
    reader.SetFileName( "%s" % (input_filename) ) 
    reader.Update()
    image_data = reader.GetOutput()
    
    # Extract the dimensions, spacing, and origin
    spacing = image_data.GetSpacing()
    
    # Extract the point values 
    field_point_data = image_data.GetPointData() 
    field_values = vtk_to_numpy(field_point_data.GetArray(0)) 
    
    # Get the coordinates of each point
    position_list = deque()
    for index in range(len(field_values)):  # do something 
        position = image_data.GetPoint(index)
        position_list.append(position)
    position_array = np.array(position_list)
    
    # Return the field distribution
    distribution_array = np.column_stack((position_array,field_values))
    distribution = pd.DataFrame(data=distribution_array[:,[0,1,3]], columns=["x", "y", "oxygen"])
    return distribution, spacing#, dimensions, origin

# Define a function to return a field and basic stats  
def get_field(paraview_data):
    
    # Change the data type
    paraview_data = paraview_data.astype('float32')
    
    # Generates pair plots
    #sns.pairplot(dataset[['x', 'y', 'conc']], diag_kind="kde")
    
    # Calculate concentration statistics
    O2_stats = paraview_data['oxygen'].describe()
    
    # Downsample data if needed
    paraview_data = paraview_data[::int(1)]
    
    # Convert dataframe into NumPy matrix
    mat = paraview_data.to_numpy()
    
    # Get the x and y axes
    x = np.unique(mat[:,0])
    y = np.unique(mat[:,1])
    
    # Create a mesh from the x and y axes 
    X,Y = np.meshgrid(x, y)
    
    # Interpolate the concentration values over the mesh
    Z = scipy.interpolate.griddata((mat[:,0], mat[:,1]), mat[:,2], (X,Y), method='nearest')

    # Return the oxygen stats
    return Z, O2_stats

# Define a function to return the hypoxic fraction of a field  
def calculate_hypoxic_fraction(oxygen_field, hypoxic_threshold):
    
    # Calculate the number of points below hypoxic threshold    
    hypoxic_points = (oxygen_field < hypoxic_threshold).sum()
    
    # Get the number of total points
    number_of_points = 1
    for dim in np.shape(oxygen_field): number_of_points *= dim
    
    # Calculate the hypoxic fraction
    hypoxic_fraction = hypoxic_points/number_of_points
    
    # Return the hypoxic fraction
    return hypoxic_fraction

# Define a function to read a file and return the hypoxic fraction
def get_hypoxic_fraction(solver_name, alpha_value, radius_threshold, hypoxic_threshold, plot=0):    
    
    # Set the file path
    file_path = '/Users/vedang/Simulations/Haematocrit Solver Comparisons/TestDichotomousNetwork4SolversAlpha1.4/' \
    + solver_name + 'Haematocrit/Alpha' + alpha_value + '/LambdaEquals4/RadiusThresholdEquals' \
    + radius_threshold + 'um/oxygen_solution_0.vti'
    
    # Import the data for the equilateral network
    field_data, field_spacing = get_vti_data(file_path)
    
    # Get the O2 distributions and stats
    O2_field, O2_stats = get_field(field_data)
    
    # Plot the O2 distribution
    if plot==1:
        fig = plt.figure()
        ax = plt.axes()
        colour_map = plt.cm.get_cmap('jet')
        plt.suptitle('O$_2$ distribution generated by the ' + solver_name + ' solver in the dichotomous vessel network with alpha = ' + alpha_value + ' and radius threshold = ' + radius_threshold + '  μm (1 unit = ' + str(field_spacing [0]) + ' μm)')
        ref_map = ax.imshow(O2_field, cmap=colour_map, origin='lower')
        fig.colorbar(ref_map, ax=ax, label='nM')
        plt.show()
    
    # Calculate the number of points below hypoxic threshold    
    hypoxic_points = (O2_field < hypoxic_threshold).sum()
    
    # Get the number of total points
    number_of_points = 1
    for dim in np.shape(O2_field): number_of_points *= dim
    
    # Calculate the hypoxic fraction
    hypoxic_fraction = hypoxic_points/number_of_points
    
    # Return the hypoxic fraction
    return hypoxic_fraction

# Define a function to return the hypoxic fraction data for a solver
def get_hypoxic_fraction_stats(solver_name, alpha_list, radius_threshold_list):
    table = np.array([])
    for alpha_value in alpha_list:
        alpha_table = np.array([])
        for radius_threshold in radius_threshold_list:
            hypoxic_fraction = get_hypoxic_fraction(solver_name, alpha_value, radius_threshold, hypoxic_threshold)
            table_entry = np.array([float(alpha_value), float(radius_threshold), hypoxic_fraction])
            alpha_table = np.vstack([alpha_table, table_entry]) if alpha_table.size else table_entry
        table = np.hstack([table, alpha_table]) if table.size else alpha_table
    return table

# Define a function to plot a 3D plot for the hypoxic fraction
def plot_solver_stats(solver_data, solver_name, hypoxic_threshold):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.gca().invert_xaxis()
    plt.suptitle(solver_name + ' haematocrit solver in the dichotomous vessel network (hypoxic threshold = ' + str(hypoxic_threshold) + ' nM )')
    ax.plot(solver_data[:,0], solver_data[:,1], solver_data[:,2])
    ax.plot(solver_data[:,3], solver_data[:,4], solver_data[:,5])
    ax.plot(solver_data[:,6], solver_data[:,7], solver_data[:,8])
    ax.plot(solver_data[:,9], solver_data[:,10], solver_data[:,11])
    ax.plot(solver_data[:,12], solver_data[:,13], solver_data[:,14])
    ax.set_xlabel('alpha')
    ax.set_ylabel('radius threshold (μm)')
    ax.set_zlabel('hypoxic fraction')
    ax.view_init(elev=40., azim=60)
    plt.show()

# =============================================================================
# HYPOXIC FRACTION DATA
# =============================================================================

# Enter details to allow looping over folders
solver_list = ['Constant', 'Pries', 'Memory', 'Betteridge']
alpha_list = ['1', '1.1', '1.2', '1.3', '1.4']
max_radius = 35  #50
radius_threshold_list = [str(x) for x in range(max_radius + 1)]

# Set the threshold for hypoxia
hypoxic_threshold = 7331.5  #O2_stats['25%'] # set hypoxia threshold (<0.1% mmHg is radioresistance, <2% is hypoxia)

# Get the hypoxic fraction data for solvers
constant_data = get_hypoxic_fraction_stats(solver_list[0], alpha_list, radius_threshold_list)
pries_data = get_hypoxic_fraction_stats(solver_list[1], alpha_list, radius_threshold_list)
memory_data = get_hypoxic_fraction_stats(solver_list[2], alpha_list, radius_threshold_list)
betteridge_data = get_hypoxic_fraction_stats(solver_list[3], alpha_list, radius_threshold_list)

# Plot the data
plot_solver_stats(constant_data, solver_list[0], hypoxic_threshold)
plot_solver_stats(pries_data, solver_list[1], hypoxic_threshold)
plot_solver_stats(memory_data, solver_list[2], hypoxic_threshold)
plot_solver_stats(betteridge_data, solver_list[3], hypoxic_threshold)

# prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))