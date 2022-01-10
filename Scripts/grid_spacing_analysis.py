#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Jan  8 18:27:10 2022

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Tested in Python 3.7.4.

Calculate metrics for unpruned forking network at different grid resolutions.

Used to generate figures for Transfer Report addendum.

"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
from collections import deque
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate
import time
import vtk
from vtk.util.numpy_support import vtk_to_numpy

# Starts stopwatch to clock execution time
start_time = time.time()

# Set LaTex-style font
from pathlib import Path
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read a .vti file and return the data
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
#    print(spacing)
    
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

# Define a function to convert the .vti data to a plottable field with some basic stats
def get_field(paraview_data):
    
    # Change the data type
    paraview_data = paraview_data.astype('float32')
    
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

# Define a function to extract the middle distribution of the field
def get_middle_portion(node_file_path, field_dataset, middle_generation_number):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(node_file_path)
    reader.Update()
    polydata = reader.GetOutput()
    points = polydata.GetPoints()
    array = points.GetData()
    point_coordinates = vtk_to_numpy(array)
    generation_coordinates = np.unique(point_coordinates[:,0])
    x_start = generation_coordinates[middle_generation_number-1]
    x_end = generation_coordinates[-middle_generation_number]
    field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end))]
    return field_dataset
#    return vessel_network, oxygen_distribution, middle_x

# Define a function to read an O2 distribution .vti file and return the basic stats
def get_distribution_stats(grid_spacing, hypoxic_threshold_list, plot=0, read=0):    
    
    # Set the file path
    folder_path = '/Users/vedang/Simulations/Grid Spacing Analysis/TestDichotomousNetwork/MemoryHaematocrit/Alpha1.4/NoPruning/GridSpacingEquals' + str(grid_spacing)
#    folder_path = '/scratch/narain/Stochastic Pruning in Forking Network with 100 Trials adjusted/TestDichotomousNetwork/' + solver_name + 'Haematocrit/Lambda4/Alpha' + alpha_value + '/Selection' + layout_selection + '/Beta' + beta
    field_path =  folder_path + '/oxygen_solution_0.vti'    
    network_path = folder_path + '/FinalHaematocrit.vtp'
    
    # Print status update
    print(field_path)
    
    # Save file or read from existing file
    if read==0:    
        field_data, field_spacing = get_vti_data(field_path)  # import the data for the network
        middle_field = get_middle_portion(network_path, field_data, 6)  # get the distribution in the middle of the field (replace with designated file)
        if plot==1:  # plot the O2 distribution
            O2_field, _ = get_field(middle_field)  # get the O2 mesh for plots and stats
            fig = plt.figure()
            ax = plt.axes()
            colour_map = plt.cm.get_cmap('jet')
            ref_map = ax.imshow(O2_field, cmap=colour_map, origin='lower')
            fig.colorbar(ref_map, ax=ax, label='nM')
            plt.show()
        flat_field = middle_field['oxygen'].to_numpy().flatten()
  
    # Get the number of total points
    number_of_points = 1
    for dim in np.shape(flat_field): number_of_points *= dim

    # Calculate the number of points below the given hypoxic thresholds
    hypoxic_fraction_list = []
    for hypoxic_threshold in hypoxic_threshold_list:
        hypoxic_points = (flat_field < hypoxic_threshold).sum()
        hypoxic_fraction = hypoxic_points/number_of_points  # calculate the hypoxic fraction
        hypoxic_fraction_list.append(hypoxic_fraction)

    # Return the stats
    return hypoxic_fraction_list, np.mean(flat_field), np.amin(flat_field), np.percentile(flat_field, 50), np.amax(flat_field), np.std(flat_field)

# Define a function to return statistics for all the heterogeneities in the data
def get_grid_stats(grid_spacing_list, hypoxic_threshold_list, plot, read):
    grid_table = np.array([])
    for grid_spacing in grid_spacing_list:
            grid_data = get_distribution_stats(grid_spacing, hypoxic_threshold_list, plot, read)
            table_entry = np.hstack([float(grid_spacing), grid_data])
            grid_table = np.vstack([grid_table, table_entry]) if grid_table.size else table_entry
    return grid_table

# =============================================================================
# COMPUTATION TIMES
# =============================================================================

# Enter the data
grid_spacing_list = [1,2,4,8,10,16,32,64,128,256,512]
computation_times = [270.43, 57.72, 15.65, 5.41, 4.43, 3.00, 2.19, 2.16, 2.03, 2.02, 1.93]
hypoxic_threshold_list = [2195, 27441] 

# Plot the computation times
fig = plt.figure(figsize=(20,12), tight_layout = {'pad': 2})
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
linecolours = ['#1f77b4','r']
linelegends = ['anoxia','hypoxia']
plt.xscale('log')
plt.axvline(10, ls='dotted', c='black')
#plt.yscale('log')
plt.plot(grid_spacing_list, computation_times, marker='o', ls = 'solid')
#for i, txt in enumerate(computation_times):
#    plt.annotate(txt, (grid_spacing_list[i], computation_times[i]),  fontsize=10)
plt.ylabel('execution time (s)')
plt.xlabel('grid spacing (μm)')
plt.grid()
plt.show()

# Save image
file_path = Path('~/Desktop/computation_times.svg').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
file_path = Path('~/Desktop/computation_times.png').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')

# =============================================================================
# DISTRIBUTION STATS & HYPOXIC FRACTIONS
# =============================================================================

# Get the stats for all solvers
grid_stats = get_grid_stats(grid_spacing_list, hypoxic_threshold_list, plot=0, read=0)

# Filter by grid spacing
mean_composite = np.array([])
min_composite = np.array([])
half_composite = np.array([])
max_composite = np.array([])
sd_composite = np.array([])
hypoxic_fraction_composite = np.array([])
for grid_spacing in grid_spacing_list:
    spacing_array = grid_stats[(grid_stats[:,0]==float(grid_spacing))]
    hypoxic_fraction_data = spacing_array[:,2:-5]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
    mean_data = spacing_array[:,-5]
    min_data = spacing_array[:,-4]
    half_data = spacing_array[:,-3]
    max_data = spacing_array[:,-2]
    sd_data = spacing_array[:,-1]
    hypoxic_fraction_composite = np.vstack([hypoxic_fraction_composite, hypoxic_fraction_data]) if hypoxic_fraction_composite.size else hypoxic_fraction_data
    mean_composite = np.vstack([mean_composite, mean_data]) if mean_composite.size else mean_data
    min_composite = np.vstack([min_composite, min_data]) if min_composite.size else min_data
    half_composite = np.vstack([half_composite, half_data]) if half_composite.size else half_data
    max_composite = np.vstack([max_composite, max_data]) if max_composite.size else max_data
    sd_composite = np.vstack([sd_composite, sd_data]) if sd_composite.size else sd_data

# Obtain the absolute difference from the finest mesh size
mean_difference_list = np.abs(mean_composite - mean_composite[0])
min_difference_list = np.abs(min_composite - min_composite[0])
max_difference_list = np.abs(max_composite - max_composite[0])
sd_difference_list = np.abs(sd_composite - sd_composite[0])

# Plot the results for all grid sizes
fig = plt.figure(figsize=(20,12), tight_layout = {'pad': 2})
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
plt.plot(grid_spacing_list, mean_difference_list, ls='dashed', label='mean', marker='o')
plt.plot(grid_spacing_list, min_difference_list, ls='dotted', label='min', marker='s')
plt.plot(grid_spacing_list, max_difference_list, ls='dashdot', label='max', marker='^')
plt.plot(grid_spacing_list, sd_difference_list, ls='solid', label='SD', marker='x')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend(loc="upper right", prop={'size': 15})
plt.axvline(10, ls='dotted', c='black')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('absolute difference from simulation with finest mesh size (nM)')
plt.xlabel('grid spacing (μm)')
plt.show()

# Save image
file_path = Path('~/Desktop/grid_size_differences.svg').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
file_path = Path('~/Desktop/grid_size_differences.png').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')

# Prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))
