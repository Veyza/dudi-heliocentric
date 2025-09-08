import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LogNorm

def read_matrix_from_file(filename):
    """Reads a matrix from a .dat file.
    
    Args:
    filename (str): The path to the .dat file containing the matrix data.
    
    Returns:
    numpy.ndarray: The matrix read from the file.
    """
    return np.loadtxt(filename)  # Reading without specifying a delimiter

def visualize_matrix(matrix, filename):
    """Visualizes a matrix using a logarithmic color map where higher values are darker and zero is white.
    The matrix is transposed and rotated 90 degrees clockwise.
    Adjusts the axes to range from -200 to 200, labeled as 'km from the asteroid center'.
    Uses a LaTeX expression in the title for the number density. Also saves the image as a PNG.
    Adds an orange arrow labeled 'Sun'.
    
    Args:
    matrix (numpy.ndarray): The matrix to visualize.
    filename (str): Filename used to save the image.
    """
    plt.figure(figsize=(10, 10))
    ax = plt.gca()
    extent = [-200, 200, -200, 200]  # Setting the extent of the axes to cover -200 to 200 km

    # Rotating the matrix 90 degrees clockwise by transposing and then flipping vertically
    rotated_matrix = np.flip(matrix.T, axis=0)

    # Masking zeros for logarithmic scale
    masked_matrix = np.ma.masked_where(rotated_matrix == 0, rotated_matrix)

    # Visualizing the matrix with logarithmic scaling and specified settings
    im = plt.imshow(masked_matrix, cmap='Blues', norm=LogNorm(), interpolation='nearest', extent=extent, origin='lower')
    cbar = plt.colorbar(im)
    cbar.set_label('Number density, $m^{-3}$', size=18)

    plt.title('Number density, $m^{-3}$', fontsize=21)
    plt.xlabel('km from the asteroid center', fontsize=18)
    plt.ylabel('km from the asteroid center', fontsize=18)

    # Drawing an orange arrow from (-75, -155) to (-170, -155) and keeping the label "Sun" at previous location
    plt.arrow(-85, -155, -85, 0, head_width=8, head_length=10, fc='orange', ec='orange')
    plt.text(-170, -180, 'Sun', color='orange', fontsize=21)
    plt.arrow(-85, -155, 27, 75, head_width=8, head_length=10, fc='grey', ec='grey')
    plt.text(-110, -80, 'apex', color='grey', fontsize=21)

    # Save the figure to a file based on the input filename with .png extension
    output_filename = os.path.splitext(filename)[0] + '.png'
    plt.savefig(output_filename)


filename = './results/result.dat'
matrix = read_matrix_from_file(filename)
visualize_matrix(matrix, filename)
