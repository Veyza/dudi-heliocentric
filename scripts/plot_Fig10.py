import numpy as np
import matplotlib.pyplot as plt
import os

def fR(R):
    Rmin = 1e-7
    Rmax = 1e-4
    q = 3.7
    return (1 - q) / (Rmax**(1 - q) - Rmin**(1 - q)) * R**(-q)

def show_image(d, imname, xcenter=0.5, ycenter=0.5, scalex=5, scaley=5):
    cols = plt.cm.Blues
    pal = cols(np.linspace(0, 1, 50))
    
    # Rotate the image 90 degrees counterclockwise
    dlog = np.rot90(np.where(d > 0, np.log10(d), np.log10(np.min(d[d > 0]) * 1e-1)))
    
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')
    
    # Display image
    im = ax.imshow(dlog, cmap='Blues_r', extent=[-1000, 1000, -1000, 1000], vmin=-12, vmax=-2)
    
    # Define contour levels
    levs = [-12, -10, -9, -8, -7, -6.5, -6, -5.5, -5, -4.5, 
            -4, -3.5, -2, -1]
    
    # Plot contour lines on top of the image, flipping the data horizontally to align
    d_flipped = np.fliplr(np.rot90(dlog, k=2))
    contour1 = ax.contour(d_flipped, levels=[lvl for lvl in levs if lvl > -7], 
                          extent=[-1000, 1000, -1000, 1000], colors='navy', linewidths=1.5, linestyles='solid')
    contour2 = ax.contour(d_flipped, levels=[lvl for lvl in levs if lvl <= -7], 
                          extent=[-1000, 1000, -1000, 1000], colors='azure', linewidths=1.5, linestyles='solid')
    
    # Add contour labels
    ax.clabel(contour1, inline=True, fontsize=14)
    ax.clabel(contour2, inline=True, fontsize=14)

    # Add the color bar with increased font size for the label
    cbar = fig.colorbar(im, ax=ax, orientation='vertical')
    cbar.set_label(r'$log_{10}$ of number density', fontsize=16)
    cbar.ax.tick_params(labelsize=12)
    
    # Set axis limits and labels with increased font size
    ax.set_xlabel("km from asteroid center", fontsize=16)
    ax.set_ylabel("km from asteroid center", fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    plt.xlim(-1000, 1000)
    plt.ylim(-1000, 1000)
    
    # Draw arrow and text above the plot area
    ax.annotate('Sun', xy=(0.05, 1.05), xycoords='axes fraction', color='orange', fontsize=14,
                ha='center', annotation_clip=False)
    ax.annotate('', xy=(0.05, 1.03), xycoords='axes fraction', xytext=(0.2, 1.03),
                arrowprops=dict(facecolor='orange', edgecolor='orange', arrowstyle="->"),
                annotation_clip=False)

    plt.savefig(imname, dpi=300)
    plt.close()

def integrate_over_beta():
    foldername = "./results"
    if not os.path.exists(foldername):
        os.makedirs(foldername)

    Rgs = np.array([1e-7, 2e-7, 3e-7, 4.2e-7, 5.5e-7, 6.7e-7, 8.5e-7, 1e-6,
                    1.2e-6, 2.5e-6, 4e-6, 6e-6, 1e-5, 99e-6])
    nrg = len(Rgs)
    imname = f"{foldername}/phaethon_1AU.png"

    # Initial file for integration
    keypart = f"{foldername}/Rg={Rgs[-1] * 1e6:.2f}micron"
    fname = f"{keypart}.dat"
    d = np.loadtxt(fname)
    mm, nn = d.shape
    integrated = np.zeros((mm, nn))

    for i in range(nrg - 2, -1, -1):
        d1 = d.copy()
        if Rgs[i] * 1e6 < 10:
            keypart = f"{foldername}/Rg= {Rgs[i] * 1e6:.2f}micron"  # Add space after "=" if less than 10
        else:
            keypart = f"{foldername}/Rg={Rgs[i] * 1e6:.2f}micron"  # No space if 10 or greater
        
        fname = f"{keypart}.dat"
        d = np.loadtxt(fname)
        integrated += 0.5 * (d1 * fR(Rgs[i + 1]) + d * fR(Rgs[i])) * (Rgs[i + 1] - Rgs[i])

    integrated /= (2.9 / 3)**2
    np.savetxt(f"{foldername}/integrated_Rmin={Rgs[0] * 1e6:.2f}.dat", integrated)
    
    # Generate the plot with contours
    show_image(integrated, f"{foldername}/integrated_Rmin={Rgs[0] * 1e6:.2f}.png")

integrate_over_beta()

