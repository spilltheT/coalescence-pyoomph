import matplotlib.pyplot as plt
import numpy as np
import os

# get all files in output directory
output_dir = "lubrication_coalescence/result"
files = sorted([f for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f))])

# open figure
fig,ax = plt.subplots()

# create folder for plots
if not os.path.exists("lubrication_coalescence/_plots"):
    os.makedirs("lubrication_coalescence/_plots")

# create list to store time and h_middle data
time_data = []
h_middle_data = []

# read data from files
for f in files:
    with open(os.path.join(output_dir, f)) as file:        
        # load the data
        data = np.loadtxt(file, delimiter="\t", skiprows=0) 
        x = data[:,0]
        h = data[:,1]

        # load only first line of data
        header = open(os.path.join(output_dir, f)).readline()
        time=float(header.split('@time=')[-1])
        time_data.append(time)
        # find h value at the x value closest to middle
        h_middle = h[np.argmin(np.abs(x-2.5))]
        h_middle_data.append(h_middle)

        # plot the data
        p,=ax.plot(x, h)

        # set the title and labels
        ax.set_title("time: {:.2f} s".format(time))
        ax.set_xlabel("x")
        ax.set_ylabel("h")

        # save the figure
        plt.savefig("lubrication_coalescence/_plots/{}.png".format(f.split(".")[0]), bbox_inches='tight')
        p.remove()

# plot the time vs h_middle data
fig,ax = plt.subplots()
ax.plot(time_data, h_middle_data)
ax.set_xlabel("time")
ax.set_ylabel("$h_{{middle}}$")
plt.savefig("lubrication_coalescence/_plots/time_vs_h_middle.png", bbox_inches='tight')

