import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 14
matplotlib.rcParams['ytick.labelsize'] = 14

folder_path = 'lubrication_coalescence/profile'

folder = 'VideoZoom'
if not os.path.exists(folder):
    os.makedirs(folder)

files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]

tt = 0
i = 0

for file in files:
    tt = 0.1*i

    file_path = os.path.join(folder_path, file)
    
    data = pd.read_csv(file_path, delimiter='\t', skiprows=1, usecols=[0, 2], header=None)
    
    data.columns = ['x', 'h']
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(data['x'], data['h'], linestyle='-', color='black') 
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$h$', fontsize=16)
    ax.set_title('$t$ = %4.2f' % tt, fontsize=20)
    # ax.axis('equal')  
    
    ax.set_xlim(-0.3, 0.3)
    ax.set_ylim(0, 0.3)
    
    output_file_path = os.path.join(folder, f'{os.path.splitext(file)[0]}.png')
    plt.savefig(output_file_path)
    plt.close()

    i = i + 1

print("All files processed and plots saved.")
