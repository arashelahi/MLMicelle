import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
import seaborn as sns

def clusters_2D(data_final):  
    fig , ax = plt.subplots(figsize=(3.0,4))
    font_size = 12
    data_final.plot(
    kind = 'scatter', 
    x = 'x', y = 'y', c = 'clusters', ax = ax, s=10,
    colormap='Set1',
    edgecolor='black',
    colorbar=False       # Disables the color bar
)
    # ax.set(title = 'clustering data')
    plt.xlabel('X(nm)',fontsize=font_size)
    plt.ylabel('Y(nm)',fontsize=font_size)
    plt.xticks([])
    plt.yticks([])
    plt.savefig('./data/output/clusters.png', dpi=600 , bbox_inches = 'tight', pad_inches = 0.01 , transparent=True)
    plt.show()
    # plt.close()

def size_2D(data_final,PO_bead):  
    # Assuming 'atom_name' is a column in data_final DataFrame
    fig , ax = plt.subplots(figsize=(4,4))
    font_size = 12
    colors = data_final['atom_name'].apply(lambda name: 'g' if name == PO_bead else 'r')
    # labels = data_final['atom_name'].apply(lambda name: 'PO' if name == PO_bead else 'EO')
    data_final.plot(
        kind='scatter',
        x='x',
        y='y',
        ax=ax,
        s=10,    # Marker size
        color=colors,
        edgecolor='black',    # Black edge color for all points
        linewidth=0.5      # Edge thickness
    )
    # plt.legend(['PO', 'EO'], fontsize=font_size)
    center_point = np.mean(data_final[data_final['atom_name']==PO_bead].loc[:,['x','y']],axis=0)
    plt.scatter(center_point[0],center_point[1],color='k',s=100)
    micelle_size = np.max(np.linalg.norm(data_final.loc[:,['x','y']]-center_point,axis=1))
    circle = Circle(center_point, micelle_size, color='black', fill=False, linewidth=1.0)  # `fill=False` for just the outline
    circle.set_linestyle('--')  # Set the line style to dashed

    plt.gca().add_patch(circle)  # Add circle to the current axis
    # plt.set(title = 'clustering data')
    plt.xlabel('X(nm)',fontsize=font_size)
    plt.ylabel('Y(nm)',fontsize=font_size)
    plt.xticks([])
    plt.yticks([])
    # plt.tick_params([], labelsize=15)
    plt.savefig('./data/output/micelle_size.png', dpi=300 , bbox_inches = 'tight', pad_inches = 0.01 , transparent=True)
    # plt.show()
    # plt.close()

def rad_dist_plot(rad_mem,clus_mem):
    ax = sns.dax = sns.distplot(rad_mem, bins = 10)
    plt.title('Hydrodynamic radius population population')
    plt.xlabel(r'$R_{H}$',fontsize=20)
    plt.ylabel('Population',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    x_mean=np.average(rad_mem,weights=clus_mem)
    x = ax.lines[0].get_xdata()
    y = ax.lines[0].get_ydata()
    x_range= np.max(x)-np.min(x)
    y_range= np.max(y)-np.min(y)
    x_text=np.max(x)-0.5*x_range;y_text=np.max(y)-0.5*y_range
    plt.text(x_text,y_text, 'mean is= %4.2f' % x_mean,fontsize=13)
    plt.savefig('Hydrodynamic',figsize=(10,10))
    plt.show()
    plt.close()
def agg_num_dist(clus_mem):
    ax = sns.distplot(clus_mem, bins = 9)
    plt.title('Aggregation number frequency')
    plt.xlabel('Aggregation number',fontsize=20)
    plt.ylabel('frequency',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    x_mean=np.average(clus_mem,weights=clus_mem)
    x = ax.lines[0].get_xdata()
    y = ax.lines[0].get_ydata()
    x_range= np.max(x)-np.min(x)
    y_range= np.max(y)-np.min(y)
    x_text=np.max(x)-0.5*x_range;y_text=np.max(y)-0.5*y_range
    plt.text(x_text,y_text, 'mean is= %4d' % x_mean,fontsize=13)
    plt.show()
    plt.savefig('Agg_number',figsize=(10,10))
    plt.close()