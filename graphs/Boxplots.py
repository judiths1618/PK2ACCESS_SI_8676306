import matplotlib.pyplot as plt
import numpy as np


# Random test data
np.random.seed(19680801)
makespan_all_data = [np.random.normal(std*(1+std)*.1, std, size=100) for std in range(1, 6)]
cost_all_data = [np.random.normal(std*(3+std)*.25, std, size=100) for std in range(1, 6)]

# sample distribution from the make-span and cost with five segments
labels = ['100', '200', '300', '400', '500']

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))

# rectangular box plot
bplot1 = axes[0].boxplot(makespan_all_data,
                         vert=True,  # vertical box alignment
                         patch_artist=True,  # fill with color
                         labels=labels)  # will be used to label x-ticks
axes[0].set_title('Make-span box plot')

# notch shape box plot
bplot2 = axes[1].boxplot(cost_all_data,
                         notch=True,  # notch shape
                         vert=True,  # vertical box alignment
                         patch_artist=True,  # fill with color
                         labels=labels)  # will be used to label x-ticks
axes[1].set_title('Cost box plot')

# fill with colors
colors = ['pink', 'lightblue', 'lightgreen','lightskyblue','yellowgreen']
for bplot in (bplot1, bplot2):
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

# adding horizontal grid lines
for ax in axes:
    ax.yaxis.grid(True)
    ax.set_xlabel('Five separate samples')
    ax.set_ylabel('Observed values')

plt.savefig('C:\\Users\\Judiths\\Desktop\\Access图表\\'+'boxplot.png')
plt.show()



