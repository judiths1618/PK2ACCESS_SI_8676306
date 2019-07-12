import numpy as np
import matplotlib.pyplot as plt

makespan= {'game-based greedy':456.961304, 'NSGA-II':305.309292, 'MOPSO':287.696525, 'DQN-based MARL':127.374440}
total_cost = {'game-based greedy':0.022345547, 'NSGA-II':0.0204265065, 'MOPSO':0.0201514714, 'DQN-based MARL':0.02129999}

m_names = list(makespan.keys())
m_values=list(makespan.values())

c_names = list(total_cost.keys())
c_values=list(total_cost.values())

fig, ax1 = plt.subplots()
width = 0.35
fig.suptitle('Comparision of different algorithms')

color = 'tab:blue'
ax1.set_xlabel('Algorithms')
ax1.set_ylabel('make-span', color=color)
ax1.bar(m_names, m_values, width, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:orange'
ax2.set_ylabel('total cost', color=color)  # we already handled the x-label with ax1
ax2.bar(c_names, c_values, width, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.show()

# fig, ax = plt.subplots()
# rects1 = ax.bar(ind - width/2, makespan_means, width, yerr=makespan_std, label='make-span')
# rects2 = ax.bar(ind + width/2, cost_means, width, yerr=cost_std,
#                 color='orange', label='cost')

# Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Merits')
# ax.set_title('Comparision of different algorithms')
# ax.set_xticks(ind)
# ax.set_xticklabels(('Game-based greedy', 'NSGA-II', 'MOPSO','DQN-based MARL'))
# ax.legend()


# def autolabel(rects, xpos='center'):
#     """
#     Attach a text label above each bar in *rects*, displaying its height.
#
#     *xpos* indicates which side to place the text w.r.t. the center of
#     the bar. It can be one of the following {'center', 'right', 'left'}.
#     """

    # xpos = xpos.lower()  # normalize the case of the parameter
    # ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    # offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off
    #
    # for rect in rects:
    #     height = rect.get_height()
    #     ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
    #             '{}'.format(height), ha=ha[xpos], va='bottom')
#
# autolabel(rects1, "left")
# autolabel(rects2, "right")
#
# plt.savefig("C:\\Users\\Judiths\\Desktop\\Access图表\\Fig6\\"+'com_algor.svg')
# plt.show()