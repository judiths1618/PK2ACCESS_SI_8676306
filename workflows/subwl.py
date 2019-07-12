import numpy as np

class subWorkflow:

    def __init__(self, init_workflow, nodes, adjacent):
        self.base = init_workflow
        self.nodes = nodes
        self.adjacent = adjacent

    def partition(self):
        bags_of_tasks = []
        flag = True
        cnt=0
        adj = self.adjacent
        while flag:
            task_bag = []
            # 按列求和
            col_sums = np.sum (adj, axis=0)
            for i, col_sum in enumerate(col_sums):
                if col_sum == 0:
                    cnt+=1
                    adj[i,:] = 0
                    adj[:,i] = -2
                    # print(i, col_sum, self.nodes[i])
                    task_bag.append(self.nodes[i])
                    if cnt == adj.shape[1]:
                        flag = False
            # print(task_bag, len(task_bag))
            bags_of_tasks.append(task_bag)
        return bags_of_tasks


