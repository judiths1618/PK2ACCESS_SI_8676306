# coding: utf-8
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
from workflows import singlewl
import time
import plotly.plotly as py
import plotly.figure_factory as ff
import datetime

# ----------------------PSO参数设置---------------------------------
class PSO():

    def __init__(self, pN, dim, max_iter):
        self.w = 0.8
        self.c1 = 2
        self.c2 = 2
        self.r1 = 0.6
        self.r2 = 0.3
        self.pN = pN  # 粒子 task
        self.dim = dim  # 搜索维度
        self.max_iter = max_iter  # 迭代次数
        self.X = np.zeros((self.pN, self.dim))  # 所有粒子的位置
        self.V = np.zeros((self.pN, self.dim))  # 所有粒子的速度
        self.pbest = np.zeros((self.pN, self.dim))  # 个体经历的最佳位置和全局最佳位置
        self.gbest = np.zeros((1, self.dim))
        self.vm_best = np.zeros((1, self.dim))  # 所有粒子对应的vm
        self.p_fit = np.zeros(self.pN)  # 每个个体的历史最佳适应值
        self.fit = 1e10  # 全局最佳适应值


    # ---------------------目标函数Sphere函数-----------------------------
    # action 对应一次匹配，(task， vm)，function对应make-span计算
    def function(self, X):
        j_keys = [j for j in range (num_workflow)]  # workflow的index
        key_count = {key: 0 for key in j_keys}      # workflow
        j_count = {key: 0 for key in j_keys}        # workflow
        m_keys = [j + 1 for j in range (num_mc)]    # vm的index
        m_count = {key: 0 for key in m_keys}        # vm的fitness value（字典）

        for index, i in enumerate(X):
            flag = TASKS[index]['id']
            gen_m = random.randint (1, 7)
            # gen_m = int (ms[i][key_count[i]])
            if flag != 'default':
                task_type_i = tTypes.index (TASKS[index]['name'])
                # single core
                # gen_t = float (ts[i][key_count[i]] * time_scp[gen_m - 1][task_type_i])
                # # multi-core
                gen_t = float(time_mcp[gen_m - 1][task_type_i])
                j_count[i] = j_count[i] + gen_t
                m_count[gen_m] = m_count[gen_m] + gen_t

                if m_count[gen_m] < j_count[i]:
                    m_count[gen_m] = j_count[i]
                elif m_count[gen_m] > j_count[i]:
                    j_count[i] = m_count[gen_m]

            key_count[i] = key_count[i] + 1

        makespan = max (j_count.values ())
        return makespan

    # ---------------------初始化种群----------------------------------
    def init_Population(self):
        for i in range(self.pN):
            for j in range(self.dim):
                self.X[i][j] = int(random.randint(0, 4))        # 粒子的位置
                self.V[i][j] = int(random.randint(0, 4))        # 粒子的速度
            self.pbest[i] = self.X[i]                           # 所有粒子的位置
            tmp = self.function(self.X[i])                      # 当前粒子的fitness function值
            self.p_fit[i] = tmp             # 存fitness
            if tmp < self.fit:
                self.fit = tmp
                self.gbest = self.X[i]

                # ----------------------更新粒子位置----------------------------------

    def iterator(self):
        fitness = []
        for t in range(self.max_iter):
            for i in range(self.pN):  # 更新gbest\pbest
                temp = self.function(self.X[i])
                if temp < self.p_fit[i]:  # 更新个体最优
                    self.p_fit[i] = temp
                    self.pbest[i] = self.X[i]
                    if self.p_fit[i] < self.fit:  # 更新全局最优
                        self.gbest = self.X[i]
                        self.fit = self.p_fit[i]
            for i in range(self.pN):
                self.V[i] = self.w * self.V[i] + self.c1 * self.r1 * (self.pbest[i] - self.X[i]) + self.c2 * self.r2 * (self.gbest - self.X[i])
                tmp= self.X[i] + self.V[i]
                self.X[i] = list(int(abs(tmp[i])%5) for i in range(len(tmp)))
            fitness.append(self.fit)
            # print(self.fit)  # 输出最优值
        return fitness

        # ----------------------程序执行-----------------------


if __name__=='__main__':
# ========================== input 参数 ===================================
    num_workflow = 5  # number of workflows
    num_mc = 7  # number of machines
    num_task = 30   # number of tasks in each workflow
    dataset_pos = "D://PycharmProjects//plan-b//data//"

    pos = 'D://PycharmProjects/plan-b/workflows/'
    files = ['CyberShake_30.xml', 'Epigenomics_24.xml', 'Inspiral_30.xml',
             'Montage_25.xml', 'Sipht_29.xml']

    w1 = singlewl.SingleWorkflow (pos + files[0])
    w2 = singlewl.SingleWorkflow (pos + files[1])
    w3 = singlewl.SingleWorkflow (pos + files[2])
    w4 = singlewl.SingleWorkflow (pos + files[3])
    w5 = singlewl.SingleWorkflow (pos + files[4])
    W = [w1, w2, w3, w4, w5]

    nodes1 = w1.tasks ()
    nodes2 = w2.tasks ()
    nodes3 = w3.tasks ()
    nodes4 = w4.tasks ()
    nodes5 = w5.tasks ()
    N = [nodes1, nodes2, nodes3, nodes4, nodes5]

    TASKS = []
    for i, ns in enumerate (N):
        for j in range (30):
            if j <= len (ns) - 1:
                TASKS.append (ns[j])
            elif j > len (ns) - 1:
                TASKS.append ({'id': 'default', 'name': 'null', 'namespace': 'null'})

    workflows = ["CyberShake", "Genome", "LIGO", "Montage", "SIPHT"]
    tTypes = ['tType1', 'tType2', 'tType3', 'tType4', 'tType5']
    vms = ['t3.medium', 't3.large', 'c5.large', 'm5.large', 'c5n.large', 'r5a.large', 'a1.4xlarge']

    # pt_tmp = pd.read_excel (dataset_pos+"WSP_dataset.xlsx", sheet_name="Processing Time", index_col=[0])
    ms_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Machines Sequence", index_col=[0])
    mp_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Machines Price", index_col=[0])
    ts_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Task Size", index_col=[0])
    scp_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Single-core Performance", index_col=[0])
    mcp_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Multi-core Performance", index_col=[0])


    ms = [list (map (int, ms_tmp.iloc[i])) for i in range (num_workflow)]
    price = [list (map (float, mp_tmp.iloc[i])) for i in range (num_mc)]
    ts = [list (map (int, ts_tmp.iloc[i])) for i in range (num_workflow)]
    scp = [list (map (float, scp_tmp.iloc[i])) for i in range (num_mc)]
    mcp = [list (map (float, mcp_tmp.iloc[i])) for i in range (num_mc)]

    time_scp = np.zeros ((7, 5))
    cost_scp = np.zeros ((7, 5))
    time_mcp = np.zeros ((7, 5))
    cost_mcp = np.zeros ((7, 5))
    for j, vm in enumerate (vms):
        for i, t in enumerate (N):
            time_scp[j][i] = round (20000 / scp[j][i], 6)
            cost_scp[j][i] = round (20000 / scp[j][i] * price[j][0], 6)
            time_mcp[j][i] = round (30000 / mcp[j][i], 6)
            cost_mcp[j][i] = round (30000 / mcp[j][i] * price[j][0], 6)


    start_time = time.time ()
    dims = num_workflow*num_task
    my_pso = PSO(pN=20, dim=dims, max_iter=500)
    my_pso.init_Population()
    fitness = my_pso.iterator()
    print(fitness[len(fitness)-1])


# -------------------画图--------------------
    plt.figure(1)
    plt.title("PSO")
    plt.xlabel("iterators", size=14)
    plt.ylabel("makespan", size=14)
    t = np.array([t for t in range(0, 500)])
    fitness = np.array(fitness)
    plt.plot(t, fitness, color='b', linewidth=3)
    plt.show()

# '''--------plot gantt chart-------'''


    #
    # m_keys = [j + 1 for j in range (num_mc)]
    # j_keys = [j for j in range (num_workflow)]
    # key_count = {key: 0 for key in j_keys}
    # j_count = {key: 0 for key in j_keys}
    # m_count = {key: 0 for key in m_keys}
    # j_record = []
    # for index, i in enumerate(best):
    #     flag=TASKS[index]['id']
    #     # gen_m = int (ms[i][key_count[i]])
    #     gen_m=random.randint(1,7)
    #     if flag!='default':
    #         task_vm_i=workflows.index(TASKS[index]['namespace'])
    #         task_type_i=tTypes.index(TASKS[index]['name'])
    #         # single core
    #         # gen_t = float (ts[i][key_count[i]] * time_scp[gen_m - 1][task_type_i])
    #         # # multi-core
    #         gen_t = float (time_mcp[gen_m - 1][task_type_i])
    #         j_count[i] = j_count[i] + gen_t
    #         m_count[gen_m] = m_count[gen_m] + gen_t
    #
    #         if m_count[gen_m] < j_count[i]:
    #             m_count[gen_m] = j_count[i]
    #         elif m_count[gen_m] > j_count[i]:
    #             j_count[i] = m_count[gen_m]
    #
    #         start_time = str (
    #             # single core
    #             # datetime.timedelta (seconds=j_count[i] - ts[i][key_count[i]] * time_scp[gen_m - 1][task_type_i]))  # convert seconds to hours, minutes and seconds
    #             # end_time = str (datetime.timedelta (seconds=j_count[i]))
    #         #        multi-core
    #             datetime.timedelta (seconds=j_count[i] - time_mcp[gen_m - 1][task_type_i]))  # convert seconds to hours, minutes and seconds
    #         end_time = str (datetime.timedelta (seconds=j_count[i]))
    #
    #         j_record.append((i, gen_m, [start_time, end_time]))
    #
    #         # j_record[(i, gen_m)] = [start_time, end_time]
    #
    #     key_count[i] = key_count[i] + 1
    #
    # df = []
    # for m in m_keys:
    #     for j in j_keys:
    #         for i in j_record:
    #             if (m,j)==(i[1], i[0]):
    #                 # df.append (dict (Task='Machine %s' % (m), Start='2018-12-22 %s' % (str (i[2][0])),
    #                 #          Finish='2018-12-22 %s' % (str (i[2][1])), Resource='Workflow %s' % (j + 1)))
    #
    #                 df.append (dict (Task=vms[m-1], Start='2018-12-22 %s' % (str (i[2][0])),
    #                              Finish='2018-12-22 %s' % (str (i[2][1])), Resource='Workflow %s' % (j + 1)))
    # fig = ff.create_gantt (df, index_col='Resource', show_colorbar=True, group_tasks=True, showgrid_x=True,
    #                    title='Workflow Schedule')
    # py.plot (fig, filename='GA_workflow_scheduling', world_readable=True)

