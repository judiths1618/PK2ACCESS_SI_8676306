import numpy as np
import pandas as pd
from workflows import singlewl, subwl


class Game:
    num_vms=7

    player1 = 'make_span'
    player2 = 'cost'
    makespan_collective = [[] for  i in range (num_vms)]     # 记录每台虚拟机上的time
    cost_collective = [[] for j in range (num_vms)]          # 记录每台虚拟机上的cost
    vms_fini = [0,0,0,0,0,0,0]  # 初始化每台虚拟机的最大完成时间
    vms_cost = [0,0,0,0,0,0,0]  # 初始化每台虚拟机的花销
    make_span = max(vms_fini)       # 记录所有虚拟机的最大完成时间
    total_cost = sum(vms_cost)      # 记录所有虚拟机的总花销
    players = {player1: make_span, player2: total_cost}
    strategy = []         # 部署结果

    def __init__(self, bag_of_tasks, tTypes, vms, time_matrix, cost_matrix, price):
        self.bag_of_tasks = bag_of_tasks        # 已经分好的task包
        self.tTypes = tTypes                    # task的任务类型集合
        self.vms = vms                          # 当前的VMs
        self.time_matrix = time_matrix          # vm-taskType 执行时间矩阵
        self.cost_matrix = cost_matrix          # vm-taskType 金钱成本矩阵
        self.price = price                      # 虚拟机的单位时间价格

    '''===（x,y）第一个元素x的选择，锁定值最小的行；第二个元素y的选择，锁定值最小的列 ==='''
    def choose_action(self, i):
        for index, task in enumerate(self.bag_of_tasks):
            task_type_i = tTypes.index (task['name'])
            w_i = workflows.index(task['namespace'])
            col_t = self.time_matrix[:, task_type_i]
            col_c = self.cost_matrix[:, task_type_i]

            if i%2 == 0:
                # 第偶数次make-span做出选择: 取列
                # print ('列', col_t)
                base_t = self.vms_fini
                tmp_new_vms_fini = np.array (col_t) + np.array (base_t)
                mini_index = tmp_new_vms_fini.argmin(axis=0)
                v_i, t_i, runTime, cost = mini_index, task_type_i, self.time_matrix[mini_index][task_type_i], self.cost_matrix[mini_index][task_type_i]

                item = (v_i, t_i, w_i, runTime, cost)
                self.makespan_collective[v_i].append(runTime)
                self.cost_collective[v_i].append(cost)
                self.vms_fini = [sum (self.makespan_collective[i]) for i in range (self.num_vms)]
                self.vms_cost = [sum (self.cost_collective[i]) for i in range (self.num_vms)]
                self.make_span = max(self.vms_fini)
                self.total_cost = sum(self.vms_cost)
                print('0 new make-span', self.make_span)
                return item, task
            if i % 2 == 1:
                # 第奇数次“礼让”cost做出选择
                base_t = self.vms_cost
                tmp_new_vms_fini = np.array (col_c) + np.array (base_t)
                mini_index = tmp_new_vms_fini.argmin (axis=0)
                v_i, t_i, runTime, cost = mini_index, task_type_i, self.time_matrix[mini_index][task_type_i], self.cost_matrix[mini_index][
                    task_type_i]
                item = (v_i, t_i, w_i, runTime, cost)
                self.makespan_collective[v_i].append (runTime)
                self.cost_collective[v_i].append (cost)
                self.vms_fini = [sum (self.makespan_collective[i]) for i in range (self.num_vms)]
                self.vms_cost = [sum (self.cost_collective[i]) for i in range (self.num_vms)]
                self.make_span = max (self.vms_fini)
                self.total_cost = sum (self.vms_cost)
                print ('1 new make-span', self.make_span)
                return item, task

    def sequential(self):
        strategy = []
        for i in range(len(self.bag_of_tasks)):
            is_done = self.choose_action(i)
            strategy.append(is_done[0])
            self.bag_of_tasks.remove(is_done[1])
        print(len(strategy))
        print ('current make-span, cost', self.make_span, self.total_cost)
        return self.makespan_collective, self.cost_collective, strategy


if __name__=='__main__':

    num_workflow = 5  # number of jobs
    num_vm = 7  # number of machines
    tTypes = ['tType1', 'tType2', 'tType3', 'tType4', 'tType5']
    vms = ['t3.medium', 't3.large', 'c5.large', 'm5.large', 'c5n.large', 'r5a.large', 'a1.4xlarge']
    workflows = ["CyberShake", "Genome", "LIGO", "Montage", "SIPHT"]

    pos = 'D://PycharmProjects/plan-b/workflows/'
    files = ['CyberShake_30.xml', 'CyberShake_50.xml', 'CyberShake_100.xml', 'CyberShake_1000.xml',
             'Epigenomics_24.xml', 'Epigenomics_46.xml', 'Epigenomics_100.xml', 'Epigenomics_997.xml',
             'Inspiral_30.xml', 'Inspiral_50.xml', 'Inspiral_100.xml', 'Inspiral_1000.xml',
             'Montage_25.xml', 'Montage_50.xml', 'Montage_100.xml', 'Montage_1000.xml',
             'Sipht_29.xml', 'Sipht_60.xml', 'Sipht_100.xml', 'Sipht_1000.xml']

    w1 = singlewl.SingleWorkflow(pos + files[2])
    w2 = singlewl.SingleWorkflow(pos + files[6])
    w3 = singlewl.SingleWorkflow(pos + files[10])
    w4 = singlewl.SingleWorkflow(pos + files[14])
    w5 = singlewl.SingleWorkflow(pos + files[18])
    W = [w1, w2, w3, w4, w5]

    nodes1 = w1.tasks ()
    nodes2 = w2.tasks ()
    nodes3 = w3.tasks ()
    nodes4 = w4.tasks ()
    nodes5 = w5.tasks ()
    N = [nodes1, nodes2, nodes3, nodes4, nodes5]
    print(sum([len(i) for i in N]))

    dataset_pos = "D://PycharmProjects//plan-b//data//"
    mp_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Machines Price", index_col=[0])
    ts_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Task Size", index_col=[0])
    scp_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Single-core Performance", index_col=[0])
    mcp_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Multi-core Performance", index_col=[0])

    # speed up the data search
    # Below code can also be  written "pt = pt_tmp.as_matrix().tolist()"
    price = [list (map (float, mp_tmp.iloc[i])) for i in range (num_vm)]
    scp = [list (map (float, scp_tmp.iloc[i])) for i in range (num_vm)]
    mcp = [list (map (float, mcp_tmp.iloc[i])) for i in range (num_vm)]

    time_scp = np.zeros ((7, 5))        # 单核的时间矩阵
    cost_scp = np.zeros ((7, 5))        # 单核的时间对应的价格矩阵
    time_mcp = np.zeros ((7, 5))        # 多核的时间矩阵
    cost_mcp = np.zeros ((7, 5))        # 多核的时间对应的价格矩阵
    for j, vm in enumerate (vms):
        for i, t in enumerate (N):
            time_scp[j][i] = round (20000 / scp[j][i], 6)
            cost_scp[j][i] = round (20000 / scp[j][i] * price[j][0], 6)
            time_mcp[j][i] = round (30000 / mcp[j][i], 6)
            cost_mcp[j][i] = round (30000 / mcp[j][i] * price[j][0], 6)


    '''------------- 拆分multiple workflows 为多个task bags ------------------'''
    adjacent1 = w1.adjacent_matrix (nodes1)
    bags1 = subwl.subWorkflow (w1, nodes1, adjacent1).partition()

    adjacent2 = w2.adjacent_matrix (nodes2)
    bags2 = subwl.subWorkflow (w2, nodes2, adjacent2).partition ()

    adjacent3 = w3.adjacent_matrix (nodes3)
    bags3 = subwl.subWorkflow (w3, nodes3, adjacent3).partition ()

    adjacent4 = w4.adjacent_matrix (nodes4)
    bags4 = subwl.subWorkflow (w4, nodes4, adjacent4).partition ()

    adjacent5 = w5.adjacent_matrix (nodes5)
    bags5 = subwl.subWorkflow (w5, nodes5, adjacent5).partition ()

    # print(len(bags1),len(bags2), len(bags3), len(bags4), len(bags5))

    import random
    # 打乱顺序
    m0=[]
    m1=[bags1[0],bags2[0],bags3[0],bags4[0],bags5[0]]
    # random.shuffle(m1)
    m2=[bags1[1],bags2[1],bags3[1],bags4[1],bags5[1]]
    # random.shuffle(m2)
    m3=[         bags2[2],bags3[2],bags4[2],bags5[2]]
    # random.shuffle(m3)
    m4=[bags1[2],bags2[3],bags3[3],bags4[3],bags5[3]]
    # random.shuffle(m4)
    m5=[bags1[3],bags2[4],bags3[4],bags4[4],bags5[4]]
    # random.shuffle(m5)
    m6=[         bags2[5],bags3[5],bags4[5]]
    # random.shuffle(m6)
    m7=[         bags2[6],         bags4[6]]
    # random.shuffle(m7)
    m8=[        bags2[7],         bags4[7],bags4[8]]
    # random.shuffle(m8)
    TS = [[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
           1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
          [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
           1.0, 1.0],
          [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
           1.0, 1.0, 1.0],
          [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
           1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
          [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
           1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]
    ts = TS[0] + TS[1] + TS[2] + TS[3] + TS[4]
    bags= m1+m2+m3+m4+m5+m6+m7+m8


    # sequential game
    strategies = []
    make_span = 0
    total_cost = 0
    for k, bag in enumerate(bags):
        g = Game(bag, tTypes, vms, time_mcp, cost_mcp, price)
        make_span_coll, cost_coll, strategy = g.sequential()
        strategies+=strategy
        if k==len(bags)-1:
            final_makespan_collective = g.makespan_collective
            final_cost_collective = g.cost_collective
            # print(final_makespan_collective, '\n', final_cost_collective)
            make_span = max (sum (final_makespan_collective[i]) for i in range (len (vms)))
            total_cost = sum (sum (final_cost_collective[i]) for i in range (len (vms)))
            # print ('final make-span, cost', round (make_span, 8), round (total_cost/3600, 8))

    print(len(strategies))

    # --------plot gantt chart-------
    import plotly.plotly as py
    import plotly.figure_factory as ff
    import datetime

    m_keys = [j + 1 for j in range (len (vms))]
    j_keys = [j for j in range (len (workflows))]
    key_count = {key: 0 for key in j_keys}
    j_count = {key: 0 for key in j_keys}
    m_count = {key: 0 for key in m_keys}
    j_record = []
    _cost = []

    for index, s in enumerate (strategies):
        # print (s)
        gen_m = s[0] + 1  # 虚拟机编号
        task_type_i = s[1]  # task类型
        i = s[2]  # workflow类型

        # print(gen_m)
        # single core
        # gen_t = float (ts[i]*time_scp[gen_m - 1][task_type_i])
        # multi-core
        gen_t = float (s[3])
        j_count[i] = j_count[i] + gen_t
        m_count[gen_m] = m_count[gen_m] + gen_t

        if m_count[gen_m] < j_count[i]:
            m_count[gen_m] = j_count[i]
        elif m_count[gen_m] > j_count[i]:
            j_count[i] = m_count[gen_m]

        _cost.append(gen_t*price[gen_m-1][0])
        start_time = str (
            # single core
            # datetime.timedelta (seconds=j_count[i] - ts[i][key_count[i]] * time_scp[gen_m - 1][task_type_i]))  # convert seconds to hours, minutes and seconds
            # end_time = str (datetime.timedelta (seconds=j_count[i]))
            # multi-core
            datetime.timedelta (seconds=j_count[i] - ts[i] * time_mcp[gen_m - 1][
                task_type_i]))  # convert seconds to hours, minutes and seconds
        end_time = str (datetime.timedelta (seconds=j_count[i]))

        j_record.append ((i, gen_m, [start_time, end_time]))

        key_count[i] = key_count[i] + 1
    makespan = max (j_count.values ())
    print (len (j_record), sum (_cost)/3600, makespan)

    df = []
    for m in m_keys:
        for j in j_keys:
            for i in j_record:
                if (m, j) == (i[1], i[0]):
                    # df.append (dict (Task='Machine %s' % (m), Start='2018-12-22 %s' % (str (i[2][0])),
                    #          Finish='2018-12-22 %s' % (str (i[2][1])), Resource='Workflow %s' % (j + 1)))

                    df.append (dict (Task=vms[m - 1], Start='2018-12-22 %s' % (str (i[2][0])),
                                     Finish='2018-12-22 %s' % (str (i[2][1])), Resource='Workflow %s' % (j + 1)))
    fig = ff.create_gantt (df, index_col='Resource', show_colorbar=True, group_tasks=True, showgrid_x=True,
                           title='Workflow Schedule')
    py.plot (fig, filename='Greedy_workflow_scheduling', world_readable=True)

