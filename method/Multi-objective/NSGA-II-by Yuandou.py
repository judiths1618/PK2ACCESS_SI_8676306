import pandas as pd
import numpy as np
import time
import copy
from workflows import singlewl
import matplotlib.pyplot as plt
import random


'''===========function==============='''
'''-------non-dominated sorting function-------'''


def non_dominated_sorting(population_size, chroms_obj_record):
    s, n = {}, {}
    front, rank = {}, {}
    front[0] = []
    for p in range (population_size * 2):
        s[p] = []
        n[p] = 0
        for q in range (population_size * 2):

            if ((chroms_obj_record[p][0] < chroms_obj_record[q][0] and chroms_obj_record[p][1] < chroms_obj_record[q][
                1]) or (chroms_obj_record[p][0] <= chroms_obj_record[q][0] and chroms_obj_record[p][1] <
                        chroms_obj_record[q][1])
                    or (chroms_obj_record[p][0] < chroms_obj_record[q][0] and chroms_obj_record[p][1] <=
                        chroms_obj_record[q][1])):
                if q not in s[p]:
                    s[p].append (q)
            elif ((chroms_obj_record[p][0] > chroms_obj_record[q][0] and chroms_obj_record[p][1] > chroms_obj_record[q][
                1]) or (chroms_obj_record[p][0] >= chroms_obj_record[q][0] and chroms_obj_record[p][1] >
                        chroms_obj_record[q][1])
                  or (chroms_obj_record[p][0] > chroms_obj_record[q][0] and chroms_obj_record[p][1] >=
                      chroms_obj_record[q][1])):
                n[p] = n[p] + 1
        if n[p] == 0:
            rank[p] = 0
            if p not in front[0]:
                front[0].append (p)

    i = 0
    while (front[i] != []):
        Q = []
        for p in front[i]:
            for q in s[p]:
                n[q] = n[q] - 1
                if n[q] == 0:
                    rank[q] = i + 1
                    if q not in Q:
                        Q.append (q)
        i = i + 1
        front[i] = Q

    del front[len (front) - 1]
    return front


'''--------calculate crowding distance function---------'''


def calculate_crowding_distance(front, chroms_obj_record):
    distance = {m: 0 for m in front}
    for o in range (2):
        obj = {m: chroms_obj_record[m][o] for m in front}
        sorted_keys = sorted (obj, key=obj.get)
        distance[sorted_keys[0]] = distance[sorted_keys[len (front) - 1]] = 999999999999
        for i in range (1, len (front) - 1):
            if len (set (obj.values ())) == 1:
                distance[sorted_keys[i]] = distance[sorted_keys[i]]
            else:
                distance[sorted_keys[i]] = distance[sorted_keys[i]] + (
                        obj[sorted_keys[i + 1]] - obj[sorted_keys[i - 1]]) / (
                                                   obj[sorted_keys[len (front) - 1]] - obj[sorted_keys[0]])

    return distance


'''----------selection----------'''


def selection(population_size, front, chroms_obj_record, total_chromosome):
    N = 0
    new_pop = []
    while N < population_size:
        for i in range (len (front)):
            N = N + len (front[i])
            if N > population_size:
                distance = calculate_crowding_distance (front[i], chroms_obj_record)
                sorted_cdf = sorted (distance, key=distance.get)
                sorted_cdf.reverse ()
                for j in sorted_cdf:
                    if len (new_pop) == population_size:
                        break
                    new_pop.append (j)
                break
            else:
                new_pop.extend (front[i])

    population_list = []
    for n in new_pop:
        population_list.append (total_chromosome[n])

    return population_list, new_pop


'''基本思路：将138个tasks以及其对应的虚拟机表示一个种群（一种调度策略）'''
if __name__ == '__main__':

    ''' ================= initialization setting ======================'''
    num_workflow = 5  # number of work-flows
    num_mc = 7  # number of virtual machines

    workflows = ["CyberShake", "Genome", "LIGO", "Montage", "SIPHT"]
    tTypes = ['tType1', 'tType2', 'tType3', 'tType4', 'tType5']
    vms = ['t3.medium','t3.large','c5.large','m5.large','c5n.large','r5a.large','a1.4xlarge']

    dataset_pos = "D://PycharmProjects//plan-b//data//"
    pos = 'D://PycharmProjects//plan-b//workflows//'

    files = ['CyberShake_30.xml', 'CyberShake_50.xml', 'CyberShake_100.xml', 'CyberShake_1000.xml',
             'Epigenomics_24.xml', 'Epigenomics_46.xml', 'Epigenomics_100.xml', 'Epigenomics_997.xml',
             'Inspiral_30.xml', 'Inspiral_50.xml', 'Inspiral_100.xml', 'Inspiral_1000.xml',
             'Montage_25.xml', 'Montage_50.xml', 'Montage_100.xml', 'Montage_1000.xml',
             'Sipht_29.xml', 'Sipht_60.xml', 'Sipht_100.xml', 'Sipht_1000.xml']


    # five workflow objects
    w1 = singlewl.SingleWorkflow (pos + files[2])
    w2 = singlewl.SingleWorkflow (pos + files[6])
    w3 = singlewl.SingleWorkflow (pos + files[10])
    w4 = singlewl.SingleWorkflow (pos + files[14])
    w5 = singlewl.SingleWorkflow (pos + files[18])
    W = [w1, w2, w3, w4, w5]

    nodes1 = w1.tasks()
    nodes2 = w2.tasks()
    nodes3 = w3.tasks()
    nodes4 = w4.tasks()
    nodes5 = w5.tasks()
    N = [nodes1, nodes2, nodes3, nodes4, nodes5]

    # tasks to be schedule for NSGA-II, and add to the default values
    NSGA_II_TASKS = []
    for i, ns in enumerate (N):
        for j in range (len(ns)):
            if j <= len (ns) - 1:
                NSGA_II_TASKS.append (ns[j])
            elif j > len (ns) - 1:
                NSGA_II_TASKS.append ({'id': 'default', 'name': 'null', 'namespace': 'null'})

    gene_num = len(NSGA_II_TASKS)
    print(gene_num)
    mp_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Machines Price", index_col=[0])
    ts_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Task Size", index_col=[0])
    scp_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Single-core Performance", index_col=[0])
    mcp_tmp = pd.read_excel (dataset_pos + "WSP_dataset.xlsx", sheet_name="Multi-core Performance", index_col=[0])

    # raw_input is used in python 3
    # population_size = int (input ('Please input the size of population: ') or 100)  # default value is 20
    # crossover_rate = float (input ('Please input the size of Crossover Rate: ') or 0.8)  # default value is 0.8
    # mutation_rate = float (input ('Please input the size of Mutation Rate: ') or 0.3)  # default value is 0.3
    # mutation_selection_rate = float (input ('Please input the mutation selection rate: ') or 0.4)
    # num_mutation_jobs = round (gene_num * mutation_selection_rate)
    # num_iteration = int (input ('Please input number of iteration: ') or 500)  # default value is 1000

    # raw_input is used in python 3
    population_size = 20  # default value is 20
    crossover_rate = 0.8  # default value is 0.8
    mutation_rate = 0.3  # default value is 0.3
    mutation_selection_rate = 0.4
    num_mutation_jobs = round (gene_num * mutation_selection_rate)
    num_iteration = 500  # default value is 1000


    # speed up the data search
    # Below code can also be  written "pt = pt_tmp.as_matrix().tolist()"
    price = [list (map (float, mp_tmp.iloc[i])) for i in range (num_mc)]
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

    '''==================== main code ==============================='''
    '''----- generate initial population -----'''
    best_list, best_obj = [], []    # 最好的list-scheduling order， 最好的obj值
    population_list = []
    best_obj1 = []  # make-span
    best_obj2 = []  # cost
    for i in range (population_size):
        nxm_random_num = list (
            np.random.permutation (gene_num))  # generate a random permutation of 0 to gene_num-1
        population_list.append (nxm_random_num)  # add to the population_list
        for j in range (gene_num):
            population_list[i][j] = population_list[i][j] % num_workflow  # convert to job number format, every job appears m times

    for n in range (num_iteration):
        '-------- two point crossover --------'
        parent_list = copy.deepcopy (population_list)
        offspring_list = []
        S = list (np.random.permutation (
            population_size))  # generate a random sequence to select the parent chromosome to crossover

        for m in range (int (population_size / 2)):
            parent_1 = population_list[S[2 * m]][:]
            parent_2 = population_list[S[2 * m + 1]][:]
            child_1 = parent_1[:]
            child_2 = parent_2[:]

            # 切片+交叉
            cutpoint = list (np.random.choice (gene_num, 2, replace=False))
            cutpoint.sort ()

            child_1[cutpoint[0]:cutpoint[1]] = parent_2[cutpoint[0]:cutpoint[1]]
            child_2[cutpoint[0]:cutpoint[1]] = parent_1[cutpoint[0]:cutpoint[1]]

            offspring_list.extend ((child_1, child_2))  # append child chromosome to offspring list

        '''----------repairment-------------'''
        for m in range (population_size):
            job_count = {}
            larger, less = [], []  # 'larger' record jobs appear in the chromosome more than m times, and 'less' records less than m times.
            for i in range (gene_num//num_workflow):
                if i in offspring_list[m]:
                    count = offspring_list[m].count (i)
                    pos = offspring_list[m].index (i)
                    job_count[i] = [count, pos]  # store the above two values to the job_count dictionary
                else:
                    count = 0
                    job_count[i] = [count, 0]
                if count > (gene_num//num_workflow):
                    larger.append (i)
                elif count < (gene_num//num_workflow):
                    less.append (i)

            for k in range (len (larger)):
                chg_job = larger[k]
                while job_count[chg_job][0] > (gene_num//num_workflow):
                    for d in range (len (less)):
                        if job_count[less[d]][0] < (gene_num//num_workflow):
                            offspring_list[m][job_count[chg_job][1]] = less[d]
                            job_count[chg_job][1] = offspring_list[m].index (chg_job)
                            job_count[chg_job][0] = job_count[chg_job][0] - 1
                            job_count[less[d]][0] = job_count[less[d]][0] + 1
                        if job_count[chg_job][0] == (gene_num//num_workflow):
                            break

        '''--------mutation--------'''
        for m in range (len (offspring_list)):
            mutation_prob = np.random.rand ()
            if mutation_rate <= mutation_prob:
                m_chg = list (np.random.choice (gene_num, num_mutation_jobs,
                                                replace=False))  # chooses the position to mutation
                t_value_last = offspring_list[m][m_chg[0]]  # save the value which is on the first mutation position
                for i in range (num_mutation_jobs - 1):
                    offspring_list[m][m_chg[i]] = offspring_list[m][m_chg[i + 1]]  # displacement

                offspring_list[m][m_chg[
                    num_mutation_jobs - 1]] = t_value_last  # move the value of the first mutation position to the last mutation position

        '''----------- 交叉之后边界检查convert to 0-4 --------------'''
        for m in range (population_size):
            for j in range (gene_num):
                offspring_list[m][j] = offspring_list[m][j] % num_workflow
                parent_list[m][j] = parent_list[m][j] % num_workflow

        '''--------fitness value(calculate  makespan and cost)-------------'''
        total_chromosome = copy.deepcopy (parent_list) + copy.deepcopy (offspring_list)  # combine parent and offspring chromosomes
        chroms_obj_record = {}  # record each chromosome objective values as chromosome_obj_record={chromosome:[cost,makespan, vm_order_list]}
        for m in range (population_size * 2):
            '''==================== make-span ===================='''
            j_keys = [j for j in range (num_workflow)]
            key_count = {key: 0 for key in j_keys}
            j_count = {key: 0 for key in j_keys}
            c_count = {key: 0 for key in j_keys}
            m_keys = [j + 1 for j in range (num_mc)]
            m_count = {key: 0 for key in m_keys}
            m_count_c = {key: 0 for key in m_keys}
            m_vm_order_record = []  # record each

            for index, i in enumerate (total_chromosome[m]):
                task_vm_i = workflows.index (NSGA_II_TASKS[index]['namespace'])
                task_type_i = tTypes.index (NSGA_II_TASKS[index]['name'])
                # gen_m = int(ms[i])      # 虚拟机的order list
                gen_m=int(random.randint(1,7))
                # single core
                # gen_t = float (ts[i] * time_scp[gen_m - 1][task_type_i])
                # # multi-core
                gen_t = float (time_mcp[gen_m - 1][task_type_i])
                gen_c = float (cost_mcp[gen_m - 1][task_type_i])
                j_count[i] = j_count[i] + gen_t         # vm finish time
                c_count[i] = c_count[i] + gen_c         # cost
                m_count[gen_m] = m_count[gen_m] + gen_t
                m_count_c[gen_m] = m_count_c[gen_m] + gen_c
                m_vm_order_record.append (gen_m)
                if m_count[gen_m] < j_count[i]:
                    m_count[gen_m] = j_count[i]
                elif m_count[gen_m] > j_count[i]:
                    j_count[i] = m_count[gen_m]

                key_count[i] = key_count[i] + 1

            print(c_count)

            makespan = max (j_count.values ())
            cost = sum (c_count.values ())
            chroms_obj_record[m] = [makespan, cost, m_vm_order_record]


        '''-------non-dominated sorting-------'''
        front = non_dominated_sorting (population_size, chroms_obj_record)

        '''----------selection----------'''
        population_list, new_pop = selection (population_size, front, chroms_obj_record, total_chromosome)
        new_pop_obj = [chroms_obj_record[k] for k in new_pop]

        '''----------comparison----------'''
        if n == 0:
            best_list = copy.deepcopy (population_list)
            best_obj = copy.deepcopy (new_pop_obj)
        else:
            total_list = copy.deepcopy (population_list) + copy.deepcopy (best_list)
            total_obj = copy.deepcopy (new_pop_obj) + copy.deepcopy (best_obj)

            now_best_front = non_dominated_sorting (population_size, total_obj)
            best_list, best_pop = selection (population_size, now_best_front, total_obj, total_list)
            best_obj = [total_obj[k] for k in best_pop]

        best_obj1.append (best_obj[0][0])
        best_obj2.append (best_obj[0][1]/3600)

    '''----------result----------'''
    # print (best_list)
    print (len (best_list), len (best_list[0]))
    print (best_obj[0][1]/3600, best_obj[0][0])
    print ('the elapsed time:%s' % (time.time () - start_time))


'''
    # %matplotlib inline
    # -------------- convergences & Pareto front
    fig, (ax0, ax1) = plt.subplots (nrows=2)

    ax0.grid(True)
    ax0.set_xlabel('episodes')
    ax0.set_ylabel('Maximum completion time(s)')
    ax0.plot ([i for i in range (len (best_obj1))], best_obj1, label='makespan')
    ax0.legend ()

    ax1.grid (True)
    ax1.set_xlabel ('episodes')
    ax1.set_ylabel ('total cost($)')
    ax1.plot ([i for i in range (len (best_obj2))], best_obj2, 'orange', label='cost')
    ax1.legend ()
    fig.subplots_adjust (hspace=0.3)
    # plt.savefig ("C:\\Users\\Judiths\\Desktop\\Access图表\\Fig4\\" + 'makespan&cost_conv_NSGA-II.svg')
    fig.tight_layout()
    plt.show ()

    plt.plot (best_obj1, best_obj2, 'bo')
    plt.xlabel ('makespan', fontsize=15)
    plt.ylabel ('cost', fontsize=15)
    # plt.savefig ("C:\\Users\\Judiths\\Desktop\\Access图表\\Fig4\\" + 'makespan&cost_NSGA-II.svg')
    # plt.show ()

    # --------plot gantt chart-------
    import plotly.plotly as py
    import plotly.figure_factory as ff
    import datetime
    m_keys = [j + 1 for j in range (num_mc)]
    j_keys = [j for j in range (num_workflow)]
    key_count = {key: 0 for key in j_keys}
    j_count = {key: 0 for key in j_keys}
    m_count = {key: 0 for key in m_keys}
    j_record = []
    _cost = []
    gen_m_list = best_obj[0][2]
    for index, i in enumerate (best_list[0]):
        gen_m = gen_m_list[index]
        task_vm_i = workflows.index (NSGA_II_TASKS[index]['namespace'])
        task_type_i = tTypes.index (NSGA_II_TASKS[index]['name'])
        # single core
        # gen_t = float (ts[i]*time_scp[gen_m - 1][task_type_i])
        # multi-core
        gen_t = float (time_mcp[gen_m - 1][task_type_i])
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
            datetime.timedelta (seconds=j_count[i] - time_mcp[gen_m - 1][
                task_type_i]))  # convert seconds to hours, minutes and seconds
        end_time = str (datetime.timedelta (seconds=j_count[i]))

        j_record.append ((i, gen_m, [start_time, end_time]))

        key_count[i] = key_count[i] + 1
    makespan = max (j_count.values ())
    print(len(j_record), sum(_cost)/3600, makespan )

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
    py.plot (fig, filename='NSGA_II_workflow_scheduling', world_readable=True)
'''