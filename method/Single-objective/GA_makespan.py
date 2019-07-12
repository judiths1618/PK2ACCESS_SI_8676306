# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 17:24:51 2018

Author: cheng-man wu
LinkedIn: www.linkedin.com/in/chengmanwu
Github: https://github.com/wurmen

"""
import matplotlib

'''==========Solving job shop scheduling problem by gentic algorithm in python======='''
# importing required modules
import pandas as pd
import numpy as np
import time
import copy
from workflows import singlewl
import random
''' ================= initialization setting ======================'''

num_workflow = 5  # number of jobs
num_mc = 7  # number of machines
num_task=30
dataset_pos = "D://PycharmProjects//plan-b//data//"

pos = 'D://PycharmProjects/plan-b/workflows/'
files = ['CyberShake_30.xml', 'Epigenomics_24.xml', 'Inspiral_30.xml',
           'Montage_25.xml', 'Sipht_29.xml']

w1 = singlewl.SingleWorkflow(pos+files[0])
w2 = singlewl.SingleWorkflow(pos+files[1])
w3 = singlewl.SingleWorkflow(pos+files[2])
w4 = singlewl.SingleWorkflow(pos+files[3])
w5 = singlewl.SingleWorkflow(pos+files[4])
W=[w1,w2,w3,w4,w5]

nodes1 = w1.tasks ()
nodes2 = w2.tasks ()
nodes3 = w3.tasks ()
nodes4 = w4.tasks ()
nodes5 = w5.tasks ()
N=[nodes1,nodes2,nodes3,nodes4,nodes5]

TASKS=[]
for i,ns in enumerate(N):
    for j in range(30):
        if j<=len(ns)-1:
            TASKS.append(ns[j])
        elif j>len(ns)-1:
            TASKS.append({'id':'default','name':'null','namespace':'null'})
# print(NSGA_II_TASKS)
# print(len(NSGA_II_TASKS))

workflows=["CyberShake","Genome","LIGO","Montage","SIPHT"]
tTypes=['tType1', 'tType2', 'tType3', 'tType4', 'tType5']
vms = ['t3.medium','t3.large','c5.large','m5.large','c5n.large','r5a.large','a1.4xlarge']


# pt_tmp = pd.read_excel (dataset_pos+"WSP_dataset.xlsx", sheet_name="Processing Time", index_col=[0])
ms_tmp = pd.read_excel (dataset_pos+"WSP_dataset.xlsx", sheet_name="Machines Sequence", index_col=[0])
mp_tmp =pd.read_excel(dataset_pos+"WSP_dataset.xlsx", sheet_name="Machines Price", index_col=[0])
ts_tmp = pd.read_excel (dataset_pos+"WSP_dataset.xlsx", sheet_name="Task Size", index_col=[0])
scp_tmp = pd.read_excel (dataset_pos+"WSP_dataset.xlsx", sheet_name="Single-core Performance", index_col=[0])
mcp_tmp = pd.read_excel (dataset_pos+"WSP_dataset.xlsx", sheet_name="Multi-core Performance", index_col=[0])

num_gene = num_workflow * num_task  # number of genes in a chromosome

# raw_input is used in python 2
population_size = int (input ('Please input the size of population: ') or 20)  # default value is 20
crossover_rate = float (input ('Please input the size of Crossover Rate: ') or 0.8)  # default value is 0.8
mutation_rate = float (input ('Please input the size of Mutation Rate: ') or 0.3)  # default value is 0.3
mutation_selection_rate = float (input ('Please input the mutation selection rate: ') or 0.4)
num_mutation_jobs = round (num_workflow * num_mc * mutation_selection_rate)
num_iteration = int (input ('Please input number of iteration: ') or 500)  # default value is 1000

# speed up the data search
# Below code can also be  written "pt = pt_tmp.as_matrix().tolist()"
# pt = [list (map (float, pt_tmp.iloc[i])) for i in range (num_workflow)]
ms = [list (map (int, ms_tmp.iloc[i])) for i in range (num_workflow)]
price = [list (map (float, mp_tmp.iloc[i])) for i in range (num_mc)]
ts = [list (map (float, ts_tmp.iloc[i])) for i in range (num_workflow)]
scp = [list (map (float, scp_tmp.iloc[i])) for i in range (num_mc)]
mcp = [list (map (float, mcp_tmp.iloc[i])) for i in range (num_mc)]

time_scp = np.zeros((7,5))
cost_scp = np.zeros((7,5))
time_mcp = np.zeros((7,5))
cost_mcp = np.zeros((7,5))
for j,vm in enumerate(vms):
    for i,t in enumerate(N):
        time_scp[j][i]=round(20000/scp[j][i],6)
        cost_scp[j][i]=round(20000/scp[j][i]*price[j][0],6)
        time_mcp[j][i]=round(30000/mcp[j][i],6)
        cost_mcp[j][i]=round(30000/mcp[j][i]*price[j][0],6)
start_time = time.time ()


'''==================== main code ==============================='''
'''----- generate initial population -----'''
Tbest = 999999999999999
best_list, best_obj = [], []
population_list = []
makespan_record = []
vm_record = []
# cost_record=[]
for i in range (population_size):
    nxm_random_num = list (np.random.permutation (num_gene))  # generate a random permutation of 0 to num_job*num_mc-1
    population_list.append (nxm_random_num)  # add to the population_list
    for j in range (num_gene):
        population_list[i][j] = population_list[i][j] % num_workflow  # convert to job number format, every job appears m times

for n in range (num_iteration):
    Tbest_now = 99999999999

    '''-------- two point crossover --------'''
    parent_list = copy.deepcopy (population_list)
    offspring_list = copy.deepcopy (population_list)
    S = list (np.random.permutation (population_size))  # generate a random sequence to select the parent chromosome to crossover

    for m in range (int (population_size / 2)):
        crossover_prob = np.random.rand ()
        if crossover_rate >= crossover_prob:
            parent_1 = population_list[S[2 * m]][:]
            parent_2 = population_list[S[2 * m + 1]][:]
            child_1 = parent_1[:]
            child_2 = parent_2[:]
            cutpoint = list (np.random.choice (num_gene, 2, replace=False))
            cutpoint.sort ()

            child_1[cutpoint[0]:cutpoint[1]] = parent_2[cutpoint[0]:cutpoint[1]]
            child_2[cutpoint[0]:cutpoint[1]] = parent_1[cutpoint[0]:cutpoint[1]]
            offspring_list[S[2 * m]] = child_1[:]
            offspring_list[S[2 * m + 1]] = child_2[:]

    '''----------repairment-------------'''
    for m in range (population_size):
        # print('lala')
        job_count = {}
        larger, less = [], []  # 'larger' record jobs appear in the chromosome more than m times, and 'less' records less than m times.
        for i in range (num_workflow):
            if i in offspring_list[m]:
                count = offspring_list[m].count (i)
                pos = offspring_list[m].index (i)
                job_count[i] = [count, pos]  # store the above two values to the job_count dictionary
            else:
                count = 0
                job_count[i] = [count, 0]
            if count > num_task:
                larger.append (i)
            elif count < num_task:
                less.append (i)

        for k in range (len (larger)):
            chg_job = larger[k]
            while job_count[chg_job][0] > num_task:
                for d in range (len (less)):
                    if job_count[less[d]][0] < num_task:
                        offspring_list[m][job_count[chg_job][1]] = less[d]
                        job_count[chg_job][1] = offspring_list[m].index (chg_job)
                        job_count[chg_job][0] = job_count[chg_job][0] - 1
                        job_count[less[d]][0] = job_count[less[d]][0] + 1
                    if job_count[chg_job][0] == num_task:
                        break

    '''--------mutation--------'''
    for m in range (len (offspring_list)):
        mutation_prob = np.random.rand ()
        if mutation_rate >= mutation_prob:
            m_chg = list (
                np.random.choice (num_gene, num_mutation_jobs, replace=False))  # chooses the position to mutation
            t_value_last = offspring_list[m][m_chg[0]]  # save the value which is on the first mutation position
            for i in range (num_mutation_jobs - 1):
                offspring_list[m][m_chg[i]] = offspring_list[m][m_chg[i + 1]]  # displacement

            offspring_list[m][m_chg[
                num_mutation_jobs - 1]] = t_value_last  # move the value of the first mutation position to the last mutation position

    '''--------fitness value(calculate makespan)-------------'''
    total_chromosome = copy.deepcopy (parent_list) + copy.deepcopy (offspring_list)  # parent and offspring chromosomes combination
    chrom_fitness, chrom_fit = [], []
    total_fitness = 0
    for m in range (population_size * 2):
        j_keys = [j for j in range (num_workflow)]
        key_count = {key: 0 for key in j_keys}
        j_count = {key: 0 for key in j_keys}
        m_keys = [j + 1 for j in range (num_mc)]
        m_count = {key: 0 for key in m_keys}

        for index, i in enumerate(total_chromosome[m]):
            # print(index, TASKS[index])
            flag = TASKS[index]['id']
            # gen_m = int (ms[i][key_count[i]])
            gen_m=random.randint(1,7)

            if flag!='default':
                task_vm_i=workflows.index(TASKS[index]['namespace'])
                task_type_i=tTypes.index(TASKS[index]['name'])
                # single core
                # gen_t = float (ts[i][key_count[i]] * time_scp[gen_m - 1][task_type_i])
                # # multi-core
                gen_t = float (ts[i][key_count[i]] * time_mcp[gen_m - 1][task_type_i])
                # print(gen_t, gen_m, price[gen_m-1][0])
                # j_count[i] = j_count[i] + gen_t* price[gen_m-1][0]    # cost
                j_count[i] = j_count[i] + gen_t         # makespan
                m_count[gen_m] = m_count[gen_m] + gen_t

                if m_count[gen_m] < j_count[i]:
                    m_count[gen_m] = j_count[i]
                elif m_count[gen_m] > j_count[i]:
                    j_count[i] = m_count[gen_m]

            key_count[i] = key_count[i] + 1

        makespan = max (j_count.values ())
        chrom_fitness.append (1 / makespan)
        chrom_fit.append (makespan)
        total_fitness = total_fitness + chrom_fitness[m]
        # cost = max (j_count.values ())
        # chrom_fitness.append (1 / cost)
        # chrom_fit.append (cost)
        # total_fitness = total_fitness + chrom_fitness[m]

    '''----------selection(roulette wheel approach)----------'''
    pk, qk = [], []

    for i in range (population_size * 2):
        pk.append (chrom_fitness[i] / total_fitness)
    for i in range (population_size * 2):
        cumulative = 0
        for j in range (0, i + 1):
            cumulative = cumulative + pk[j]
        qk.append (cumulative)

    selection_rand = [np.random.rand () for i in range (population_size)]

    for i in range (population_size):
        if selection_rand[i] <= qk[0]:
            population_list[i] = copy.deepcopy (total_chromosome[0])
        else:
            for j in range (0, population_size * 2 - 1):
                if selection_rand[i] > qk[j] and selection_rand[i] <= qk[j + 1]:
                    population_list[i] = copy.deepcopy (total_chromosome[j + 1])
                    break
    '''----------comparison----------'''
    for i in range (population_size * 2):
        if chrom_fit[i] < Tbest_now:
            Tbest_now = chrom_fit[i]
            sequence_now = copy.deepcopy (total_chromosome[i])
    if Tbest_now <= Tbest:
        Tbest = Tbest_now
        sequence_best = copy.deepcopy (sequence_now)

    makespan_record.append (Tbest)
'''----------result----------'''
print ("optimal sequence", sequence_best)
print("makespan", makespan_record)
print ("optimal value:%f" % Tbest)
print ('the elapsed time:%s' % (time.time () - start_time))

import matplotlib.pyplot as plt

# %matplotlib inline
plt.title('GA')
plt.plot ([i for i in range (len (makespan_record))], makespan_record, 'b')
plt.ylabel ('makespan', fontsize=15)
plt.xlabel ('generation', fontsize=15)
plt.show ()

# '''--------plot gantt chart-------'''
# import pandas as pd
# import plotly.plotly as py
# import plotly.figure_factory as ff
# import datetime
#
# m_keys = [j + 1 for j in range (num_mc)]
# j_keys = [j for j in range (num_workflow)]
# key_count = {key: 0 for key in j_keys}
# j_count = {key: 0 for key in j_keys}
# m_count = {key: 0 for key in m_keys}
# j_record = []
# for index, i in enumerate(sequence_best):
#     flag=TASKS[index]['id']
#     gen_m = int (ms[i][key_count[i]])
#     # gen_m=random.randint(1,7)
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
#             #     multi-core
#             datetime.timedelta (seconds=j_count[i] - ts[i][key_count[i]] * time_mcp[gen_m - 1][task_type_i]))  # convert seconds to hours, minutes and seconds
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
#                                  Finish='2018-12-22 %s' % (str (i[2][1])), Resource='Workflow %s' % (j + 1)))
# fig = ff.create_gantt (df, index_col='Resource', show_colorbar=True, group_tasks=True, showgrid_x=True,
#                        title='Workflow Schedule')
# py.plot (fig, filename='GA_workflow_scheduling', world_readable=True)
