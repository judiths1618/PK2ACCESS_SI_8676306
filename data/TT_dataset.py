import pandas as pd
from workflows import singlewl
import xlsxwriter
import numpy as np


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

NSGA_II_TASKS=[]
for i,ns in enumerate(N):
    for j in range(30):
        if j<=len(ns)-1:
            NSGA_II_TASKS.append(ns[j])
        elif j>len(ns)-1:
            NSGA_II_TASKS.append({'id':'default','name':'null','namespace':'null'})
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

# speed up the data search
# Below code can also be  written "pt = pt_tmp.as_matrix().tolist()"
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

print(time_scp.tolist())
print(cost_scp.tolist())
print('\n')
print(time_mcp.tolist())
print(cost_mcp.tolist())



#
# def createPT(pt):
#     for i,w in enumerate(W):
#         for q,ns in enumerate(N):
#             for k, t in enumerate(ns):
#                 l=tTypes.index(t['name'])
#                 # print(k, t['id'], t['name'], tTypes.index(t['name']), t['namespace'], workflows.index(t['namespace']))
#                 # print(i,k,l, ts[i][k], ms[i][k])
#                 print(pt[i][k])
#                 j = ms[i][k]-1
#                 pt[i][k]=round(tus[j][l]*ts[i][k],6)
#                 # print(i,k,l,j, vms[j], tus[j][l], pt[i][k])
#     print(pt)
#     processing_time=np.array(pt)
#     input = {}
#     for i in range (5):
#         for j in range (30):
#             input['T' + str (j)] = processing_time[:, j]
#     print (input)
#     df = pd.DataFrame (input)
#     writer = pd.ExcelWriter ("WSP_dataset2.xlsx", engine='xlsxwriter')
#     df.to_excel (writer, sheet_name="Processing Time-x")
#
#     writer.save ()
#     return pt


