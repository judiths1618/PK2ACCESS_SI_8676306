from xml.etree.ElementTree import ElementTree
from workflows.DAG import Graph
import numpy as np
import re


class SingleWorkflow:
    # '分解单个的scientific workflow；构造新的job类型'
    job_tag = '{http://pegasus.isi.edu/schema/DAX}job'
    child_tag = '{http://pegasus.isi.edu/schema/DAX}child'
    parent_tag = '{http://pegasus.isi.edu/schema/DAX}parent'

    tTypes = ['tType1', 'tType2', 'tType3', 'tType4', 'tType5']
    vms = ['m5a.xlarge', 'm5.large', 'm5.2xlarge', 't3.2xlarge', 't2.xlarge', 'm4.large', 't2.2xlarge']

    def __init__(self, xml_file):
        self.xml_file = xml_file

    def jobs(self):
        # '解析DAG节点'
        tree = ElementTree(file=self.xml_file)
        root = tree.getroot()
        simple_jobs = []
        for elem in root.iter(tag=self.job_tag):
            simple_job = {'id':elem.attrib['id'], 'name':elem.attrib['name'],'namespace':elem.attrib['namespace']}
            simple_jobs.append(simple_job)
        return simple_jobs

    # 解析Job类型
    def types(self):
        types = []
        res = []
        for job in self.jobs():
            types.append(job['name'])
        for i, type in enumerate({}.fromkeys (types).keys()):
            res.append(type)
        return res

    # 将Job原来的类型转化为我们匹配的任务类型
    def types_trans(self):
        indexes = []
        types_inter = []
        for p in range(len(self.types())):
            tmp = p % len(self.tTypes)
            q = {self.types()[p]: self.tTypes[tmp]}
            indexes.append(tmp)
            types_inter.append(q)
        return types_inter

    # 将变身后的任务存起来
    def tasks(self):
        tasks = []
        for item in self.jobs():
            for j in self.types_trans():
                if item['name']==list(j.keys())[0]:
                    item['name'] = j[item['name']]
                    tasks.append(item)
        # print(tasks)
        return tasks

    # 将工作流的拓扑结构转化为有向图的邻接矩阵
    def adjacent_matrix(self, tasks):
        init = np.zeros((len(tasks), len(tasks)))
        tree = ElementTree (file=self.xml_file)
        for child in tree.iter(tag=self.child_tag):
            child_id = child.get('ref')
            for parent in child.findall(self.parent_tag):
                parent_id = parent.get('ref')
                row_i = int(re.findall('\d+', parent_id)[0])
                col_i = int(re.findall('\d+', child_id)[0])
                # print(parent_id, row_i, child_id, col_i)
                # 有向图的邻接矩阵，parent_id -> child_id 设置为1
                init[row_i][col_i] = 1
                # init[col_i][row_i] = 1
        # print(init)
        return init
