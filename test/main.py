from workflows import singlewl, subwl
from method import *


if __name__=='__main__':
    pos = 'D://PycharmProjects/plan-b/workflows/'
    set = ['CyberShake_30.xml', 'Epigenomics_24.xml', 'Inspiral_30.xml',
           'Montage_25.xml', 'Sipht_29.xml']

    w1 = singlewl.SingleWorkflow(pos+set[0])
    w2 = singlewl.SingleWorkflow(pos+set[1])
    w3 = singlewl.SingleWorkflow(pos+set[2])
    w4 = singlewl.SingleWorkflow(pos+set[3])
    w5 = singlewl.SingleWorkflow(pos+set[4])

    tTypes = ['tType1', 'tType2', 'tType3', 'tType4', 'tType5']
    vms = ['m5a.xlarge', 'm5.large', 'm5.2xlarge', 't3.2xlarge', 't2.xlarge', 'm4.large', 't2.2xlarge']

    price = {'m5a.xlarge': 0.172, 'm5.large': 0.096, 'm5.2xlarge': 0.384, 't3.2xlarge': 0.334,
             't2.xlarge': 0.1856, 'm4.large': 0.100, 't2.2xlarge': 0.3712}


    # -------------拆分workflow------------------
    nodes1 = w1.tasks ()
    adjacent1 = w1.adjacent_matrix (nodes1)
    bags1 = subwl.subWorkflow (w1, nodes1, adjacent1).partition()

    nodes2 = w2.tasks ()
    adjacent2 = w2.adjacent_matrix (nodes2)
    bags2 = subwl.subWorkflow (w2, nodes2, adjacent2).partition ()

    nodes3 = w3.tasks ()
    adjacent3 = w3.adjacent_matrix (nodes3)
    bags3 = subwl.subWorkflow (w3, nodes3, adjacent3).partition ()

    nodes4 = w4.tasks ()
    adjacent4 = w4.adjacent_matrix (nodes4)
    bags4 = subwl.subWorkflow (w4, nodes4, adjacent4).partition ()

    nodes5 = w5.tasks ()
    adjacent5 = w5.adjacent_matrix (nodes5)
    bags5 = subwl.subWorkflow (w5, nodes5, adjacent5).partition ()


    bags=bags1+bags2+bags3+bags4+bags5
    counter = 0
    MAKESPAN_COLL=[[] for i in range(7)]
    COST_COLL=[[] for i in range(7)]

    for k, bag in enumerate(bags):
        g = game.Game(bag)
        e = g.e_table()
        p = g.payoff_table(e)
        # print(p)
        make_span_coll, cost_coll = g.sequential()
        for i in range(7):
            MAKESPAN_COLL[i]+=make_span_coll[i]
            COST_COLL[i]+=cost_coll[i]
    # print(MAKESPAN_COLL, COST_COLL)
    make_span = max(sum (MAKESPAN_COLL[i]) for i in range (len (vms)))
    cost = sum(sum(COST_COLL[i]) for i in range (len (vms)))

    print('Make-span', round(make_span, 8))
    print('Cost', round(cost, 8))

