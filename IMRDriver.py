# encoding: utf-8
'''
@company:Xidian University
@Author: sylvester.
@Email: cczxsong@163.com
@Software: PyCharm
@File: IMRDriver.py
@Time: 2021/12/10 20:48
@Desc:
'''
__author__ = 'sylvester.'

import json
import operator
import os

import numpy as np
import copy
import time
import networkx as nx
import pandas as pd
from tqdm import tqdm

def generate_matrix(G, mode="in_degree"):
    """
    :param G: 输入的有向图
    :param mode:
            in_degree: 通过入度的方式生成matrix.
            default = “in_degree”
    :return: matrix
    """
    if mode == "in_degree":
        A = np.array(nx.adjacency_matrix(G).todense())
        A[np.abs(A) != 0] = 1
        s = (len(G.nodes), len(G.nodes))
        matrix = np.zeros(s)
        for i, n in enumerate(G.nodes):
            regulators = {n1: data for n1, n2, data in G.in_edges(nbunch=[n], data=True)}
            if not regulators:
                continue
            matrix[:, i] = 1 / len(regulators)  # 传播概率
        matrixres = matrix * A
        return matrixres
    return None


def generate_Mr(G, alpha=1):
    """

    :param G:
    :return: 生成Mr
    """
    Mr = []
    for n in G.nodes:
        # Mr.append(G.nodes[n]["nodeValue"])
        if G.nodes[n]['Type'] != 'miRNA':
            Mr.append(G.nodes[n]["nodeValue"])
        else:
            Mr.append(G.nodes[n]["nodeValue"] * alpha)
    return Mr

def LFA(G,matrix, r, alpha, t):
    n = len(matrix)
    # Mr = [1 for i in range(n)]
    Mr = generate_Mr(G=G, alpha=alpha)  # use node value
    for i_ in tqdm(range(1, n), desc="node"):
        i = n - i_
        for j in range(0, i + 1):
            Mr[r[j]] = Mr[r[j]] + matrix[r[j]][r[i]] * Mr[r[i]]  # i -> j 的
            Mr[r[i]] = (1 - matrix[r[j]][r[i]]) * Mr[r[i]]
    return Mr


def IMRDriver(G,matrix, alpha, max_iteration=10):
    print("alpha:{}  max_iteration:{}".format(alpha, max_iteration))
    t = 0
    r0 = [i for i in range(len(matrix))]
    r = r0
    while (True):
        t = t + 1
        Mr = LFA(G,matrix, r, alpha, t)
        r = np.argsort(-np.array(Mr))
        if operator.eq(list(r0), list(r)) or t == max_iteration:
            break
        r0 = copy.copy(r)
    print(r)
    print('迭代次数 ： {}'.format(t))
    return r

def gold():
    gold_standard = pd.read_csv("Data/Gold/Census_allMon Jun  7 06_21_37 2021.csv", sep=',')
    return gold_standard["Gene Symbol"].tolist()

def TOPk_(res, k=500):
    Gold = gold()
    cnt = 0
    tp = res[:k]
    for n in tp:
        if n in Gold:
            cnt += 1
    return cnt


def top1_n(res, n):
    topklst = []
    for i in range(1, n):
        topklst.append(TOPk_(res, i))
    return topklst

def _run(gml, alpha=0.1):
    mode = ["in_degree"]
    for m in mode:
        G = nx.read_gml(gml)
        nodeslst = [n for n in G.nodes]
        matrix = generate_matrix(G, mode=m)
        r = IMRDriver(G=G, matrix=matrix, alpha=alpha, max_iteration=100)
        res_rank = []
        for x in r:
            res_rank.append(nodeslst[x])
        print(res_rank)
        res = top1_n(res_rank, n=500)
        print(res)
        return res, res_rank

if __name__ == "__main__":

    cancerdir = ["BRCA"]
    # cancerdir = ["BRCA", "LUAD", "PRAD", "BLCA", "HNSC", "LUSC", "SKCM", "UCEC"]
    ppi = "Default"
    m = "PCA"
    for c in cancerdir:
        # if not os.path.exists("new_res_rank/{}".format(c)):
        #    os.makedirs("new_res_rank/{}".format(c))

        gml = "Data/{}/newgraph/{}_miRNA_mRNA_TF_v3--others_regnetwork--ppi_from{}.gml".format(c, m, ppi)

        res, res_rank = _run(gml, alpha=1)
