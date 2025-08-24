#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 18:55:40 2023

@author: tezcan.b

Experiments
"""
import src 
from src import Network, Prob
import timeit 
import numpy as np
from itertools import product, repeat
import networkx as nx 
import gurobipy as gp
from gurobipy import GRB
import pandas as pd 

onlyPath = True
NETWORKS = ['grid_3_3_1']
CUTOFFS = [0.5]
#CUTOFFS = [0]
BUDGETS = [1, 2]

result_df = pd.DataFrame()
#scenario clustering analysis 
for n_id in NETWORKS:
    for budget in BUDGETS:
        for cutoff in CUTOFFS:
            print(f'network: {n_id}, budget: {budget}, cutoff: {cutoff}')
            params = {'cutoff': cutoff, 'budget': budget}
            n = Network(n_id)
            n.create_scen_tree_v2(params)
            #iteration_df = n.store_output(params)
            #P = Prob(n, params)
            #m = P.path_based_formulation()
            #iteration_df = P.store_output(params, iteration_df).to_frame()
            #iteration_df = iteration_df.T
            # Append the DataFrame to the result_df
            #result_df = pd.concat([result_df, iteration_df], axis = 0, ignore_index=True)    
            #result_df.columns = iteration_df.columns
            #with pd.ExcelWriter('current_exp.xlsx', engine='openpyxl', mode='a', if_sheet_exists ='overlay') as writer:   
            #    sheet = writer.book['PathExperiments']
            #    iteration_df.to_excel(writer, sheet_name='PathExperiments', index=False, header=False, startrow=sheet.max_row)
            if cutoff == 0 and onlyPath == False:
                iteration_df = n.store_output(params)
                P = Prob(n, params)

                m = P.arc_based_formulation()
                iteration_df = P.store_output(params, iteration_df).to_frame()
                iteration_df = iteration_df.T
                result_df = pd.concat([result_df, iteration_df], axis = 0, ignore_index=True)    
                result_df.columns = iteration_df.columns
                with pd.ExcelWriter('current_exp.xlsx', engine='openpyxl', mode='a', if_sheet_exists = 'overlay') as writer:   
                    sheet = writer.book['PathExperiments_newcomptest']
                    iteration_df.to_excel(writer, sheet_name='PathExperiments_newcomptest', index=False, header=False, startrow=sheet.max_row)            
#P.get_lower_bound()
#print(f'lb: {P.lb}')
#%%Experiments  
# n_id = 'grid_3_3'

# CUTOFFS = [0, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50]
# BUDGETS = [1, 2]
# result_df = pd.DataFrame()

# #scenario clustering analysis 
# for budget in BUDGETS:
#     for cutoff in CUTOFFS:
#         print(f'budget: {budget}, cutoff: {cutoff}')
#         params = {'cutoff': cutoff, 'budget': budget}
#         n = Network(n_id)
#         n.create_scen_tree(params)
#         iteration_df = n.store_output(params)
#         P = Prob(n, params)
#         m = P.path_based_formulation()
#         iteration_df = P.store_output(params, iteration_df).to_frame()
#         iteration_df = iteration_df.T
#         # Append the DataFrame to the result_df
#         result_df = pd.concat([result_df, iteration_df], axis = 0, ignore_index=True)    
#         result_df.columns = iteration_df.columns
#         with pd.ExcelWriter('current_exp.xlsx', engine='openpyxl', mode='a', if_sheet_exists ='overlay') as writer:   
#             sheet = writer.book['PathExperiments']
#             iteration_df.to_excel(writer, sheet_name='PathExperiments', index=False, header=False, startrow=sheet.max_row)


#         if cutoff == 0:
#             iteration_df = n.store_output(params)
#             P = Prob(n, params)
#             m = P.arc_based_formulation()
#             iteration_df = P.store_output(params, iteration_df).to_frame()
#             iteration_df = iteration_df.T
#             result_df = pd.concat([result_df, iteration_df], axis = 0, ignore_index=True)    
#             result_df.columns = iteration_df.columns
#             with pd.ExcelWriter('current_exp.xlsx', engine='openpyxl', mode='a', if_sheet_exists = 'overlay') as writer:   
#                 sheet = writer.book['ArcExperiments']
#                 iteration_df.to_excel(writer, sheet_name='ArcExperiments', index=False, header=False, startrow=sheet.max_row)
      

