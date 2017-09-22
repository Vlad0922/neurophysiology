# -*- coding: utf-8 -*-

import os
import sys

# from PyQt4.QtGui import *

import matplotlib.pyplot as plt

import xlrd
import xlwt
from xlutils.copy import copy

def get_workbook(fname, indices, data_keys):
    if not os.path.isfile(fname):
        wb = xlwt.Workbook()
        ws = wb.add_sheet('all_results')

        headers = list(indices)
        headers.extend(data_keys)

        for idx, h in enumerate(headers):
            ws.write(0, idx, h)
        wb.save(fname)
    
    wb = xlrd.open_workbook(fname)
    return wb



def write_to_excel(fname, sheet, df, indices):
    data_keys = sorted(list(set(df.keys()) - set(indices)))
    rb = get_workbook(fname, indices, data_keys)
    
    row = [df[ind] for ind in indices]
    row.extend([df[ind] for ind in data_keys])
    
    nrows = rb.sheet_by_index(0).nrows

    wb = copy(rb) 
    ws = wb.get_sheet(0) 
    for idx, v in enumerate(row):
        ws.write(nrows, idx, v)

    wb.save(fname)    
    
def rem_key(data, keys):
    return {k:data[k] for k in data.keys() if k not in keys}


# class MyTable(QTableWidget):
#     def __init__(self, data, index_header):
#         QTableWidget.__init__(self, len(data.values()[0]), len(data))
#         self.index_header = index_header
#         self.data = data
#         self.setmydata()
#         self.resizeColumnsToContents()
#         self.resizeRowsToContents()
        
 
#     def setmydata(self):
#         horHeaders = []
        
#         self.init_index(horHeaders)
        
#         for n, key in enumerate(sorted(rem_key(self.data, self.index_header).keys())):
#             horHeaders.append(key)
#             for m, item in enumerate(self.data[key]):
#                 newitem = QTableWidgetItem(str(item))
#                 self.setItem(m, n+1, newitem)
                
#         self.setHorizontalHeaderLabels(horHeaders)
  
    
#     def init_index(self, headers):
#         headers.append(self.index_header)
#         for m, item in enumerate(self.data[self.index_header]):
#             newitem = QTableWidgetItem(str(item))
#             self.setItem(m, 0, newitem)
        

# def draw_table(data):   
#     a = QApplication(sys.argv)
 
#     table = MyTable(data, 'data_name')
#     table.show()
#     return a.exec_()