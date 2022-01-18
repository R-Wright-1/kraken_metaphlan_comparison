#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 31 13:05:32 2021

@author: robynwright
"""

import os
import pandas as pd
import pickle

with open('ncbi_taxid.dict', 'rb') as f:
    taxid_dict = pickle.load(f)

print(len(taxid_dict))

