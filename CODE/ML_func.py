# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 12:16:19 2020
@author: liamo
"""
import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder,MinMaxScaler
from sklearn.feature_selection import SelectKBest
from sklearn.metrics import make_scorer,matthews_corrcoef
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

def prepare_targets(Y):
    y_test_enc = LabelEncoder().fit_transform(Y)
    return y_test_enc

def feature_select(k,X,Y):
    fs = SelectKBest(k=k)
    fs.fit(X, Y)
    X_fs = fs.transform(X)
    return X_fs

def ClassificationCV(X,Y,cv=5,rep=1):
    results=[]
    mcc_scorer = make_scorer(matthews_corrcoef)

    for i in range(rep):
        model = RandomForestClassifier()

        cv_results = cross_val_score(model,X,Y,cv=cv,scoring=mcc_scorer)
        results.append(cv_results)
        # print(f"{i+1}: {round(cv_results.mean(),4)}") # ({cv_results.std()})

    return results
