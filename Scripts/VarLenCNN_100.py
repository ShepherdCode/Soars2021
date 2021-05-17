#!/usr/bin/env python
# coding: utf-8

# # Variable Length CNN
# In theory, a CNN can operate on inputs of variable dimensions.
# In practice, the model is compiled to a specific input dimension.
# 
# Here, work out how to feed the CNN on sequences of different length.
# Center the real bases within padding bases.
# For example, 'A'=[1,0,0,0] and 'N'=[0,0,0,0].
# 
# Organize training sequences so batches are of similar length.
# May need to implement a per-batch callback mechanism.

# In[ ]:




