#!/usr/bin/env python
# coding: utf-8

# # Graph
# Listing all paths from source to sink in an undirected weighted graph.  
# Jason Miller.  

# In[1]:


import numpy as np


# ## Define the Graph
# Represent the graph on page 377, Fig 5.3, of Kurose & Ross 7th ed.

# In[2]:


G=np.zeros( (6,6), dtype=np.int8)
NODES=['u','v','w','x','y','z']
nodecount = len(NODES)
u,v,w,x,y,z=0,1,2,3,4,5


# In[3]:


def add_edge(graph,node1,node2,weight):
    graph[node1,node2]=graph[node2,node1]=weight
add_edge(G,u,v,2)
add_edge(G,u,x,1)
add_edge(G,u,w,5)
add_edge(G,x,v,2)
add_edge(G,x,w,3)
add_edge(G,x,y,1)
add_edge(G,w,v,3)
add_edge(G,w,y,1)
add_edge(G,w,z,5)
add_edge(G,y,z,2)


# In[4]:


print("Node count:",len(G))
print("Edge count:",np.count_nonzero(G)//2)
print(G)


# ## Graph utility functions

# In[11]:


# Recursive DFS
def list_paths_from(graph,this_path,all_paths):
    this_node=this_path[-1]
    for i in range(0,len(graph)):
        if i not in this_path: # avoid cycles
            if graph[this_node,i]>0: # positive weight => adjacency
                next_path = this_path.copy()
                next_path.append(i)
                #print("this_path",this_path)
                #print("next_path",next_path)
                all_paths.append(next_path)
                all_paths = list_paths_from(graph,next_path,all_paths)
    return all_paths
# Convert node names from ints to letters
def get_path_name(one_path):
    pname=''
    for i in one_path:
        pname = pname + NODES[i]
    return pname
# Print a list of all paths
def show_all_paths(paths):
    c = 0
    for p in paths:
        c = c+1
        print(c,get_path_name(p))
# Restrict a list based on last node
def list_paths_a_to_b(graph,a,b):
    paths_from_a = [[a]]
    paths_from_a = list_paths_from(G,[a],paths_from_a)
    paths_a_to_b = []
    for p in paths_from_a:
        if p[-1]==b:
            paths_a_to_b.append(p)
    show_all_paths(paths_a_to_b)


# ## Graph exploration 

# In[12]:


# Long output, was used for debugging
paths_from_Y = [[y]]
paths_from_Y = list_paths_from(G,[y],paths_from_Y)
# show_all_paths(paths_from_Y)


# In[7]:


list_paths_a_to_b(G,y,u)


# In[8]:


list_paths_a_to_b(G,x,z)


# In[9]:


list_paths_a_to_b(G,z,u)


# In[10]:


list_paths_a_to_b(G,z,w)


# In[ ]:




