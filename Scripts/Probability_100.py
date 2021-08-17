#!/usr/bin/env python
# coding: utf-8

# In[32]:


def to_string(num,base=2):
    if num==0:
        return '0'
    q = num//base
    r = num%base
    return to_string(q,base) + str(r)
print("16 in binary:",to_string(16))
print("11 base 3:",to_string(11,3))


# In[45]:


def has_consec(string,symbol,min):
    consec=symbol*min
    has=0
    for i in range(0,len(string)-min+1):
        if string[i:i+min] == consec:
            has = 1
            print("YES! {:>10}".format(string))
            break
    return has
print ("15 has 3 consecutive ones?",has_consec(to_string(15),'1',3))
print ("16 has 3 consecutive ones?",has_consec(to_string(16),'1',3))


# In[46]:


def sum_consec_one(n,m,base):
    symbol='1'
    c = 0
    for i in range(0,base**n):
        s = to_string(i,base)
        c = c + has_consec(s,symbol,m)
    return c
base = 2
n = 8
m = 7
c = sum_consec_one(n,m,base)
print ("Count strings with at least %d consec ones in 2^%d bits = %d"%(m,n,c) )


# In[49]:


base = 2
n = 8
for m in range(8,0,-1):
    c = sum_consec_one(n,m,base)
    print ("Count strings with at least %d consec ones in 2^%d bits = %d"%(m,n,c) )


# In[ ]:




