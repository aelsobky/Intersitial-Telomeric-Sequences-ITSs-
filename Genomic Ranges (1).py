#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pyranges as pr


# In[47]:


import os
path = os.path.abspath("chr_all.bed")
chr_all = pr.read_bed(path, False, nrows=1258)
#print(chr_all)
path = os.path.abspath("chr_published.bed")
chr_published = pr.read_bed(path, False, nrows=125)
#print(chr_published)
print(chr_all.intersect(chr_published))


# In[ ]:




