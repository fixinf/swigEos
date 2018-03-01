import numpy as np
import Models2
# coding: utf-8

# # Solving the branches from contours

# ## Separating minima and maxima

# In[1]:


import fixedNF


# In[2]:


folder = '/home/const/MEGA/Physics/RhoCond/second_sol/data/contours_HDp/'


# In[3]:


#List of a
alist = np.array([0.7, 0.6, 0.5, 0.4, 0.3])

# List for b and f
plist = np.array([[0.16911015273196256, 0.20000000000000004], 
                  [0.2056222713669307, 0.30000000000000004], 
                  [0.2801188525986087, 0.4], 
                  [0.4037903815564332, 0.5], 
                  [0.6316556702588175, 0.6],
                  [0., 100]])


# In[4]:


#create models with varying parameters
models = []
for a, p in zip(alist, plist.tolist()):
    b, f = p
    wr = Models2.MKVOR_phi_rho(a, 0.01, f, b)
    wr.setDeltaPotential(-50.)
    m = wr.rcond_delta_phi
    models.append(m)
    m.loadEos()


# In[5]:


m = models[-3]
a = alist[-3]
print(a)


# In[6]:


#paths = np.load(join(folder, 'pathes_%.2f.npy'%a))

#branches = fixedNF.getBranches(m, paths, iterations=500)

#br1 = fixedNF.refineBranches(m, branches, iterations=100)

#n, init, f = m.reset_from_branches(br1, nmax=3.8)


# In[6]:


m.loadEos()


# In[7]:


m.processMaxw()


# In[8]:


any(np.diff(m.nrange_maxw) < 0)


# In[ ]:


mass = m.dumpMassesCrust(write=1)


# In[ ]:


m

