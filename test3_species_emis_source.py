# Calcualte values for CH2O or other species:
# Developer : Chi-Tsan Wang
# Date 2018.1.12
# Version :  0.0.2
# Requirment package : PERMM, PseudoNetCDF, matplotlib2.0.2, pandas
# check point : self.LIM, self.RMlist
#
import sys
import numpy as np
from permm import Mechanism
from netCDF4 import Dataset
import pylab
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch

np.set_printoptions(suppress=True)

tspcs = sys.argv[1]
#layers = sys.argv[2]

# Setup
############################################################################
pypa_path = './'
#pypa_path = '/Users/chi-tsanwang/Dropbox/PycharmProjects/MBP2018a/'
file_list = ['region_OG_EPA2014_newWRF','region_FHsouth_EPA2014_newWRF','region_FHnorth_EPA2014_newWRF','region_City_EPA2014_newWRF'] #list the file tags
#file_list = ['region_FHsouth_NEI2011_adj','region_FHnorth_NEI2011_adj','region_City_NEI2011_adj','region_OG_NEI2011_adj'] #list the file tags
#file_list = ['region_City_NEI2011_adj']
mosico = pypa_path+'/full_mosaic_4bin_mozart_web.yaml'
base_mech = Mechanism(mosico)
base_mech.globalize(globals())
#############################################################################
class SPCOP:

    def __init__(self,tspc):
        self.spc = tspc
        self.ydt = []
        self.pdt = []
        self.others_o = []
        self.leftpdt =[]
        self.RMList = ['NO','NO2','HNO3','M','OH', 'O3', 'NO3','HO2','H2O','HOCH2OO','PAN','MPAN','CH3OOH'] # should consider rm "CH3OOH" or not
        self.RMList_cycle = [-HOCH2OO,-PAN]## Removed cycle reaction species, -OOH can be emission, don't remove
#        self.LIM = 0.007 * mech_RunA_downtown.make_net_rxn([], CH2O).get(CH2O).sum()*1000 #cutoff point, used CH2O as a standard species
        self.LIM = 0.01#(ppb) setup the range for combining small species #
        self.LIM1 = 0.05
        self.LIM2 = 0.1
        self.make_data()

    def make_data(self):
        if  self.spc in ['CH4', 'C2H4','C2H2','C3H8','C2H6','ISOP','BIGENE','BIGALK', 'APIN', 'BPIN', 'MYRC','LIMON' ,'BCARY', 'TOLUENE', 'XYLENES','C3H6','MBO']: #from emission
            self.source = pd.DataFrame(list([['-',0.0],]))
            self.sink = pd.DataFrame(list([[self.spc, np.round(mech_RunA_downtown.make_net_rxn(self.spc, []).get(self.spc).sum() * 1000, 4)],]))
            self.sink = np.round(mech_RunA_downtown.make_net_rxn(self.spc, []).get(self.spc).sum() * 1000, 4)
            self.ratio = 1
        elif self.spc in ['-']:
            self.source = pd.DataFrame(list([['-', 0.0], ]))
            self.ratio = 0
            self.sink = 0
        else :
            for xpc in mech_RunA_downtown.make_net_rxn([], self.spc).sum().reactants():
                self.pdt.append((xpc, round(mech_RunA_downtown.make_net_rxn(xpc, self.spc).get(self.spc).sum() * 1000, 4),))
            for ypc in mech_RunA_downtown.make_net_rxn(self.spc,[] ).sum().products():
                self.ydt.append((ypc, round(mech_RunA_downtown.make_net_rxn(self.spc, ypc ).get(ypc).sum() * 1000, 4),))
            self.sink = np.round(mech_RunA_downtown.make_net_rxn(self.spc, self.RMList_cycle).get(self.spc).sum() * 1000, 4)
            # print self.pdf
            self.pdt.sort(key=lambda x: x[1], reverse=True)
            self.ydt.sort(key=lambda x: x[1], reverse=True)
            self.source = pd.DataFrame(self.pdt)
            self.product = pd.DataFrame(self.ydt)
            for x in self.RMList:
                self.source = self.source[self.source[0] != x]
                self.product = self.product[self.product[0] != x]
#            if level > 2:
#                for y in self.product[0]:           # try to remove the self cycle species
#                    self.source = self.source[self.source[0] != y]          # try to remove the self cycle species
#            else:
#                self.source = self.source[self.source[0] != 'CH2O']
            if level < 3 :
                self.others = self.source[self.source[1] < self.LIM]
                self.source = self.source[self.source[1] > self.LIM]
                self.source = self.source.append([['-', self.others[1].sum()], ], ignore_index=True)
            elif level > 2 & level < 5 :
                self.others = self.source[self.source[1] < self.LIM1]
                self.source = self.source[self.source[1] > self.LIM1]
                self.source = self.source.append([['-', self.others[1].sum()], ], ignore_index=True)
            elif level > 4 :
                self.others = self.source[self.source[1] < self.LIM2]
                self.source = self.source[self.source[1] > self.LIM2]
                self.source = self.source.append([['-', self.others[1].sum()], ], ignore_index=True)
            self.total_ppb = self.source[1].sum()
            self.source.index = range(len(self.source))
            if (self.source[1].sum()) == 0:
                self.ratio = 0
            else:
                self.r = (self.source[1].sum())/(self.sink * -1)
                if self.r < 1:  #(usage is bigger than generate)
                    self.ratio = self.r
                else:   #(usage is smaller than generate)
                    self.ratio = 1
######################################################################################################################

def combiname(namelist, pctlist):
    df = pd.DataFrame({'name':namelist,'pct':pctlist})
    return list(df.name +" "+  df.pct.map(str) + "%")


####################################################################################################################
def emis_contribution(vr,espcs):
    pf = pd.DataFrame(list([['-', 0.0], ]))
    pll = pd.DataFrame(list([['-', 0.0], ]))
    level = 0
    if vr[1,0].ratio == 1:
        value_all=vr[1,0].source[1].sum()
    elif vr[1, 0].ratio != 1:
        value_all = vr[1, 0].sink * -1
    label = ('%s \n %s ppb' % (vr[1,0].spc, np.round(value_all,1)))
    print label
    level = 1
    if vr[1, 0].ratio  == 1:
        value_ppb = vr[1, 0].source[1].sum()
    elif vr[1, 0].ratio != 1:
        value_ppb = vr[1, 0].sink * -1
    print vr[1, 0].source[1], vr[1, 0].sink
    pf = pf.append(vr[1, 0].source)
    for m in range(len(vr[level,0].source[1])-1):
        level = 2
        if vr[level,m].source[1].sum() == 0:
            ratio0 = 0
        else:
            ratio0 = vr[level-1,0].source[1][m]/vr[level,m].source[1].sum()  #(The ratio of use/gen of species)
            if ratio0 > 1:  #used> gen : The total used species is more than the formation
                r= 1 * vr[level,m].ratio  ## all gen species will jump into used part, but only chemical, so use vr[level,m].ratio to corrent the  emission and physics process
            else:           #gen> used
                r = ratio0 * 1
        temp_pf = pd.concat([vr[level, m].source[0],vr[level, m].source[1]*r], axis=1)
        pf = pf.append(temp_pf, ignore_index=True)
        for n in range(len(vr[level, m].source[1]) - 1):
            level = 3
            r1=0
            if vr[level, m, n].source[1].sum() == 0:
                ratio1 = 0
            else:
                ratio1 = vr[level-1, m].source[1][n] / vr[level, m, n].source[1].sum()  # (use/gen)
                if ratio1 > 1:  # (use>gen)
                    r1 = 1 * r * vr[level, m, n].ratio
                else:  #GEN > use
                    r1 = ratio1 * r
            temp_pf1 = pd.concat([vr[level, m, n].source[0], vr[level, m, n].source[1] * r1], axis=1)
            pf = pf.append(temp_pf1, ignore_index=True)
            for o in range(len(vr[level, m, n].source[1]) - 1):
                level = 4
                r2 = 0
                if vr[level, m, n, o].source[1].sum() == 0:
                    ratio2 = 0
                elif vr[level, m, n, o].source[1].sum() != 0:
                    ratio2 =  vr[level - 1, m, n].source[1][o]/vr[level, m, n, o].source[1].sum()  #(use/gen)
                    if ratio2 > 1: #( used>gen )
                        r2 = 1 * r1 * vr[level, m, n, o].ratio
                    else:
                        r2 = ratio2 * r1
                temp_pf2 = pd.concat([vr[level, m, n, o].source[0], vr[level, m, n, o].source[1] * r2], axis=1)
                pf = pf.append(temp_pf2, ignore_index=True)
                for p in range(len(vr[level, m, n, o].source[1]) - 1):
                    level = 5
                    r3 = 0
                    if vr[level, m, n, o, p].source[1].sum() == 0:
                        ratio3 = 0
                    elif vr[level, m, n, o, p].source[1].sum() != 0:
                        ratio3 = vr[level - 1, m, n, o].source[1][p] / vr[level, m, n, o, p].source[1].sum()  # (use/gen)
                        if ratio3 > 1:  # ( used>gen )
                            r3 = 1 * r2 * vr[level, m, n, o, p].ratio
                        else:
                            r3 = ratio3 * r2
                    temp_pf3 = pd.concat([vr[level, m, n, o, p].source[0], vr[level, m, n, o, p].source[1] * r3], axis=1)
                    pf = pf.append(temp_pf3, ignore_index=True)
                    for q in range(len(vr[level, m, n, o, p].source[1]) - 1):
                        level = 6
                        r4 = 0
                        if vr[level, m, n, o, p, q].source[1].sum() == 0:
                            ratio4 = 0
                        elif vr[level, m, n, o, p, q].source[1].sum() != 0:
                            ratio4 = vr[level - 1, m, n, o, p].source[1][q] / vr[level, m, n, o, p, q].source[1].sum()  # (use/gen)
                            if ratio4 > 1:  # ( used>gen )
                                r4 = 1 * r3 * vr[level, m, n, o, p, q].ratio
                            else:
                                r4 = ratio4 * r3
                        temp_pf4 = pd.concat([vr[level, m, n, o, p, q].source[0], vr[level, m, n, o, p, q].source[1] * r4], axis=1)
                        pf = pf.append(temp_pf4, ignore_index=True)
                        for s in range(len(vr[level, m, n, o, p, q].source[1]) - 1):
                            level = 7
                            r5 = 0
                            if vr[level, m, n, o, p, q, s].source[1].sum() == 0:
                                ratio5 = 0
                            elif vr[level, m, n, o, p, q, s].source[1].sum() != 0:
                                ratio5 = vr[level - 1, m, n, o, p, q].source[1][s] / vr[level, m, n, o, p, q, s].source[1].sum()  # (use/gen)
                                if ratio5 > 1:  # ( used>gen )
                                    r5 = 1 * r4 * vr[level, m, n, o, p, q, s].ratio
                                else:
                                    r5 = ratio5 * r4
                            temp_pf5 = pd.concat([vr[level, m, n, o, p, q, s].source[0], vr[level, m, n, o, p, q, s].source[1] * r5], axis=1)
                            pf = pf.append(temp_pf5, ignore_index=True)
                            for t in range(len(vr[level, m, n, o, p, q, s].source[1]) - 1):
                                level = 8
                                r6 = 0
                                if vr[level, m, n, o, p, q, s, t].source[1].sum() == 0:
                                    ratio6 = 0
                                elif vr[level, m, n, o, p, q, s, t].source[1].sum() != 0:
                                    ratio6 = vr[level - 1, m, n, o, p, q, s].source[1][t] / vr[level, m, n, o, p, q, s, t].source[1].sum()  # (use/gen)
                                    if ratio6 > 1:  # ( used>gen )
                                        r6 = 1 * r5 * vr[level, m, n, o, p, q, s, t].ratio
                                    else:
                                        r6 = ratio6 * r5
                                temp_pf6 = pd.concat([vr[level, m, n, o, p, q, s, t].source[0], vr[level, m, n, o, p, q, s, t].source[1] * r6],
                                    axis=1)
                                pf = pf.append(temp_pf6, ignore_index=True)
                                pll = pll.append(temp_pf6, ignore_index=True)
                                #print level, m,n,p,p,q,s,t, vr[level, m, n, o, p, q, s, t].source[1]

    pff = pf.groupby(pf[0]).sum()/value_ppb #combine the spc from all layers, for emission species
    pllf = pll.groupby([pll[0]]).sum()/value_ppb #combine the spc from last layers, for emission
    return pff, pllf
#    return pf, pll
##################################################################################################
# Operate Code start
##################################################################################################
#tspcs = CH2O
emission_spc = ['CH4','CH3OH','C2H4','C2H6','C2H5OH','C3H6','C3H8','BIGALK','BIGENE','ISOP','XYLENES','TOLUENE','APIN', 'BPIN', 'MYRC','LIMON','BCARY','MBO' ]
for loc in file_list:
#for p in range(1):
    nc_irr_file = pypa_path + 'IRR_out_'+loc+'.nc'
    print(nc_irr_file)
    pypa_RunA_downtown = Dataset(nc_irr_file)
    mech_RunA_downtown = Mechanism(mosico)
    mech_RunA_downtown.set_mrg(pypa_RunA_downtown)
    vr = {}
    level = 1
    vr[1, 0] = SPCOP(tspcs)
    for i in range(len(vr[1, 0].source[0])):
        level = 2
        print vr[1,0].source
        vr[level, i] = SPCOP(vr[1, 0].source[0][i])
        for k in range(len(vr[2, i].source[0])):
            level = 3
            vr[level, i, k] = SPCOP(vr[2, i].source[0][k])
            for l in range(len(vr[3, i, k].source[0])):
                level = 4
                vr[level, i, k, l] = SPCOP(vr[3, i, k].source[0][l])
                for m in range(len(vr[4, i, k, l].source[0])):
                    level = 5
                    vr[level, i, k, l, m] = SPCOP(vr[4, i, k, l].source[0][m])
                    for p in range(len(vr[5, i, k, l, m].source[0])):
                        level = 6
                        vr[level, i, k, l, m, p] = SPCOP(vr[5, i, k, l, m].source[0][p])
                        for q in range(len(vr[6, i, k, l, m, p].source[0])):
                            level = 7
                            vr[level, i, k, l, m, p, q] = SPCOP(vr[6, i, k, l, m, p].source[0][q])
                            for r in range(len(vr[7, i, k, l, m, p, q].source[0])):
                                level = 8
                                vr[level, i, k, l, m, p, q, r] = SPCOP(vr[7, i, k, l, m, p, q].source[0][r])
                                #print i,k,l,m,p,q,r, vr[level, i, k, l, m, p, q, r].spc
    (pfftable, pllftable) = emis_contribution(vr, emission_spc)
    emistable = pd.DataFrame(list([['-', 0.0], ])) # creat a new pandas table for emission
    for emis in emission_spc: # input the select emission species and make the table
        if emis in pfftable[1] :
            emistable = emistable.append(pd.DataFrame(list([[emis, pfftable[1][emis]],])),ignore_index=True)
        else :
            print 'There is no %s here.' %(emis)
    print 'outupt file %s_%s_emission_with_total_l8_list.csv' % (loc,tspcs)
    print 'outupt file %s_%s_lastlayer_with_total_l8_list.csv' % (loc,tspcs)
    emistable.to_csv('%s_%s_emission_with_total_l8_list.csv' % (loc,tspcs))
    pllftable.to_csv('%s_%s_lastlayer_with_total_l8_list.csv' % (loc,tspcs))
    pfftable = []
    pllftable = [] #clean the table after output