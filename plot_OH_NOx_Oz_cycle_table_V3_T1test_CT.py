### ggp (Nov 2018): took out HONO to work with T1

import sys
import numpy as np
import numpy.ma as ma
import pylab
from permm import Mechanism
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import string
import netCDF4
import matplotlib
import code 
#mosico = './full_mosaic_4bin_mozart_CT.yaml'
#mosico = './full_mosaic_4bin_mozart_CT_mods.yaml'
mosico = 'MZ148_T1mozcart_20170323_CT.yaml'
base_mech = Mechanism(mosico)
base_mech.globalize(globals())

### Downtown PA domain first
i1=0


trajnumber = ['region_City_EPA2014_newWRF_T1','region_OG_EPA2014_newWRF_T1']

i=0
#pypa_path = '../../'
pypa_path = '/home/pfister/IRR_analysis/'
#nc_f = nc_irr_file
#nc_fid = Dataset(nc_f, 'r')  # Dataset is the class behavior to open the file
#ydata = nc_fid.variables['ALD_ppb'][:]
#i2=len(ydata)
#nc_fid.close()
#irange = i2

#print(irange)
mech_RunA_downtown = Mechanism(mosico)

########################################################################
############################class and methods###########################
class IRROP:

    def __init__(self,irrdata):
        self.irrdata = irrdata
#        self.spcop = SPCOP()
        self.ydt = []
        self.pdt = []
        self.OH_VOC = pd.DataFrame([])
        self.HOx_src = pd.DataFrame([])
        self.chem_bal = pd.DataFrame([])
#        self.chem_loss = pd.DataFrame([])
        self.inirxn()
        self.hoxcyc()
        self.spcbal()
        self.noxcyc()

    def inirxn(self):

        voc_source = self.irrdata.make_net_rxn([OH,VOC], []).sum().reactants()
        for x in voc_source:
            if x == 'OH':
                self.OH = pd.DataFrame({x: np.round(self.irrdata.make_net_rxn([OH, x], [])[OH.reactant()][i1:i2] * -1000, 4)})
            else:
                self.OH_VOC = pd.concat(
                [self.OH_VOC, pd.DataFrame({x: np.round(self.irrdata.make_net_rxn([OH, x], [])[OH.reactant()][i1:i2] * -1000, 4)})],
                axis=1)
                
                
        self.OH_VOC = pd.concat([pd.DataFrame({'VOC': np.round(self.irrdata.make_net_rxn([OH, VOC], [])[OH.reactant()][i1:i2] * -1000, 4)}), self.OH_VOC],
                axis=1)
        self.VOCrank = self.OH_VOC.sum().sort_values(ascending=False)
        
    def hoxcyc(self):
        hox_source = ['ALD','CH2O','CH3COCHO','GLYOXAL','CH3CHO']
        for x in hox_source:
            self.HOx_src = pd.concat(
               [self.HOx_src, pd.DataFrame({x: np.round(self.irrdata.make_net_rxn([x,-Radical], [HOx])[HOx.product()][i1:i2] * 1000, 4)})],
                axis=1)
        self.HOx_src = pd.concat([pd.DataFrame({'O1D': np.round(self.irrdata.make_net_rxn([O1D], [OH])[OH.product()][i1:i2] * 1000, 4)}), self.HOx_src],
                axis=1)
        self.HOx_src = pd.concat([self.HOx_src,pd.DataFrame({'all_new_HOx' : self.HOx_src['ALD']+self.HOx_src['O1D']})],axis=1)
        self.HOx_term = pd.DataFrame({'HOx_term': np.round(self.irrdata.make_net_rxn(HOx, [HOz, -Radical])[HOx.reactant()][i1:i2] * -1000, 4),
                                      'HOx_return': np.round(self.irrdata.make_net_rxn([HOz, -Radical], HOx )[HOx.product()][i1:i2] * 1000, 4),
                                      'Bal_termHOx' : np.round(self.irrdata.make_net_rxn(HOx, [HOz, -Radical])[HOx.reactant()][i1:i2] * -1000, 4)+
                                                      np.round(self.irrdata.make_net_rxn([HOz, -Radical], HOx)[HOx.product()][i1:i2] * 1000, 4)})
        self.ohcyc = self.OH_VOC['VOC']/self.HOx_src['all_new_HOx']

    def noxcyc(self):
        self.NOtoNO2 = pd.DataFrame({'NOtoNO2': np.round(self.irrdata.make_net_rxn([NO], [NO2])[NO2.product()][i1:i2] * 1000, 4)})
        Ox_source = ['O3','TRO2','RO2','XO2','HO2']
        for x in Ox_source:
            self.NOtoNO2 = pd.concat(
               [self.NOtoNO2, pd.DataFrame({x: np.round(self.irrdata.make_net_rxn([NO,x], [NO2])[NO2.product()][i1:i2] * 1000, 4)})],
                axis=1)

        self.bal_NOz = pd.DataFrame({'NOz': np.round(self.irrdata.make_net_rxn([], [])[NOz][i1:i2] * 1000, 4)})
        NOz_list = ['HNO3','PAN','MPAN','ONITR']
        for x in NOz_list:
            self.bal_NOz = pd.concat(
               [self.bal_NOz, pd.DataFrame({x: np.round(self.irrdata.make_net_rxn([], [])[x][i1:i2] * 1000, 4)})],
                axis=1)

    def spcbal(self):
        spclist=['O3','CH3CHO','CH2O']
        for x in spclist:
            self.chem_bal = pd.concat(
                [self.chem_bal, pd.DataFrame(
                {x: np.round(self.irrdata.make_net_rxn([],[])[x][i1:i2]*1000, 4)})],
                axis=1)

#            self.chem_gen = pd.concat([self.chem_gen, pd.DataFrame(
#                {x: np.round(self.irrdata.make_net_rxn([], x)[x.product()][i1:i2]*1000, 4)})],
#            axis=1)
#            self.chem_loss = pd.concat([self.chem_loss, pd.DataFrame(
#                {x: np.round(self.irrdata.make_net_rxn(x, [])[x.reactant()][i1:i2]*1000, 4)})],
#            axis=1)
#        self.chem_bal = self.chem_gen + self.chem_loss
#HNO3 + PAN + MPAN + ONITR

#class SPCOP:

#    def __init__(self,tspc):
#        self.spc = tspc
#        self.ydt = []
#        self.pdt = []
#        self.spc_balance()

#    def spc_balance(self):
#        for xpc in self.irrdata.make_net_rxn([], self.spc).sum().reactants():
#            self.pdt.append((xpc, round(self.irrdata.make_net_rxn(xpc, self.spc).get(self.spc).sum() * 1000, 4),))
#        for ypc in self.irrdata.make_net_rxn(self.irrdata.spc,[] ).sum().products():
#            self.ydt.append((ypc, round(self.irrdata.make_net_rxn(self.spc, ypc ).get(self.spc).sum() * 1000, 4),))
            # print self.pdf
#        self.pdt.sort(key=lambda x: x[1], reverse=True)
#        self.ydt.sort(key=lambda x: x[1], reverse=True)
#        self.source = pd.DataFrame(self.pdt)
#        self.sink = pd.DataFrame(self.ydt)

########################################################################################################################
#plot function section
def barplot(data):
    fig = pylab.figure()
    ax  = fig.add_subplot(111)
    ax.set_title(' traj'+trajnumber[f]+' OH+VOCs')
    ax.set_ylabel('OH react with VOCs (ppb/hr)', size=14)
    ax.set_xlabel('Time (hour)', size=14)
    ind = np.arange(len(data.OH_VOC))
    width = 0.5
##    colors_list = list(colors._colors_full_map.values())
    colors_list = list(matplotlib.colors.cnames)
    print colors_list
   
#    code.interact(local=locals())

#    ax.bar(ind, data.OH_VOC['VOC'], color='gray', alpha=0.5)
#    ax.bar(ind, data.OH_VOC['CH2O'],  color='red', alpha=1)
    bottoms = 0
    icol=np.append([0,107,103,1,116,92,71],[31]*len(data.OH_VOC.columns))
    for i in  range(1,len(data.OH_VOC.columns)):
#        ax.bar(ind, data.OH_VOC[data.VOCrank.index[i]],bottom = bottoms ,color=colors_list[i*3], alpha=1)
        ax.bar(ind, data.OH_VOC[data.VOCrank.index[i]],bottom = bottoms ,color=colors_list[icol[i]], alpha=1)
        bottoms = data.OH_VOC[data.VOCrank.index[i]] + bottoms
#        print  data.VOCrank.index[i],bottoms[7:18]
    ax.set_ylim(0, 4, 1)
    ax.set_xticks(range(len(data.OH_VOC)))
    ax.set_xticklabels(np.arange(len(data.OH_VOC))
        )
    leg = ax.legend((data.VOCrank.index[1:7]), ncol=1, loc=1)
    pylab.setp(leg.get_texts(), fontsize=12)
    fig.savefig('./t' + trajnumber[f] + '_OHwVOC_bar_for_Mazcart_chitsan_T1.pdf')
    plt.close()

def spc_bal_barplot(data):
    for spc in data.chem_bal.columns:
#        print spc
        fig2 = pylab.figure()
        ax = fig2.add_subplot(111)
        ax.set_title(' traj' + trajnumber[f] + ' Net_'+ spc +'_Production')
        ax.set_ylabel('ppb/hour', size=14)
        ax.set_xlabel('Time (hour)', size=14)
        ind = np.arange(len(data.chem_bal))
        width = 0.5
        ax.bar(ind, data.chem_bal[spc], color='g', alpha=0.5)
#        print(min(data.chem_bal[spc]) * 0.9)
#        print(max(data.chem_bal[spc]) * 1.1)
        ax.set_ylim(min(data.chem_bal[spc]) * 0.9, max(data.chem_bal[spc]) * 1.1, 1)
        ax.grid(color='lightgray', linestyle=':', linewidth=1)
        ax.set_xticks(range(len(data.chem_bal)))
        ax.set_xticklabels(np.arange(len(data.OH_VOC)))
#        print('./t' + trajnumber[f] + '_'+spc+'_production.pdf')
        fig2.savefig('./t' + trajnumber[f] + '_'+spc+'_production_T1.pdf')
        plt.close()

def HOx_source_plot(data):
    fig = pylab.figure()
    linetypelist = ['r-','b-','g-','c-','p-']
    HOx_list = ['all_new_HOx','ALD','O1D']
    for n in range(0,len(HOx_list)):
        plt.plot(data.HOx_src[HOx_list[n]],linetypelist[n],label='HOx from %s'%(HOx_list[n]))
    plt.plot(data.OH_VOC['VOC'],'k-',label = 'OH+VOC')
    plt.title(' traj'+trajnumber[f])
    plt.ylabel('ppb')
    plt.xlabel('Hours')
    plt.legend()
    fig.savefig('./t'+trajnumber[f]+'_newOH_CT.pdf')
    plt.close()

############################################
###Operate section#########################
for f in range(len(trajnumber)): # running all files by for loop
    nc_irr_file = pypa_path + 'IRR_out_' + trajnumber[f] + '_short.nc'
    pypa_RunA_downtown = Dataset(nc_irr_file)
    mech_RunA_downtown.set_mrg(pypa_RunA_downtown)
    i2 = len(mech_RunA_downtown.make_net_rxn([OH],[])[OH.reactant()])
    i1=0
    irrdata = mech_RunA_downtown
    IRD=IRROP(irrdata) # generate a data include all tables that needed, such as OH+VOC, new OH source, NOx cycle....
#   print IRD.OH_VOC  # Show OH_VOC table
#   print IRD.HOx_src # Show HOx source table
#   All tables list generate after IRROP , and you can echo the data by using:
#   IRD.OH_VOC    IRD.HOx_src    IRD.HOx_term     IRD.ohcyc     IRD.NOtoNO2     IRD.bal_NOz     IRD.chem_bal
#   For example:
#    print IRD.OH_VOC
    print IRD.HOx_src
#    print IRD.ohcyc
#    print IRD.NOtoNO2
#    print IRD.bal_NOz
#    print IRD.chem_bal
#######   plot the figures  ###############################
    barplot(IRD)  # apply barplot function for IRD OH_VOC data.
    spc_bal_barplot(IRD)   #apply spc_bal_barplot function for spcies balance table in IRD.
    HOx_source_plot(IRD)
#############################################


