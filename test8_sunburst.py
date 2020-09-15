# Calcualate values for CH2O or other species:
# Developer : Chi-Tsan Wang
# Date 2018.12.05
# Version :  0.0.5
# Requirment package : PERMM, PseudoNetCDF, matplotlib2.0.2, pandas
# Changes in 0.0.4 : Can do color change, fontsize change, fix the algorithm for the ratio problem for non-chemical
# processes (etc CH3OOH),can setup the cutoff point (self.LIM). Fixed L4 location bug. Adding output layers selections
# check point : at line 47  self.LIM  , self.RMlist
# python test7_Sunburst.py CH2O 4
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
layers = sys.argv[2]

# Setup
############################################################################
pypa_path = './'
#pypa_path = '/Users/chi-tsanwang/Dropbox/PycharmProjects/MBP2018a/'
#file_list = ['Traj_270','Traj_341','Traj_196','Traj_115','Traj_71','region_FHsouth','region_FHnorth','region_City'] #list the file tags
file_list = ['region_City_NEI2011_adj']
mosico = pypa_path+'/full_mosaic_4bin_mozart_web.yaml'
base_mech = Mechanism(mosico)
base_mech.globalize(globals())
#layers = 2
##################################################################################################
#Function and Class method section
#####################################################################
class SPCOP:

    def __init__(self,tspc):
        self.spc = tspc
        self.ydt = []
        self.pdt = []
        self.others_o = []
        self.leftpdt =[]
        self.RMList = ['NO','NO2','HNO3','M','OH', 'O3', 'NO3','HOCH2OO','H2O','PAN']
        self.RMList_cycle = [-HOCH2OO,-PAN]## Removed cycle reaction species, -OOH can be emission, don't remove
#        self.LIM = 0.01 * mech_RunA_downtown.make_net_rxn([], 'CH2O').get('CH2O').sum()*1000 #cutoff point, used CH2O as a standard species
        self.LIM = 0.1#(ppb) setup the range for combining small species #
        self.make_data()

    def make_data(self):
        if  self.spc in ['CH4', 'C2H4','C2H2','C3H8','C2H6','ISOP','BIGENE', 'APIN', 'BPIN', 'MYRC','LIMON' ,'BCARY', 'TOLUENE', 'XYLENES','C3H6']: #from emission
            self.source = pd.DataFrame(list([['-',0.0],]))
#            self.sink = pd.DataFrame(list([[self.spc, np.round(mech_RunA_downtown.make_net_rxn(self.spc, []).get(self.spc).sum() * 1000, 4)],]))
            self.sink = np.round(mech_RunA_downtown.make_net_rxn(self.spc, []).get(self.spc).sum() * 1000, 4)
            self.ratio = 1
        elif self.spc in ['-']:
            self.source = pd.DataFrame(list([['-', 0.0], ]))
            self.ratio = 0
            self.sink = 0
#        elif self.spc in :
#            self.source = pd.DataFrame(list([['-', 0.0], ]))
#            self.sink = 0
#            self.ratio = 0
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
            for y in self.product[0]:           # try to remove the self cycle species
                self.source = self.source[self.source[0] != y]
            self.others = self.source[self.source[1] < self.LIM]
            self.source = self.source[self.source[1] > self.LIM]
            self.source = self.source.append([['-', self.others[1].sum()], ],ignore_index=True)

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

#######################################################################################################################
#def cutoff(pdlist, limit, ra):
#    pdlist = pdlist * ra
#    vr[level, 0].others = vr[level, 0].source[vr[level, 0].source[1] < vr[level, 0].LIM]
#    vr[level, 0].source = vr[level, 0].source[vr[level, 0].source[1] > vr[level, 0].LIM]
#    vr[level, 0].source = vr[level, 0].source.append([['-', vr[level, 0].others[1].sum()], ], ignore_index=True)
#    labels = combiname(list(vr[1,0].source[0]),list(np.round(vr[1,0].source[1]/value_ppb*100,1)))
########################################################################################################################
#######################Plot function####################################################################################

def sunburst(vr,trajnumber,tspcs,layers):
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='polar')
    level = 0
    if vr[1,0].ratio == 1:
        value_all=vr[1,0].source[1].sum()
    elif vr[1, 0].ratio != 1:
        value_all = vr[1, 0].sink * -1
    label = ('%s \n %s ppb' % (vr[1,0].spc, np.round(value_all,1)))
    value = 100
    ax.bar([0], [1.5], [np.pi * 2],color='red',edgecolor='red')
    ax.text(0, 0, label, ha='center', va='center', fontsize=5)  # adjust Fontsize
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_axis_off()
    #    return tspcs, 'source',SPCOP(CH2O).source[1].sum(), 'sink', SPCOP(CH2O).sink[1].sum()
    ax = plt.subplot(111, projection='polar')
    #   total=np.pi * 2
    level = 1
    offset = 0
    d = np.pi * 2 / 100
    labels = []
    widths = []
    values = []
    rotations=[]
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_axis_off()
    if vr[1, 0].ratio  == 1:
        value_ppb = vr[1, 0].source[1].sum()
    elif vr[1, 0].ratio != 1:
        value_ppb = vr[1, 0].sink * -1
#    vr[level, 0].source[1] = vr[level, 0].source[1] * 1
#    vr[level, 0].others_o = vr[level, 0].source[vr[level, 0].source[1] < vr[level, 0].LIM]
#    vr[level, 0].source = vr[level, 0].source[vr[level, 0].source[1] > vr[level, 0].LIM]
#    vr[level, 0].source = vr[level, 0].source[:-1].append([['-',vr[level, 0].source[1][-1:].sum()+vr[level, 0].others_o[1].sum()],], ignore_index=True)
    labels = combiname(list(vr[1,0].source[0]),list(np.round(vr[1,0].source[1]/value_ppb*100,1)))
    labels = labels[:-1]
    widths = list(vr[1,0].source[1]/value_ppb*100*d)
    values = np.cumsum([offset * d] + widths[:-1])
    position2 =  np.cumsum((vr[1,0].source[1]/value_ppb*100))
    heights = np.array([1] * len(vr[1,0].source[1]))+ 1
    bottoms = np.zeros(len(vr[1,0].source[1])) + level + 0.5
    icol1 = np.append(['lightgreen']*(len(values)-1),['black'])   # adjust the second layer color to orange and the "others" species color to black
    rects2 = ax.bar(values, heights, widths, bottoms, linewidth=1, edgecolor='white',align='edge',color = icol1)
    for rect2, label in zip(rects2, labels):
        x = rect2.get_x() + rect2.get_width() / 2
        y = rect2.get_y() + rect2.get_height()/2
        rotation = (90 + (360 - np.degrees(x) % 180)) % 360
        ax.text(x, y, label, rotation=rotation, ha='center', va='center',fontsize=3.7)
#####################  LAYER1  ######################
    if layers == "1":
        plt.savefig('sunburst_%s_%s_ALKR05_l1_debug_test01.png' % (trajnumber,tspcs), dpi = 300 )
        plt.close()
        exit()
    else:
#####################  LAYER2  ######################
#    for m in range(6): #simple select for first 6 species
        for m in range(len(vr[level,0].source[1])-1):
            total = 100
            level = 2
            if m == 0:
                offset = 0
            else:
                offset = position2[m-1]
            if vr[level,m].source[1].sum() == 0:
                ratio0 = 0
            else:
                ratio0 = vr[level-1,0].source[1][m]/vr[level,m].source[1].sum()  # (The ratio of use/gen of species)
                if ratio0 > 1:  # used> gen : The total used species is more than the formation
                    r= 1 * vr[level,m].ratio  # all gen species will jump into used part, but only chemical, so use vr[level,m].ratio to corrent the  emission and physics process
                else:           # gen> used
                    r = ratio0 * 1
    #        print 'spc_from'  vr[level-1,0].source[m] 'spc_gen' vr[level,m].source[1].sum()
            d = np.pi * 2 / total
            print m, r, vr[level, m].spc
            labels = []
            widths = []
            values = []
            local_offset = offset
            #labels = combiname(list(vr[level,m].source[0]), list(np.round(vr[level,m].source[1]*r/value_ppb*100,1))) with precentage
    #        vr[level, m].source[1] = vr[level, m].source[1]* r
    #        vr[level, m].others_o = vr[level, m].source[vr[level, m].source[1] < vr[level, m].LIM]
    #        vr[level, m].source = vr[level, m].source[vr[level, m].source[1] > vr[level, m].LIM]
    #        vr[level, m].source = vr[level, m].source[:-1].append(
    #            [['-', vr[level, m].source[1][-1:].sum() + vr[level, m].others_o[1].sum()], ], ignore_index=True)
            labels = list(vr[level, m].source[0])
            labels = labels[:-1]
            widths = list(vr[level,m].source[1]*r/value_ppb*100 * d)
            local_offset = offset
            values = np.cumsum([offset * d] + widths[:-1])
            position3 = np.cumsum((vr[level, m].source[1]*r / value_ppb * 100))
            heights = np.array([1] * len(vr[level,m].source[1])) + 1
            bottoms = np.zeros(len(vr[level,m].source[1])) + level + 1.5
            icol2 = np.append(['lightblue'] * (len(values) - 1), ['black'])
            rects2 = ax.bar(values, heights, widths, bottoms, linewidth=1, edgecolor='white', align='edge',color=icol2)
            for rect2, label in zip(rects2, labels):
                x2 = rect2.get_x() + rect2.get_width() / 2
                y2 = rect2.get_y() + rect2.get_height() / 2
                rotation = (90 + (360 - np.degrees(x2) % 180)) % 360
                ax.text(x2, y2, label, rotation=rotation, ha='center', va='center', fontsize=4)  #adjust Fontsize
#####################  LAYER2  ######################
    if layers == "2":
        plt.savefig('sunburst_%s_%s_ALKR05_l2_debug_test01.png' % (trajnumber, tspcs), dpi=300)
        plt.close()
        exit()
    else:
#####################  LAYER3  ######################
#    for n in range(2):
        for n in range(len(vr[level, m].source[1]) - 1):
            total = 100
            level = 3
            if n == 0 and m == 0:
                offset = 0
            elif n == 0 and m != 0:
                offset = position2[m-1]
            elif n != 0 and m == 0:
                offset = position3[n-1]
            else:
                offset = position2[m-1] + position3[n-1] # This step defines the pie locations for each species
            r1=0
            if vr[level, m, n].source[1].sum() == 0:
                ratio1 = 0
            else:
                ratio1 = vr[level-1, m].source[1][n] / vr[level, m, n].source[1].sum()  # (use/gen)
                if ratio1 > 1:  # (use>gen)
                    r1 = 1 * r * vr[level, m, n].ratio
                else:  #GEN > use
                    r1 = ratio1 * r
            d = np.pi * 2 / total
            print m, n,  r, r1, vr[level, m, n].spc
            labels = []
            widths = []
            values = []
            local_offset = offset
#            vr[level, m,n].source[1] = vr[level, m,n].source[1] * r1
#            vr[level, m,n].others_o = vr[level, m,n].source[vr[level, m,n].source[1] < vr[level, m,n].LIM]
#            vr[level, m,n].source = vr[level, m,n].source[vr[level, m,n].source[1] > vr[level, m,n].LIM]
#            vr[level, m,n].source = vr[level, m,n].source[:-1].append(
#                [['-', vr[level, m,n].source[1][-1:].sum() + vr[level, m,n].others_o[1].sum()], ], ignore_index=True)
            position4 = np.cumsum((vr[level, m, n].source[1]*r1/value_ppb * 100))
#           labels = combiname(list(vr[level,m,n].source[0]), list(np.round(vr[level,m,n].source[1]*r1/value_ppb*100,1))) # with precentage
            labels = list(vr[level, m, n].source[0])
            labels = labels[:-1]
            widths = list(vr[level,m,n].source[1]*r1/value_ppb*100 * d)
            local_offset = offset
            values = np.cumsum([offset * d] + widths[:-1])
            heights = np.array([1] * len(vr[level,m,n].source[1])) + 1
            bottoms = np.zeros(len(vr[level,m,n].source[1])) + level + 2.5
            icol3 = np.append(['c'] * (len(values) - 1), ['black'])
            rects3 = ax.bar(values, heights, widths, bottoms, linewidth=1, edgecolor='white', align='edge',color=icol3)
            for rect3, label in zip(rects3, labels):
                x3 = rect3.get_x() + rect3.get_width() / 2
                y3 = rect3.get_y() + rect3.get_height() / 2
                rotation = (90 + (360 - np.degrees(x3) % 180)) % 360
                ax.text(x3, y3, label, rotation=rotation, ha='center', va='center', fontsize=4)  #adjust Fontsize
#####################  LAYER3  ######################
        if layers == "3":
            plt.savefig('sunburst_%s_%s_ALKR05_l3_debug_test01.png' % (trajnumber, tspcs), dpi=300)
            plt.close()
            exit()
        elif layers == "4":
    #####################  LAYER4  ######################
            for o in range(len(vr[level, m, n].source[1]) - 1):
                total = 100
                level = 4
                # m = 3rd layer position, n= 4th layer postion, o= 5th layer postion m!=0, n=0, o!=0
                if o == 0 :
                    offset2 = offset  # zero point
                else :
                    offset2 = position4[o-1] +offset
                r2 = 0
    #                else:
    #                    offset2 = offset + position4[o - 1]                r2 = 0
    #   debug          print 'level3_loc',offset, 'level4_loc',offset2,vr[level, m, n, o].spc, 'm=',m,'n=',n, 'o=',o,'r2=',r2, 'r1=',r1,'r=',r
    #   debug          print  vr[level, m, n, o].source
                if vr[level, m, n, o].source[1].sum() == 0:
                    ratio2 = 0
                elif vr[level, m, n, o].source[1].sum() != 0:
                    ratio2 =  vr[level - 1, m, n].source[1][o]/vr[level, m, n, o].source[1].sum()  #(use/gen)
                    if ratio2 > 1: #( used>gen )
                        r2 = 1 * r1 * vr[level, m, n, o].ratio
                    else:
                        r2 = ratio2 * r1
                d = np.pi * 2 / total
                print m, n, o, r, r1, r2, vr[level, m,n,o].spc
                labels = []
                widths = []
                values = []
                vr[level, m,n,o].source[1] = vr[level, m,n,o].source[1] * r2
                vr[level, m,n,o].others_o = vr[level, m,n,o].source[vr[level, m,n,o].source[1] < vr[1,0].LIM/5]
                vr[level, m,n,o].source = vr[level, m,n,o].source[vr[level, m,n,o].source[1] > vr[1,0].LIM/5]
                vr[level, m,n,o].source = vr[level, m,n,o].source[:-1].append(
                    [['-', vr[level, m,n,o].source[1][-1:].sum() + vr[level, m,n,o].others_o[1].sum()], ], ignore_index=True)
                labels = list(vr[level, m, n, o].source[0])
                labels = labels[:-1]
                widths = list(vr[level, m, n, o].source[1] / value_ppb * 100 * d)
                local_offset =  offset2
                values = np.cumsum([offset2 * d] + widths[:-1])
                heights = np.array([1] * len(vr[level, m, n, o].source[1])) + 1
                bottoms = np.zeros(len(vr[level, m, n, o].source[1])) + level + 3.5
                icol4 = np.append(['thistle'] * (len(values) - 1), ['black'])
    #                alpha4 = np.append([0.5] * (len(values) - 1), [1])
                rects4 = ax.bar(values, heights, widths, bottoms, linewidth=1, edgecolor='white', align='edge',
                                color=icol4)
                for rect4, label in zip(rects4, labels):
                    x4 = rect4.get_x() + rect4.get_width() / 2
                    y4 = rect4.get_y() + rect4.get_height() / 2
                    rotation = (90 + (360 - np.degrees(x4) % 180)) % 360
                    ax.text(x4, y4, label, rotation=rotation, ha='center', va='center', fontsize=3.5) #adjust Fontsize
            legend_elements = [Patch(facecolor='black', edgecolor='white',label='Other small reactants for %s' % (tspcs))]
            ax.legend(handles=legend_elements, fontsize = 'xx-small',loc='lower left')
            plt.savefig('Sunburst_plot_%s_%s_Ly4_t1_r2.png' % (trajnumber,tspcs), dpi = 300 )
            plt.close()

######################################################################################################################
################################################################################################
######################################################################################################################
################################################################################################


##################################################################################################
# Operate Code start
##################################################################################################
tspcs = "HO2"
for p in range(len(file_list)):
#for p in range(1):
    nc_irr_file = pypa_path + 'IRR_out_'+file_list[p]+'_short.nc'
    print(nc_irr_file)
    print layers
    pypa_RunA_downtown = Dataset(nc_irr_file)
    mech_RunA_downtown = Mechanism(mosico)
    mech_RunA_downtown.set_mrg(pypa_RunA_downtown)
    vr = {}
    vr[1, 0] = SPCOP(tspcs)
    for i in range(len(vr[1, 0].source[0])):
        level = 2
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
    sunburst(vr, file_list[p], tspcs, layers)
    # if layers == '1':
    #     sunburst1(vr, file_list[p],tspcs)
    #     print "*"
    # elif layers == '2':
    #     sunburst2(vr, file_list[p],tspcs)
    #     print "**"
    # elif layers == '3':
    #     sunburst3(vr, file_list[p],tspcs)
    #     print "***"
    # elif layers == '4':
    #     sunburst4(vr, file_list[p],tspcs)
    #     print "****"
    # elif layers == 'all':
    #     sunburst1(vr, file_list[p],tspcs)
    #     sunburst2(vr, file_list[p],tspcs)
    #     sunburst3(vr, file_list[p],tspcs)
    #     sunburst4(vr, file_list[p],tspcs)
    #     print "all"
