# Figure - Simulation of fedbatch reactor with yeast
#          with functions added to facilitate explorative simulation work
#
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2022-09-17 - Updated for FMU-explore 0.9.5
# 2023-02-09 - Updated to FMU-explore 0.9.6e
# 2023-02-13 - Consolidate FMU-explore to 0.9.6 and means parCheck and par() udpate and simu() with opts as arg
# 2023-02-20 - Updatation to updated BPL.Control and block VarLimPID used
# 2023-02-22 - Adjusted parDict, parLocation and simu('cont')
# 2023-02-23 - Added Kla_O2 and Kla_CO2 to the derived parameters that can be reached by describe()
# 2023-02-28 - Update FMU-explore for FMPy 0.9.6 in one leap and added list key_variables for logging
# 2023-03-20 - Update FMU-explore for FMPy 0.9.7b includding prel simu('cont') ans for LimPID
# 2023-03-21 - Correcting the script by including logging of states in a pedestrian way
# 2023-03-23 - Update FMU-explore 0.9.7c
# 2023-03-28 - Update FMU-explore 0.9.7
# 2023-04-21 - Compiled for Ubuntu 20.04 and changed BPL_version
# 2023-05-31 - Adjusted to from importlib.meetadata import version
# 2024-01-22 - Update FMU-explore 0.9.9 icluding function process_diagram() although not used since no GUI gasphase
# 2024-03-08 - Update to finalize transition to FMU-explore 0.9.9
# 2024-05-08 - Look through early part and call it all FMU-explore 1.0.0
# 2024-05-20 - Updated the OpenModelica version to 1.23.0-dev
# 2024-06-01 - Corrected model_get() to handle string values as well - improvement very small and keep ver 1.0.0
# 2024-07-03 - Updated with adjustments of names (connections to FixValue) due to new FMU with BPL 2.2.1 
# 2024-08-13 - Corrected model_get() to handle calculatedParameters - call it ver 1.0.1
# 2024-10-13 - Updated information aboout BPL 2.2.2 - GUI
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Framework
#------------------------------------------------------------------------------------------------------------------

# Setup framework
import sys
import platform
import locale
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.image as img
import zipfile

from fmpy import simulate_fmu
from fmpy import read_model_description
import fmpy as fmpy

from itertools import cycle
from importlib.metadata import version  

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
global fmu_model, model
if platform.system() == 'Windows':
   print('Windows - run FMU pre-compiled JModelica 2.14')
   fmu_model ='BPL_YEAST_AIR_Fedbatch_DOcontrol_windows_jm_cs.fmu'        
   model_description = read_model_description(fmu_model)  
   flag_vendor = 'JM'
   flag_type = 'CS'
elif platform.system() == 'Linux':  
   flag_vendor = 'OM'
   flag_type = 'ME'
   if flag_vendor in ['OM','om']:
      print('Linux - run FMU pre-comiled OpenModelica') 
      if flag_type in ['CS','cs']:         
         fmu_model ='BPL_YEAST_AIR_Fedbatch_DOcontrol_om_cs.fmu'    
         model_description = read_model_description(fmu_model) 
      if flag_type in ['ME','me']:         
         fmu_model ='BPL_YEAST_AIR_Fedbatch_DOcontrol_linux_om_me.fmu'    
         model_description = read_model_description(fmu_model) 
   else:    
      print('There is no FMU for this platform')

# Provide various opts-profiles
if flag_type in ['CS', 'cs']:
   opts_std = {'ncp': 500}
elif flag_type in ['ME', 'me']:
   opts_std = {'ncp': 500}
else:    
   print('There is no FMU for this platform')
  
# Provide various MSL and BPL versions
if flag_vendor in ['JM', 'jm']:
   constants = [v for v in model_description.modelVariables if v.causality == 'local'] 
   MSL_usage = [x[1] for x in [(constants[k].name, constants[k].start) \
                     for k in range(len(constants))] if 'MSL.usage' in x[0]][0]   
   MSL_version = [x[1] for x in [(constants[k].name, constants[k].start) \
                       for k in range(len(constants))] if 'MSL.version' in x[0]][0]
   BPL_version = [x[1] for x in [(constants[k].name, constants[k].start) \
                       for k in range(len(constants))] if 'BPL.version' in x[0]][0] 
elif flag_vendor in ['OM', 'om']:
   MSL_usage = '3.2.3 - used components: RealInput, RealOutput, LimPID-components' 
   MSL_version = '3.2.3'
   BPL_version = 'Bioprocess Library version 2.2.2 - GUI' 
else:    
   print('There is no FMU for this platform')

# Simulation time
global simulationTime; simulationTime = 20.0
global prevFinalTime; prevFinalTime = 0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture', 'bioreactor.gas_liquid_transfer']

# Provide process diagram on disk
fmu_process_diagram ='BPL_GUI_TEST2_Fedbatch_process_diagram_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateDict, parDict, diagrams, newplot(), describe()
#------------------------------------------------------------------------------------------------------------------

# Create stateDict that later will be used to store final state and used for initialization in 'cont':
global stateDict; stateDict =  {}
stateDict = {variable.derivative.name:None for variable in model_description.modelVariables \
                                            if variable.derivative is not None}
stateDict.update(timeDiscreteStates) 

global stateDictInitial; stateDictInitial = {}
for key in stateDict.keys():
    if not key[-1] == ']':
         if key[-3:] == 'I.y':
            stateDictInitial[key] = key[:-10]+'I_start'
         elif key[-3:] == 'D.x':
            stateDictInitial[key] = key[:-10]+'D_start'
         else:
            stateDictInitial[key] = key+'_start'
    elif key[-3] == '[':
        stateDictInitial[key] = key[:-3]+'_start'+key[-3:]
    elif key[-4] == '[':
        stateDictInitial[key] = key[:-4]+'_start'+key[-4:]
    elif key[-5] == '[':
        stateDictInitial[key] = key[:-5]+'_start'+key[-5:] 
    else:
        print('The state vector has more than 1000 states')
        break

global stateDictInitialLoc; stateDictInitialLoc = {}
for value in stateDictInitial.values(): stateDictInitialLoc[value] = value

# Create dictionaries parDictLocation[] and parLocation[]
global parDict; parDict = {}
parDict['V_start'] = 4.5
parDict['VX_start'] = 4.5*1.0
parDict['VG_start'] = 4.5*5.0
parDict['VE_start'] = 0.0

parDict['V_tot'] = 8.0
parDict['V_diss_O2_start'] = 0
parDict['V_diss_CO2_start'] = 0
parDict['V_gas_N2_start'] = 2.4
parDict['V_gas_O2_start'] = 0.6
parDict['V_gas_CO2_start'] = 0

parDict['qGmax'] = 20.0e-3
parDict['Ks'] = 10.0e-3
parDict['qO2lim'] = 6.9e-3

parDict['alpha_O2'] = 1.0

parDict['feedtank_V_start'] = 50.0
parDict['G_in'] = 500.0
parDict['F_start'] = 0.0
parDict['mu_feed'] = 0.10
parDict['t_startExp'] = 3.0
parDict['F_startExp'] = 0.00133
parDict['F_max'] = 0.3

parDict['airFlow_setpoint'] = 120.0

parDict['DO_setpoint'] = 40.0
parDict['DO_sensor_x_start'] = 87.0
parDict['K'] = 10.0
parDict['Ti'] = 0.5
parDict['Td'] = 0.0
parDict['I_start'] = 0
parDict['D_start'] = 0
parDict['N_low'] = 500
parDict['N_high'] = 2000

global parLocation; parLocation = {}
parLocation['V_start'] = 'bioreactor.V_start'
parLocation['VX_start'] = 'bioreactor.m_start[1]' 
parLocation['VG_start'] = 'bioreactor.m_start[2]' 
parLocation['VE_start'] = 'bioreactor.m_start[3]' 

parLocation['V_tot'] = 'bioreactor.V_tot'
parLocation['V_diss_O2_start'] = 'bioreactor.m_start[4]'
parLocation['V_diss_CO2_start'] = 'bioreactor.m_start[5]'
parLocation['V_gas_N2_start'] = 'bioreactor.V_gas_start[1]'
parLocation['V_gas_O2_start'] = 'bioreactor.V_gas_start[2]'
parLocation['V_gas_CO2_start'] = 'bioreactor.V_gas_start[3]'

parLocation['qGmax'] = 'bioreactor.culture.qGmax' 
parLocation['Ks'] = 'bioreactor.culture.Ks' 
parLocation['qO2lim'] = 'bioreactor.culture.qO2lim' 

parLocation['alpha_O2'] = 'bioreactor.gas_liquid_transfer.alpha_O2'

parLocation['feedtank_V_start'] = 'feedtank.V_start'
parLocation['G_in'] = 'feedtank.c_in[2]'
parLocation['F_start'] = 'dosagescheme.F_start'
parLocation['mu_feed'] = 'dosagescheme.mu_feed'
parLocation['t_startExp'] = 'dosagescheme.t_startExp'
parLocation['F_startExp'] = 'dosagescheme.F_startExp'
parLocation['F_max'] = 'dosagescheme.F_max'

parLocation['airFlow_setpoint'] = 'airFlow_setpoint.val'

parLocation['DO_setpoint'] = 'DO_setpoint.val'
parLocation['DO_sensor_x_start'] = 'DOsensor.x_start'
parLocation['K'] = 'PIreg.K'
parLocation['Ti'] = 'PIreg.Ti'
parLocation['Td'] = 'PIreg.Td'
parLocation['I_start'] = 'PIreg.I_start'
parLocation['D_start'] = 'PIreg.D_start'
parLocation['N_low'] = 'N_low.val'
parLocation['N_high'] = 'N_high.val'

# Extended list of parameters and variables only for describe and not change
global key_variables; key_variables = []
parLocation['mu'] = 'bioreactor.culture.mu'; key_variables.append(parLocation['mu'])
parLocation['qO2'] = 'bioreactor.culture.qO2'; key_variables.append(parLocation['qO2'])

parLocation['Kla_O2'] = 'bioreactor.gas_liquid_transfer.Kla_O2'; key_variables.append(parLocation['Kla_O2'])
parLocation['Kla_CO2'] = 'bioreactor.gas_liquid_transfer.Kla_CO2'; key_variables.append(parLocation['Kla_CO2'])

parLocation['V_tot'] = 'bioreactor.V_tot'; key_variables.append(parLocation['V_tot'])
parLocation['bioreactor.V'] = 'bioreactor.V'; key_variables.append(parLocation['bioreactor.V'])
parLocation['bioreactor.m[1]'] = 'bioreactor.m[1]'; key_variables.append(parLocation['bioreactor.m[1]'])
parLocation['bioreactor.m[2]'] = 'bioreactor.m[2]'; key_variables.append(parLocation['bioreactor.m[2]'])
parLocation['bioreactor.m[3]'] = 'bioreactor.m[3]'; key_variables.append(parLocation['bioreactor.m[3]'])
parLocation['bioreactor.m[4]'] = 'bioreactor.m[4]'; key_variables.append(parLocation['bioreactor.m[4]'])
parLocation['bioreactor.m[5]'] = 'bioreactor.m[5]'; key_variables.append(parLocation['bioreactor.m[5]'])

parLocation['bioreactor.V_gasphase'] = 'bioreactor.V_gasphase'; key_variables.append(parLocation['bioreactor.V_gasphase'])
parLocation['bioreactor.V_gas[1]'] = 'bioreactor.V_gas[1]'; key_variables.append(parLocation['bioreactor.V_gas[1]'])
parLocation['bioreactor.V_gas[2]'] = 'bioreactor.V_gas[2]'; key_variables.append(parLocation['bioreactor.V_gas[2]'])
parLocation['bioreactor.V_gas[3]'] = 'bioreactor.V_gas[3]'; key_variables.append(parLocation['bioreactor.V_gas[3]'])
parLocation['bioreactor.V_gas[4]'] = 'bioreactor.V_gas[4]'; key_variables.append(parLocation['bioreactor.V_gas[4]'])

parLocation['DO_setpoint.out'] = 'DO_setpoint.out'; key_variables.append(parLocation['DO_setpoint.out'])
parLocation['DOsensor.x'] = 'DOsensor.x'; key_variables.append(parLocation['DOsensor.x'])
parLocation['PIreg.limPID.D.x'] = 'PIreg.limPID.D.x'; key_variables.append(parLocation['PIreg.limPID.D.x'])
parLocation['PIreg.limPID.I.y'] = 'PIreg.limPID.I.y'; key_variables.append(parLocation['PIreg.limPID.I.y'])

parLocation['airtube.V'] = 'airtube.V'; key_variables.append(parLocation['airtube.V'])
parLocation['atmosphere.V'] = 'atmosphere.V'; key_variables.append(parLocation['atmosphere.V'])
parLocation['atmosphere.V_gas[1]'] = 'atmosphere.V_gas[1]'; key_variables.append(parLocation['atmosphere.V_gas[1]'])
parLocation['atmosphere.V_gas[2]'] = 'atmosphere.V_gas[2]'; key_variables.append(parLocation['atmosphere.V_gas[2]'])
parLocation['atmosphere.V_gas[3]'] = 'atmosphere.V_gas[3]'; key_variables.append(parLocation['atmosphere.V_gas[3]'])
parLocation['atmosphere.V_gas[4]'] = 'atmosphere.V_gas[4]'; key_variables.append(parLocation['atmosphere.V_gas[4]'])

parLocation['feedtank.V'] = 'feedtank.V'; key_variables.append(parLocation['feedtank.V'])

# Parameter value check - especially for hysteresis to avoid runtime error
global parCheck; parCheck = []
parCheck.append("parDict['V_start'] > 0")
parCheck.append("parDict['VX_start'] >= 0")
parCheck.append("parDict['VG_start'] >= 0")

# Create list of diagrams to be plotted by simu()
global diagrams
diagrams = []

# Define standard diagrams
def newplot(title='Yeast fedbatch cultivation', plotType='TimeSeries'):
   """ Standard plot window
        title = ''
       two possible diagrams
          diagram = 'TimeSeries' default
          diagram = 'TimeSeriesExtended', 'Extended' """ 
          
   # Reset pens
   setLines()

   # Transfer of global axes to simu()
   global ax1, ax2, ax3, ax4
   global ax11, ax21, ax31, ax41, ax51, ax61, ax71
   global ax12, ax22, ax32, ax42, ax52, ax62, ax72     
   
   if plotType in ['Overview']:

      plt.figure()       
      ax11 = plt.subplot(7,2,1); ax12 = plt.subplot(7,2,2)   
      ax21 = plt.subplot(7,2,3); ax22 = plt.subplot(7,2,4)    
      ax31 = plt.subplot(7,2,5); ax32 = plt.subplot(7,2,6)   
      ax41 = plt.subplot(7,2,7); ax42 = plt.subplot(7,2,8) 
      ax51 = plt.subplot(7,2,9); ax52 = plt.subplot(7,2,10) 
      ax61 = plt.subplot(7,2,11); ax62 = plt.subplot(7,2,12) 
      ax71 = plt.subplot(7,2,13);  

      ax11.set_title(title)

      ax11.grid()
      ax11.set_ylabel('G [g/L]')

      ax21.grid()
      ax21.set_ylabel('E [g/L]')
      
      ax31.grid()
      ax31.set_ylabel('X [g/L]')

      ax41.grid()
      ax41.set_ylabel('DO [%]')

      ax51.grid()
      ax51.set_ylabel('N [rpm]')

      ax61.grid()
      ax61.set_ylabel('F [L/h]')

      ax71.grid()
      ax71.set_ylabel('V [L]')      
      ax71.set_xlabel('Time [h]')

      ax12.grid()
      ax12.set_ylabel('qG [mole/(h*g)]')

      ax22.grid()
      ax22.set_ylabel('qE [mole/(h*g)]')

      ax32.grid()
      ax32.set_ylabel('mu [1/h]')

      ax42.grid()
      ax42.set_ylabel('qO2 [mole/h,g]')

      ax52.grid()
      ax52.set_ylabel('OUR [mole/h]')    

      ax62.grid()
      ax62.set_ylabel('Q [W]')        
      ax62.set_xlabel('Time [h]')     
      
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax11.plot(sim_res['time'],sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")
      diagrams.append("ax21.plot(sim_res['time'],sim_res['bioreactor.c[3]'],color='b',linestyle=linetype)")
      diagrams.append("ax31.plot(sim_res['time'],sim_res['bioreactor.c[1]'],color='b',linestyle=linetype)")
      diagrams.append("ax41.plot(sim_res['time'],sim_res['DOsensor.out'],color='b',linestyle=linetype)")
      diagrams.append("ax41.plot(sim_res['time'],sim_res['DO_setpoint.out'],color='y',linestyle='--')")
      diagrams.append("ax51.step(sim_res['time'],sim_res['bioreactor.N'],color='c',linestyle=linetype)")
      diagrams.append("ax51.set_ylim([0,1700])")
      diagrams.append("ax61.plot(sim_res['time'],sim_res['bioreactor.inlet[1].F'],color='c',linestyle=linetype)")
      diagrams.append("ax71.plot(sim_res['time'],sim_res['bioreactor.V'],color='b',linestyle=linetype)")
      diagrams.append("ax71.plot(sim_res['time'],sim_res['bioreactor.V_tot'],color='y', linestyle='--')")
      diagrams.append("ax71.set_ylim([0, 9.0])")

      diagrams.append("ax12.plot(sim_res['time'],sim_res['bioreactor.culture.qGm'], color='r', linestyle=linetype)")
      diagrams.append("ax12.plot(sim_res['time'],sim_res['bioreactor.culture.qGr'], color='b', linestyle=linetype)")  
      diagrams.append("ax22.plot(sim_res['time'],-sim_res['bioreactor.culture.qEm'], color='r', linestyle=linetype)")
      diagrams.append("ax22.plot(sim_res['time'],sim_res['bioreactor.culture.qEr'], color='b', linestyle=linetype)")
      diagrams.append("ax32.plot(sim_res['time'],sim_res['bioreactor.culture.q[1]'],color='b',linestyle=linetype)")
      
      diagrams.append("ax42.plot(sim_res['time'],sim_res['bioreactor.culture.qO2'],color='b',linestyle=linetype)")
      diagrams.append("ax52.plot(sim_res['time'],sim_res['bioreactor.m[1]']*sim_res['bioreactor.culture.qO2'],color='b',linestyle=linetype)")
      diagrams.append("ax62.plot(sim_res['time'],sim_res['bioreactor.m[1]']*sim_res['bioreactor.culture.Qspec'],color='b',linestyle=linetype)")

   elif plotType in ['Focus DO-control']:
   
      plt.figure()
      ax1 = plt.subplot(5,1,1)
      ax2 = plt.subplot(5,1,2)
      ax3 = plt.subplot(5,1,3)
      ax4 = plt.subplot(5,1,4)

      ax1.set_title(title)    
      ax1.grid()
      ax1.set_ylabel('DO [%]')    
    
      ax2.grid()
      ax2.set_ylabel('N [rpm]') 

      ax3.grid()
      ax3.set_ylabel('OUR [mole/h]') 

      ax4.grid()
      ax4.set_ylabel('F [L/h]') 
      ax4.set_xlabel('Time [h]') 

      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax1.plot(sim_res['time'],sim_res['DOsensor.out'],color='b',linestyle=linetype)")
      diagrams.append("ax1.plot(sim_res['time'],sim_res['DO_setpoint.out'],color='r',linestyle='--')")
      diagrams.append("ax2.step(sim_res['time'],sim_res['bioreactor.N'],color='b',linestyle=linetype)")
      diagrams.append("ax2.set_ylim([0,2500])")
      diagrams.append("ax3.plot(sim_res['time'],sim_res['bioreactor.m[1]']*sim_res['bioreactor.culture.qO2'],color='b',linestyle=linetype)")
      diagrams.append("ax4.plot(sim_res['time'],sim_res['bioreactor.inlet[1].F'],color='b',linestyle=linetype)")
           
   else:
      print("Plot window type not correct")

#def eigValReactor(model):
#   """Calculate from the model the eigenvalues for the reactor"""
#   (A,B,C,D) = model.get_state_space_representation(
#                   A=True, B=False, C=False, D=False,
#                   use_structure_info=False)
#   (eigValues, eigVectors) = linalg.eig(A[0:4,0:4])
#   return eigValues

def describe(name, decimals=3):
   """Look up description of culture, media, as well as parameters and variables in the model code"""
           
   if name == 'culture':
      print('Saccharomyces cerevisae - default parameters for strain H1022')        
        
   elif name in ['broth', 'liquidphase', 'liquid-phase','media']:
      X = model_get('liquidphase.X'); 
      X_description = model_get_variable_description('liquidphase.X'); 
      X_mw = model_get('liquidphase.mw[1]')
        
      G = model_get('liquidphase.G'); 
      G_description = model_get_variable_description('liquidphase.G'); 
      G_mw = model_get('liquidphase.mw[2]')
        
      E = model_get('liquidphase.E'); 
      E_description = model_get_variable_description('liquidphase.E'); 
      E_mw = model_get('liquidphase.mw[3]')
        
      Diss_O2 = model_get('liquidphase.O2'); 
      Diss_O2_description = model_get_variable_description('liquidphase.O2'); 
      O2_mw = model_get('liquidphase.mw[4]')
        
      Diss_CO2 = model_get('liquidphase.CO2'); 
      Diss_CO2_description = model_get_variable_description('liquidphase.CO2'); 
      CO2_mw = model_get('liquidphase.mw[5]')

      print('Reactor broth substances included in the model')
      print()
      print(X_description, '  index       = ', X, '- molecular weight = ', X_mw, 'Da')
      print(G_description, 'index       = ', G, '- molecular weight = ', G_mw, 'Da')
      print(E_description, 'index       = ', E, '- molecular weight = ', E_mw, 'Da')
      print(Diss_O2_description, 'index  = ', Diss_O2, '- molecular weight = ', O2_mw, 'Da')
      print(Diss_CO2_description, 'index = ', Diss_CO2, '- molecular weight = ', CO2_mw, 'Da')

   elif name in ['gasphase', 'gas-phase']:
      N2 = model_get('gasphase.N2'); 
      N2_description = model_get_variable_description('gasphase.N2'); 
      N2_mw = model_get('gasphase.mw[1]')
        
      O2 = model_get('gasphase.O2'); 
      O2_description = model_get_variable_description('gasphase.O2'); 
      O2_mw = model_get('gasphase.mw[2]')
      
      CO2 = model_get('gasphase.CO2'); 
      CO2_description = model_get_variable_description('gasphase.CO2'); 
      CO2_mw = model_get('gasphase.mw[3]')
        
      E = model_get('gasphase.E'); 
      E_description = model_get_variable_description('gasphase.E'); 
      E_mw = model_get('gasphase.mw[4]')
         
      print('Reactor gasphase substances included in the model')
      print()
      print(N2_description, 'index  = ',N2, '- molecular weight = ', N2_mw, 'Da')
      print(O2_description, 'index      = ',O2, '- molecular weight = ', O2_mw, 'Da')
      print(CO2_description, 'index     = ',CO2, '- molecular weight = ', CO2_mw, 'Da')
      print(E_description, 'index = ',E, '- molecular weight = ', E_mw, 'Da') 
 
   elif name in ['parts']:
      describe_parts(component_list_minimum)

   elif name in ['MSL']:
      describe_MSL()
     
   else:
      describe_general(name, decimals)

#------------------------------------------------------------------------------------------------------------------
#  General code 
FMU_explore = 'FMU-explore for FMPy version 1.0.1'
#------------------------------------------------------------------------------------------------------------------

# Define function par() for parameter update
def par(parDict=parDict, parCheck=parCheck, parLocation=parLocation, *x, **x_kwarg):
   """ Set parameter values if available in the predefined dictionaryt parDict. """
   x_kwarg.update(*x)
   x_temp = {}
   for key in x_kwarg.keys():
      if key in parDict.keys():
         x_temp.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an accessible parameter - check the spelling')
   parDict.update(x_temp)
   
   parErrors = [requirement for requirement in parCheck if not(eval(requirement))]
   if not parErrors == []:
      print('Error - the following requirements do not hold:')
      for index, item in enumerate(parErrors): print(item)

# Define function init() for initial values update
def init(parDict=parDict, *x, **x_kwarg):
   """ Set initial values and the name should contain string '_start' to be accepted.
       The function can handle general parameter string location names if entered as a dictionary. """
   x_kwarg.update(*x)
   x_init={}
   for key in x_kwarg.keys():
      if '_start' in key: 
         x_init.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an initial value, use par() instead - check the spelling')
   parDict.update(x_init)

# Define fuctions similar to pyfmi model.get(), model.get_variable_descirption(), model.get_variable_unit()
def model_get(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get() but returns just a value and not a list"""
   par_var = model_description.modelVariables
   for k in range(len(par_var)):
      if par_var[k].name == parLoc:
         try:
            if (par_var[k].causality in ['local']) & (par_var[k].variability in ['constant']):
               value = float(par_var[k].start)                 
            elif par_var[k].causality in ['parameter']: 
               value = float(par_var[k].start)  
            elif par_var[k].causality in ['calculatedParameter']: 
               value = float(sim_res[par_var[k].name][0]) 
            elif par_var[k].name in start_values.keys():
               value = start_values[par_var[k].name]   
            elif par_var[k].variability == 'continuous':
               try:
                  timeSeries = sim_res[par_var[k].name]
                  value = float(timeSeries[-1])
               except (AttributeError, ValueError):
                  value = None
                  print('Variable not logged')
            else:
               value = None
         except NameError:
            print('Error: Information available after first simulation')
            value = None          
   return value
   
def model_get_variable_description(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get_variable_description() but returns just a value and not a list"""
   par_var = model_description.modelVariables
#   value = [x[1] for x in [(par_var[k].name, par_var[k].description) for k in range(len(par_var))] if parLoc in x[0]]
   value = [x.description for x in par_var if parLoc in x.name]   
   return value[0]
   
def model_get_variable_unit(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get_variable_unit() but returns just a value and not a list"""
   par_var = model_description.modelVariables
#   value = [x[1] for x in [(par_var[k].name, par_var[k].unit) for k in range(len(par_var))] if parLoc in x[0]]
   value = [x.unit for x in par_var if parLoc in x.name]
   return value[0]
      
# Define function disp() for display of initial values and parameters
def disp(name='', decimals=3, mode='short'):
   """ Display intial values and parameters in the model that include "name" and is in parLocation list.
       Note, it does not take the value from the dictionary par but from the model. """
   
   def dict_reverser(d):
      seen = set()
      return {v: k for k, v in d.items() if v not in seen or seen.add(v)}
   
   if mode in ['short']:
      k = 0
      for Location in [parLocation[k] for k in parDict.keys()]:
         if name in Location:
            if type(model_get(Location)) != np.bool_:
               print(dict_reverser(parLocation)[Location] , ':', np.round(model_get(Location),decimals))
            else:
               print(dict_reverser(parLocation)[Location] , ':', model_get(Location))               
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parDict.keys():
            if name in parName:
               if type(model_get(Location)) != np.bool_:
                  print(parName,':', np.round(model_get(parLocation[parName]),decimals))
               else: 
                  print(parName,':', model_get(parLocation[parName])[0])

   if mode in ['long','location']:
      k = 0
      for Location in [parLocation[k] for k in parDict.keys()]:
         if name in Location:
            if type(model_get(Location)) != np.bool_:       
               print(Location,':', dict_reverser(parLocation)[Location] , ':', np.round(model_get(Location),decimals))
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parDict.keys():
            if name in parName:
               if type(model_get(Location)) != np.bool_:
                  print(parLocation[parName], ':', dict_reverser(parLocation)[Location], ':', parName,':', 
                     np.round(model_get(parLocation[parName]),decimals))

# Line types
def setLines(lines=['-','--',':','-.']):
   """Set list of linetypes used in plots"""
   global linecycler
   linecycler = cycle(lines)

# Show plots from sim_res, just that
def show(diagrams=diagrams):
   """Show diagrams chosen by newplot()"""
   # Plot pen
   linetype = next(linecycler)    
   # Plot diagrams 
   for command in diagrams: eval(command)

# Define simulation
def simu(simulationTime=simulationTime, mode='Initial', options=opts_std, diagrams=diagrams):
   """Model loaded and given intial values and parameter before, and plot window also setup before."""   
   
   # Global variables
   global sim_res, prevFinalTime, stateDict, stateDictInitial, stateDictInitialLoc, start_values
   
   # Simulation flag
   simulationDone = False
   
   # Internal help function to extract variables to be stored
   def extract_variables(diagrams):
       output = []
       variables = [v for v in model_description.modelVariables if v.causality == 'local']
       for j in range(len(diagrams)):
           for k in range(len(variables)):
               if variables[k].name in diagrams[j]:
                   output.append(variables[k].name)
       return output

   # Run simulation
   if mode in ['Initial', 'initial', 'init']: 
      
      start_values = {parLocation[k]:parDict[k] for k in parDict.keys()}
      
      # Simulate
      sim_res = simulate_fmu(
         filename = fmu_model,
         validate = False,
         start_time = 0,
         stop_time = simulationTime,
         output_interval = simulationTime/options['ncp'],
         record_events = True,
         start_values = start_values,
         fmi_call_logger = None,
         output = list(set(extract_variables(diagrams) + list(stateDict.keys()) + key_variables))
      )
      
      simulationDone = True
      
   elif mode in ['Continued', 'continued', 'cont']:
      
      if prevFinalTime == 0: 
         print("Error: Simulation is first done with default mode = init'")
         
      else:         
         # Update parDictMod and create parLocationMod
         parDictRed = parDict.copy()
         parLocationRed = parLocation.copy()
         for key in parDict.keys():
            if parLocation[key] in stateDictInitial.values(): 
               del parDictRed[key]  
               del parLocationRed[key]
         parLocationMod = dict(list(parLocationRed.items()) + list(stateDictInitialLoc.items()))
   
         # Create parDictMod and parLocationMod
         parDictMod = dict(list(parDictRed.items()) + 
            [(stateDictInitial[key], stateDict[key]) for key in stateDict.keys()])      

         start_values = {parLocationMod[k]:parDictMod[k] for k in parDictMod.keys()}
  
         # Simulate
         sim_res = simulate_fmu(
            filename = fmu_model,
            validate = False,
            start_time = prevFinalTime,
            stop_time = prevFinalTime + simulationTime,
            output_interval = simulationTime/options['ncp'],
            record_events = True,
            start_values = start_values,
            fmi_call_logger = None,
            output = list(set(extract_variables(diagrams) + list(stateDict.keys()) + key_variables))
         )
      
         simulationDone = True
   else:
      
      print("Error: Simulation mode not correct")

   if simulationDone:
      
      # Plot diagrams from simulation
      linetype = next(linecycler)    
      for command in diagrams: eval(command)
   
      # Store final state values in stateDict:        
      for key in stateDict.keys(): stateDict[key] = model_get(key)  
         
      # Store time from where simulation will start next time
      prevFinalTime = sim_res['time'][-1]
      
   else:
      print('Error: No simulation done')
            
# Describe model parts of the combined system
def describe_parts(component_list=[]):
   """List all parts of the model""" 
       
   def model_component(variable_name):
      i = 0
      name = ''
      finished = False
      if not variable_name[0] == '_':
         while not finished:
            name = name + variable_name[i]
            if i == len(variable_name)-1:
                finished = True 
            elif variable_name[i+1] in ['.', '(']: 
                finished = True
            else: 
                i=i+1
      if name in ['der', 'temp_1', 'temp_2', 'temp_3', 'temp_4', 'temp_5', 'temp_6', 'temp_7']: name = ''
      return name
    
#   variables = list(model.get_model_variables().keys())
   variables = [v.name for v in model_description.modelVariables]
        
   for i in range(len(variables)):
      component = model_component(variables[i])
      if (component not in component_list) \
      & (component not in ['','BPL', 'Customer', 'today[1]', 'today[2]', 'today[3]', 'temp_2', 'temp_3']):
         component_list.append(component)
      
   print(sorted(component_list, key=str.casefold))

# Describe MSL   
def describe_MSL(flag_vendor=flag_vendor):
   """List MSL version and components used"""
   print('MSL:', MSL_usage)
 
# Describe parameters and variables in the Modelica code
def describe_general(name, decimals):
  
   if name == 'time':
      description = 'Time'
      unit = 'h'
      print(description,'[',unit,']')
      
   elif name in parLocation.keys():
      description = model_get_variable_description(parLocation[name])
      value = model_get(parLocation[name])
      try:
         unit = model_get_variable_unit(parLocation[name])
      except FMUException:
         unit =''
      if unit =='':
         if type(value) != np.bool_:
            print(description, ':', np.round(value, decimals))
         else:
            print(description, ':', value)            
      else:
        print(description, ':', np.round(value, decimals), '[',unit,']')
                  
   else:
      description = model_get_variable_description(name)
      value = model_get(name)
      try:
         unit = model_get_variable_unit(name)
      except FMUException:
         unit =''
      if unit =='':
         if type(value) != np.bool_:
            print(description, ':', np.round(value, decimals))
         else:
            print(description, ':', value)     
      else:
         print(description, ':', np.round(value, decimals), '[',unit,']')

# Plot process diagram
def process_diagram(fmu_model=fmu_model, fmu_process_diagram=fmu_process_diagram):   
   try:
       process_diagram = zipfile.ZipFile(fmu_model, 'r').open('documentation/processDiagram.png')
   except KeyError:
       print('No processDiagram.png file in the FMU, but try the file on disk.')
       process_diagram = fmu_process_diagram
   try:
       plt.imshow(img.imread(process_diagram))
       plt.axis('off')
       plt.show()
   except FileNotFoundError:
       print('And no such file on disk either')
         
# Describe framework
def BPL_info():
   print()
   print('Model for bioreactor has been setup. Key commands:')
   print(' - par()       - change of parameters and initial values')
   print(' - init()      - change initial values only')
   print(' - simu()      - simulate and plot')
   print(' - newplot()   - make a new plot')
   print(' - show()      - show plot from previous simulation')
   print(' - disp()      - display parameters and initial values from the last simulation')
   print(' - describe()  - describe culture, broth, parameters, variables with values/units')
   print()
   print('Note that both disp() and describe() takes values from the last simulation')
   print('and the command process_diagram() brings up the main configuration')
   print()
   print('Brief information about a command by help(), eg help(simu)') 
   print('Key system information is listed with the command system_info()')

def system_info():
   """Print system information"""
#   FMU_type = model.__class__.__name__
   constants = [v for v in model_description.modelVariables if v.causality == 'local']
   
   print()
   print('System information')
   print(' -OS:', platform.system())
   print(' -Python:', platform.python_version())
   try:
       scipy_ver = scipy.__version__
       print(' -Scipy:',scipy_ver)
   except NameError:
       print(' -Scipy: not installed in the notebook')
   print(' -FMPy:', version('fmpy'))
   print(' -FMU by:', read_model_description(fmu_model).generationTool)
   print(' -FMI:', read_model_description(fmu_model).fmiVersion)
   if model_description.modelExchange is None:
      print(' -Type: CS')
   else:
      print(' -Type: ME')
   print(' -Name:', read_model_description(fmu_model).modelName)
   print(' -Generated:', read_model_description(fmu_model).generationDateAndTime)
   print(' -MSL:', MSL_version)    
   print(' -Description:', BPL_version)   
   print(' -Interaction:', FMU_explore)
   
#------------------------------------------------------------------------------------------------------------------
#  Startup
#------------------------------------------------------------------------------------------------------------------

BPL_info()