# Figure - Simulation of fedbatch reactor with yeast
#          with functions added to facilitate explorative simulation work
#
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2020-03-06 - Python 2 script for compilation 
#            - Import platform and locale (for later use with OpenModelica FMU)
#            - Added system print of system information
#            - Improved check of platform to adapt code for Windows/Linux in dialog
#            - Move plt.show() from newplot() to simu()
#            - Skip initialization of stateDict  
#            - Change newplot and simu using objectoriented diagrams
#            - Simplify handling of simulation data and eliminate Trajectory
# 2020-03-06 - Tested with JModelica 2.14 and seems ok
# 2020-03-06 - Tested with OCT 1.14.1 and changed names of FMUs
# 2020-03-16 - Indluced in system_info() information if FMU is ME or CS
#------------------------------------------------------------------------------------------------------------------
# 2020-07-20 - Adapted for BP6a
# 2020-07-20 - Corrected the handling of stateDict and simu('cont') 
# 2020-07-21 - Further improved handling of stateDict, use model.get_states_list
# 2020-07-21 - Impoved safety around stateDict and simu('cont')
# 2020-07-27 - Introduce choice of Linux FMU - JModelica or OpenModelica
# 2020-10-02 - Updated for new BP6a with development from BP6c
# 2020-10-04 - Adapted for DO-control
# 2020-10-08 - Polished description of the broth and gas phase
# 2020-10-08 - Extended newplot() and simu() to include DO and N
# 2020-10-10 - Harmonize PIDregType and PIregDiscreteType for interchange
# 2020-10-12 - Decided now to keep the discrete tine controller and samplePeriod
# 2020-11-22 - Adapted to ReactorType with n_inlets, n_outlets and n_ports
# 2021-02-04 - Adjust describe() for change to liquidphase
#------------------------------------------------------------------------------------------------------------------
# 2021-02-10 - Adapt for BPL_v2
# 2021-03-20 - Adapt for BPL ver 2.0.3
#------------------------------------------------------------------------------------------------------------------
# 2021-05-12 - Adapted for BPL ver 2.0.5 and updated BPL interface
# 2021-05-19 - Default value given for disp()
# 2021-05-22 - Improved decribe() to include short names in par[]
# 2021-05-25 - Change to parDict and parLocation and use par() for parDict.update()
# 2021-06-06 - Polish and include diagram "extended"
# 2021-06-24 - Moved to BPL_dev and test on integration of newplot() and simu()
# 2021-07-29 - Modifed describe() and lifted out the general part and put that in FMU-explore part instead
# 2021-07-30 - Introduced describe_parts()
# 2021-07-31 - Corrected funcetion disp() to handle number of displayed decimals and improved describe_parts()
# 2021-08-05 - Improved describe() and describe_general() to also handle decimals
# 2021-08-07 - Adapted for teaching purpose
# 2021-10-28 - Adapted for FMU-explore 0.8.5
# 2022-01-19 - Updated with FMU-explore 0.8.7
# 2022-01-28 - Updated with FMU-explore 0.8.8
# 2022-03-26 - Updated to FMU-explore 0.9.0 - model.reset(), and par(), init()
# 2022-08-21 - Update with FMU-explore 0.9.2 and also include a slightly modified plot
# 2022-08-26 - Updated newplot() with new diagram for DO-control
# 2022-08-26 - Test with Linux-FMU
# 2022-08-26 - Take away not necessary plot-types fron newplot() amd only keep waht used for Colab-demo
# 2022-09-01 - Included parameter alpha_O2 and also took away mw[O2] from parLocation list
# 2022-09-13 - Updated for FMU-explore 0.9.3
# 2022-09-17 - Updated for FMU-explore 0.9.5
# 2023-02-09 - Updated to FMU-explore 0.9.6e
# 2023-02-13 - Consolidate FMU-explore to 0.9.6 and means parCheck and par() udpate and simu() with opts as arg
# 2023-02-20 - Updatation to updated BPL.Control and block VarLimPID used
# 2023-02-22 - Adjusted parDict, parLocation and simu('cont')
# 2023-02-23 - Added Kla_O2 and Kla_CO2 to the derived parameters that can be reached by describe()
# 2023-03-28 - Update FMU-explore 0.9.7
# 2023-04-21 - Compiled for Ubuntu 20.04 and changed BPL_version
# 2023-05-31 - Adjusted to from importlib.meetadata import version
# 2024-01-22 - Update FMU-explore 0.9.9 icluding function process_diagram() although not used since no GUI gasphase
# 2024-05-08 - Look through early part and call it all FMU-explore 1.0.0
# 2024-05-20 - Updated the OpenModelica version to 1.23.0-dev
# 2024-07-03 - Updated with adjustments of names (connections to FixValue) due to new FMU with BPL 2.2.1
# 2024-10-13 - Updated information aboout BPL 2.2.2 - GUI
# 2024-11-06 - Updated for BPL 2.3.0
# 2025-02-14 - Change from qO2lim to qO2max and use KsO2 and also introduce VO2_start and VCO2_start
# 2025-06-28 - Updated for BPL 2.3.1 Linux MSL 4.1.0 and bring back parameters Td and D_start
# 2025-07-02 - Change PIreg to PIDreg
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

from pyfmi import load_fmu
from pyfmi.fmi import FMUException

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
   flag_vendor = 'JM'
   flag_type = 'CS'
   fmu_model ='BPL_YEAST_AIR_Fedbatch_DOcontrol_windows_jm_cs.fmu'        
   model = load_fmu(fmu_model, log_level=0)  
elif platform.system() == 'Linux':
   flag_vendor = 'OM'
   flag_type = 'ME'
   if flag_vendor in ['OM','om']:
      print('Linux - run FMU pre-compiled OpenModelica') 
      if flag_type in ['CS','cs']:         
         fmu_model ='BPL_YEAST_AIR_Fedbatch_DOcontrol_linux_om_cs.fmu'    
         model = load_fmu(fmu_model, log_level=0) 
      if flag_type in ['ME','me']:         
         fmu_model ='BPL_YEAST_AIR_Fedbatch_DOcontrol_linux_om_me.fmu'    
         model = load_fmu(fmu_model, log_level=0)
   else:    
      print('There is no FMU for this platform')

# Provide various opts-profiles
if flag_type in ['CS', 'cs']:
   opts_std = model.simulate_options()
   opts_std['silent_mode'] = True
   opts_std['ncp'] = 500 
   opts_std['result_handling'] = 'binary'     
elif flag_type in ['ME', 'me']:
   opts_std = model.simulate_options()
   opts_std["CVode_options"]["verbosity"] = 50 
   opts_std['ncp'] = 500 
   opts_std['result_handling'] = 'binary'  
else:    
   print('There is no FMU for this platform')
  
# Provide various MSL and BPL versions
if flag_vendor in ['JM', 'jm']:
   MSL_usage = model.get('MSL.usage')[0]
   MSL_version = model.get('MSL.version')[0]
   BPL_version = model.get('BPL.version')[0]
elif flag_vendor in ['OM', 'om']:
   MSL_usage = '4.1.0 - used components: RealInput, RealOutput, LimPID-components' 
   MSL_version = '4.1.0'
   BPL_version = 'Bioprocess Library version 2.3.1' 
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
fmu_process_diagram ='BPL_YEAST_AIR_Fedbatch_DOcontrol_process_diagram_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateDict, parDict, diagrams, newplot(), describe()
#------------------------------------------------------------------------------------------------------------------

# Create stateDict that later will be used to store final state and used for initialization in 'cont':
global stateDict; stateDict =  {}
stateDict = model.get_states_list()
stateDict.update(timeDiscreteStates)

# Create dictionaries parDictLocation[] and parLocation[]
global parDict; parDict = {}
parDict['V_start'] = 4.5
parDict['VX_start'] = 4.5*1.0
parDict['VG_start'] = 4.5*5.0
parDict['VE_start'] = 0.0
parDict['V_diss_O2_start'] = 0.0067
parDict['V_diss_CO2_start'] = 1.25

parDict['V_tot'] = 8.0
parDict['V_gas_N2_start'] = 2.4
parDict['V_gas_O2_start'] = 0.6
parDict['V_gas_CO2_start'] = 0

parDict['qGmax'] = 20.0e-3
parDict['Ks'] = 10.0e-3
parDict['qO2max'] = 6.9e-3
parDict['KsO2'] = 1.0e-5

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
#parDict['t_regStart'] = 0.0
#parDict['samplePeriod'] = 0.1
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
parLocation['V_diss_O2_start'] = 'bioreactor.m_start[4]'
parLocation['V_diss_CO2_start'] = 'bioreactor.m_start[5]'

parLocation['V_tot'] = 'bioreactor.V_tot'
parLocation['V_gas_N2_start'] = 'bioreactor.V_gas_start[1]'
parLocation['V_gas_O2_start'] = 'bioreactor.V_gas_start[2]'
parLocation['V_gas_CO2_start'] = 'bioreactor.V_gas_start[3]'

parLocation['qGmax'] = 'bioreactor.culture.qGmax' 
parLocation['Ks'] = 'bioreactor.culture.Ks' 
parLocation['qO2max'] = 'bioreactor.culture.qO2max' 
parLocation['KsO2'] = 'bioreactor.culture.KsO2'

parLocation['alpha_O2'] = 'bioreactor.gas_liquid_transfer.alpha_O2'

parLocation['feedtank_V_start'] = 'feedtank.V_start'
parLocation['G_in'] = 'feedtank.c_in[2]'
parLocation['F_start'] = 'dosagescheme.F_start'
parLocation['mu_feed'] = 'dosagescheme.mu_feed'
parLocation['t_startExp'] = 'dosagescheme.t_startExp'
parLocation['F_startExp'] = 'dosagescheme.F_startExp'
parLocation['F_max'] = 'dosagescheme.F_max'

parLocation['airFlow_setpoint'] = 'airFlow_setpoint.value'

parLocation['DO_setpoint'] = 'DO_setpoint.value'
parLocation['DO_sensor_x_start'] = 'DOsensor.x_start'
#parLocation['t_regStart'] = 'PIDreg.t_regStart'
#parLocation['samplePeriod'] = 'PIDreg.samplePeriod'
parLocation['K'] = 'PIDreg.K'
parLocation['Ti'] = 'PIDreg.Ti'
parLocation['Td'] = 'PIDreg.Td'
parLocation['I_start'] = 'PIDreg.I_start'
parLocation['D_start'] = 'PIDreg.D_start'
parLocation['N_low'] = 'N_low.value'
parLocation['N_high'] = 'N_high.value'

# Extended list of parameters and variables only for display and not change
parLocation['mu'] = 'bioreactor.culture.mu'
parLocation['Kla_O2'] = 'bioreactor.gas_liquid_transfer.Kla_O2'
parLocation['Kla_CO2'] = 'bioreactor.gas_liquid_transfer.Kla_CO2'
parLocation['qO2lim'] = 'bioreactor.culture.qO2lim'

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
      diagrams.append("ax11.plot(t,sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")
      diagrams.append("ax21.plot(t,sim_res['bioreactor.c[3]'],color='b',linestyle=linetype)")
      diagrams.append("ax31.plot(t,sim_res['bioreactor.c[1]'],color='b',linestyle=linetype)")
      diagrams.append("ax41.plot(t,sim_res['DOsensor.out'],color='b',linestyle=linetype)")
      diagrams.append("ax41.plot(t,sim_res['DO_setpoint.out'],color='y',linestyle='--')")
      diagrams.append("ax51.step(t,sim_res['bioreactor.N'],color='c',linestyle=linetype)")
      diagrams.append("ax51.set_ylim([0,1700])")
      diagrams.append("ax61.plot(t,sim_res['bioreactor.inlet[1].F'],color='c',linestyle=linetype)")
      diagrams.append("ax71.plot(t,sim_res['bioreactor.V'],color='b',linestyle=linetype)")
      diagrams.append("ax71.plot(t,sim_res['bioreactor.V_tot'],color='y', linestyle='--')")
      diagrams.append("ax71.set_ylim([0, 9.0])")

      diagrams.append("ax12.plot(t,sim_res['bioreactor.culture.qGm'], color='r', linestyle=linetype)")
      diagrams.append("ax12.plot(t,sim_res['bioreactor.culture.qGr'], color='b', linestyle=linetype)")  
      diagrams.append("ax22.plot(t,-sim_res['bioreactor.culture.qEm'], color='r', linestyle=linetype)")
      diagrams.append("ax22.plot(t,sim_res['bioreactor.culture.qEr'], color='b', linestyle=linetype)")
      diagrams.append("ax32.plot(t,sim_res['bioreactor.culture.q[1]'],color='b',linestyle=linetype)")
      
      diagrams.append("ax42.plot(t,sim_res['bioreactor.culture.qO2'],color='b',linestyle=linetype)")
      diagrams.append("ax52.plot(t,sim_res['bioreactor.m[1]']*sim_res['bioreactor.culture.qO2'],color='b',linestyle=linetype)")
      diagrams.append("ax62.plot(t,sim_res['bioreactor.m[1]']*sim_res['bioreactor.culture.Qspec'],color='b',linestyle=linetype)")

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
      diagrams.append("ax1.plot(t,sim_res['DOsensor.out'],color='b',linestyle=linetype)")
      diagrams.append("ax1.plot(t,sim_res['DO_setpoint.out'],color='r',linestyle='--')")
      diagrams.append("ax2.step(t,sim_res['bioreactor.N'],color='b',linestyle=linetype)")
      diagrams.append("ax2.set_ylim([0,2500])")
      diagrams.append("ax3.plot(t,sim_res['bioreactor.m[1]']*sim_res['bioreactor.culture.qO2'],color='b',linestyle=linetype)")
      diagrams.append("ax4.plot(t,sim_res['bioreactor.inlet[1].F'],color='b',linestyle=linetype)")
           
   else:
      print("Plot window type not correct")

def eigValReactor(model):
   """Calculate from the model the eigenvalues for the reactor"""
   (A,B,C,D) = model.get_state_space_representation(
                   A=True, B=False, C=False, D=False,
                   use_structure_info=False)
   (eigValues, eigVectors) = linalg.eig(A[0:4,0:4])
   return eigValues

def describe(name, decimals=3):
   """Look up description of culture, media, as well as parameters and variables in the model code"""
           
   if name == 'culture':
      print('Saccharomyces cerevisae - default parameters for strain H1022')        
        
   elif name in ['broth', 'liquidphase', 'liquid-phase','media']:
      X = model.get('liquidphase.X')[0]; 
      X_description = model.get_variable_description('liquidphase.X'); 
      X_mw = model.get('liquidphase.mw[1]')[0]
        
      G = model.get('liquidphase.G')[0]; 
      G_description = model.get_variable_description('liquidphase.G'); 
      G_mw = model.get('liquidphase.mw[2]')[0]
        
      E = model.get('liquidphase.E')[0]; 
      E_description = model.get_variable_description('liquidphase.E'); 
      E_mw = model.get('liquidphase.mw[3]')[0]
        
      Diss_O2 = model.get('liquidphase.O2')[0]; 
      Diss_O2_description = model.get_variable_description('liquidphase.O2'); 
      O2_mw = model.get('liquidphase.mw[4]')[0]
        
      Diss_CO2 = model.get('liquidphase.CO2')[0]; 
      Diss_CO2_description = model.get_variable_description('liquidphase.CO2'); 
      CO2_mw = model.get('liquidphase.mw[5]')[0]

      print('Reactor broth substances included in the model')
      print()
      print(X_description, '  index       = ', X, '- molecular weight = ', X_mw, 'Da')
      print(G_description, 'index       = ', G, '- molecular weight = ', G_mw, 'Da')
      print(E_description, 'index       = ', E, '- molecular weight = ', E_mw, 'Da')
      print(Diss_O2_description, 'index  = ', Diss_O2, '- molecular weight = ', O2_mw, 'Da')
      print(Diss_CO2_description, 'index = ', Diss_CO2, '- molecular weight = ', CO2_mw, 'Da')

   elif name in ['gasphase', 'gas-phase']:
      N2 = model.get('gasphase.N2')[0]; 
      N2_description = model.get_variable_description('gasphase.N2'); 
      N2_mw = model.get('gasphase.mw[1]')[0]
        
      O2 = model.get('gasphase.O2')[0]; 
      O2_description = model.get_variable_description('gasphase.O2'); 
      O2_mw = model.get('gasphase.mw[2]')[0]
      
      CO2 = model.get('gasphase.CO2')[0]; 
      CO2_description = model.get_variable_description('gasphase.CO2'); 
      CO2_mw = model.get('gasphase.mw[3]')[0]
        
      E = model.get('gasphase.E')[0]; 
      E_description = model.get_variable_description('gasphase.E'); 
      E_mw = model.get('gasphase.mw[4]')[0]
         
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
FMU_explore = 'FMU-explore version 1.0.0'
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
   
# Define function disp() for display of initial values and parameters
def dict_reverser(d):
   seen = set()
   return {v: k for k, v in d.items() if v not in seen or seen.add(v)}
   
def disp(name='', decimals=3, mode='short'):
   """ Display intial values and parameters in the model that include "name" and is in parLocation list.
       Note, it does not take the value from the dictionary par but from the model. """
   global parLocation, model
   
   if mode in ['short']:
      k = 0
      for Location in [parLocation[k] for k in parDict.keys()]:
         if name in Location:
            if type(model.get(Location)[0]) != np.bool_:
               print(dict_reverser(parLocation)[Location] , ':', np.round(model.get(Location)[0],decimals))
            else:
               print(dict_reverser(parLocation)[Location] , ':', model.get(Location)[0])               
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parDict.keys():
            if name in parName:
               if type(model.get(Location)[0]) != np.bool_:
                  print(parName,':', np.round(model.get(parLocation[parName])[0],decimals))
               else: 
                  print(parName,':', model.get(parLocation[parName])[0])
   if mode in ['long','location']:
      k = 0
      for Location in [parLocation[k] for k in parDict.keys()]:
         if name in Location:
            if type(model.get(Location)[0]) != np.bool_:       
               print(Location,':', dict_reverser(parLocation)[Location] , ':', np.round(model.get(Location)[0],decimals))
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parDict.keys():
            if name in parName:
               if type(model.get(Location)[0]) != np.bool_:
                  print(parLocation[parName], ':', dict_reverser(parLocation)[Location], ':', parName,':', 
                     np.round(model.get(parLocation[parName])[0],decimals))

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

# Simulation
def simu(simulationTimeLocal=simulationTime, mode='Initial', options=opts_std, \
         diagrams=diagrams,timeDiscreteStates=timeDiscreteStates):         
   """Model loaded and given intial values and parameter before,
      and plot window also setup before."""
    
   # Global variables
   global model, parDict, stateDict, prevFinalTime, simulationTime, sim_res, t
   
   # Simulation flag
   simulationDone = False
   
   # Transfer of argument to global variable
   simulationTime = simulationTimeLocal 
      
   # Check parDict
   value_missing = 0
   for key in parDict.keys():
      if parDict[key] in [np.nan, None, '']:
         print('Value missing:', key)
         value_missing =+1
   if value_missing>0: return
         
   # Load model
   if model is None:
      model = load_fmu(fmu_model) 
   model.reset()
      
   # Run simulation
   if mode in ['Initial', 'initial', 'init']:
      # Set parameters and intial state values:
      for key in parDict.keys():
         model.set(parLocation[key],parDict[key])   
      # Simulate
      sim_res = model.simulate(final_time=simulationTime, options=options)  
      simulationDone = True
   elif mode in ['Continued', 'continued', 'cont']:

      if prevFinalTime == 0: 
         print("Error: Simulation is first done with default mode = init'")      
      else:
         
         # Set parameters and intial state values:
         for key in parDict.keys():
            model.set(parLocation[key],parDict[key])                

         for key in stateDict.keys():
            if not key[-1] == ']':
               if key[-3:] == 'I.y': 
                  model.set(key[:-10]+'I_start', stateDict[key]) 
               elif key[-3:] == 'D.x': 
                  model.set(key[:-10]+'D_start', stateDict[key]) 
               else:
                  model.set(key+'_start', stateDict[key])
            elif key[-3] == '[':
               model.set(key[:-3]+'_start'+key[-3:], stateDict[key]) 
            elif key[-4] == '[':
               model.set(key[:-4]+'_start'+key[-4:], stateDict[key]) 
            elif key[-5] == '[':
               model.set(key[:-5]+'_start'+key[-5:], stateDict[key]) 
            else:
               print('The state vecotr has more than 1000 states')
               break

         # Simulate
         sim_res = model.simulate(start_time=prevFinalTime,
                                 final_time=prevFinalTime + simulationTime,
                                 options=options) 
         simulationDone = True             
   else:
      print("Simulation mode not correct")

   if simulationDone:
    
      # Extract data
      t = sim_res['time']
 
      # Plot diagrams
      linetype = next(linecycler)    
      for command in diagrams: eval(command)
            
      # Store final state values stateDict:
      for key in list(stateDict.keys()): stateDict[key] = model.get(key)[0]        

      # Store time from where simulation will start next time
      prevFinalTime = model.time
   
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
    
   variables = list(model.get_model_variables().keys())
        
   for i in range(len(variables)):
      component = model_component(variables[i])
      if (component not in component_list) \
      & (component not in ['','BPL', 'Customer', 'today[1]', 'today[2]', 'today[3]', 'temp_2', 'temp_3']):
         component_list.append(component)
      
   print(sorted(component_list, key=str.casefold))
   
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
      description = model.get_variable_description(parLocation[name])
      value = model.get(parLocation[name])[0]
      try:
         unit = model.get_variable_unit(parLocation[name])
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
      description = model.get_variable_description(name)
      value = model.get(name)[0]
      try:
         unit = model.get_variable_unit(name)
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
   print('Model for the process has been setup. Key commands:')
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
   FMU_type = model.__class__.__name__
   print()
   print('System information')
   print(' -OS:', platform.system())
   print(' -Python:', platform.python_version())
   try:
       scipy_ver = scipy.__version__
       print(' -Scipy:',scipy_ver)
   except NameError:
       print(' -Scipy: not installed in the notebook')
   print(' -PyFMI:', version('pyfmi'))
   print(' -FMU by:', model.get_generation_tool())
   print(' -FMI:', model.get_version())
   print(' -Type:', FMU_type)
   print(' -Name:', model.get_name())
   print(' -Generated:', model.get_generation_date_and_time())
   print(' -MSL:', MSL_version)    
   print(' -Description:', BPL_version)   
   print(' -Interaction:', FMU_explore)
   
#------------------------------------------------------------------------------------------------------------------
#  Startup
#------------------------------------------------------------------------------------------------------------------

BPL_info()