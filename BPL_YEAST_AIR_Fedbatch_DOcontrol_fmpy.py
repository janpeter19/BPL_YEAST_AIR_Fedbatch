# setup data YEAST_AIR_Fedbatch_fmpy 
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2026-09-08 - Created
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
from fmpy import simulate_fmu
from fmpy import read_model_description

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
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
      print('Linux - run FMU pre-compiled OpenModelica') 
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
   opts_std = {'NCP': 500}
elif flag_type in ['ME', 'me']:
   opts_std = {'NCP': 500}
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
   MSL_usage = '4.1.0 - used components: RealInput, RealOutput, LimPID-components' 
   MSL_version = '4.1.0'
   BPL_version = 'Bioprocess Library version 2.3.2' 
else:    
   print('There is no FMU for this platform')
   
#------------------------------------------------------------------------------------------------------------------

# Simulation time
simulationTime = 20.0
prevFinalTime = 0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Create stateValue that later will be used to store final state and used for initialization in 'cont':
stateValue =  {}
stateValue = {variable.derivative.name:None for variable in model_description.modelVariables \
                                            if variable.derivative is not None}
stateValue.update(timeDiscreteStates) 

stateValueInitial = {}
for key in stateValue.keys():
    if not key[-1] == ']':
         if key[-3:] == 'I.y':
            stateValueInitial[key] = key[:-10]+'I_start'
         elif key[-3:] == 'D.x':
            stateValueInitial[key] = key[:-10]+'D_start'
         else:
            stateValueInitial[key] = key+'_start'
    elif key[-3] == '[':
        stateValueInitial[key] = key[:-3]+'_start'+key[-3:]
    elif key[-4] == '[':
        stateValueInitial[key] = key[:-4]+'_start'+key[-4:]
    elif key[-5] == '[':
        stateValueInitial[key] = key[:-5]+'_start'+key[-5:] 
    else:
        print('The state vector has more than 1000 states')
        break

stateValueInitialLoc = {}
for value in stateValueInitial.values(): stateValueInitialLoc[value] = value

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture', 'bioreactor.gas_liquid_transfer']

# Provide process diagram on disk
fmu_process_diagram ='BPL_YEAST_AIR_Fedbatch_DOcontrol_process_diagram_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateValue, parValue, parLocation, parCheck, diagrams, ax
#------------------------------------------------------------------------------------------------------------------

# Create dictionaries parValue[] and parLocation[]
parValue = {}
parValue['V_start'] = 4.5
parValue['VX_start'] = 4.5*1.0
parValue['VG_start'] = 4.5*5.0
parValue['VE_start'] = 0.0
parValue['V_diss_O2_start'] = 0.0067
parValue['V_diss_CO2_start'] = 1.25

parValue['V_tot'] = 8.0
parValue['V_gas_N2_start'] = 2.4
parValue['V_gas_O2_start'] = 0.6
parValue['V_gas_CO2_start'] = 0

parValue['qGmax'] = 20.0e-3
parValue['Ks'] = 10.0e-3
parValue['qO2max'] = 6.9e-3
parValue['KsO2'] = 1.0e-5

parValue['alpha_O2'] = 1.0

parValue['feedtank_V_start'] = 50.0
parValue['G_in'] = 500.0
parValue['F_start'] = 0.0
parValue['mu_feed'] = 0.10
parValue['t_startExp'] = 3.0
parValue['F_startExp'] = 0.00133
parValue['F_max'] = 0.3

parValue['airFlow_setpoint'] = 120.0

parValue['DO_setpoint'] = 40.0
parValue['DO_sensor_x_start'] = 87.0
parValue['K'] = 10.0
parValue['Ti'] = 0.5
parValue['Td'] = 0.0
parValue['Nd'] = 3.0
parValue['I_start'] = 0
parValue['D_start'] = 0.0
parValue['N_low'] = 500
parValue['N_high'] = 2000

parLocation = {}
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
parLocation['K'] = 'PIDreg.K'
parLocation['Ti'] = 'PIDreg.Ti'
parLocation['Td'] = 'PIDreg.Td'
parLocation['Nd'] = 'PIDreg.Nd'
parLocation['I_start'] = 'PIDreg.I_start'
parLocation['D_start'] = 'PIDreg.D_start'
parLocation['N_low'] = 'N_low.value'
parLocation['N_high'] = 'N_high.value'

# Extended list of parameters and variables only for describe and not change
keyVariables = []
parLocation['mu'] = 'bioreactor.culture.mu'; keyVariables.append(parLocation['mu'])
parLocation['qO2'] = 'bioreactor.culture.qO2'; keyVariables.append(parLocation['qO2'])
parLocation['qO2lim'] = 'bioreactor.culture.qO2lim'; keyVariables.append(parLocation['qO2lim'])

parLocation['Kla_O2'] = 'bioreactor.gas_liquid_transfer.Kla_O2'; keyVariables.append(parLocation['Kla_O2'])
parLocation['Kla_CO2'] = 'bioreactor.gas_liquid_transfer.Kla_CO2'; keyVariables.append(parLocation['Kla_CO2'])

parLocation['V_tot'] = 'bioreactor.V_tot'; keyVariables.append(parLocation['V_tot'])
parLocation['bioreactor.V'] = 'bioreactor.V'; keyVariables.append(parLocation['bioreactor.V'])
parLocation['bioreactor.m[1]'] = 'bioreactor.m[1]'; keyVariables.append(parLocation['bioreactor.m[1]'])
parLocation['bioreactor.m[2]'] = 'bioreactor.m[2]'; keyVariables.append(parLocation['bioreactor.m[2]'])
parLocation['bioreactor.m[3]'] = 'bioreactor.m[3]'; keyVariables.append(parLocation['bioreactor.m[3]'])
parLocation['bioreactor.m[4]'] = 'bioreactor.m[4]'; keyVariables.append(parLocation['bioreactor.m[4]'])
parLocation['bioreactor.m[5]'] = 'bioreactor.m[5]'; keyVariables.append(parLocation['bioreactor.m[5]'])

parLocation['bioreactor.V_gasphase'] = 'bioreactor.V_gasphase'; keyVariables.append(parLocation['bioreactor.V_gasphase'])
parLocation['bioreactor.V_gas[1]'] = 'bioreactor.V_gas[1]'; keyVariables.append(parLocation['bioreactor.V_gas[1]'])
parLocation['bioreactor.V_gas[2]'] = 'bioreactor.V_gas[2]'; keyVariables.append(parLocation['bioreactor.V_gas[2]'])
parLocation['bioreactor.V_gas[3]'] = 'bioreactor.V_gas[3]'; keyVariables.append(parLocation['bioreactor.V_gas[3]'])
parLocation['bioreactor.V_gas[4]'] = 'bioreactor.V_gas[4]'; keyVariables.append(parLocation['bioreactor.V_gas[4]'])

parLocation['DO_setpoint.out'] = 'DO_setpoint.out'; keyVariables.append(parLocation['DO_setpoint.out'])
parLocation['DOsensor.x'] = 'DOsensor.x'; keyVariables.append(parLocation['DOsensor.x'])
parLocation['PIDreg.limPID.D.x'] = 'PIDreg.limPID.D.x'; keyVariables.append(parLocation['PIDreg.limPID.D.x'])
parLocation['PIDreg.limPID.I.y'] = 'PIDreg.limPID.I.y'; keyVariables.append(parLocation['PIDreg.limPID.I.y'])

parLocation['airtube.V'] = 'airtube.V'; keyVariables.append(parLocation['airtube.V'])
parLocation['atmosphere.V'] = 'atmosphere.V'; keyVariables.append(parLocation['atmosphere.V'])
parLocation['atmosphere.V_gas[1]'] = 'atmosphere.V_gas[1]'; keyVariables.append(parLocation['atmosphere.V_gas[1]'])
parLocation['atmosphere.V_gas[2]'] = 'atmosphere.V_gas[2]'; keyVariables.append(parLocation['atmosphere.V_gas[2]'])
parLocation['atmosphere.V_gas[3]'] = 'atmosphere.V_gas[3]'; keyVariables.append(parLocation['atmosphere.V_gas[3]'])
parLocation['atmosphere.V_gas[4]'] = 'atmosphere.V_gas[4]'; keyVariables.append(parLocation['atmosphere.V_gas[4]'])

parLocation['feedtank.V'] = 'feedtank.V'; keyVariables.append(parLocation['feedtank.V'])

# Parameter value check - especially for hysteresis to avoid runtime error
parCheck = []
parCheck.append("parValue['V_start'] > 0")
parCheck.append("parValue['VX_start'] >= 0")
parCheck.append("parValue['VG_start'] >= 0")

# Create list of diagrams to be plotted by simu()
diagrams = []

# Create an empty list axes to be defined in newplot() and plotted by simu() or show()
ax = []

# Create list of pens for the diagrams
lines = ['-','--',':','-.']