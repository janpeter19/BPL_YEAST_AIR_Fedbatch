# setup applicateion data BPL_YEAST_AIR_Fedbatch 
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2026-08-24 - Created
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Framework
#------------------------------------------------------------------------------------------------------------------

# Setup framework
import sys
import platform
import locale
import matplotlib.pyplot as plt
from pyfmi import load_fmu

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
      
#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
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
   BPL_version = 'Bioprocess Library version 2.3.2' 
else:    
   print('There is no FMU for this platform')

# Simulation time
simulationTime = 20.0
prevFinalTime = 0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Create stateValue that later will be used to store final state and used for initialization in 'cont':
stateValue =  {}
stateValue = model.get_states_list()
stateValue.update(timeDiscreteStates)

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture', 'bioreactor.gas_liquid_transfer']

# Provide process diagram on disk
fmu_process_diagram ='BPL_YEAST_AIR_Fedbatch_DOcontrol_process_diagram_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateValue, parValue, parLocation, parCheck, diagrams, ax, lines
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
#parValue['t_regStart'] = 0.0
#parValue['samplePeriod'] = 0.1
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
#parLocation['t_regStart'] = 'PIDreg.t_regStart'
#parLocation['samplePeriod'] = 'PIDreg.samplePeriod'
parLocation['K'] = 'PIDreg.K'
parLocation['Ti'] = 'PIDreg.Ti'
parLocation['Td'] = 'PIDreg.Td'
parLocation['Nd'] = 'PIDreg.Nd'
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
