# setup application functions BPL_YEAST_AIR_Fedbatch, dependent on previous import of functions from fmu_explore 
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2026-08-24 - Created
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Specific application functions: newplot(), describe()
#------------------------------------------------------------------------------------------------------------------

# Define standard diagrams
def newplot(title='Yeast fedbatch cultivation', plotType='TimeSeries'):
   """ Standard plot window
        title = ''
       two possible diagrams
          diagram = 'TimeSeries' default
          diagram = 'TimeSeriesExtended', 'Extended' """ 
          
   # Reset pens
   resetPen()
   
   if plotType in ['Overview']:
   
      ax11 = plt.subplot(7,2,1); ax12 = plt.subplot(7,2,2)   
      ax21 = plt.subplot(7,2,3); ax22 = plt.subplot(7,2,4)    
      ax31 = plt.subplot(7,2,5); ax32 = plt.subplot(7,2,6)   
      ax41 = plt.subplot(7,2,7); ax42 = plt.subplot(7,2,8) 
      ax51 = plt.subplot(7,2,9); ax52 = plt.subplot(7,2,10) 
      ax61 = plt.subplot(7,2,11); ax62 = plt.subplot(7,2,12) 
      ax71 = plt.subplot(7,2,13);  
      
      ax.clear()
      ax.append(ax11) #  0
      ax.append(ax12) #  1
      ax.append(ax21) #  2
      ax.append(ax22) #  3
      ax.append(ax31) #  4
      ax.append(ax32) #  5
      ax.append(ax41) #  6
      ax.append(ax42) #  7
      ax.append(ax51) #  8
      ax.append(ax52) #  9
      ax.append(ax61) # 10
      ax.append(ax62) # 11
      ax.append(ax71) # 12

      ax[0].set_title(title)

      ax[0].grid()
      ax[0].set_ylabel('G [g/L]')

      ax[2].grid()
      ax[2].set_ylabel('E [g/L]')
      
      ax[4].grid()
      ax[4].set_ylabel('X [g/L]')

      ax[6].grid()
      ax[6].set_ylabel('DO [%]')

      ax[8].grid()
      ax[8].set_ylabel('N [rpm]')

      ax[10].grid()
      ax[10].set_ylabel('F [L/h]')

      ax[12].grid()
      ax[12].set_ylabel('V [L]')      
      ax[12].set_xlabel('Time [h]')

      ax[1].grid()
      ax[1].set_ylabel('qG [mole/(h*g)]')

      ax[3].grid()
      ax[3].set_ylabel('qE [mole/(h*g)]')


      ax[5].grid()
      ax[5].set_ylabel('mu [1/h]')

      ax[7].grid()
      ax[7].set_ylabel('qO2 [mole/h,g]')

      ax[9].grid()
      ax[9].set_ylabel('OUR [mole/h]')    

      ax[11].grid()
      ax[11].set_ylabel('Q [W]')        
      ax[11].set_xlabel('Time [h]')     
      
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax[0].plot(t,sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")
      diagrams.append("ax[2].plot(t,sim_res['bioreactor.c[3]'],color='b',linestyle=linetype)")
      diagrams.append("ax[4].plot(t,sim_res['bioreactor.c[1]'],color='b',linestyle=linetype)")
      diagrams.append("ax[6].plot(t,sim_res['DOsensor.out'],color='b',linestyle=linetype)")
      diagrams.append("ax[6].plot(t,sim_res['DO_setpoint.out'],color='y',linestyle='--')")
      diagrams.append("ax[8].step(t,sim_res['bioreactor.N'],color='c',linestyle=linetype)")
      diagrams.append("ax[8].set_ylim([0,1700])")
      diagrams.append("ax[10].plot(t,sim_res['bioreactor.inlet[1].F'],color='c',linestyle=linetype)")
      diagrams.append("ax[12].plot(t,sim_res['bioreactor.V'],color='b',linestyle=linetype)")
      diagrams.append("ax[12].plot(t,sim_res['bioreactor.V_tot'],color='y', linestyle='--')")
      diagrams.append("ax[12].set_ylim([0, 9.0])")

      diagrams.append("ax[1].plot(t,sim_res['bioreactor.culture.qGm'], color='r', linestyle=linetype)")
      diagrams.append("ax[1].plot(t,sim_res['bioreactor.culture.qGr'], color='b', linestyle=linetype)")  
      diagrams.append("ax[3].plot(t,-sim_res['bioreactor.culture.qEm'], color='r', linestyle=linetype)")
      diagrams.append("ax[3].plot(t,sim_res['bioreactor.culture.qEr'], color='b', linestyle=linetype)")
      diagrams.append("ax[5].plot(t,sim_res['bioreactor.culture.q[1]'],color='b',linestyle=linetype)")
      
      diagrams.append("ax[7].plot(t,sim_res['bioreactor.culture.qO2'],color='b',linestyle=linetype)")
      diagrams.append("ax[9].plot(t,sim_res['bioreactor.m[1]']*sim_res['bioreactor.culture.qO2'],color='b',linestyle=linetype)")
      diagrams.append("ax[11].plot(t,sim_res['bioreactor.m[1]']*sim_res['bioreactor.culture.Qspec'],color='b',linestyle=linetype)")

   elif plotType in ['Focus DO-control']:
   
      ax1 = plt.subplot(5,1,1)
      ax2 = plt.subplot(5,1,2)
      ax3 = plt.subplot(5,1,3)
      ax4 = plt.subplot(5,1,4)
      
      ax.clear()
      ax.append(ax1)
      ax.append(ax2)
      ax.append(ax3)
      ax.append(ax4)

      ax[0].set_title(title)    
      ax[0].grid()
      ax[0].set_ylabel('DO [%]')    
    
      ax[1].grid()
      ax[1].set_ylabel('N [rpm]') 

      ax[2].grid()
      ax[2].set_ylabel('OUR [mole/h]') 

      ax[3].grid()
      ax[3].set_ylabel('F [L/h]') 
      ax[3].set_xlabel('Time [h]') 

      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax[0].plot(t,sim_res['DOsensor.out'],color='b',linestyle=linetype)")
      diagrams.append("ax[0].plot(t,sim_res['DO_setpoint.out'],color='r',linestyle='--')")
      diagrams.append("ax[1].step(t,sim_res['bioreactor.N'],color='b',linestyle=linetype)")
      diagrams.append("ax[1].set_ylim([0,2500])")
      diagrams.append("ax[2].plot(t,sim_res['bioreactor.m[1]']*sim_res['bioreactor.culture.qO2'],color='b',linestyle=linetype)")
      diagrams.append("ax[3].plot(t,sim_res['bioreactor.inlet[1].F'],color='b',linestyle=linetype)")
           
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
#  Startup
#------------------------------------------------------------------------------------------------------------------

FMU_explore_info()