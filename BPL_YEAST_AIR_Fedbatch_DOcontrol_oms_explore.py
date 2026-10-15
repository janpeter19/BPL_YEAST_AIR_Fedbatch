# Setup application functions BPL_YEAST_AIR_Fedbatch, dependent on previous import of functions from fmu_explore
# Author: Jan Peter Axelsson
# ------------------------------------------------------------------------------------------------------------------
# 2026-08-24 - Created
# 2026-09-18 - Brought in matplotlib.pylot since used here
# 2026-09-18 - Brought numpy back here since used in several notebooks
# 2026-09-22 - Change indentaiton from 3 spaces to 4 using black
# 2026-10-15 - Import for oms related handling using DyMat
# ------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
#  Framework
# --------------------------------------------------------------------------------------------------

import numpy as np
import scipy.io
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------------------------------
#  Specific application functions: newplot(), describe()
# --------------------------------------------------------------------------------------------------


def newplot(title="Yeast fedbatch cultivation", plotType="Overview"):
    """Standard plot window
     title = ''
    two possible diagrams
       diagram = 'Overview' default
       diagram = 'TimeSeriesExtended', 'Extended'"""

    # Reset pens
    resetPen()

    if plotType in ["Overview"]:

        ax11 = plt.subplot(7, 2, 1)
        ax12 = plt.subplot(7, 2, 2)
        ax21 = plt.subplot(7, 2, 3)
        ax22 = plt.subplot(7, 2, 4)
        ax31 = plt.subplot(7, 2, 5)
        ax32 = plt.subplot(7, 2, 6)
        ax41 = plt.subplot(7, 2, 7)
        ax42 = plt.subplot(7, 2, 8)
        ax51 = plt.subplot(7, 2, 9)
        ax52 = plt.subplot(7, 2, 10)
        ax61 = plt.subplot(7, 2, 11)
        ax62 = plt.subplot(7, 2, 12)
        ax71 = plt.subplot(7, 2, 13)

        ax.clear()
        ax.append(ax11)  #  0
        ax.append(ax12)  #  1
        ax.append(ax21)  #  2
        ax.append(ax22)  #  3
        ax.append(ax31)  #  4
        ax.append(ax32)  #  5
        ax.append(ax41)  #  6
        ax.append(ax42)  #  7
        ax.append(ax51)  #  8
        ax.append(ax52)  #  9
        ax.append(ax61)  # 10
        ax.append(ax62)  # 11
        ax.append(ax71)  # 12

        ax[0].set_title(title)

        ax[0].grid()
        ax[0].set_ylabel("G [g/L]")

        ax[2].grid()
        ax[2].set_ylabel("E [g/L]")

        ax[4].grid()
        ax[4].set_ylabel("X [g/L]")

        ax[6].grid()
        ax[6].set_ylabel("DO [%]")

        ax[8].grid()
        ax[8].set_ylabel("N [rpm]")

        ax[10].grid()
        ax[10].set_ylabel("F [L/h]")

        ax[12].grid()
        ax[12].set_ylabel("V [L]")
        ax[12].set_xlabel("Time [h]")

        ax[1].grid()
        ax[1].set_ylabel("qG [mole/(h*g)]")

        ax[3].grid()
        ax[3].set_ylabel("qE [mole/(h*g)]")

        ax[5].grid()
        ax[5].set_ylabel("mu [1/h]")

        ax[7].grid()
        ax[7].set_ylabel("qO2 [mole/h,g]")

        ax[9].grid()
        ax[9].set_ylabel("OUR [mole/h]")

        ax[11].grid()
        ax[11].set_ylabel("Q [W]")
        ax[11].set_xlabel("Time [h]")

        # List of commands to be executed by simu() after a simulation
        diagrams.clear()
        diagrams.append(
            "ax[0].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.c[2]'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[2].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.c[3]'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[4].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.c[1]'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[6].plot(t,sim_res['Fedbatch_DOcontrol.DOsensor.out'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[6].plot(t,sim_res['Fedbatch_DOcontrol.DO_setpoint.out'],color='y',linestyle='--')"
        )
        diagrams.append(
            "ax[8].step(t,sim_res['Fedbatch_DOcontrol.bioreactor.N'],color='c',linestyle=linetype)"
        )
        diagrams.append("ax[8].set_ylim([0,1700])")
        diagrams.append(
            "ax[10].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.inlet[1].F'],color='c',linestyle=linetype)"
        )
        diagrams.append(
            "ax[12].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.V'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[12].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.V_tot'],color='y', linestyle='--')"
        )
        diagrams.append("ax[12].set_ylim([0, 9.0])")

        diagrams.append(
            "ax[1].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.culture.qGm'], color='r', linestyle=linetype)"
        )
        diagrams.append(
            "ax[1].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.culture.qGr'], color='b', linestyle=linetype)"
        )
        diagrams.append(
            "ax[3].plot(t,-sim_res['Fedbatch_DOcontrol.bioreactor.culture.qEm'], color='r', linestyle=linetype)"
        )
        diagrams.append(
            "ax[3].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.culture.qEr'], color='b', linestyle=linetype)"
        )
        diagrams.append(
            "ax[5].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.culture.q[1]'],color='b',linestyle=linetype)"
        )

        diagrams.append(
            "ax[7].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.culture.qO2'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[9].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.m[1]']*sim_res['Fedbatch_DOcontrol.bioreactor.culture.qO2'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[11].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.m[1]']*sim_res['Fedbatch_DOcontrol.bioreactor.culture.Qspec'],color='b',linestyle=linetype)"
        )

    elif plotType in ["Focus DO-control"]:

        ax1 = plt.subplot(5, 1, 1)
        ax2 = plt.subplot(5, 1, 2)
        ax3 = plt.subplot(5, 1, 3)
        ax4 = plt.subplot(5, 1, 4)

        ax.clear()
        ax.append(ax1)
        ax.append(ax2)
        ax.append(ax3)
        ax.append(ax4)

        ax[0].set_title(title)
        ax[0].grid()
        ax[0].set_ylabel("DO [%]")

        ax[1].grid()
        ax[1].set_ylabel("N [rpm]")

        ax[2].grid()
        ax[2].set_ylabel("OUR [mole/h]")

        ax[3].grid()
        ax[3].set_ylabel("F [L/h]")
        ax[3].set_xlabel("Time [h]")

        # List of commands to be executed by simu() after a simulation
        diagrams.clear()
        diagrams.append(
            "ax[0].plot(t,sim_res['Fedbatch_DOcontrol.DOsensor.out'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[0].plot(t,sim_res['Fedbatch_DOcontrol.DO_setpoint.out'],color='r',linestyle='--')"
        )
        diagrams.append(
            "ax[1].step(t,sim_res['Fedbatch_DOcontrol.bioreactor.N'],color='b',linestyle=linetype)"
        )
        diagrams.append("ax[1].set_ylim([0,2500])")
        diagrams.append(
            "ax[2].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.m[1]']*sim_res['Fedbatch_DOcontrol.bioreactor.culture.qO2'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[3].plot(t,sim_res['Fedbatch_DOcontrol.bioreactor.inlet[1].F'],color='b',linestyle=linetype)"
        )

    else:
        print("Plot window type not correct")


def describe(name, decimals=3):
    """Look up description of culture, media, as well as parameters and variables in the model code"""

    if name == "culture":
        print("Saccharomyces cerevisae - default parameters for strain H1022")

    elif name in ["MSL"]:
        describe_MSL()

    else:
        describe_general(name, decimals)


# --------------------------------------------------------------------------------------------------
#  Startup
# --------------------------------------------------------------------------------------------------

FMU_explore_info()
