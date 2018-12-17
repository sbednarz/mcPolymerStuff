
# include Interface
source experiment_setup.tcl
source userPlugin.tcl
source mcPolymerInterface.tcl

set addPlugin(PLP,DT) $DT
set addPlugin(PLP,PRR) $PRR
set addPlugin(PLP,NP) $NP
set addPlugin(PLP,BD) $BD
set addPlugin(PLP,BN) $BN
set addPlugin(PLP,NP) $NP
set addPlugin(PLP,R0) $R0

usePlugin



# define low molecular species
Monomer M $MW
Species I_

Concentration M $M0
Concentration I_ $I_0

# define macromolecules
SpeciesMacro P_
SpeciesMacro D
disableLogDistribution

# reaction coefficients
RateConstant ki $ki
RateConstant kp $kp
RateConstant ktd $ktd
RateConstant ktc $ktc
RateConstant ktr $ktr
RateConstant k1 $k1


Initiation I_ + M --> P_ ki
Propagation P_ + M --> P_ kp
TerminationDisp P_ + P_ --> D + D ktd
Termination P_ + P_ --> D ktc
Transfer2monomer P_ + M --> D + P_ ktr
Transfer_PL-P P_ + I_ --> D k1

set addDataProcessing(PolymerAnalysis) "PolymerAnalysis.cfg"



# debug function
ListAllSpecies

# debug function
ListAllRateConstants

# Interface function
ListAllReactions

# set frequency for writing intermediate results
Setdt $DT

# define number of molecules for simulation 1e8..1e10
InitSimulation $SNM

# start simulation to overall reaction time
Simulation $ST
