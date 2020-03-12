# mcPolymer - Monte Carlo simulation program of radical polymerization
# Copyright (C) 2008 - 2012, 2014 Marco Drache, Georg Drache, Benjamin Hosemann
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# __________________________________________________________________
# 0. global settings
set tcl_precision 10
set nSpecies 0
set nReactions 0
set nRateConstants 0
set simulationControl(cleanupresults) 1
set simulationControl(modemacro) 0
set simulationControl(modemacroLCB) 0
set simulationControl(usePlugin) 0

# __________________________________________________________________
# 1. Species definition

# find species and doubly definitions
proc checkDefinedSpecies { sname } {
	global nSpecies
	global Species
	for {set i 1} {$i<=$nSpecies} {incr i} {
		if {$Species($i,name)==$sname} {
			return $i
		}
	}
	return 0
}

# 1.1 low molecular species

proc Species { sname } {
	global nSpecies
	global Species
	if {[checkDefinedSpecies $sname]!=0} {
		puts "Error - Species - doubly definition:\tSpecies $sname"
		exit 1
	}
	incr nSpecies
	set Species($nSpecies,name) $sname
	set Species($nSpecies,type) "Species"
	set Species($nSpecies,macromolecular) 0
	set Species($nSpecies,concentration) 0.0
	set Species($nSpecies,TransferSpecies) 0
}

proc TransferSpecies { sname } {
	global nSpecies
	global Species
	
	set sindex [checkDefinedSpecies $sname]
	if {$sindex==0} {
		puts "Error - TransferSpecies definition - undefined species:\t$sname"
		exit
	}
	set Species($sindex,TransferSpecies) 1
}

# Monomer Name
proc Monomer { sname MM} {
	global nSpecies
	global Species
	if {[checkDefinedSpecies $sname]!=0} {
		puts "Error - Species - doubly definition:\tMonomer $sname"
		exit 1
	}
	incr nSpecies
	set Species($nSpecies,name) $sname
	set Species($nSpecies,type) "Monomer"
	set Species($nSpecies,macromolecular) 0
	set Species($nSpecies,concentration) 0.0
	set Species($nSpecies,TransferSpecies) 0
	set Species($nSpecies,mass) $MM
}

# Initiator Name f
proc Initiator { sname sf } {
	global nSpecies
	global Species
	if {[checkDefinedSpecies $sname]!=0} {
		puts "Error - Species - doubly definition:\tInitiator $sname"
		exit 1
	}
	incr nSpecies
	set Species($nSpecies,name) $sname
	set Species($nSpecies,type) "Initiator"
	set Species($nSpecies,f) $sf
	set Species($nSpecies,macromolecular) 0
	set Species($nSpecies,concentration) 0.0
	set Species($nSpecies,TransferSpecies) 0
}


# 1.2 macromolecular species

# SpeciesMacro Name
proc SpeciesMacro { sname } {
	global nSpecies
	global Species
	if {[checkDefinedSpecies $sname]!=0} {
		puts "Error - Species - doubly definition:\tSpeciesMacro $sname"
		exit 1
	}
	incr nSpecies
	set Species($nSpecies,name) $sname
	set Species($nSpecies,type) "Macro"
	set Species($nSpecies,macromolecular) 1
	set Species($nSpecies,concentration) 0.0
	set Species($nSpecies,excludeLogDistribution) 0
}

proc disableLogDistribution { } {
	global nSpecies
	global Species
	
	for {set i 1} {$i<=$nSpecies} {incr i} {
		if {$Species($i,macromolecular) == 1} {
			set Species($i,excludeLogDistribution) 1
		}
	}
	return
}

# __________________________________________________________________
# 2. set concentration of defined species / mol/L

# Concentration SpeciesName Concentration
proc Concentration { sname sc} {
	global nSpecies
	global Species
	set sindex [checkDefinedSpecies $sname]
	if {$sindex!=0} {
		set Species($sindex,concentration) $sc
		return
	}
	set efile [open "errormsg.txt" w]
	puts "Error - Concentration - species not found:\t$sname"
	puts $efile "Error - Concentration - species not found:\t$sname"
	close $efile
	exit 1
}

# debug function
# ListAllSpecies
proc ListAllSpecies { } {
	global nSpecies
	global Species	
	puts "$nSpecies species defined:"
	puts "Name\tType\tConcentration / mol/L"
	for {set i 1} {$i<=$nSpecies} {incr i} {
		if { $Species($i,type) == "Macro"} {
			if { $Species($i,excludeLogDistribution) == 1} {
				puts [format "$Species($i,name)\t$Species($i,type)\t$Species($i,concentration)\texcludeLogDistribution"]
				continue
			}
			puts [format "$Species($i,name)\t$Species($i,type)\t$Species($i,concentration)"]
			continue
		} else {
			if { $Species($i,TransferSpecies) == 1} {
				puts [format "$Species($i,name)\t$Species($i,type)\t$Species($i,concentration)\tTransferSpecies"]
				continue
			}
		}
		puts [format "$Species($i,name)\t$Species($i,type)\t$Species($i,concentration)"]
	}
	puts ""
}

# __________________________________________________________________
# 3. reaction rate constants definition

# find reaction rates and doubly definitions
proc checkDefinedRateConstant { cname } {
	global nRateConstants
	global Constant
	for {set i 1} {$i<=$nRateConstants} {incr i} {
		if {$Constant($i,name)==$cname} {
			return $i
		}
	}
	return 0
}

# temperature independent
# RateConstant RateConstantName Value (s,L,mol)
proc RateConstant { cname cvalue} {
	global nRateConstants
	global Constant
	if {[checkDefinedRateConstant $cname]!=0} {
		set efile [open "errormsg.txt" w]
		puts "Error - RateConstant - doubly definition:\t$cname"
		puts $efile "Error - RateConstant - doubly definition:\t$cname"
		close $efile	
		exit 1
	}
	incr nRateConstants
	set Constant($nRateConstants,name) $cname
	set Constant($nRateConstants,type) "tempIndep"
	set Constant($nRateConstants,value) $cvalue
	set Constant($nRateConstants,k0) 0.0
	set Constant($nRateConstants,ea) 0.0
	set Constant($nRateConstants,Va) 0.0
}

# temperature dependent
# RateConstantArrhenius RateConstantName k0 (s,L,mol) Ea (kJ/mol)
proc RateConstantArrhenius { cname ck0 cea} {
	global nRateConstants
	global Constant
	if {[checkDefinedRateConstant $cname]!=0} {
		set efile [open "errormsg.txt" w]
		puts "Error - RateConstantArrhenius - doubly definition:\t$cname"
		puts $efile "Error - RateConstantArrhenius - doubly definition:\t$cname"
		close $efile	
		exit 1
	}
	
	incr nRateConstants
	set Constant($nRateConstants,name) $cname
	set Constant($nRateConstants,type) "tempDep"
	set Constant($nRateConstants,k0) $ck0
	set Constant($nRateConstants,ea) [expr $cea * 1000.0 / 8.314]
	set Constant($nRateConstants,Va) 0.0
	set Constant($nRateConstants,value) 0.0
}

# temperature dependent
# activation Volume --> pressure dependent
# RateConstantArrheniusV# RateConstantName k0 (s,L,mol) Ea (kJ/mol) V# (mL/mol)
proc RateConstantArrheniusV# { cname ck0 cea Va} {
	global nRateConstants
	global Constant
	if {[checkDefinedRateConstant $cname]!=0} {
		set efile [open "errormsg.txt" w]
		puts "Error - RateConstantArrheniusV# - doubly definition:\t$cname"
		puts $efile "Error - RateConstantArrheniusV# - doubly definition:\t$cname"
		close $efile	
		exit 1
	}
	incr nRateConstants
	set Constant($nRateConstants,name) $cname
	set Constant($nRateConstants,type) "tempDep_pressureDep"
	set Constant($nRateConstants,k0) $ck0
	set Constant($nRateConstants,ea) [expr $cea * 1000.0 / 8.314]
	set Constant($nRateConstants,Va) [expr $Va / 10 / 8.314]
	set Constant($nRateConstants,value) 0.0
}


# debug function
# ListAllRateConstants
proc ListAllRateConstants { } {
	global nRateConstants
	global Constant	
	puts "$nRateConstants rate constants defined:"
	puts "Name\tValue\t\[Value\]\t\[Value\]"
	for {set i 1} {$i<=$nRateConstants} {incr i} {
		puts -nonewline [format "$Constant($i,name)\t"]
		if {$Constant($i,type)=="tempIndep"} {
			puts [format "$Constant($i,value)"]
		} else {
			if {$Constant($i,type)=="tempDep"} {
				puts [format "$Constant($i,k0)\t$Constant($i,ea)"]
			} else {
				if {$Constant($i,type)=="tempDep_pressureDep"} {
					puts [format "$Constant($i,k0)\t$Constant($i,ea)\t$Constant($i,Va)"]
				}
			}
		}
	}
	puts ""
}

# __________________________________________________________________
# 4. Reaction definition

# find species and check species type
proc checkDefinedSpeciesType {sourcename sname stype} {
	global nSpecies
	global Species
	set sindex [checkDefinedSpecies $sname]
	if {$sindex==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition $sourcename - undefined species:\t$sname"
		puts $efile "Error - reaction definition $sourcename - undefined species:\t$sname"
		close $efile	
		exit 1
	}
	if {$stype=="anySpecies"} {		
		return 1
	}
	if {$stype=="TransferSpecies"} {
		if {$Species($sindex,TransferSpecies) != 1} {
			set efile [open "errormsg.txt" w]
			puts "Error - $sourcename: $Species($sindex,name) $Species($sindex,type) is not TransferSpecies"
			puts $efile "Error - $sourcename: $Species($sindex,name) $Species($sindex,type) is not TransferSpecies"
			close $efile	
			exit 1
		} else {
			return 1
		}	
	}
	if {$stype=="Macro"} {
		if {[string compare "Macro" $Species($sindex,type)]>0} {
			set efile [open "errormsg.txt" w]
			puts "Error - $sourcename: $Species($sindex,name) is not type Macro"
			puts $efile "Error - $sourcename: $Species($sindex,name) is not type Macro"
			close $efile	
			exit 1
		}
	} 
	return 1
}

# debug function
# ListAllReactions
proc ListAllReactions { } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	
	puts "$nReactions reactions defined:"
	puts "Name\tEductlist\t-->\tProductlist\trate constant"
	for {set i 1} {$i<=$nReactions} {incr i} {
		puts -nonewline [format "$Reaction($i,type) "]
		for {set j 1} {$j<=$Reaction($i,nEducts)} {incr j} {
			if {$j>1} {puts -nonewline "+ "}
			puts -nonewline [format "$Reaction($i,educt$j) "]
		}
		puts -nonewline "--> "
		for {set j 1} {$j<=$Reaction($i,nProducts)} {incr j} {
			if {$j>1} {puts -nonewline "+ "}
			puts -nonewline [format "$Reaction($i,product$j) "]
		}
		puts [format "$Reaction($i,constant)"]
	}
	puts ""
}

# 4.1 elemental reactions

# elemental reaction: 1 educt, 1 product
# Elemental_1E1P SpeciesName1 --> SpeciesName2 RateConstantName
proc Elemental_1E1P { anyA arrow anyB rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Elemental_1E1P $anyA anySpecies]==0} {
		exit
	}
	if {[checkDefinedSpeciesType Elemental_1E1P $anyB anySpecies]==0} {
		exit
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_1E1P - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Elemental_1E1P - undefined reaction rate:\t$rconstant"
		close $efile		
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_1E1P - \"-->\" missing"
		puts $efile "Error - reaction definition Elemental_1E1P - \"-->\" missing"
		close $efile		
		exit 1
	}
	incr nReactions
	set Reaction($nReactions,educt1) $anyA 
	set Reaction($nReactions,nEducts) 1
	set Reaction($nReactions,product1) $anyB 
	set Reaction($nReactions,nProducts) 1	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Elemental_1E1P"	
}

# elemental reaction: 2 educts, 1 product
# Elemental_2E1P SpeciesName1 + SpeciesName2 --> SpeciesName3 RateConstantName
proc Elemental_2E1P { anyA plus anyB arrow anyC rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Elemental_2E1P $anyA anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_2E1P $anyB anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_2E1P $anyC anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_2E1P - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Elemental_2E1P - undefined reaction rate:\t$rconstant"
		close $efile
		exit
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_2E1P - \"+\" missing"
		puts $efile "Error - reaction definition Elemental_2E1P - \"+\" missing"
		close $efile
		exit
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts  "Error - reaction definition Elemental_2E1P - \"-->\" missing"
		puts  $efile "Error - reaction definition Elemental_2E1P - \"-->\" missing"
		close $efile
		exit
	}
	incr nReactions
	set Reaction($nReactions,educt1) $anyA 
	set Reaction($nReactions,educt2) $anyB 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $anyC 
	set Reaction($nReactions,nProducts) 1	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Elemental_2E1P"	
}

# elemental reaction: 1 educt, 2 products
# Elemental_1E2P SpeciesName1 --> SpeciesName2 + SpeciesName3 RateConstantName
proc Elemental_1E2P { anyA arrow anyB plus anyC rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Elemental_1E2P $anyA anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_1E2P $anyB anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_1E2P $anyC anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_1E2P - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Elemental_1E2P - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_1E2P - \"+\" missing"
		puts $efile "Error - reaction definition Elemental_1E2P - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts  "Error - reaction definition Elemental_1E2P - \"-->\" missing"
		puts  $efile "Error - reaction definition Elemental_1E2P - \"-->\" missing"
		close $efile
		exit 1
	}
	
	incr nReactions
	set Reaction($nReactions,educt1) $anyA 
	set Reaction($nReactions,nEducts) 1
	set Reaction($nReactions,product1) $anyB 
	set Reaction($nReactions,product2) $anyC 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Elemental_1E2P"	
}

# elemental reaction: 2 educts, 2 products
# Elemental_2E2P SpeciesName1 + SpeciesName2 --> SpeciesName3 + SpeciesName4 RateConstantName
proc Elemental_2E2P { anyA plus1 anyB arrow anyC plus2 anyD  rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Elemental_2E2P $anyA anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_2E2P $anyB anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_2E2P $anyC anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_2E2P $anyD anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_2E2P - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Elemental_2E2P - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus1 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_2E2P - \"+\" missing"
		puts $efile "Error - reaction definition Elemental_2E2P - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus2 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_2E2P - \"+\" missing"
		puts $efile "Error - reaction definition Elemental_2E2P - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts  "Error - reaction definition Elemental_2E2P - \"-->\" missing"
		puts  $efile "Error - reaction definition Elemental_2E2P - \"-->\" missing"
		close $efile
		exit 1
	}
	incr nReactions
	set Reaction($nReactions,educt1) $anyA 
	set Reaction($nReactions,educt2) $anyB 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $anyC 
	set Reaction($nReactions,product2) $anyD 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Elemental_2E2P"	
}

# elemental reaction: 2 educts, 2 products
# Elemental_2E2P SpeciesName1 + SpeciesName2 --> SpeciesName3 + SpeciesName4 + SpeciesName5 RateConstantName
proc Elemental_2E3P { anyA plus1 anyB arrow anyC plus2 anyD plus3 anyE rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Elemental_2E3P $anyA anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_2E3P $anyB anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_2E3P $anyC anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_2E3P $anyC anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_2E3P $anyC anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Elemental_2E3P $anyC anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_2E3P - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Elemental_2E3P - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus1 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_2E3P - \"+\" missing"
		puts $efile "Error - reaction definition Elemental_2E3P - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus2 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_2E3P - \"+\" missing"
		puts $efile "Error - reaction definition Elemental_2E3P - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus3 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Elemental_2E3P - \"+\" missing"
		puts $efile "Error - reaction definition Elemental_2E3P - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts  "Error - reaction definition Elemental_2E3P - \"-->\" missing"
		puts  $efile "Error - reaction definition Elemental_2E3P - \"-->\" missing"
		close $efile
		exit 1
	}
	
	incr nReactions
	set Reaction($nReactions,educt1) $anyA 
	set Reaction($nReactions,educt2) $anyB 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $anyC 
	set Reaction($nReactions,product2) $anyD 
	set Reaction($nReactions,product3) $anyE 
	set Reaction($nReactions,nProducts) 3	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Elemental_2E3P"	
}

# initiator decomposition
# InitiatorDecomposition InitiatorName  --> SpeciesName1 + SpeciesName2 RateConstantName
proc InitiatorDecomposition { rinitiator arrow rradical1 plus rradical2 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType InitiatorDecomposition $rinitiator Initiator]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType InitiatorDecomposition $rradical1 anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType InitiatorDecomposition $rradical2 anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition InitiatorDecomposition - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition InitiatorDecomposition - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition InitiatorDecomposition - \"+\" missing"
		puts $efile "Error - reaction definition InitiatorDecomposition - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts  "Error - reaction definition InitiatorDecomposition - \"-->\" missing"
		puts  $efile "Error - reaction definition InitiatorDecomposition - \"-->\" missing"
		close $efile
		exit 1
	}
	
	incr nReactions
	set Reaction($nReactions,educt1) $rinitiator 
	set Reaction($nReactions,nEducts) 1
	set Reaction($nReactions,product1) $rradical1 
	set Reaction($nReactions,product2) $rradical1 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "InitiatorDecomposition"	
}

# 4.2 macromolecular reactions

# ThermalInitiation
# ThermalInitiation MonomerName + MonomerName + MonomerName  --> SpeciesMacroName + SpeciesMacroName RateConstantName
proc ThermalInitiation { m1 plus1 m2 plus2 m3 arrow p1 plus3 p2 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType ThermalInitiation $m1 Monomer]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType ThermalInitiation $p1 Macro]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition ThermalInitiation - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition ThermalInitiation - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus1 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition ThermalInitiation - \"+\" missing"
		puts $efile "Error - reaction definition ThermalInitiation - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus2 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition ThermalInitiation - \"+\" missing"
		puts $efile "Error - reaction definition ThermalInitiation - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus3 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition ThermalInitiation - \"+\" missing"
		puts $efile "Error - reaction definition ThermalInitiation - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition ThermalInitiation - \"-->\" missing"
		puts $efile "Error - reaction definition ThermalInitiation - \"-->\" missing"
		close $efile
		exit 1
	}
	incr nReactions
	set Reaction($nReactions,educt1) $m1 
	set Reaction($nReactions,educt2) $m1
	set Reaction($nReactions,educt3) $m1
	set Reaction($nReactions,nEducts) 3
	set Reaction($nReactions,product1) $p1 
	set Reaction($nReactions,product2) $p2 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "ThermalInitiation"	
}

# Propagation
# Propagation SpeciesMacroName + MonomerName --> SpeciesMacroName RateConstantName
proc Propagation { p1 plus m arrow p11 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Propagation $p1 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Propagation $m Monomer]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Propagation $p11 Macro]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Propagation - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Propagation - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Propagation - \"+\" missing"
		puts $efile "Error - reaction definition Propagation - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Propagation - \"-->\" missing"
		puts $efile "Error - reaction definition Propagation - \"-->\" missing"
		close $efile
		exit 1
	}
	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,educt2) $m 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p11 
	set Reaction($nReactions,nProducts) 1	 
	set Reaction($nReactions,constant) $rconstant	
	if {$p1 == $p11} {
		set Reaction($nReactions,type) "Propagation"	
	} else {
		set Reaction($nReactions,type) "CrossPropagation"	
	}
}

# Termination
# Termination SpeciesMacroName1 + SpeciesMacroName1 --> SpeciesMacroName2 RateConstantName
proc Termination { p1 plus p2 arrow p12 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Termination $p1 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Termination $p2 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Termination $p12 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Termination - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Termination - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Termination - \"+\" missing"
		puts $efile "Error - reaction definition Termination - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Termination - \"-->\" missing"
		puts $efile "Error - reaction definition Termination - \"-->\" missing"
		close $efile
		exit 1
	}

	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,educt2) $p2 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p12 
	set Reaction($nReactions,nProducts) 1	 
	set Reaction($nReactions,constant) $rconstant
	if {$p1 == $p2} {
		set Reaction($nReactions,type) "Termination"	
	} else {
		set Reaction($nReactions,type) "CrossTermination"	
	}
}

# TerminationDisp
# TerminationDisp SpeciesMacroName1 + SpeciesMacroName1 --> SpeciesMacroName2 + SpeciesMacroName3 RateConstantName
proc TerminationDisp { p1 plus p2 arrow p3 plus1 p4 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Termination $p1 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Termination $p2 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Termination $p3 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Termination $p4 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TerminationDisp - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition TerminationDisp - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TerminationDisp - \"+\" missing"
		puts $efile "Error - reaction definition TerminationDisp - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus1 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TerminationDisp - \"+\" missing"
		puts $efile "Error - reaction definition TerminationDisp - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TerminationDisp - \"-->\" missing"
		puts $efile "Error - reaction definition TerminationDisp - \"-->\" missing"
		close $efile
		exit 1
	}

	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,educt2) $p2 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p3 
	set Reaction($nReactions,product2) $p4 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant
	if {$p1 == $p2} {
		set Reaction($nReactions,type) "TerminationDisp"	
	} else {
		set Reaction($nReactions,type) "CrossTerminationDisp"	
	}
}
# TerminationLCB
# TerminationLCB SpeciesMacroName1 + SpeciesMacroName1 --> SpeciesMacroName2 RateConstantName
proc TerminationLCB { p1 plus p2 arrow p12 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType TerminationLCB $p1 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType TerminationLCB $p2 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType TerminationLCB $p12 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TerminationLCB - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition TerminationLCB - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TerminationLCB - \"+\" missing"
		puts $efile "Error - reaction definition TerminationLCB - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TerminationLCB - \"-->\" missing"
		puts $efile "Error - reaction definition TerminationLCB - \"-->\" missing"
		close $efile
		exit 1
	}

	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,educt2) $p2 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p12 
	set Reaction($nReactions,nProducts) 1	 
	set Reaction($nReactions,constant) $rconstant
	if {$p1 == $p2} {
		set Reaction($nReactions,type) "TerminationLCB"	
	} else {
		set Reaction($nReactions,type) "CrossTerminationLCB"	
	}
}
# Transfer2monomer
# Transfer2monomer SpeciesMacroName1 + MonomerName --> SpeciesMacroName1 + SpeciesMacroName1 RateConstantName
proc Transfer2monomer { p1 plus m arrow p1D plus1 p2 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Transfer2monomer $p1 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer2monomer $m Monomer]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer2monomer $p1D MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer2monomer $p2 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2monomer - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Transfer2monomer - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2monomer - \"+\" missing"
		puts $efile "Error - reaction definition Transfer2monomer - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus1 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2monomer - \"+\" missing"
		puts $efile "Error - reaction definition Transfer2monomer - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2monomer - \"-->\" missing"
		puts $efile "Error - reaction definition Transfer2monomer - \"-->\" missing"
		close $efile
		exit 1
	}
	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,educt2) $m 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p1D 
	set Reaction($nReactions,product2) $p2 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Transfer2monomer"	
}


# Initiation
# Initiation SpeciesName1 + SpeciesName2 --> SpeciesMacroName RateConstantName
proc Initiation { r plus m arrow p1 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Initiation $r anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Initiation $m anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Initiation $p1 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Initiation - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Initiation - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Initiation - \"+\" missing"
		puts $efile "Error - reaction definition Initiation - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Initiation - \"-->\" missing"
		puts $efile "Error - reaction definition Initiation - \"-->\" missing"
		close $efile
		exit 1
	}

	incr nReactions
	set Reaction($nReactions,educt1) $r 
	set Reaction($nReactions,educt2) $m 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p1 
	set Reaction($nReactions,nProducts) 1	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Initiation"	
}

# InitiationB
# InitiationB SpeciesName1 + SpeciesName2 --> SpeciesMacroName RateConstantName
proc InitiationB { r plus m arrow p1 plus1 r1 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType InitiationB $r anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType InitiationB $m anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType InitiationB $p1 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType InitiationB $m anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition InitiationB - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition InitiationB - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition InitiationB - \"+\" missing"
		puts $efile "Error - reaction definition InitiationB - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition InitiationB - \"-->\" missing"
		puts $efile "Error - reaction definition InitiationB - \"-->\" missing"
		close $efile
		exit 1
	}
	if {$plus1 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition InitiationB - \"+\" missing"
		puts $efile "Error - reaction definition InitiationB - \"+\" missing"
		close $efile
		exit 1
	}
	incr nReactions
	set Reaction($nReactions,educt1) $r 
	set Reaction($nReactions,educt2) $m 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p1 
	set Reaction($nReactions,product2) $r1 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Initiation"	
}

# Transfer_P-P
# Transfer_P-P SpeciesMacroName1 --> SpeciesMacroName2 RateConstantName
proc Transfer_P-P { p arrow pa rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Transfer_P-P $p Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer_P-P $pa Macro]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_P-P - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Transfer_P-P - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_P-P - \"-->\" missing"
		puts $efile "Error - reaction definition Transfer_P-P - \"-->\" missing"
		close $efile
		exit 1
	}

	incr nReactions
	set Reaction($nReactions,educt1) $p 
	set Reaction($nReactions,nEducts) 1
	set Reaction($nReactions,product1) $pa 
	set Reaction($nReactions,nProducts) 1	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Transfer_P-P"	
}

# Transfer_PL-P
# Transfer_PL-P SpeciesMacroName1 + SpeciesName --> SpeciesMacroName2 RateConstantName
proc Transfer_PL-P { p plus a arrow pa rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Transfer_PL-P $p Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer_PL-P $a anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer_PL-P $pa Macro]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_PL-P - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Transfer_PL-P - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_PL-P - \"+\" missing"
		puts $efile "Error - reaction definition Transfer_PL-P - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_PL-P - \"-->\" missing"
		puts $efile "Error - reaction definition Transfer_PL-P - \"-->\" missing"
		close $efile
		exit 1
	}

	incr nReactions
	set Reaction($nReactions,educt1) $p 
	set Reaction($nReactions,educt2) $a 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $pa 
	set Reaction($nReactions,nProducts) 1	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Transfer_PL-P"	
}

# Transfer_P-PL
# Transfer_P-PL SpeciesMacroName1 --> SpeciesMacroName2 + SpeciesName RateConstantName
proc Transfer_P-PL { pa arrow p plus a rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Transfer_P-PL $pa Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer_P-PL $p Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer_P-PL $p anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_P-PL - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Transfer_P-PL - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_P-PL - \"+\" missing"
		puts $efile "Error - reaction definition Transfer_P-PL - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_P-PL - \"-->\" missing"
		puts $efile "Error - reaction definition Transfer_P-PL - \"-->\" missing"
		close $efile
		exit 1
	}

	incr nReactions
	set Reaction($nReactions,educt1) $pa 
	set Reaction($nReactions,nEducts) 1
	set Reaction($nReactions,product1) $p 
	set Reaction($nReactions,product2) $a 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Transfer_P-PL"	
}

# Transfer_PL-PL
# Transfer_PL-PL SpeciesMacroName1 + SpeciesName1 --> SpeciesMacroName2 + SpeciesName2 RateConstantName
proc Transfer_PL-PL { p1 plus1 a arrow p2 plus2 b rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Transfer_PL-PL $p1 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer_PL-PL $a anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer_PL-PL $p2 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer_PL-PL $b anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_PL-PL - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Transfer_PL-PL - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus1 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_PL-PL - \"+\" missing"
		puts $efile "Error - reaction definition Transfer_PL-PL - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus2 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_PL-PL - \"+\" missing"
		puts $efile "Error - reaction definition Transfer_PL-PL - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_PL-PL - \"-->\" missing"
		puts $efile "Error - reaction definition Transfer_PL-PL - \"-->\" missing"
		close $efile
		exit 1
	}

	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,educt2) $a 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p2 
	set Reaction($nReactions,product2) $b 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Transfer_PL-PL"	
}


# Transfer_PP-PP
# Transfer_PP-PP SpeciesMacroName1 + SpeciesMacroName1 --> SpeciesMacroName2 + SpeciesMacroName2 RateConstantName
proc Transfer_PP-PP { p1 plus1 p2 arrow p3 plus2 p4 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Transfer_PP-PP $p1 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer_PP-PP $p2 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer_PP-PP $p3 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer_PP-PP $p4 Macro]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_PP-PP - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Transfer_PP-PP - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus1 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_PP-PP - \"+\" missing"
		puts $efile "Error - reaction definition Transfer_PP-PP - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus2 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_PP-PP - \"+\" missing"
		puts $efile "Error - reaction definition Transfer_PP-PP - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer_PP-PP - \"-->\" missing"
		puts $efile "Error - reaction definition Transfer_PP-PP - \"-->\" missing"
		close $efile
		exit 1
	}

	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,educt2) $p2 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p3 
	set Reaction($nReactions,product2) $p4 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Transfer_PP-PP"	
}



# TransferPropagation_PL-PL_SCB
# TransferPropagation_PL-PL_SCB SpeciesMacroName1 + SpeciesName1 --> SpeciesMacroName2 + SpeciesName2 RateConstantName
proc TransferPropagation_PL-PL_SCB { p1 plus1 a arrow p2 plus2 b rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType TransferPropagation_PL-PL_SCB $p1 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType TransferPropagation_PL-PL_SCB $a anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType TransferPropagation_PL-PL_SCB $p2 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType TransferPropagation_PL-PL_SCB $b anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TransferPropagation_PL-PL_SCB - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition TransferPropagation_PL-PL_SCB - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus1 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TransferPropagation_PL-PL_SCB - \"+\" missing"
		puts $efile "Error - reaction definition TransferPropagation_PL-PL_SCB - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus2 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TransferPropagation_PL-PL_SCB - \"+\" missing"
		puts $efile "Error - reaction definition TransferPropagation_PL-PL_SCB - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TransferPropagation_PL-PL_SCB - \"-->\" missing"
		puts $efile "Error - reaction definition TransferPropagation_PL-PL_SCB - \"-->\" missing"
		close $efile
		exit 1
	}

	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,educt2) $a 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p2 
	set Reaction($nReactions,product2) $b 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "TransferPropagation_PL-PL_SCB"	
}
# TransferPropagation_PL-PL_LCB
# TransferPropagation_PL-PL_LCB SpeciesMacroName1 + SpeciesName1 --> SpeciesMacroName2 + SpeciesName2 RateConstantName
proc TransferPropagation_PL-PL_LCB { p1 plus1 a arrow p2 plus2 b rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType TransferPropagation_PL-PL_LCB $p1 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType TransferPropagation_PL-PL_LCB $a anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType TransferPropagation_PL-PL_LCB $p2 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType TransferPropagation_PL-PL_LCB $b anySpecies]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TransferPropagation_PL-PL_LCB - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition TransferPropagation_PL-PL_LCB - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus1 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TransferPropagation_PL-PL_LCB - \"+\" missing"
		puts $efile "Error - reaction definition TransferPropagation_PL-PL_LCB - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus2 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TransferPropagation_PL-PL_LCB - \"+\" missing"
		puts $efile "Error - reaction definition TransferPropagation_PL-PL_LCB - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition TransferPropagation_PL-PL_LCB - \"-->\" missing"
		puts $efile "Error - reaction definition TransferPropagation_PL-PL_LCB - \"-->\" missing"
		close $efile
		exit 1
	}
	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,educt2) $a 
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p2 
	set Reaction($nReactions,product2) $b 
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "TransferPropagation_PL-PL_LCB"	
}
# Fragmentation
# Fragmentation SpeciesMacro(nChains)Name1 --> SpeciesMacro(n-1Chains)Name2 + SpeciesMacroName3 RateConstantName
proc Fragmentation { p1 arrow p2 plus p3 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Fragmentation $p1 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Fragmentation $p2 Macro]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Fragmentation $p3 Macro]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Fragmentation - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Fragmentation - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Fragmentation - \"+\" missing"
		puts $efile "Error - reaction definition Fragmentation - \"+\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Fragmentation - \"-->\" missing"
		puts $efile "Error - reaction definition Fragmentation - \"-->\" missing"
		close $efile
		exit 1
	}
	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,nEducts) 1
	set Reaction($nReactions,product1) $p2 
	set Reaction($nReactions,product2) $p3
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Fragmentation"	
}
# Transfer2polymer
# Transfer2polymer SpeciesMacroName1 + SpeciesMacroName2 TransferSpecies nameTransferSpecies --> SpeciesMacroName3 + SpeciesMacroName4 RateConstantName
proc Transfer2polymer { p1 plus p2 TransferSpecies trs arrow p3 plus1 p4 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Transfer2polymer $p1 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer2polymer $p2 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer2polymer $p3 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer2polymer $p4 MacroExact]==0} {
		exit 1
	}
	if {[checkDefinedSpeciesType Transfer2polymer $trs TransferSpecies]==0} {
		exit 1
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2polymer - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Transfer2polymer - undefined reaction rate:\t$rconstant"
		close $efile
		exit 1
	}
	if {$plus != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2polymer - \"+\" missing"
		puts $efile "Error - reaction definition Transfer2polymer - \"+\" missing"
		close $efile
		exit 1
	}
	if {$plus1 != "+"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2polymer - \"+\" missing"
		puts $efile "Error - reaction definition Transfer2polymer - \"+\" missing"
		close $efile
		exit 1
	}
	if {$TransferSpecies != "TransferSpecies"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2polymer - \"TransferSpecies\" missing"
		puts $efile "Error - reaction definition Transfer2polymer - \"TransferSpecies\" missing"
		close $efile
		exit 1
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2polymer - \"-->\" missing"
		puts $efile "Error - reaction definition Transfer2polymer - \"-->\" missing"
		close $efile
		exit 1
	}
	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,educt2) $p2
	set Reaction($nReactions,educt3) $trs
	set Reaction($nReactions,nEducts) 3
	set Reaction($nReactions,product1) $p3 
	set Reaction($nReactions,product2) $p4
	set Reaction($nReactions,nProducts) 2	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Transfer2polymer"	
}

# Transfer2polymerIntramolecular
# Transfer2polymerIntramolecular SpeciesMacroName1 TransferSpecies nameTransferSpecies --> SpeciesMacroName2 RateConstantName
proc Transfer2polymerIntramolecular { p1 TransferSpecies trs arrow p2 rconstant } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	

	if {[checkDefinedSpeciesType Transfer2polymerIntramolecular $p1 MacroExact]==0} {
		exit
	}
	if {[checkDefinedSpeciesType Transfer2polymerIntramolecular $p2 MacroExact]==0} {
		exit
	}
	if {[checkDefinedSpeciesType Transfer2polymerIntramolecular $trs TransferSpecies]==0} {
		exit
	}
	if {[checkDefinedRateConstant $rconstant]==0} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2polymerIntramolecular - undefined reaction rate:\t$rconstant"
		puts $efile "Error - reaction definition Transfer2polymerIntramolecular - undefined reaction rate:\t$rconstant"
		close $efile
		exit
	}
	if {$TransferSpecies != "TransferSpecies"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2polymerIntramolecular - \"TransferSpecies\" missing"
		puts $efile "Error - reaction definition Transfer2polymerIntramolecular - \"TransferSpecies\" missing"
		close $efile
		exit
	}
	if {$arrow != "-->"} {
		set efile [open "errormsg.txt" w]
		puts "Error - reaction definition Transfer2polymerIntramolecular - \"-->\" missing"
		puts $efile "Error - reaction definition Transfer2polymerIntramolecular - \"-->\" missing"
		close $efile
		exit
	}

	incr nReactions
	set Reaction($nReactions,educt1) $p1 
	set Reaction($nReactions,educt2) $trs
	set Reaction($nReactions,nEducts) 2
	set Reaction($nReactions,product1) $p2 
	set Reaction($nReactions,nProducts) 1	 
	set Reaction($nReactions,constant) $rconstant	
	set Reaction($nReactions,type) "Transfer2polymerIntramolecular"	
}


# __________________________________________________________________
# 5. Simulation  settings

# Temperature
# Temperature temperature (°C)   
proc Temperature { ttemp } {
	global Temperature
	set Temperature $ttemp
}

# Pressure
# Pressure pressure (bar)   
proc Pressure { pp } {
	global Pressure
	set Pressure $pp
}

# Setdt
# Setdt TimeInterval (s)
proc Setdt { sdt } {
	global simulationControl
	set simulationControl(dt) $sdt	
}

# CleanupResults
# CleanupResults on / off
proc CleanupResults { cr } {
	global simulationControl
	if {$cr=="off"} {
		set simulationControl(cleanupresults) 0
	} else {
		set simulationControl(cleanupresults) 1
	} 
}

# setCoefficient
proc setCoefficient { name k } {
	global nRateConstants
	global Constant
	
	for {set i 1} {$i<=$nRateConstants} {incr i} {
		if {$Constant($i,name)==$name} {
			if {$Constant($i,used)!=0} {
				#puts "--> $Constant($i,index)"
				setk $Constant($i,index) $k
				return
			} else {
				puts "Error - setCoefficient $name not used"
				exit
			}
		}
	}
	puts "Error - setCoefficient $name not found"
	exit
}

# setReactant
proc setReactant { name mole } {
	global nSpecies
	global Species
	
	for {set i 1} {$i<=$nSpecies} {incr i} {
		if {$Species($i,name)==$name} {
			setc $i $mole
			return
		}
	}
	puts "Error - Reactant $name not found"
	exit
}

# InitSimulation
# InitSimulation molecules
proc InitSimulation { nMolecules } {
	global nReactions
	global Reaction
	global nSpecies
	global Species	
	global nRateConstants
	global Constant
	global server1
	global Temperature
	global Pressure
	global time
	global timesteps
	global ratioCnMolecules 
	global simulationControl
	global initialnsx
	
	set initialnsx $nMolecules
	set debugm 0
	if {$debugm==1} {
		set cmdf [open "init.sim-doc" w]
	} else {
		set cmdf [open "init.sim" w]
	}
	set modeLCB 0
	set modeMacro 0
	for {set i 1} {$i<=$nReactions} {incr i} {
		set rtype 0
		switch -exact $Reaction($i,type) {
			"InitiatorDecomposition" { set rtype 1 }
			"Propagation" { set rtype 10 }
			"CrossPropagation" { set rtype 11 }
			"Termination" { set rtype 12 }
			"CrossTermination" { set rtype 13 }
			"Transfer2monomer" { set rtype 14 }
			"Initiation" { set rtype 15 }
			"Transfer_P-P" { set rtype 16 }
			"Transfer_PL-P" { set rtype 16 }
			"Transfer_PL-PL" { set rtype 16 }
			"Transfer_P-PL" { set rtype 16 }
			"Transfer_PP-PP" { set rtype 24 }
			"TransferPropagation_PL-PL_SCB" { set rtype 17 }
			"TransferPropagation_PL-PL_LCB" { set rtype 18 }
			"Transfer2polymer" { 
				set rtype 40 
				set modeLCB 1
			}
			"Fragmentation" { set rtype 20 }
			"ThermalInitiation" { set rtype 21 }
			"TerminationLCB" { set rtype 22 }
			"CrossTerminationLCB" { set rtype 23 }
			"TerminationDisp" { set rtype 24 }
			"CrossTerminationDisp" { set rtype 24 }
			
			"Transfer2polymerIntramolecular" { 
				set rtype 42 
				set modeLCB 1
			}
			
			
		}
		if {$rtype >= 10} {
			set modeMacro 1
		}
	}
	
	if {$modeMacro == 0} {
		puts $cmdf "calculateMacromolecules\t0"
	} else {
		if {$modeMacro == 1} {
			if {$modeLCB == 1} {
				puts $cmdf "calculateMacromolecules\t2"
			} else {
				puts $cmdf "calculateMacromolecules\t1"
			}
		}
	}
	
	puts $cmdf [format "Species\t%u" $nSpecies]

	set csum 0.0
	for {set i 1} {$i<=$nSpecies} {incr i} {
		set csum [expr $csum + $Species($i,concentration)]
	}

	set ratioCnMolecules [expr $csum/$nMolecules]
	
	for {set i 1} {$i<=$nSpecies} {incr i} {
		set mole($i) [expr $nMolecules*$Species($i,concentration)/$csum]
		set Species($i,mole0) $mole($i)
	}


	puts -nonewline "register species "
	set resf [open "result.0" w]
	puts $resf 0.0
	for {set i 1} {$i<=$nSpecies} {incr i} {
		
		if {$Species($i,macromolecular)>0} {
			if {$Species($i,excludeLogDistribution)>0} {
				set stype [expr 1000 + $Species($i,macromolecular)]
			} else {
				set stype [expr 100 + $Species($i,macromolecular)]
			}
		} else {
			if {$Species($i,type)=="Monomer"} {
				set stype 2
			} else {
				
				if {$Species($i,type)=="Initiator"} {
					set stype 1
				} else {
					set stype 0
				}
			}
			if {$Species($i,TransferSpecies) == 1} {
				set stype [expr 50 + $stype]
			}
		}
		if {$debugm==1} {
			if {($stype==2)||($stype==3)} {
				puts $cmdf [format "%s:\t%u\t%u\t%.8e\t%.1f" $Species($i,name) $i $stype $mole($i)]
			} else {
				puts $cmdf [format "%s:\t%u\t%u" $Species($i,name) $i $stype $mole($i)]
			}
		} else {
			if {($stype==2)||($stype==3)} {
				puts $cmdf [format "%u\t%u\t%.8e" $i $stype $mole($i)]
			} else {
				if {$stype==1} {
					puts $cmdf [format "%u\t%u\t%.8e\t%.3f" $i $stype $mole($i) $Species($i,f)]
				} else {
					puts $cmdf [format "%u\t%u\t%.8e\t%.3f" $i $stype $mole($i) 0.0]
				}
			}
		}
		#puts $resf [format "%u\t%u" $i [expr int($mole($i))]]
		puts $resf [format "%u\t%8.14e" $i $mole($i)]
		
	}
	#close $resf

	puts "ok"
	
	puts -nonewline "register rate constants "
	for {set i 1} {$i<=$nRateConstants} {incr i} {
		set Constant($i,used) 0
		set Constant($i,order) 0
		set Constant($i,cf) 1.0
	}
	for {set i 1} {$i<=$nReactions} {incr i} {
		set ik [checkDefinedRateConstant $Reaction($i,constant)]
		set Constant($ik,used) 1
		set rorder $Reaction($i,nEducts)
		set rtype 0
		switch -exact $Reaction($i,type) {
			"Termination2P" { set rtype 3 }
			"TerminationD" { set rtype 9 }
			"Transfer2polymer" { set rtype 40 }
			"Transfer2PolymerIntramolecular" { set rtype 42 }
		}
		set cfaktor 1.0
		if {$rorder==2} {
			if {$Reaction($i,educt1) == $Reaction($i,educt2)} {
				set cfaktor 2.0
			} 
		} 
		if {$rtype == 40} {
			set rorder 2
		}
		if {$rtype == 41} {
			set rorder 2
		}
		if {$rtype == 42} {
			set rorder 2
		}
		if {$rtype == 43} {
			set rorder 2
		}
		if {$rtype == 44} {
			set rorder 2
		}
		set Constant($ik,cf) $cfaktor 
		set Constant($ik,order) $rorder
	}
	set usedConstants 0
	for {set i 1} {$i<=$nRateConstants} {incr i} {
		if {$Constant($i,used)!=0} {
			if {$Constant($i,type)=="tempIndep"} {
				set k $Constant($i,value)
				for {set j 2} {$j<=$Constant($i,order)} {incr j} {
					set k [expr $k*$ratioCnMolecules]
				}
				set Constant($i,MCvalue) [expr $k*$Constant($i,cf)]
			} else {
				if {$Constant($i,type)=="tempDep"} {
					set k0 $Constant($i,k0)
					set ea $Constant($i,ea)
					set k [expr $k0*exp(-1.0*$ea/(273.15+$Temperature))]
					for {set j 2} {$j<=$Constant($i,order)} {incr j} {
						set k [expr $k*$ratioCnMolecules]
					}
					set Constant($i,MCvalue) [expr $k*$Constant($i,cf)]
				} else {
					if {$Constant($i,type)=="tempDep_pressureDep"} {
						set k0 $Constant($i,k0)
						set ea $Constant($i,ea)
						set Va $Constant($i,Va)
						set k [expr exp(log($k0)-$ea/(273.15+$Temperature)-$Va*$Pressure/(273.15+$Temperature))]
						for {set j 2} {$j<=$Constant($i,order)} {incr j} {
							set k [expr $k*$ratioCnMolecules]
						}
						set Constant($i,MCvalue) [expr $k*$Constant($i,cf)]
					}
				}
			}
			incr usedConstants

		}
	}
	puts $cmdf [format "MCRates\t%u" $usedConstants]	
	for {set i 1} {$i<=$nRateConstants} {incr i} {
		if {$Constant($i,used)!=0} {
			if {$debugm==1} {
				puts $cmdf [format "%s\t%u\t%.12e" $Constant($i,name) $i $Constant($i,MCvalue)]
			} else {
				puts $cmdf [format "%u\t%.12e" $i $Constant($i,MCvalue)]
			}
		}
		set Constant($i,index) $i
	}
	close $resf
	puts "ok"
	puts -nonewline "register reactions "
	puts $cmdf [format "Reactions\t%u" $nReactions]	
	for {set i 1} {$i<=$nReactions} {incr i} {
		if {$debugm==1} {
			puts -nonewline $cmdf [format "%s\t%u" $Reaction($i,type) $i]
		} else {
			puts -nonewline $cmdf [format "%u" $i]
		}
		set rtype 0
		switch -exact $Reaction($i,type) {
			"InitiatorDecomposition" { set rtype 1 }
			"Propagation" { set rtype 10 }
			"CrossPropagation" { set rtype 11 }
			"Termination" { set rtype 12 }
			"CrossTermination" { set rtype 13 }
			"Transfer2monomer" { set rtype 14 }
			"Initiation" { set rtype 15 }
			"Transfer_P-P" { set rtype 16 }
			"Transfer_PL-P" { set rtype 16 }
			"Transfer_P-PL" { set rtype 16 }
			"Transfer_PL-PL" { set rtype 16 }
			"Transfer_PP-PP" { set rtype 24 }
			"TransferPropagation_PL-PL_SCB" { set rtype 17 }
			"TransferPropagation_PL-PL_LCB" { set rtype 18 }
			"Fragmentation" { set rtype 20 }
			"ThermalInitiation" { set rtype 21 }
			"TerminationLCB" { set rtype 22 }
			"CrossTerminationLCB" { set rtype 23 }
			"TerminationDisp" { set rtype 24 }
			"CrossTerminationDisp" { set rtype 24 }
			"Transfer2polymer" { set rtype 40 }
			"Transfer2polymerIntramolecular" { set rtype 42 }
			
		}
		puts -nonewline $cmdf [format "\t%u" $rtype]
		set ik [checkDefinedRateConstant $Reaction($i,constant)]
		puts -nonewline $cmdf [format "\t%u" $ik]
		puts -nonewline $cmdf [format "\t%u" $Reaction($i,nEducts)]
		for {set j 1} {$j<=$Reaction($i,nEducts)} {incr j} {
			set side [checkDefinedSpecies $Reaction($i,educt$j)]
			puts -nonewline $cmdf [format "\t%u" $side]
		}
		puts -nonewline $cmdf [format "\t%u" $Reaction($i,nProducts)]	
		for {set j 1} {$j<=$Reaction($i,nProducts)} {incr j} {
			set sidp [checkDefinedSpecies $Reaction($i,product$j)]
			puts -nonewline $cmdf [format "\t%u" $sidp]
		}
		puts $cmdf ""
	}
	close $cmdf
	puts "ok" 
	set time 0.0
	set timesteps 0
	set rf [open "result.init" w]
	puts $rf [format "%.8e" $ratioCnMolecules]
	for {set i 1} {$i<=$nSpecies} {incr i} {
		puts $rf [format "%u\t%s" $i $Species($i,name)]
	}
	close $rf
	initSim
	initDataProcessing
}

# __________________________________________________________________
# 6. Simulation  control

# internal getConversion
proc getConversion {SpeciesIndex timestep} {
	global Species
	global mole
	set conversion [expr ($Species($SpeciesIndex,mole0) - $mole($timestep,$SpeciesIndex)) / $Species($SpeciesIndex,mole0)]
	return $conversion
}

proc SpeciesName2Index { _name } {
	global nSpecies
	global Species
	
	for {set i 1} {$i<=$nSpecies} {incr i} {
		if {$Species($i,name) == $_name} {
			return $i
		}
	}
	return -1
}

proc usePlugin {} {
	global simulationControl
	set simulationControl(usePlugin) 1
	initPlugin
}

# Simulation
# Simulation time (s)
proc Simulation {endtime} {
	global simulationControl
	global timestep
	set commandcount 0
	set fcount 0
	set time 0.0
	
	resultFileProcessing $fcount
	puts "time\t $time s\t[clock seconds]"
	while {$time<$endtime} {
		set time [expr $time + $simulationControl(dt)]
		if {$time>$endtime} {set time $endtime} 
				
		incr commandcount
		gotoTime $commandcount $time
		puts "time\t $time s\t[clock seconds]"
		incr fcount
		resultFileProcessing $fcount
		runDataProcessing $fcount $timestep($fcount)
		if {$simulationControl(usePlugin)==1} {
			runPlugin $fcount $timestep($fcount)
		}
	} 
	finalizeDataProcessing $fcount
}
# include Data Processing Plugin
source dataProcessing.tcl
