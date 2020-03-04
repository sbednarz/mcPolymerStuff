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

proc initDataProcessing {} {
	global addDataProcessing
	global nPolymerAnalysis
	global PolymerAnalysisSubdirs
	
	set PolymerAnalysisSubdirs "false"
	set nPolymerAnalysis 0
	if {[info exists addDataProcessing(PolymerAnalysis)]} {
		if {[readPolymerAnalysis $addDataProcessing(PolymerAnalysis)] == 0} {
			puts "Error - initDataProcessing $addDataProcessing(PolymerAnalysis) not found"
			exit
		}
	}
	if {$PolymerAnalysisSubdirs == "true"} {
		file mkdir cld
		file mkdir HlogM
		file mkdir analysis
	}
	puts "Data Processing initialized"
}
proc runDataProcessing {times time} {
	global PolymerAnalysisSubdirs
	cldProcessing $times $time
	doPolymerAnalysis $times $time
	if {$PolymerAnalysisSubdirs == "true"} {
		cleanupfiles
	}
}
proc finalizeDataProcessing {times} {
    #sbed
    PnPwProcessing $times
}
proc cleanupfiles {} {
	catch {
		foreach i [glob *.cld] {
			file rename -force $i "cld/$i"
	}	}

	catch {
		foreach i [glob *.HlogM.dat] {
			file rename -force $i "HlogM/$i"
		}
	}

	catch {
		foreach i [glob *.*analysis] {
			file rename -force $i "analysis/$i"
		}
	}
}

proc resultFileProcessing {times} {
	global mole
	global timestep
	global simulationControl
	global nSpecies
	global Species
	global ratioCnMolecules 
	global cf
	
	if {$times == 0} {
		for {set i 1} {$i<=$nSpecies} {incr i} {
			set cf($i) [open "c$Species($i,name).dat" w]
		}
	}
	
	set resf [open "result.$times" r]
	gets $resf zeile1
	set lzeile1 [split $zeile1]
	set timestep($times) [lindex $lzeile1 0]
	for {set i 1} {$i<=$nSpecies} {incr i} {
		gets $resf zeile1
		set lzeile1 [split $zeile1]
		set mole($times,$i) [lindex $lzeile1 1]
		puts $cf($i) [format "%f\t%.8e\t%.8e" $timestep($times) [expr $mole($times,$i)*$ratioCnMolecules] $mole($times,$i)]
		flush $cf($i)
	}
	close $resf
	if {$simulationControl(cleanupresults) == 1} {
		file delete -force "result.$times"
	}	
}
proc readPolymerAnalysis {filename} {
	global PolymerAnalysis
	global nPolymerAnalysis
	global PolymerAnalysisSubdirs
	global nSpecies
	global Species
	if {[file exists $filename]!=1} {
		return 0
	}
	set PolymerAnalysisSubdirs "false"
	set cfgf [open "$filename" r]
	set ndata 0
	while {[gets $cfgf zeile]!= -1} {
		if {[string index $zeile 0] == "#"} { continue }
		
		set lzeile [split $zeile]
		set timestep [lindex $lzeile 0]
		set Speciesname [lindex $lzeile 1]
		if {[lindex $lzeile 0] == "Subdirsbytype"} {
			set PolymerAnalysisSubdirs "true"
			continue
		}
		if {[checkDefinedSpecies $Speciesname]  == 0} {
			puts "Error: readPolymerAnalysis - Species $Speciesname not found"
			exit
		}
		for {set i 2} {$i<=[llength $lzeile]} {incr i} {
			if {[lindex $lzeile $i] == "cld" } {
				incr ndata
				set dataq(timestep,$ndata) $timestep
				set dataq(Speciesname,$ndata) $Speciesname	
				set dataq(type,$ndata) "cld"	
				set dataq(delete,$ndata) "false"
			}
			if {[lindex $lzeile $i] == "HlogM" } {
				incr ndata
				set dataq(timestep,$ndata) $timestep
				set dataq(Speciesname,$ndata) $Speciesname	
				set dataq(type,$ndata) "HlogM"	
				set dataq(delete,$ndata) "false"
			}
		}
		if {[lindex $lzeile 2] == "ConnectType" } {
			incr ndata
			set dataq(timestep,$ndata) $timestep
			set dataq(Speciesname,$ndata) $Speciesname	
			set dataq(type,$ndata) "ConnectType"
			set dataq(cType,$ndata) [lindex $lzeile 3]
			set dataq(equation,$ndata) [lindex $lzeile 4]
			set dataq(value,$ndata) [lindex $lzeile 5]
		}
		if {[lindex $lzeile 2] == "Fragments" } {
			incr ndata
			set dataq(timestep,$ndata) $timestep
			set dataq(Speciesname,$ndata) $Speciesname	
			set dataq(type,$ndata) "Fragments"
			set dataq(equation,$ndata) [lindex $lzeile 3]
			set dataq(value,$ndata) [lindex $lzeile 4]
		}
		if {[lindex $lzeile 2] == "SCB" } {
			incr ndata
			set dataq(timestep,$ndata) $timestep
			set dataq(Speciesname,$ndata) $Speciesname	
			set dataq(type,$ndata) "SCB"
			set dataq(equation,$ndata) [lindex $lzeile 3]
			set dataq(value,$ndata) [lindex $lzeile 4]
		}
		if {[lindex $lzeile 2] == "FragmentLength" } {
			incr ndata
			set dataq(timestep,$ndata) $timestep
			set dataq(Speciesname,$ndata) $Speciesname	
			set dataq(type,$ndata) "FragmentLength"
			set dataq(equation,$ndata) [lindex $lzeile 3]
			set dataq(value,$ndata) [lindex $lzeile 4]
		}
		if {[lindex $lzeile 2] == "Monomer" } {
			incr ndata
			set dataq(timestep,$ndata) $timestep
			set dataq(Speciesname,$ndata) $Speciesname	
			set dataq(type,$ndata) "Monomer"
			if {[checkDefinedSpecies [lindex $lzeile 3]]  == 0} {
				puts "Error: readPolymerAnalysis - Monomer [lindex $lzeile 3] not found"
				exit
			}
			set dataq(Monomername,$ndata) [lindex $lzeile 3]
			set dataq(equation,$ndata) [lindex $lzeile 4]
			set dataq(value,$ndata) [lindex $lzeile 5]
			#puts "$dataq(timestep,$ndata)\t$dataq(Speciesname,$ndata)\t$dataq(type,$ndata)\t$dataq(Monomername,$ndata)\t$dataq(equation,$ndata)\t$dataq(value,$ndata)"
		}
	}
	close $cfgf
	
	set nPolymerAnalysis 0
	for {set i 1} {$i<=$ndata} {incr i} {
		if {$dataq(type,$i) == "cld"} {
			incr nPolymerAnalysis
			set PolymerAnalysis(timestep,$nPolymerAnalysis) $dataq(timestep,$i)	
			set PolymerAnalysis(Speciesname,$nPolymerAnalysis) $dataq(Speciesname,$i)
			set PolymerAnalysis(type,$nPolymerAnalysis) $dataq(type,$i)
		}
		if {$dataq(type,$i) == "HlogM"} {
			incr nPolymerAnalysis
			set PolymerAnalysis(timestep,$nPolymerAnalysis) $dataq(timestep,$i)	
			set PolymerAnalysis(Speciesname,$nPolymerAnalysis) $dataq(Speciesname,$i)
			set PolymerAnalysis(type,$nPolymerAnalysis) $dataq(type,$i)
		}
		if {$dataq(type,$i) == "ConnectType"} {
			incr nPolymerAnalysis
			set PolymerAnalysis(timestep,$nPolymerAnalysis) $dataq(timestep,$i)	
			set PolymerAnalysis(Speciesname,$nPolymerAnalysis) $dataq(Speciesname,$i)
			set PolymerAnalysis(type,$nPolymerAnalysis) $dataq(type,$i)
			set PolymerAnalysis(cType,$nPolymerAnalysis) $dataq(cType,$i)
			set PolymerAnalysis(equation,$nPolymerAnalysis) $dataq(equation,$i)
			set PolymerAnalysis(value,$nPolymerAnalysis) $dataq(value,$i)
		}
		if {$dataq(type,$i) == "Fragments"} {
			incr nPolymerAnalysis
			set PolymerAnalysis(timestep,$nPolymerAnalysis) $dataq(timestep,$i)	
			set PolymerAnalysis(Speciesname,$nPolymerAnalysis) $dataq(Speciesname,$i)
			set PolymerAnalysis(type,$nPolymerAnalysis) $dataq(type,$i)
			set PolymerAnalysis(equation,$nPolymerAnalysis) $dataq(equation,$i)
			set PolymerAnalysis(value,$nPolymerAnalysis) $dataq(value,$i)
		}
		if {$dataq(type,$i) == "SCB"} {
			incr nPolymerAnalysis
			set PolymerAnalysis(timestep,$nPolymerAnalysis) $dataq(timestep,$i)	
			set PolymerAnalysis(Speciesname,$nPolymerAnalysis) $dataq(Speciesname,$i)
			set PolymerAnalysis(type,$nPolymerAnalysis) $dataq(type,$i)
			set PolymerAnalysis(equation,$nPolymerAnalysis) $dataq(equation,$i)
			set PolymerAnalysis(value,$nPolymerAnalysis) $dataq(value,$i)
		}
		if {$dataq(type,$i) == "FragmentLength"} {
			incr nPolymerAnalysis
			set PolymerAnalysis(timestep,$nPolymerAnalysis) $dataq(timestep,$i)	
			set PolymerAnalysis(Speciesname,$nPolymerAnalysis) $dataq(Speciesname,$i)
			set PolymerAnalysis(type,$nPolymerAnalysis) $dataq(type,$i)
			set PolymerAnalysis(equation,$nPolymerAnalysis) $dataq(equation,$i)
			set PolymerAnalysis(value,$nPolymerAnalysis) $dataq(value,$i)
		}
		if {$dataq(type,$i) == "Monomer"} {
			incr nPolymerAnalysis
			set PolymerAnalysis(timestep,$nPolymerAnalysis) $dataq(timestep,$i)	
			set PolymerAnalysis(Speciesname,$nPolymerAnalysis) $dataq(Speciesname,$i)
			set PolymerAnalysis(type,$nPolymerAnalysis) $dataq(type,$i)
			set PolymerAnalysis(Monomername,$nPolymerAnalysis) $dataq(Monomername,$i)
			set PolymerAnalysis(equation,$nPolymerAnalysis) $dataq(equation,$i)
			set PolymerAnalysis(value,$nPolymerAnalysis) $dataq(value,$i)
		}
		
	}
	puts "init polymer analysis"
	for {set i 1} {$i<=$nPolymerAnalysis} {incr i} {
		puts "$PolymerAnalysis(timestep,$i)\t$PolymerAnalysis(Speciesname,$i)\t$PolymerAnalysis(type,$i)"
	}	
}

proc cldProcessing {times time} {
	global PolymerAnalysisSubdirs
	global nSpecies
	global Species
	global Pn
	global Pw
	set nMonomers 0
	for {set j 1} {$j<=$nSpecies} {incr j} {
		if {$Species($j,type)=="Monomer"} {
			incr nMonomers
		}
	}
	for {set i 1} {$i<=$nSpecies} {incr i} {
		if {$Species($i,macromolecular)>0} {
			if {$Species($i,excludeLogDistribution)==0} {
				set pf [open "polymer.$times.$i" r]
				gets $pf  zeile1
				set lzeile1 [split $zeile1]
				if {[lindex $lzeile1 1]<100} {
					set Pn($times,$i) 0.0
					set Pw($times,$i) 0.0
				} else {
					gets $pf  zeile1
					set lzeile1 [split $zeile1] 
					set Pn($times,$i) [lindex $lzeile1 1]
					set Pw($times,$i) [lindex $lzeile1 3]
				}
				close $pf
				file rename -force "polymer.$times.$i" "$Species($i,name).$time.cld"
				set MonomerCount 0
				for {set j 1} {$j<=$nSpecies} {incr j} {
					if {$Species($j,type)=="Monomer"} {
						if {$nMonomers > 1} {
							#file rename -force "polymer.$times.$i.$j" "$Species($i,name).$Species($j,name).$time.cldM"
						}
						incr MonomerCount
						set massMonomer $Species($j,mass)
						set index $j
					}
				}
				if {$MonomerCount == 1} {
					read_cld "$Species($i,name).$time.cld" nk ck Pk
					if {$nk > 10} {
						cld2HlogM $massMonomer logM H Hraw $nk ck Pk
						writeHlogM "$Species($i,name).$time.HlogM.dat" logM H Hraw
					}
				}
				set ccountM 0
				if {$MonomerCount > 1} {
						read_cld "$Species($i,name).$time.cld" nk ck Pk
					for {set j 1} {$j<=$nk} {incr j} {
						set cldM($j) 0.0
					}
					for {set j 1} {$j<=$nSpecies} {incr j} {
						if {$Species($j,type)=="Monomer"} {	
							incr ccountM
							read_cldM $ccountM "$Species($i,name).$time.cld" $Species($j,mass) cldM ck Pk
						}
					}
					if {$nk > 10} {
						cld2HlogMmulti cldM logM H Hraw $nk ck Pk
						writeHlogM "$Species($i,name).$time.HlogM.dat" logM H Hraw
					}
				}
				
			}
		} 
	}
}
proc doPolymerAnalysis {times time} {
	global nSpecies
	global Species
	global PolymerAnalysis
	global nPolymerAnalysis
	set nMonomers 0
	for {set j 1} {$j<=$nSpecies} {incr j} {
		if {$Species($j,type)=="Monomer"} {
			incr nMonomers
		}
	}
	for {set i 1} {$i<=$nPolymerAnalysis} {incr i} {
		if {$time == $PolymerAnalysis(timestep,$i) || $PolymerAnalysis(timestep,$i) == "dt"} {
			set index [SpeciesName2Index $PolymerAnalysis(Speciesname,$i)]
			if {$index != -1} {
				if {$Species($index,macromolecular)	!= 1} {
					puts "Error doPolymerAnalysis: $PolymerAnalysis(Speciesname,$i) not macromolecular" 
					exit
				}
				if {$PolymerAnalysis(type,$i) == "Monomer"} {
					set monomerindex [SpeciesName2Index $PolymerAnalysis(Monomername,$i)]
					if {$PolymerAnalysis(equation,$i) == "<"} {
						set operator 1
					}
					if {$PolymerAnalysis(equation,$i) == "="} {
						set operator 2
					}
					if {$PolymerAnalysis(equation,$i) == ">"} {
						set operator 3
					}
					polymerAnalysis $index	2	$monomerindex	$operator	$PolymerAnalysis(value,$i)
					set analysisfile "$PolymerAnalysis(Speciesname,$i).$PolymerAnalysis(Monomername,$i).$operator.$PolymerAnalysis(value,$i).$time.monomeranalysis"
					file rename -force "polymerMonomerAnalysis" "$analysisfile"
					set HlogMname "$analysisfile.HlogM.dat"
				}
				if {$PolymerAnalysis(type,$i) == "ConnectType"} {
					if {$PolymerAnalysis(equation,$i) == "<"} {
						set operator 1
					}
					if {$PolymerAnalysis(equation,$i) == "="} {
						set operator 2
					}
					if {$PolymerAnalysis(equation,$i) == ">"} {
						set operator 3
					}
					
					polymerAnalysis $index	1	$PolymerAnalysis(cType,$i)	$operator	$PolymerAnalysis(value,$i)
					set analysisfile "$PolymerAnalysis(Speciesname,$i).connect.$PolymerAnalysis(cType,$i).$operator.$PolymerAnalysis(value,$i).$time.analysis"
					file rename -force "polymerConnectAnalysis" "$analysisfile"
					set HlogMname "$analysisfile.HlogM.dat"
					
				}
				if {$PolymerAnalysis(type,$i) == "Fragments"} {
					if {$PolymerAnalysis(equation,$i) == "<"} {
						set operator 1
					}
					if {$PolymerAnalysis(equation,$i) == "="} {
						set operator 2
					}
					if {$PolymerAnalysis(equation,$i) == ">"} {
						set operator 3
					}
					# connect type 100 --> case number of fragments
					polymerAnalysis $index	1	100	$operator	$PolymerAnalysis(value,$i)
					set analysisfile "$PolymerAnalysis(Speciesname,$i).fragments.$operator.$PolymerAnalysis(value,$i).$time.analysis"
					file rename -force "polymerConnectAnalysis" "$analysisfile"
					set HlogMname "$analysisfile.HlogM.dat"
					
				}
				if {$PolymerAnalysis(type,$i) == "SCB"} {
					if {$PolymerAnalysis(equation,$i) == "<"} {
						set operator 1
					}
					if {$PolymerAnalysis(equation,$i) == "="} {
						set operator 2
					}
					if {$PolymerAnalysis(equation,$i) == ">"} {
						set operator 3
					}
					polymerAnalysis $index	3	0	$operator	$PolymerAnalysis(value,$i)
					set analysisfile "$PolymerAnalysis(Speciesname,$i).SCB.$operator.$PolymerAnalysis(value,$i).$time.analysis"
					file rename -force "polymerSCBAnalysis" "$analysisfile"
					set HlogMname "$analysisfile.HlogM.dat"
					
				}
				if {$PolymerAnalysis(type,$i) == "FragmentLength"} {
					if {$PolymerAnalysis(equation,$i) == "<"} {
						set operator 1
					}
					if {$PolymerAnalysis(equation,$i) == "="} {
						set operator 2
					}
					if {$PolymerAnalysis(equation,$i) == ">"} {
						set operator 3
					}
					polymerAnalysis $index	4	0	$operator	$PolymerAnalysis(value,$i)
					set analysisfile "$PolymerAnalysis(Speciesname,$i).FragmentLength.$operator.$PolymerAnalysis(value,$i).$time.analysis"
					file rename -force "polymerFragmentLengthAnalysis" "$analysisfile"
					set HlogMname "$analysisfile.HlogM.dat"
					
				}
				if {$PolymerAnalysis(type,$i) == "HlogM"} {
					writePolymerSpecies $index
					set analysisfile "$PolymerAnalysis(Speciesname,$i).$time"
					file rename -force "polymer.0.0" "$analysisfile.cld"
					set HlogMname "$analysisfile.HlogM.dat"
					set analysisfile "$analysisfile.cld"
				}
				if {$PolymerAnalysis(type,$i) == "cld"} {
					writePolymerSpecies $index
					set analysisfile "$PolymerAnalysis(Speciesname,$i).$time"
					file rename -force "polymer.0.0" "$analysisfile.cld"
					continue
				}
				read_cld "$analysisfile" nk ck Pk
				if {$nk > 10} {
					if {$nMonomers == 1} {
						set MonomerCount 0
						for {set j 1} {$j<=$nSpecies} {incr j} {
							if {$Species($j,type)=="Monomer"} {
								set massMonomer $Species($j,mass)
								break
							}
						}
						cld2HlogM $massMonomer logM H Hraw $nk ck Pk
						writeHlogM "$HlogMname" logM H Hraw
					}
									
					if {$nMonomers > 1} {
						#read_cld "$PolymerAnalysis(Speciesname,$i).$time.cld" nk ck Pk
						for {set j 1} {$j<=$nk} {incr j} {
							set cldM($j) 0.0
						}
						set ccountM 0
						for {set j 1} {$j<=$nSpecies} {incr j} {
							if {$Species($j,type)=="Monomer"} {	
								incr ccountM
								read_cldM $ccountM "$analysisfile" $Species($j,mass) cldM ck Pk
							}
						}
						#if {$nk > 10} {}
						cld2HlogMmulti cldM logM H Hraw $nk ck Pk
						writeHlogM "$HlogMname" logM H Hraw					
					}
				}
				
			}
		}
	}
}
proc mol2concentrationProcessing {times} {
	global mole
	global nSpecies
	global Species
	global ratioCnMolecules 
	global timestep
	
	for {set i 1} {$i<=$nSpecies} {incr i} {
		set cf [open "c$Species($i,name).dat" w]
		puts "write c$Species($i,name)"
		for {set j 0} {$j<=$times} {incr j} {
			puts $cf [format "%f\t%.8e\t%.8e" $timestep($j) [expr $mole($j,$i)*$ratioCnMolecules] $mole($j,$i)]
		}
		close $cf
	}
	
}
proc PnPwProcessing {times} {
	global Pn
	global Pw
	global nSpecies
	global Species
	global timestep
	for {set i 1} {$i<=$nSpecies} {incr i} {
		if {$Species($i,macromolecular)>0} {
			if {$Species($i,excludeLogDistribution)==0} {
				set cf [open "$Species($i,name).polymer.dat" w]
				for {set j 1} {$j<=$times} {incr j} {
					puts $cf [format "%f\t%.1f\t%.1f" $timestep($j) $Pn($j,$i) $Pw($j,$i)]
				}
				close $cf
			}
		}
	}
}
# analysis functions
proc read_cld {filename nk k Pk} {
	upvar $k kk
	upvar $Pk pkk
	upvar $nk nkk
	set fin [open $filename r]
	gets $fin zeile
	gets $fin zeile
	gets $fin zeile
	set nkk 0
	while {[gets $fin zeile]!= -1} {
		set lzeile [split $zeile]
		incr nkk
		set kk($nkk) [lindex $lzeile 0]
		set pkk($nkk) [lindex $lzeile 1]	
	}
	close $fin
}
proc read_cldM {indexM filename MMonomer cM k Pk} {
	upvar $cM cMM
	upvar $k kk
	upvar $Pk pkk
	set fin [open $filename r]
	gets $fin zeile
	gets $fin zeile
	gets $fin zeile
	set nkk 0
	#puts $MMonomer
	while {[gets $fin zeile]!= -1} {
		set lzeile [split $zeile]
		incr nkk
		set nM [lindex $lzeile [expr 3 + $indexM]]	
		set cMM($nkk) [expr $cMM($nkk)+$nM*$MMonomer/$kk($nkk)/$pkk($nkk)]
	}
	close $fin
}


proc interpolation {logM H nk logMk Hk} {
	upvar $logM lM
	upvar $H hh
	upvar $logMk lMk
	upvar $Hk hhk
	set iv1 1
	set iv2 2
	for {set i 1} {$i <=1600} {incr i} {
		set lM($i) [expr $i/200.0]
		set hh($i) 0
		if {$lM($i)<$lMk(1)} continue
		if {$lM($i)>$lMk($nk)} continue
		set ok 0
		while {$ok==0} {
			if { ($lM($i)>=$lMk($iv1))&&($lM($i)<=$lMk($iv2)) } {
				set ok 1
			} else {
				incr iv1
				incr iv2		
			}
			if {$iv2>$nk} {
				set ok 1
			}
			if {$iv2>$nk} continue
			set dx [expr $lMk($iv2) - $lMk($iv1)]
			#puts "$lMk($iv2)\t$lMk($iv1)\t[expr pow(10,$lMk($iv1))]"
			set dy [expr $hhk($iv2) - $hhk($iv1)]
			set dxx [expr $lM($i) - $lMk($iv1)]
			set dyy [expr $dxx/$dx*$dy ]
			set hh($i) [expr $hhk($iv1)+$dyy]
		}
	}
} 
proc integral {n x y} {
	upvar $x xx
	upvar $y yy
	set integ 0
	set h [expr $xx(2)-$xx(1)]
	set SummeRand [expr $yy(1)+$yy($n)]
	set SummeGerade 0
	set SummeUngerade 0
	set gerade 1
	for {set i 2} {$i <= $n} {incr i} {
		if {$gerade == 1} {
			set SummeGerade [expr $SummeGerade + $yy($i)]
			set gerade 0
		} else {
			set SummeUngerade [expr $SummeUngerade + $yy($i)]
			set gerade 1
		}
	}
	set integ [expr ($SummeRand + 2 * $SummeGerade + 4 * $SummeUngerade) * $h / 3]
	return $integ
}
proc cld2HlogM { MMonomer logM H Hraw nk k Pk} {
	upvar $logM lM
	upvar $H hh
	upvar $k kk
	upvar $Pk pkk
	upvar $Hraw hhraw
	
	for {set i 1} {$i <= $nk} {incr i} {
		set lmk($i) [expr log10($kk($i)*$MMonomer)]
		set hhk($i) [expr $pkk($i)*$kk($i)*$kk($i)*$MMonomer*$MMonomer]
	}
	interpolation lM hh $nk lmk hhk
	set integ [integral 1600 lM hh]
	for {set i 1} {$i <=1600} {incr i} {
		set hhraw($i) $hh($i)
		set hh($i) [expr $hh($i)/$integ]
	}
}
proc cld2HlogMmulti { Mcld logM H Hraw nk k Pk} {
	upvar $Mcld MM
	upvar $logM lM
	upvar $H hh
	upvar $Hraw hhraw
	upvar $k kk
	upvar $Pk pkk
	
	for {set i 1} {$i <= $nk} {incr i} {
		set lmk($i) [expr log10($kk($i)*$MM($i))]
		set hhk($i) [expr $pkk($i)*$kk($i)*$kk($i)*$MM($i)*$MM($i)]
	}
	interpolation lM hh $nk lmk hhk
	set integ [integral 1600 lM hh]
	for {set i 1} {$i <=1600} {incr i} {
		set hhraw($i) $hh($i)
		set hh($i) [expr $hh($i)/$integ]
		
	}
}
proc writeHlogM {filename logM H Hraw} {
	upvar $logM lM
	upvar $H hh
	upvar $Hraw hhraw
	
	set fout [open $filename w]
	for {set i 1} {$i <= 1600} {incr i} {
		puts $fout [format "%.4f %.8f %.4e" $lM($i) $hh($i) $hhraw($i)]
	}
	close $fout
}