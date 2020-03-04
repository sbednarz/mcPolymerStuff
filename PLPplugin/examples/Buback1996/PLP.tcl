# PLP plugin for mcPolymer
# 13/3/2018 v1.0 (c) Szczepan Bednarz



# 1 = idle, 2 = pause, 3 = dark, 4 = shot 
proc logPLP {chanel time mode} {
	puts $chanel [format "%.6f $mode" $time]
}

proc initPLP {} {
    puts "PLP plugin initialized ----*"
	global addPlugin
	global I_0

	set zero 1e-13
	set R0 $addPlugin(PLP,R0)

	set addPlugin(PLP,pulsesCounter) 1
	set addPlugin(PLP,burstsCounter) 1
	set addPlugin(PLP,pauseCounter) 0

	set lfile [open "PLPlog.dat" w]
	puts $lfile "time\tstatus"

	# MODE: 1 BurstDelay, 2 Shots, 3 Idle
	# BD = Burst Delay, BN = Burst number
	#
	if {$addPlugin(PLP,BD) == 0 && $addPlugin(PLP,BN) == 1} {
		# Shots mode
		set I_0 $R0
		set addPlugin(PLP,MODE) 2
		logPLP $lfile 0 4

	} elseif {$addPlugin(PLP,BD) > 0 && $addPlugin(PLP,BN) >= 1} {
		# Burst mode
		set I_0 $zero
		set addPlugin(PLP,MODE) 1

	} else {
		puts "!!! PLP plugin error: unknown mode. Check BD & BN"
		exit 1
	}
	

	set addPlugin(PLP,DARK) [expr int(1.0/($addPlugin(PLP,DT)*$addPlugin(PLP,PRR)))]

	set test [expr 1.0/($addPlugin(PLP,PRR)*10)]
	if {$addPlugin(PLP,DT) > $test} {
		puts "!!! PLP plugin error: DT=$addPlugin(PLP,DT) is too big (prr=$addPlugin(PLP,PRR), testval=$test)"
		exit 1
	}
	

	close $lfile
}



# Hz = 1/ (dt * dark)

# BurstDelay -> Shots (PRR, NP) -> BurstRepetition - 1 
proc runPLP {timestep time} {
	global addPlugin
	global mole
	
	set MODE $addPlugin(PLP,MODE)
    
	set lfile [open "PLPlog.dat" a]
	
	set idle 1
	set pause 2
	set dark 3
	set pulse 4
	set burst 5

	if {$MODE == 2} {
		set NP $addPlugin(PLP,NP)
		set R0 $addPlugin(PLP,R0)
		set pulsesCounter $addPlugin(PLP,pulsesCounter)
		global ratioCnMolecules
		set idx [SpeciesName2Index I_]
		if {[expr $timestep % $addPlugin(PLP,DARK)] == 0} {
			if {$pulsesCounter < $NP } {
				incr addPlugin(PLP,pulsesCounter)
				set II [expr $mole($timestep,$idx) ]
				set dI [expr $R0 / $ratioCnMolecules]
				setReactant I_ [expr $II + $dI]
				puts "----* laser pulse: $pulsesCounter of $NP (burst $addPlugin(PLP,burstsCounter) of $addPlugin(PLP,BN))"
				logPLP $lfile $time $pulse
			} else {
				if {$addPlugin(PLP,BD) == 0 && $addPlugin(PLP,BN) == 1} {
					set addPlugin(PLP,MODE) 3
					logPLP $lfile $time $idle
				} else {
					set addPlugin(PLP,MODE) 1
					incr addPlugin(PLP,burstsCounter)
					set addPlugin(PLP,pulsesCounter) 0
					set addPlugin(PLP,pauseCounter) 0
					logPLP $lfile $time $pause
				}
			}
		} else {
				logPLP $lfile $time $dark
		}
	} elseif {$MODE == 1} {
		logPLP $lfile $time $pause
		incr addPlugin(PLP,pauseCounter)
		if {$addPlugin(PLP,pauseCounter) == [expr $addPlugin(PLP,BD)/$addPlugin(PLP,DT)] - 1} {
			if {$addPlugin(PLP,burstsCounter) < [expr $addPlugin(PLP,BN) + 1]} {
				set addPlugin(PLP,MODE) 2
				set addPlugin(PLP,pulsesCounter) 0
				set addPlugin(PLP,pauseCounter) 0
			} else {
				set addPlugin(PLP,MODE) 3
			}
		}
	} elseif {$MODE == 3} {
		logPLP $lfile $time $idle
		#IDLE time
	}
	close $lfile
}

