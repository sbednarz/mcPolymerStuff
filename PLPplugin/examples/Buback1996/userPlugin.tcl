proc initPlugin {} {
    source PLP.tcl
    initPLP
}

proc runPlugin {times time} {
    runPLP $times $time
}

proc finalizePlugin {times} {
}
