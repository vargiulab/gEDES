puts "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
puts "!FUNZIONA SOLO SE L'INDEX DI REF E COMP SONO IDENTICI!!!!!!!!!!!!!!"
puts "!RICORDA DI CREARE DEI FILES PDBQT-LIKE SENZA COLONNE!!!!!!!!!!!!!!"
puts "!DOPO LE COORDINATE E CON LO STESSO ORDINE DI Hs PER IL LIGANDO !!!"
puts "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
set cut1 2.2; set cut2 1.8; set cut3 1.4; set cut 1.0
#set LIG XXX ; set p YYY ; set bs ZZZ
set LIG ADP ; set p AK ; set bs adp_exp
set outfile [open "rmsd-traj-energy-sorted-ali-bs.dat" w]
puts $outfile "#RMSD:FRAME\tBS\tLIG\tBM\tprotein-alpha" 
set logfile [open "fit-and-deviation-traj.log" w]
set refprot [atomselect 0 "name CA and protein"]
set reference [atomselect 0 "protein and noh and BS33_${bs}"]
set reference_plus [atomselect 0 "protein and noh and same residue as within 2 of BS33_${bs}"]
set ref_plus_list [$reference_plus get resid]
set reference_bb_cb [atomselect 0 "protein and (backbone or name CB) and BS33_${bs}"]
set reference_bb_cb_plus [atomselect 0 "protein and (backbone or name CB) and (same residue as within 2 of BS33_${bs})"]
set ref_bb_cb_plus_list [$reference_bb_cb_plus get resid]
#
set comprot [atomselect 1 "name CA and protein"]
set compare [atomselect 1 "protein and noh and BS33_${bs}"]
set compare_plus [atomselect 1 "protein and noh and resid $ref_plus_list"]
set compare_bb_cb [atomselect 1 "protein and (backbone or name CB) and BS33_${bs}"]
set compare_bb_cb_plus [atomselect 1 "protein and (backbone or name CB) and (resid $ref_bb_cb_plus_list)"]
set referencelig [atomselect 0 "noh and resname ${LIG}"]
set comparelig [atomselect 1 "noh and resname ${LIG}"]
#
set referencebmode [atomselect 0 "(protein and noh and BS33_${bs}) or (noh and resname ${LIG})"]
set referencebmode_bb_cb [atomselect 0 "(protein and (backbone or name CB) and BS33_${bs}) or (noh and resname ${LIG})"]
set referencebmode_bb_cb_plus [atomselect 0 "((protein and (backbone or name CB)) and (same residue as within 2 of BS33_${bs})) or (noh and resname ${LIG})"]
set refbmode_bb_cb_plus_list [$referencebmode_bb_cb_plus get resid]
#
set comparebmode [atomselect 1 "(protein and noh and BS33_${bs}) or (noh and resname ${LIG})"]
set comparebmode_bb_cb [atomselect 1 "(protein and (backbone or name CB) and BS33_${bs}) or (noh and resname ${LIG})"]
set comparebmode_bb_cb_plus [atomselect 1 "((protein and (backbone or name CB) and (resid $refbmode_bb_cb_plus_list))) or (noh and resname ${LIG})"]
set combmode_bb_cb_plus_list [$comparebmode_bb_cb_plus get resid]
#
set pdbselection [atomselect 1 "all"]
array set rmssel {} ; foreach at [lsort -integer [$reference get index]] { set rmssel($at) 0.00 }
array set rmssel_bb_cb {} ; foreach at [lsort -integer [$reference_bb_cb get index]] { set rmssel_bb_cb($at) 0.00 }
array set rmssel_bb_cb_plus {} ; foreach at [lsort -integer [$reference_bb_cb_plus get index]] { set rmssel_bb_cb_plus($at) 0.00 }
array set rmsbmsel {} ; foreach at [lsort -integer [$referencebmode get index]] { set rmsbmsel($at) 0.00 } 
array set rmsbmsel_bb_cb {} ; foreach at [lsort -integer [$referencebmode_bb_cb get index]] { set rmsbmsel_bb_cb($at) 0.00 }
array set rmsbmsel_bb_cb_plus {} ; foreach at [lsort -integer [$referencebmode_bb_cb_plus get index]] { set rmsbmsel_bb_cb_plus($at) 0.00 }
array set oldrms [array get rmssel] ; array set oldrms_bb_cb [array get rmssel_bb_cb] ; array set oldrms_bb_cb_plus [array get rmssel_bb_cb_plus]
array set oldbmrms [array get rmsbmsel] ; array set oldbmrms_bb_cb [array get rmsbmsel_bb_cb] ; array set oldbmrms_bb_cb_plus [array get rmsbmsel_bb_cb_plus]
set num_steps [molinfo 1 get numframes]
puts $logfile "NUMSTEPS: $num_steps"
for {set frame 0} {$frame < $num_steps} {incr frame} {
    $compare_bb_cb_plus frame $frame
    $comprot frame $frame
    animate goto $frame
    set trans_mat [measure fit $compare_bb_cb_plus $reference_bb_cb_plus]
    $pdbselection move $trans_mat
    unset trans_mat
    set rmsd [measure rmsd $comparebmode $referencebmode]
    set rmsd_lig [measure rmsd $comparelig $referencelig]
    set rmsd_bs [measure rmsd $compare $reference]
    set trans_mat [measure fit $comprot $refprot]
    $pdbselection move $trans_mat
    unset trans_mat
    set rmsd_prot [measure rmsd $comprot $refprot]
    puts $outfile [ format "\t%6d\t%.2f\t%.2f\t%.2f\t%.2f" ${frame} ${rmsd_bs} ${rmsd_lig} ${rmsd} ${rmsd_prot} ]
    unset rmsd_lig rmsd_bs rmsd rmsd_prot
}
array unset rmssel
array unset rmssel_bb_cb
array unset rmssel_bb_cb_plus
array unset rmsbmsel
array unset rmsbmsel_bb_cb
array unset rmsbmsel_bb_cb_plus
array unset oldrms
array unset oldrms_bb_cb
array unset oldrms_bb_cb_plus
array unset oldbmrms
array unset oldbmrms_bb_cb
array unset oldbmrms_bb_cb_plus
close $outfile
close $logfile
quit
