source /usr/local/amber20/amber.sh
module load vmd193 amber20
source $AMBERHOME/amber.sh
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
file=$1
cutwat=16
nions1=68
#nions1=280 #put number of K or Cl ions to have 0.15 Mi
#nions2=28
rm leap.log

#first check charge of complex
cat<<EOF>charge.leap
    apo=loadpdb ${file}
    charge apo
    saveoff apo apo.off
    saveamberparm apo apo.parm7 apo.rst7
    quit
EOF
tleap -f leaprc.protein.ff14SB -f leaprc.water.tip3p -f leaprc.gaff2 -f charge.leap
ambpdb -p apo.parm7 < apo.rst7 > apo.pdb

#check charge and load second leap file
charge=$(grep "Total unperturbed charge" leap.log | awk '{printf"%5.0f",$4}')
echo $charge
charge_abs=$(echo "sqrt($charge*$charge)" | bc)
nions2=$(($charge_abs+$nions1))

if [ $charge -lt 0 ] 
then
    cat<<EOF>solvate_TIP3P.leap
    loadamberparams /usr/local/amber14/dat/leap/parm/frcmod.ionsjc_tip3p
    loadoff apo.off
    apo_solv=copy apo
    alignaxes apo_solv
    solvatebox apo_solv TIP3PBOX $cutwat
    addionsrand apo_solv K+ ${nions2}
    addionsrand apo_solv Cl- ${nions1}
    saveoff apo_solv apo_solv.off
    saveamberparm apo_solv apo_solv.parm7 apo_solv.rst7
    quit
EOF
    tleap -f solvate_TIP3P.leap
    
else
    cat<<EOF>solvate_TIP3P.leap
    loadamberparams /usr/local/amber20/dat/leap/parm/frcmod.ionsjc_tip3p
    loadoff apo.off
    apo_solv=copy apo
    alignaxes apo_solv
    solvatebox apo_solv TIP3PBOX $cutwat
    addionsrand apo_solv K+ ${nions1}
    addionsrand apo_solv Cl- ${nions2}
    saveoff apo_solv apo_solv.off
    saveamberparm apo_solv apo_solv.parm7 apo_solv.rst7
    quit
EOF
    tleap -f leaprc.protein.ff14SB -f leaprc.water.tip3p -f leaprc.gaff2 -f solvate_TIP3P.leap 
fi

ambpdb -p apo_solv.parm7 < apo_solv.rst7 > apo_solv.pdb
