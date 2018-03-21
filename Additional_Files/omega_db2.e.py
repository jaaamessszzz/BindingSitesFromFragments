#!/usr/bin/env python3

import os
import sys
import logging
from openeye.oechem import *
from openeye.oeomega import *

def write_molecule(outfile, mol):
    ofs = oemolostream(outfile)
    OEWriteMolecule(ofs, mol)
    ofs.close()

def set_defaults(omega, limitConfs=1200, energyWindow=40):
    # File Parameters
    omega.SetCommentEnergy(True)
    omega.SetIncludeInput(False)
    omega.SetRotorOffset(False) #Enkhee
    omega.SetSDEnergy(False)
    omega.SetWarts(True)
    # 3D Construction Parameters
    omega.SetBuildForceField('mmff94s')
    omega.SetCanonOrder(True)
    omega.SetFixDeleteH(True)
    omega.SetDielectric(1.0)
    omega.SetExponent(1.0)
    omega.SetFixRMS(0.15)
    omega.SetFromCT(False)
    omega.SetFixMaxMatch(1)
    omega.SetFixUniqueMatch(True)
    # Structure Enumeration Parameters
    omega.SetEnumNitrogen(False)
    omega.SetEnumRing(False)#Enkhee  # TS & MK 20160524 (from False) (Improves ring pucker sampling)
    # Torsion Driving Parameters
    omega.SetEnergyWindow(energyWindow)  # JJI 20160329  # TS & MK 20160524 (from 6)
    omega.SetMaxConfs(limitConfs)
    omega.SetMaxRotors(-1)
    omega.SetMaxSearchTime(120.0)
    omega.SetRangeIncrement(5)
    omega.SetRMSThreshold(0.50)  # JJI 20160329 # TS & MK 20160524 (from .5)
    omega.SetSearchForceField('mmff94s')
    omega.SetTorsionDrive(True)
    # Stereochemsitry
    #omega.SetStrictStereo(False)
    omega.SetStrictAtomTypes(False)

# Parse arguments here
infile = sys.argv[1]

if os.environ.get('OMEGA_ENERGY_WINDOW', '').strip() != '':
    energyWindow = int(os.environ['OMEGA_ENERGY_WINDOW'])
else:
    energyWindow = 40

if os.environ.get('OMEGA_MAX_CONFS', '').strip() != '':
    limitConfs = int(os.environ['OMEGA_MAX_CONFS'])
else:
    limitConfs = 1200

if len(sys.argv) > 2:
    if sys.argv[2].isdigit():
        numHs = int(sys.argv[2])
        if numHs > 5:
            print("Refusing to build conformations with > 5 rotatable hydrogens")
            sys.exit(-1)
        if numHs >= 4:
            logging.warn("4-5 Rotatable hydrogens  reported. Reducing confs by a factor of 30")
            limitConfs = limitConfs // 30
        elif numHs >= 2:
            logging.warn("2-3 Rotatable hydrogens  reported. Reducing confs by a factor of 3")
            limitConfs = limitConfs // 3


logging.warning('Setting energy window to %d and max confs to %d' % (energyWindow, limitConfs))

omegaOpts = OEOmegaOptions(OEOmegaSampling_Dense)
omega = OEOmega(omegaOpts)
set_defaults(omega, limitConfs=limitConfs, energyWindow=energyWindow)

mol = OEMol()
inroot, inext = os.path.splitext(infile)
ifs = oemolistream(infile)
OEReadMolecule(ifs, mol)
ifs.close()
OEDetermineComponents(mol)
count, ringlist = OEDetermineRingSystems(mol)
rcf = open('ring_count', 'w')
rcf.write('%d\n' % count)
rcf.close()
outfiles = []
fail_count = 0
print('energy: ', omega.GetEnergyWindow())
print('conf: ', omega.GetMaxConfs() )
if count == 0:
    omega.SetMaxConfs(30)
    outfile = "%s.%d.db2in.mol2" % (inroot, count)
    outfiles.append(outfile)
    if omega(mol):
        write_molecule(outfile, mol)
    else:
        fail_count +=1
else:
    pred = OEPartPredAtom(ringlist)
    for i in range(1, count+1):
        pred.SelectPart(i)
        outfile = "%s.%d.db2in.mol2" % (inroot, i)
        outfiles.append(outfile)
        molcopy = OEMol(mol)
        if omega(molcopy, pred):
            write_molecule(outfile, molcopy)
        else:
            fail_count +=1

print(outfiles)

sys.exit(fail_count)
