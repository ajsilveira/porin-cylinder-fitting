import dataPoints
import mdtraj as md
import parmed as pmd
import cylinderFitting
from simtk.openmm import unit
import simtk.openmm.app as app

from dataPoints import DataPoints
from cylinderFitting import CylinderFitting
from simtk.openmm.app import ForceField, PDBxFile

pdbx = PDBxFile('final-atomistic-system.pdbx')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdbx.topology, rigidWater=False, flexibleConstraints=False, nonbondedMethod=app.PME, nonbondedCutoff=1*unit.nanometer)
structure = pmd.openmm.load_topology(pdbx.topology, system=system, xyz= pdbx.positions)
topo = md.Topology.from_openmm(pdbx.topology)
data =  DataPoints(structure, topo)
data.writeCoordinates('cylinder.xyz')
cylinder = CylinderFitting(data.coordinates)
print(cylinder.vmdCommands())
bottom_atoms, top_atoms = cylinder.atomsInExtremes(data.coordinates,3)
cylinder.writeExtremesCoords(data.coordinates, bottom_atoms, top_atoms, 'extremes.xyz')
