import h5py
import numpy as np
import sys

inputName = "particle_dump.dat"
outputName = "particle_dump.hdf5"

boxSize = 1.

if len(sys.argv) > 1:
    inputName = sys.argv[1]
    outputName = sys.argv[2]

print "Reading", inputName
data = np.loadtxt(inputName)

IDs = data[:,0]
mass = data[:,1]
pos = data[:,2:5]

numPart = np.size(IDs)

print "Found", numPart, "particles"

# Create file
print "Writing", outputName
f = h5py.File("particle_dump.hdf5", "w")


# Write particles
grp = f.create_group("PartType1")
grp.create_dataset("Coordinates", (numPart, 3), data=pos)
grp.create_dataset("Velocities", (numPart, 3), data=np.zeros(np.shape(pos)))
grp.create_dataset("Masses", (numPart, 1), data=mass)
grp.create_dataset("ParticleIDs", (numPart, 1), data=IDs, dtype='i')



#Write header
grp = f.create_group("Header")
grp.attrs.create("BoxSize", boxSize, dtype='d')
grp.attrs.create("Flag_Entropy_ICs", np.array([0,0,0,0,0,0]), dtype='i')
grp.attrs.create("MassTable", np.array([0.,0.,0.,0.,0.,0.]), dtype='d')
grp.attrs.create("NumFilesPerSnapshot", 1., dtype='i')
grp.attrs.create("NumPart_ThisFile", np.array([0,numPart,0,0,0,0]), dtype='i')
grp.attrs.create("NumPart_Total", np.array([0,numPart,0,0,0,0]), dtype='i')
grp.attrs.create("NumPart_Total_HighWord", np.array([0,0,0,0,0,0]), dtype='i')
grp.attrs.create("Time", 0., dtype='d')
