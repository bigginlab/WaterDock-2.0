import numpy as np
import MDAnalysis
import os
import sys
import scipy.cluster


def writewaterfile(filename, watercoods, finalwaterscores):
    numwater = watercoods.shape[0]
    f1 = open(filename, 'w')

    for j in range(0, numwater):
        header = 'HETATM'
        serial = j+1
        name = 'OW'
        resname = 'SOL'
        chainID = 'A'
        resSeq = j+1
        icode = ' '
        occupancy = 1.0
        tempfactor = np.abs(finalwaterscores[j, 0])
        x = watercoods[j, 0]
        y = watercoods[j, 1]
        z = watercoods[j, 2]
        f1.write("%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                 % (header, serial, name, icode, resname, chainID, resSeq, icode, x, y, z, occupancy, tempfactor))

    f1.close()


def main(proteinfilename, ligandfilename):
    U1 = MDAnalysis.Universe(proteinfilename)
    proteins = U1.select_atoms('protein and not type HD')
    proteincoods = proteins.positions

    U2 = MDAnalysis.Universe(ligandfilename)
    polaratoms = U2.select_atoms('type N or type O or type F or type Cl or type Br')

    numpolaratoms = polaratoms.n_atoms

    Z = int(os.path.isfile('placedwaters.pdb'))
    if Z == 0:
        f1 = open('dockedwaters.pdb', 'w')
        f1.close()
        sys.exit('Empty file')

    U3 = MDAnalysis.Universe('placedwaters.pdb')
    trialwaters = U3.select_atoms('resname SOL and name OW')
    trialwatercoods = trialwaters.positions

    numtrialwaters = trialwatercoods.shape[0]
    waterscores = np.zeros((numtrialwaters), dtype=float)

    tempdist = MDAnalysis.lib.distances.distance_array(trialwatercoods, proteincoods)
    watprodist = np.amin(tempdist, axis=1)

    for i in range(0, numtrialwaters):

        if watprodist[i] < 3.6 and watprodist[i] > 2.00:

            comd = 'vina --receptor ' + proteinfilename + ' --num_modes 1 --exhaustiveness 20 \
                   --ligand water.pdbqt --size_x 0.5 --size_y 0.5 --size_z 0.5 --out waterout.pdbqt --center_x '\
                        + str(trialwatercoods[i, 0]) + ' --center_y ' + str(trialwatercoods[i, 1]) + ' --center_z '\
                             + str(trialwatercoods[i, 2])
            os.system(comd)
            os.system("grep 'RESULT' waterout.pdbqt > water.txt")
            A = np.genfromtxt('water.txt', usecols=3, dtype=float)
            waterscores[i] = A
            os.remove('water.txt')
            os.remove('waterout.pdbqt')

    predictedwatercoods = np.compress(waterscores <= -0.6, trialwatercoods, axis=0)
    predictedwatercoods = np.float32(predictedwatercoods)
    numpredictedwaters = predictedwatercoods.shape[0]

    waterdata = np.genfromtxt('waterdetails.txt', dtype=int)
    predictedwaterscores1 = np.compress(waterscores <= -0.6, waterscores, axis=0)
    predictedwaterscores2 = np.reshape(predictedwaterscores1, (numpredictedwaters, 1))

    ##############################################################################################################
    ##############################################################################################################
    if numpredictedwaters > 1:
        fit = scipy.cluster.hierarchy.fclusterdata(predictedwatercoods, 2.0, criterion='distance', metric='euclidean')
        fit = fit.astype(int)
        numclust = np.max(fit)

        temppredictedwatercoods = np.zeros((numclust, 3), dtype=float)
        temppredictedwatercoods = np.float32(temppredictedwatercoods)
        temppredictedwaterscores = np.zeros((numclust, 1), dtype=float)

        for i in range(1, numclust+1):
            clusttemp = np.compress(fit == i, predictedwatercoods, axis=0)
            tempavg = np.mean(clusttemp, axis=0)
            temppredictedwatercoods[i-1, :] = tempavg

            clusttemp2 = np.compress(fit == i, predictedwaterscores2, axis=0)
            tempavg2 = np.mean(clusttemp2, axis=0)
            temppredictedwaterscores[i-1, 0] = tempavg2

    elif numpredictedwaters <= 1:
        temppredictedwatercoods = predictedwatercoods.copy()
        temppredictedwaterscores = predictedwaterscores2.copy()

    ##############################################################################################################
    ##############################################################################################################

    allligand = U2.select_atoms('all')
    allligandcoods = allligand.positions
    numpredictedwaters = temppredictedwatercoods.shape[0]

    discardindex = np.zeros((numpredictedwaters, 1), dtype=float)

    count = 0
    for i in range(0, numpolaratoms):

        if waterdata.size > 2:
            atomindex = waterdata[i, 0]
            allowedwaters = waterdata[i, 1]

        elif waterdata.size == 2:
            atomindex = waterdata[0]
            allowedwaters = waterdata[1]

        atomcoods = np.zeros((1, 3), dtype=float)
        atomcoods[0, :] = allligandcoods[atomindex, :].copy()
        atomcoods = np.float32(atomcoods)

        atwatdist = MDAnalysis.lib.distances.distance_array(temppredictedwatercoods, atomcoods)
        B = np.where(atwatdist < 3.1)
        mates = np.ravel_multi_index(B, atwatdist.shape)
        nummates = np.size(mates)

        matescores = temppredictedwaterscores[mates]

        if nummates > allowedwaters:

            numdiscardedwaters = nummates - allowedwaters
            for j in range(0, numdiscardedwaters):

                high = np.argmax(matescores)
                removedindex = mates[high]
                matescores = np.delete(matescores, high)
                mates = np.delete(mates, high)
                discardindex[count, 0] = removedindex
                count = count + 1

    trimmeddiscardindex = discardindex[0:count, :].copy()
    trimmeddiscardindex = np.transpose(trimmeddiscardindex)
    trimmeddiscardindex = np.ndarray.astype(trimmeddiscardindex, dtype=int)

    clusteredwatercoods = np.delete(temppredictedwatercoods, trimmeddiscardindex, axis=0)
    finalwaterscores = np.delete(temppredictedwaterscores, trimmeddiscardindex, axis=0)

    writewaterfile('predictedwaters.pdb', clusteredwatercoods, finalwaterscores)

####################################################################################################


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
