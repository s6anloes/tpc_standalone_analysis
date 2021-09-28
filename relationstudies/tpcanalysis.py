#!/usr/bin/env python4
# -*- coding: utf-8 -*-

##########################################################################
# basf2 (Belle II Analysis Software Framework)                           #
# Author: The Belle II Collaboration                                     #
#                                                                        #
# See git log for contributors and copyright holders.                    #
# This file is licensed under LGPL-3.0, see LICENSE.md.                  #
##########################################################################

from scipy import stats
import matplotlib.pyplot as plt
import basf2 as b2
from array import array
import numpy as np
import root_numpy
from simulation import add_simulation
from ROOT import Belle2
import ROOT
import glob
ROOT.gROOT.SetBatch(False)

b2.set_random_seed(3)

createMC = False

#: Use particle gun or EvtGen (for Y(4S) events), you can do both at the same time
useParticleGun = False
useEvtGen = True

# When you want to use MC, you can use a previously taken file as input to
# save some time. To do so, delete the 'inputFilesData = None' statement
inputFilesMC = glob.glob('../../1000Muons_1to4GeV_seed11111.root')


class getTPCSimHits(b2.Module):

    def initialize(self):
        """initialize collections to store data if needed"""
        self.digitcount = []
        self.mcpartcount = []
        self.xdiffnormed = []
        self.ydiffnormed = []

    def event(self):
        #: retrieve PXDTrueHits and event meta data from data store
        tpcSimHits = Belle2.PyStoreArray('TPCSimHits')
        # similar for TPCDigits:
        # tpcDigits = Belle2.PyStoreArray('TPCDigits')
        eventdata = Belle2.PyStoreObj('EventMetaData')
        #: get the current event
        currentevent = eventdata.getEvent()

        #: This works on every (Py)StoreArray to get the size = number of entries
        nTPCSimHits = tpcSimHits.getEntries()

        #: loop over all PXDTrueHits
        for simhit in tpcSimHits:
            # position is now a TVector3, so you can use it like a TVector3, e.g. position.Theta(), position.X(), ....
            position = simhit.getPos()
            # # here you can do whatever you want with the position
            endplatez = -83.12

            driftlength = position.Z() - endplatez

            relatedTPCDigits = simhit.getRelationsFrom("TPCDigits")

            if not relatedTPCDigits:
                self.digitcount.append(0)
                continue

            for digit in relatedTPCDigits:
                xdiff = position.X() - digit.getRecoX()
                ydiff = position.Y() - digit.getRecoY()

                self.xdiffnormed.append(xdiff * 10000 * np.sqrt(2 / driftlength))
                self.ydiffnormed.append(ydiff * 10000 * np.sqrt(2 / driftlength))

            self.digitcount.append(relatedTPCDigits.size())

            relatedMCParticles = simhit.getRelationsFrom("MCParticles")

            if not relatedMCParticles:
                continue

            self.mcpartcount.append(relatedMCParticles.size())

    def terminate(self):
        output_root_file = ROOT.TFile('TPCSimHits.root', 'recreate')

        def hist_diffusion(data, bins, color, direction, linecolor):
            fig, ax = plt.subplots()
            (mu, sigma) = stats.norm.fit(data)
            n, axbins, patches = ax.hist(data, bins=bins, color=color, density=True)
            y = stats.norm.pdf(axbins, mu, sigma)
            theplot = plt.plot(axbins, y, color=linecolor, ls='--', linewidth=2)
            ax.set_xlabel(r'diffusion coefficient / $\mu$m/$\sqrt{cm}$', fontsize=15)
            ax.set_ylabel('entries', fontsize=15)
            ax.tick_params(labelsize=15)
            plt.title('Diffusion of TPCDigits in ' + direction + '-Direction', fontsize=15)
            plt.figtext(.74, .82, r'$\sigma$ = $%.2f$' % (sigma), fontsize=12)
            fig.tight_layout()
            plt.savefig(direction + 'diffusion.pdf')

        hist_diffusion(self.xdiffnormed, 50, 'darkred', 'x', linecolor='#410200')
        hist_diffusion(self.ydiffnormed, 50, 'darkred', 'y', linecolor='#410200')

        def plt_hist(data, color, dataname):
            fig, ax = plt.subplots()
            ax.hist(data, color=color, bins=[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
            ax.set_xlabel('number of' + dataname, fontsize=15)
            ax.set_ylabel('entries', fontsize=15)
            ax.set_yscale("log")
            ax.tick_params(labelsize=15)
            # plt.figtext(.612, .7, r'mean = $%.2f$' % (sum(data) / len(data)), fontsize=12, color='darkred')
            plt.title('Number of ' + dataname + ' per TPCSimHit', fontsize=15)
            fig.tight_layout()
            plt.savefig('SimHitsTo' + dataname + '.pdf')

        # plt_hist(self.mcpartcount, 'black', 'MCParticles', '#410200')
        plt_hist(self.digitcount, 'darkred', 'TPCDigits')

        output_root_file.Write()
        output_root_file.Close()


class MCParticleTruthInformation(b2.Module):
    def initialize(self):
        """initialize collections to store data if needed"""
        self.tree = ROOT.TTree('tree', 'tree')
        self.simhitcount = []
        self.simhitarrays = np.array(self.simhitcount, dtype=np.int32)
        self.digitcount = []

    def event(self):
        #: retrieve MCParticle and event meta data from data store
        mcparticles = Belle2.PyStoreArray('MCParticles')  # c.f. mdst/dataobjects/include/MCParticle.h
        eventdata = Belle2.PyStoreObj('EventMetaData')
        currentevent = eventdata.getEvent()
        nTrackableMCTracks = 0
        for particle in mcparticles:
            #: only continue if particles are trackable, which means that they are charged and stable
            if (particle.hasStatus(Belle2.MCParticle.c_PrimaryParticle) and
                particle.hasStatus(Belle2.MCParticle.c_StableInGenerator) and
                    particle.getCharge() != 0):
                #: this track is trackable, add to number of MC tracks in event
                nTrackableMCTracks += 1
                #: get momentum information
                momentum = particle.getMomentum()  # TVector3, in GeV/c
                absmomentum = momentum.Mag()  # in GeV/c
                pt = momentum.Perp()          # in GeV/c
                phi = momentum.Phi()          # in rad
                theta = momentum.Theta()      # in rad
                #: Store data if needed

                # print(f"{pt}  {phi}  {theta}")

                #: get the TPCSimHits related to this MCParticle
                relatedTPCSimHits = particle.getRelationsTo('TPCSimHits')
                relatedTPCDigits = particle.getRelationsTo('TPCDigits')
                count = 0
                # print(relatedTPCSimHits.size())

                if (not relatedTPCSimHits):
                    continue

                self.simhitcount.append(relatedTPCSimHits.size())
                self.simhitarrays = np.append(self.simhitarrays, relatedTPCSimHits.size())
                self.digitcount.append(relatedTPCDigits.size())

    def terminate(self):
        output_root_file = ROOT.TFile('MCparticles_TPC.root', 'recreate')
        simhitarray = np.asarray(self.simhitcount, dtype=np.int32)
        digitarray = np.array(self.digitcount)
        N = len(simhitarray)

        # for i in range(len(simhitarray)):
        #     print('Array = ', simhitarray[i], '\t List = ', self.simhitcount[i])

        # print(simhitarray)
        # simhitarray.dtype = [('SimHitCounts', 'int32')]
        # simhitarray.dtype.names = ['SimHitCounts']

        def plt_hist(axis, data, bins, color, hatch, label):
            counts, edges = np.histogram(data, bins=bins)
            edges = np.repeat(edges, 2)
            hist = np.hstack((0, np.repeat(counts, 2), 0))

            outline, = ax.plot(edges, hist, linewidth=1, color=color)
            axis.fill_between(edges, hist, 0,
                              edgecolor=outline.get_color(), hatch=hatch, label=label,
                              facecolor='none')  # < removes facecolor
            axis.set_ylim(0, None, auto=True)
            ax.set_xlim(-2000, 42000)

        fig, ax = plt.subplots(1)

        plt_hist(ax, self.simhitcount, 500, 'blue', r'\ \ \ \ ', 'TPCSimHits')
        plt_hist(ax, self.digitcount, 500, 'darkred', '////', 'TPCDigits')

        ax.set_xlabel('number of hits', fontsize=15)
        ax.set_ylabel('entries', fontsize=15)
        ax.tick_params(labelsize=15)
        plt.title('Number of SimHits and Digits per MCParticle', fontsize=15)
        plt.legend()
        plt.figtext(.595, .75, r'SimHits: mean = $%.2f$' % (simhitarray.mean()), fontsize=12, color='blue')
        plt.figtext(.612, .7, r' Digits: mean = $%.2f$' % (digitarray.mean()), fontsize=12, color='darkred')
        fig.tight_layout()
        plt.savefig('simhitdistribution.pdf')

        # root_numpy.array2root(self.simhitarrays, 'MC_output.root')

        self.tree.Branch('simhitcount', self.simhitarrays, 'simhitcount/I')  # [' + str(N) + ']/I')
        self.tree.Fill()
        self.tree.Write()

        output_root_file.Write()
        output_root_file.Close()


#: Do all the actual simulation and reconstruction stuff
path = b2.create_path()


if createMC:
    #: specify number of events to be generated
    path.add_module('EventInfoSetter', evtNumList=[1000], runList=[0], expList=[0])

    path.add_module('Gearbox')
    path.add_module('Geometry', useDB=False, excludedComponents=['PXD', 'SVD', 'CDC'],
                    additionalComponents=['VTX-CMOS-5layer', 'TPC'])
    if useEvtGen:
        #: generate BBbar events
        path.add_module('EvtGenInput')
    if useParticleGun:
        particlegun = b2.register_module('ParticleGun')
        #: number of primaries per event
        particlegun.param('nTracks', 1)
        #: useful PDG codes: +/- 11: electron, +/- 13: muon, +/- 211: pion, 22: photon
        particlegun.param('pdgCodes', [-11, -13, 11, 13, 211, -211])
        particlegun.param('momentumGeneration', 'uniform')
        particlegun.param('momentumParams', [0.4, 4])
        #: Direction
        particlegun.param('thetaGeneration', 'uniform')
        particlegun.param('thetaParams', [20, 150])
        particlegun.param('phiGeneration', 'uniform')
        particlegun.param('phiParams', [0, 360])
        #: Vertex
        particlegun.param('vertexGeneration', 'normal')
        particlegun.param('xVertexParams', [0, 0.0])
        particlegun.param('yVertexParams', [0.0, 0])
        particlegun.param('zVertexParams', [0.0, 0])
        particlegun.param('independentVertices', True)

        path.add_module(particlegun)

    # for more information check simulation/scripts/simulation.py
    add_simulation(
        path,
        components=['TPC'],
        bkgfiles=None,
        bkgOverlay=False,
        forceSetPXDDataReduction=True,
        usePXDDataReduction=False,
        cleanupPXDDataReduction=False,
        useVTX=True,
        useTPC=True)

else:
    path.add_module('RootInput', inputFileNames=inputFilesMC)
    path.add_module('Gearbox')
    path.add_module('Geometry', useDB=False, excludedComponents=['PXD', 'SVD', 'CDC'],
                    additionalComponents=['VTX-CMOS-5layer', 'TPC'])

#: every module written in python can be called like this when it doesn't have parameters
path.add_module(getTPCSimHits())
path.add_module(MCParticleTruthInformation())

#: Show progress, i.e. number of events processed so far
path.add_module('Progress')

#: do the actual stuff, without this line nothing will work
b2.process(path)

#: Print statistics at the end of the run
print(b2.statistics)
