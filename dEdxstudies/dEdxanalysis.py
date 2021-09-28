#!/usr/bin/env python4
# -*- coding: utf-8 -*-

##########################################################################
# basf2 (Belle II Analysis Software Framework)                           #
# Author: The Belle II Collaboration                                     #
#                                                                        #
# See git log for contributors and copyright holders.                    #
# This file is licensed under LGPL-3.0, see LICENSE.md.                  #
##########################################################################

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
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


class getTrackInfos(b2.Module):

    def initialize(self):
        """initialize collections to store data if needed"""
        self.muondEdx = []
        self.muondNdx = []
        self.muonmomentum = []
        self.piondEdx = []
        self.piondNdx = []
        self.pionmomentum = []
        self.kaondEdx = []
        self.kaondNdx = []
        self.kaonmomentum = []
        self.protdEdx = []
        self.protdNdx = []
        self.protmomentum = []
        self.elecdEdx = []
        self.elecdNdx = []
        self.elecmomentum = []

    def event(self):
        #: retrieve PXDTrueHits and event meta data from data store
        tpcTrackInfos = Belle2.PyStoreArray('TPCTrackInfos')
        # similar for TPCDigits:
        # tpcDigits = Belle2.PyStoreArray('TPCDigits')
        eventdata = Belle2.PyStoreObj('EventMetaData')
        #: get the current event
        currentevent = eventdata.getEvent()

        #: This works on every (Py)StoreArray to get the size = number of entries
        nTrackInfos = tpcTrackInfos.getEntries()
        print('nTrackInfos = ', nTrackInfos)

        eleccount = 0
        muoncount = 0
        pioncount = 0

        #: loop over all PXDTrueHits
        for trackinfo in tpcTrackInfos:
            # position is now a TVector3, so you can use it like a TVector3, e.g. position.Theta(), position.X(), ....
            edep = trackinfo.getEnergyDeposit() * 1000.0
            # # here you can do whatever you want with the position
            tracklength = trackinfo.getTrackLength() / 10.0

            relatedMCParticles = trackinfo.getRelationsFrom("MCParticles")

            if not relatedMCParticles:
                continue

            for particle in relatedMCParticles:
                momentum = particle.getMomentum().Mag()

                if (abs(particle.getPDG()) == 13) and (abs(trackinfo.getPDGCode()) ==
                                                       13) and (trackinfo.getTrackID() < 51):
                    self.muondEdx.append(edep / tracklength)
                    self.muondNdx.append(edep / (tracklength * 26.174e-3))
                    self.muonmomentum.append(momentum)
                    muoncount += 1
                    break
                elif (abs(particle.getPDG()) == 211) and (abs(trackinfo.getPDGCode()) == 211) and \
                        (particle.getSecondaryPhysicsProcess() == 0):
                    self.piondEdx.append(edep / tracklength)
                    self.piondNdx.append(edep / (tracklength * 26.174e-3))
                    self.pionmomentum.append(momentum)
                    pioncount += 1
                    break
                elif (abs(particle.getPDG()) == 321) and (abs(trackinfo.getPDGCode()) == 321) and \
                        (particle.getSecondaryPhysicsProcess() == 0):
                    self.kaondEdx.append(edep / tracklength)
                    self.kaondNdx.append(edep / (tracklength * 26.174e-3))
                    self.kaonmomentum.append(momentum)
                    break
                elif (abs(particle.getPDG()) == 2212) and (abs(trackinfo.getPDGCode()) == 2212) and \
                        (momentum > 0.2) and (particle.getSecondaryPhysicsProcess() == 0):
                    self.protdEdx.append(edep / tracklength)
                    self.protdNdx.append(edep / (tracklength * 26.174e-3))
                    self.protmomentum.append(momentum)
                    break
                elif (abs(particle.getPDG()) == 11) and (abs(trackinfo.getPDGCode()) == 11) and \
                        (momentum > 0.1) and (trackinfo.getTrackID() < 51):
                    if (edep / tracklength) > 4.0:
                        print(trackinfo.getTrackID())
                    self.elecdEdx.append(edep / tracklength)
                    self.elecdNdx.append(edep / (tracklength * 26.174e-3))
                    self.elecmomentum.append(momentum)
                    eleccount += 1
                    break

        print('eleccount = ', eleccount, ' ; muoncount = ', muoncount, ' ; pioncount = ', pioncount)

    def terminate(self):
        output_root_file = ROOT.TFile('TPCTrackInfos.root', 'recreate')

        fig, ax = plt.subplots()

        elecs = ax.scatter(self.elecmomentum, self.elecdEdx, marker='.', s=1, color='blue', label=r'$e$')
        pions = ax.scatter(self.pionmomentum, self.piondEdx, marker='.', s=1, color='black', label=r'$\pi$')
        muons = ax.scatter(self.muonmomentum, self.muondEdx, marker='.', s=1, color='cyan', label=r'$\mu$')
        kaons = ax.scatter(self.kaonmomentum, self.kaondEdx, marker='.', s=1, color='darkred', label=r'$K$')
        prots = ax.scatter(self.protmomentum, self.protdEdx, marker='.', s=1, color='yellow', label=r'$p$')

        ax.set_xlabel('$p$ / GeV', fontsize=25)
        ax.xaxis.set_major_locator(plt.MaxNLocator(7))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.set_ylabel('d$E$/d$x$ / keV/cm', fontsize=25)
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.set_ylim(0, 10)
        ax.tick_params(labelsize=25)

        fig.set_size_inches(8, 6.9)
        plt.title('dE/dx for different particle species', fontsize=25)
        lgnd = plt.legend(handles=[elecs, pions, muons, kaons, prots], loc=1, prop={'size': 15})
        lgnd.legendHandles[0]._sizes = [35]
        lgnd.legendHandles[1]._sizes = [35]
        lgnd.legendHandles[2]._sizes = [35]
        lgnd.legendHandles[3]._sizes = [35]
        lgnd.legendHandles[4]._sizes = [35]

        fig.tight_layout()
        plt.savefig('dEdxstudies_MCmomentum_25000.pdf')
        plt.show()
        plt.close()

        fig, ax = plt.subplots()

        elecs = ax.scatter(self.elecmomentum, self.elecdNdx, marker='.', s=1, color='blue', label=r'$e$')
        pions = ax.scatter(self.pionmomentum, self.piondNdx, marker='.', s=1, color='black', label=r'$\pi$')
        muons = ax.scatter(self.muonmomentum, self.muondNdx, marker='.', s=1, color='cyan', label=r'$\mu$')
        kaons = ax.scatter(self.kaonmomentum, self.kaondNdx, marker='.', s=1, color='darkred', label=r'$K$')
        prots = ax.scatter(self.protmomentum, self.protdNdx, marker='.', s=1, color='yellow', label=r'$p$')

        ax.set_xlabel('$p$ / GeV', fontsize=25)
        ax.xaxis.set_major_locator(plt.MaxNLocator(7))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.set_ylabel('d$N$/d$x$ / 1/cm', fontsize=25)
        ax.yaxis.set_major_locator(plt.MaxNLocator(8))
        ax.yaxis.set_minor_locator(MultipleLocator(10))
        ax.set_ylim(0, 10 / (26.174e-3))
        ax.tick_params(labelsize=25)

        fig.set_size_inches(8, 6.9)
        plt.title('dN/dx for different particle species', fontsize=25)
        lgnd = plt.legend(handles=[elecs, pions, muons, kaons, prots], loc=1, prop={'size': 15})
        lgnd.legendHandles[0]._sizes = [35]
        lgnd.legendHandles[1]._sizes = [35]
        lgnd.legendHandles[2]._sizes = [35]
        lgnd.legendHandles[3]._sizes = [35]
        lgnd.legendHandles[4]._sizes = [35]

        fig.tight_layout()
        plt.savefig('dNdxstudies_MCmomentum_25000.pdf')
        plt.show()
        plt.close()

        def plt_hist(axis, data, bins, color, hatch, label):
            counts, edges = np.histogram(data, bins=bins)
            edges = np.repeat(edges, 2)
            hist = np.hstack((0, np.repeat(counts, 2), 0))

            outline, = ax.plot(edges, hist, linewidth=1, color=color)
            axis.fill_between(edges, hist, 0,
                              edgecolor=outline.get_color(), hatch=hatch, label=label,
                              facecolor='none')  # < removes facecolor
            axis.set_ylim(0, None, auto=True)

        fig, ax = plt.subplots()

        labels = [r'$e$', r'$\pi$', r'$\mu$', r'$K$']
        colors = ['blue', 'black', 'cyan', 'darkred']
        bins = []
        for i in range(150):
            bins.append(i - 0.5)

        plt_hist(ax, self.elecdNdx, bins, 'blue', r'|||', r'$e$')
        plt_hist(ax, self.piondNdx, bins, 'black', '---', r'$\pi$')
        plt_hist(ax, self.muondNdx, bins, 'cyan', r'\ \ \ ', r'$\mu$')
        plt_hist(ax, self.kaondNdx, bins, 'darkred', '////', r'$K$')
        plt_hist(ax, self.protdNdx, bins, 'yellow', r'++', r'$p$')

        ax.set_xlabel('Number of hits per cm', fontsize=15)

        ax.set_ylabel('entries', fontsize=15)

        ax.set_xlim(0, 150)
        ax.tick_params(labelsize=15)

        handles, labels = ax.get_legend_handles_labels()
        lgnd = ax.legend(handles, labels)
        lgnd.legendHandles[0]._sizes = [35]
        lgnd.legendHandles[1]._sizes = [35]
        lgnd.legendHandles[2]._sizes = [35]
        lgnd.legendHandles[3]._sizes = [35]
        lgnd.legendHandles[4]._sizes = [35]
        plt.title('Number of hits per cm for different particle species', fontsize=15)

        fig.tight_layout()
        plt.savefig('dNdxhisto.pdf')
        plt.show()
        plt.close()

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

                #: get the TPCTrackInfos related to this MCParticle
                relatedTPCTrackInfos = particle.getRelationsTo('TPCTrackInfos')
                relatedTPCDigits = particle.getRelationsTo('TPCDigits')
                count = 0
                # print(relatedTPCTrackInfos.size())

                if (not relatedTPCTrackInfos):
                    continue

                self.simhitcount.append(relatedTPCTrackInfos.size())
                self.simhitarrays = np.append(self.simhitarrays, relatedTPCTrackInfos.size())
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

        fig, ax = plt.subplots(1)

        plt_hist(ax, self.simhitcount, 50, 'blue', r'\ \ \ \ ', 'TPCTrackInfos')
        plt_hist(ax, self.digitcount, 50, 'darkred', '////', 'TPCDigits')

        ax.set_xlabel('number of hits', fontsize=15)
        ax.set_ylabel('entries', fontsize=15)
        ax.tick_params(labelsize=15)
        plt.title('Number of TrackInfos and Digits per MCParticle', fontsize=15)
        plt.legend()
        plt.figtext(.56, .75, r'TrackInfos: mean = $%f$' % (simhitarray.mean()), fontsize=12)
        plt.figtext(.577, .7, r' Digits: mean = $%f$' % (digitarray.mean()), fontsize=12)
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
path.add_module(getTrackInfos())
# path.add_module(MCParticleTruthInformation())

#: Show progress, i.e. number of events processed so far
path.add_module('Progress')

#: do the actual stuff, without this line nothing will work
b2.process(path)

#: Print statistics at the end of the run
print(b2.statistics)
