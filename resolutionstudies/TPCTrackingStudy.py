#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################################################
# basf2 (Belle II Analysis Software Framework)                           #
# Author: The Belle II Collaboration                                     #
#                                                                        #
# See git log for contributors and copyright holders.                    #
# This file is licensed under LGPL-3.0, see LICENSE.md.                  #
##########################################################################

#####################################################################
# Just to create some ParticleGun events and fill the TPC containers
# to see if everything works
#####################################################################


import basf2 as b2
import ROOT as r
from beamparameters import add_beamparameters
from tracking.harvesting_validation.combined_module import CombinedTrackingValidationModule
import simulation as sim

import os
import sys
import random

reco_tracks = 'RecoTracks'
mc_reco_tracks = 'MCRecoTracks'
useVTX = True

# Set Random Seed for reproducable simulation. 0 means really random.
rndseed = 11111
# assume the first argument is the random seed
# if(len(sys.argv) > 1):
# rndseed = sys.argv[1]
# pitch = '50'
# see = '95'



pT = 0.0
theta = 0.0
# assume argument is pT
if(len(sys.argv) > 2):
    pT = sys.argv[1]
    theta = sys.argv[2]

# outputDir = '/gpfs/group/belle2/users/loeschca/Resolution/VTXtended/Muon/1000/'
if useVTX:
    outputDir = './pTtheta/VTX/'
else:
    outputDir = './pTtheta/'

b2.set_random_seed(rndseed)
# print(pT)

outFileName = 'pT'+str(pT)+'theta'+str(theta)+'.root'
outFileName = outputDir + outFileName

pT = float(pT)
pT = pT/1000.0
theta = int(theta)
# print(pT)

# Set log level. Can be overridden with the "-l LEVEL" flag for basf2.
b2.set_log_level(b2.LogLevel.INFO)

# ---------------------------------------------------------------------------------------
main = b2.create_path()

eventinfosetter = b2.register_module('EventInfoSetter')
# default phase3 geometry:
exp_number = 0
eventinfosetter.param("expList", [exp_number])
main.add_module(eventinfosetter)

# main.add_module('EvtGenInput')


particlegun = b2.register_module('ParticleGun')
particlegun.logging.log_level = b2.LogLevel.INFO
param_pGun = {
    'pdgCodes': [13],
    'nTracks': 1,
    'momentumGeneration': 'uniformPt',
    'momentumParams': [pT, pT],
    'vertexGeneration': 'uniform',
    'phiGeneration': 'uniform',
    'phiParams': [0, 360],
    'thetaGeneration': 'uniform',
    'thetaParams': [theta, theta],
    'xVertexParams': [0, 0],            # in cm...
    'yVertexParams': [0, 0],
    'zVertexParams': [0, 0]
}

particlegun.param(param_pGun)
main.add_module(particlegun)


# Geometry parameter loader
gearbox = b2.register_module('Gearbox')
# geometryparinit = b2.register_module('TPCGeometryParInitializer', useDB=False)

# Geometry builder
geometry = b2.register_module('Geometry')
geometry.param('useDB', False)
geometry.param('excludedComponents', ['PXD', 'SVD', 'CDC'])
if useVTX:
    geometry.param('additionalComponents', ['VTX-CMOS-7layer-plus3inCDCforTPC', 'TPC'])
else:
    geometry.param('additionalComponents', ['TPC'])
geometry.logging.log_level = b2.LogLevel.ERROR

main.add_module(gearbox)
# main.add_module(geometryparinit)

main.add_module(geometry)

sim.add_simulation(
    main,
    bkgOverlay=False,
    forceSetPXDDataReduction=True,
    usePXDDataReduction=False,
    cleanupPXDDataReduction=False,
    useVTX=useVTX,
    useTPC=True)

if useVTX:
    main.add_module('VTXClusterizer')

# needed for fitting
main.add_module('SetupGenfitExtrapolation')

main.add_module('TrackFinderMCTruthRecoTracks',
                RecoTracksStoreArrayName=reco_tracks,
                WhichParticles=['primary'],
                UseSecondCDCHits=False,
                UsePXDHits=False,
                UseSVDHits=False,
                UseCDCHits=False,
                UseVTXHits=useVTX,
                UseTPCHits=True,
                UseAllTPCHits=True,
                UseOnlyBeforeTOP=True,
                UseNLoops=0.5,
                SplitAfterDeltaT=-1.0)

main.add_module("DAFRecoFitter", recoTracksStoreArrayName=reco_tracks)  # , tpcHitsStoreArrayName='TPCDigits')

# , trackColName='Tracks', trackFitResultColName='TrackFitResults')
main.add_module('TrackCreator', recoTrackColName=reco_tracks, pdgCodes=[13])

main.add_module('TrackFinderMCTruthRecoTracks',
                RecoTracksStoreArrayName=mc_reco_tracks,
                WhichParticles=['primary'],
                UseSecondCDCHits=False,
                UsePXDHits=False,
                UseSVDHits=False,
                UseCDCHits=False,
                UseVTXHits=useVTX,
                UseTPCHits=True,
                UseAllTPCHits=True,
                UseOnlyBeforeTOP=True,
                UseNLoops=0.5,
                SplitAfterDeltaT=-1.0)

main.add_module('MCRecoTracksMatcher',
                mcRecoTracksStoreArrayName=mc_reco_tracks,
                prRecoTracksStoreArrayName=reco_tracks,
                UsePXDHits=False,
                UseSVDHits=False,
                UseCDCHits=False,
                UseVTXHits=useVTX,
                UseTPCHits=True)

main.add_module(
    CombinedTrackingValidationModule(
        name='',
        contact='',
        reco_tracks_name="RecoTracks",
        output_file_name=outFileName,
        expert_level=200))


# Diese beiden Module können hinzugefügt werden, aber da sie soweit ich sehen kann keine Hilfe sind bei den Resolution-studies,
# kannst du sie auch weglassen, spart dir aber nicht wirklich Zeit.
# main.add_module('TrackingPerformanceEvaluation', ParticleHypothesis=13)
# main.add_module('StandardTrackingPerformance', recoTracksStoreArrayName='RecoTracks')


# build the name of the output file
# #outputFileName = outputDir + './' + "SimEvts_Belle2_" + str(rndseed) + '.root'
# outputFileName = 'testScriptForAndreas.root'

# Root output. Default filename can be overriden with '-o' basf2 option.
# rootOutput = b2.register_module('RootOutput')
# rootOutput.param('outputFileName', outputFileName)
# to save some space exclude everything except stuff needed for tracking
# rootOutput.param('excludeBranchNames', ["ARICHAeroHits",
# "ARICHDigits",
# "ARICHSimHits",
# "BKLMDigits",
# "BKLMSimHitPositions",
# "BKLMSimHits",
# "CDCHits",
# "CDCHits4Trg",
# "CDCSimHits",
# "CDCSimHitsToCDCHits4Trg",
# "ECLDigits",
# "ECLDsps",
# "ECLHits",
# "ECLSimHits",
# "ECLTrigs",
# "ECLDiodeHits",
# "EKLMDigits",
# "EKLMSimHits",
# "TOPBarHits",
# "TOPDigits",
# "TOPRawDigits",
# "TOPSimHits"
# ])
# main.add_module(rootOutput)

b2.print_path(main)

main.add_module('Progress')

b2.process(main)
print(b2.statistics)
