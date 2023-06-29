import FWCore.ParameterSet.Config as cms

###################################################################
## ID MODULES
###################################################################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
from EgammaAnalysis.TnPTreeProducer.cmssw_version import isReleaseAbove

def setIDs(process, options):

    switchOnVIDElectronIdProducer(process, DataFormat.AOD if options['useAOD'] else DataFormat.MiniAOD)

    # define which IDs we want to produce
    my_id_modules = [
        'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
        'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
        'ZprimeTo4l.ModifiedHEEP.Identification.modifiedHeepElectronID_cff',
        'ZprimeTo4l.ModifiedHEEP.Identification.modifiedCutBasedElectronID_Fall17_94X_V2_cff'
       ]

    ### add only miniAOD supported IDs
    if not options['useAOD'] :
        my_id_modules.append( 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff' )
        my_id_modules.append( 'EgammaAnalysis.TnPTreeProducer.Identification.cutBasedDoubleElectronHLTPreselecition_Summer16_V1_cff')

    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

    process.egmGsfElectronIDs.physicsObjectSrc            = cms.InputTag(options['ELECTRON_COLL'])
    process.electronMVAValueMapProducer.src               = cms.InputTag(options['ELECTRON_COLL'])

    if not isReleaseAbove(10, 6): # only for CMSSW_10_2
      process.electronRegressionValueMapProducer.srcMiniAOD = cms.InputTag(options['ELECTRON_COLL'])

    #
    # One tag module --> cut based tight 94X V2
    #
    process.tagEleCutBasedTight = cms.EDProducer('GsfElectronSelectorByValueMap' if options['useAOD'] else 'PatElectronSelectorByValueMap',
                                                     input     = cms.InputTag("goodElectrons"),
                                                     cut       = cms.string(options['ELECTRON_TAG_CUTS']),
                                                     selection = cms.InputTag("egmGsfElectronIDs:modifiedCutBasedElectronID-Fall17-94X-V2-tight"),
                                                     id_cut    = cms.bool(True)
                                                )


    #
    # Add many probe modules, use the PatElectronNm1Selector in case we want to check the effect of one cut
    #
    def addNewProbeModule(sequence, name, inputTag, cutNamesToMask=None):
      if cutNamesToMask:
        temp = cms.EDProducer('PatElectronNm1Selector',
                              input          = cms.InputTag("goodElectrons"),
                              cut            = cms.string(options['ELECTRON_CUTS']),
                              selection      = cms.InputTag(inputTag),
                              cutNamesToMask = cutNamesToMask,
                              )
      else:
        temp = cms.EDProducer('GsfElectronSelectorByValueMap' if options['useAOD'] else 'PatElectronSelectorByValueMap',
                              input     = cms.InputTag("goodElectrons"),
                              cut       = cms.string(options['ELECTRON_CUTS']),
                              selection = cms.InputTag(inputTag),
                              id_cut    = cms.bool(True)
                             )
      setattr(process, 'probeEle%s' % name, temp)
      sequence += temp


    probeSequence = cms.Sequence()
    if not options['useAOD']:
      addNewProbeModule(probeSequence, 'HLTsafe',          'egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1')
      addNewProbeModule(probeSequence, 'DoubleEleHLTsafe', 'egmGsfElectronIDs:cutBasedDoubleElectronHLTPreselection-Summer16-V1')

    for wp in ['Veto', 'Loose', 'Medium', 'Tight']:
      addNewProbeModule(probeSequence, 'CutBased%s94XV2' % wp, 'egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-%s' % wp.lower())

    for wp in ['wp80', 'wp90', 'wpLoose']:
      addNewProbeModule(probeSequence, 'MVA94X%snoiso' %wp,   'egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-%s' % wp)
      addNewProbeModule(probeSequence, 'MVA94X%siso' %wp,     'egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-%s' % wp)
      addNewProbeModule(probeSequence, 'MVA94X%snoisoV2' %wp, 'egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-%s' % wp)
      addNewProbeModule(probeSequence, 'MVA94X%sisoV2' %wp,   'egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-%s' % wp)

    addNewProbeModule(probeSequence, 'MVA94XwpHZZisoV2', 'egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wpHZZ')
    addNewProbeModule(probeSequence, 'CutBasedHEEPV70', 'egmGsfElectronIDs:heepElectronID-HEEPV70')
    addNewProbeModule(probeSequence, 'CutBasedModifiedHEEP', 'egmGsfElectronIDs:modifiedHeepElectronID')
    addNewProbeModule(probeSequence, 'CutBasedModified94XTightV2', 'egmGsfElectronIDs:modifiedCutBasedElectronID-Fall17-94X-V2-tight')

    modifiedHeepCuts = ["MinPt", "GsfEleModifiedDEtaInSeed", "GsfEleDPhiIn", "GsfEleFull5x5SigmaIEtaIEtaWithSat", "GsfEleModifiedFull5x5E2x5OverE5x5WithSat",
                        "GsfEleHadronicOverEMLinear", "GsfEleValueMapIsoRho", "GsfEleModifiedEmHadD1IsoRho", "GsfEleDxy", "GsfEleMissingHits", "GsfEleEcalDriven"]

    for cut in modifiedHeepCuts:
      addNewProbeModule(probeSequence, 'CutBasedModifiedHEEPNm1%sCut' % (cut), 'egmGsfElectronIDs:modifiedHeepElectronID', cutNamesToMask=cms.vstring(cut+"Cut_0"))

    #
    # Optional: SUSY variables (broken?)
    #
    if options['addSUSY'] :
      from EgammaAnalysis.TnPTreeProducer.electronsExtrasSUSY_cff  import workingPoints
      for wp in workingPoints: addNewProbeModule(probeSequence, wp, 'susyEleVarHelper:pass' + wp)

    return probeSequence
