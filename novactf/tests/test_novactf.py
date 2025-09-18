# *
# * Authors:     Scipion Team
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************
from typing import Union

from cistem.protocols import CistemProtTsCtffind
from imod.constants import OUTPUT_TILTSERIES_NAME
from imod.protocols import ProtImodImportTransformationMatrix, ProtImodTsNormalization, ProtImodExcludeViews
from imod.protocols.protocol_base_preprocess import FLOAT_DENSITIES_CHOICES
from novactf.protocols import ProtNovaCtf
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import cyanStr, magentaStr
from tomo.objects import SetOfTiltSeries, SetOfCTFTomoSeries
from tomo.protocols import ProtImportTs, ProtImportTsCTF
from tomo.protocols.protocol_import_ctf import ImportChoice
from tomo.tests import RE4_STA_TUTO, TS_03, TS_54, DataSetRe4STATuto
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


class TestNovaCtf(TestBaseCentralizedLayer):
    importedTs = None
    importedTomos = None
    unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value
    expectedTsSetSize = 2
    UNMODIFIED = 'unmodified'
    EXC_VIEWS = 'exc. views'
    RE_STACKED = 're-stacked'
    ctfExcludedViewsDict = {
        TS_03: [0, 1, 38, 39],
        TS_54: [0, 39, 40]
    }
    ctfAnglesCountDictExcluded = {
        TS_03: 36,
        TS_54: 38,
    }
    intersectAnglesCountDictExcluded = {
        TS_03: 36,
        TS_54: 36,
    }

    testAcqObjDict = {
        TS_03: DataSetRe4STATuto.testAcq03.value,
        TS_54: DataSetRe4STATuto.testAcq54.value,
    }
    #
    # testInterpAcqObjDict = {
    #     TS_03: DataSetRe4STATuto.testAcq03Interp.value,
    #     TS_54: DataSetRe4STATuto.testAcq54Interp.value,
    # }
    #
    # anglesCountDict = {
    #     TS_03: 40,
    #     TS_54: 41,
    # }
    # Excluded views stuff
    excludedViewsDict = {
        TS_03: [0, 38, 39],
        TS_54: [0, 1, 38, 39, 40]
    }
    # anglesCountDictExcluded = {
    #     TS_03: 37,
    #     TS_54: 36,
    # }

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls.runPrevProtocols()

    @classmethod
    def runPrevProtocols(cls):
        print(cyanStr('--------------------------------- RUNNING PREVIOUS PROTOCOLS ---------------------------------'))
        cls.importedTs = cls._runImportTs()
        print(
            cyanStr('\n-------------------------------- PREVIOUS PROTOCOLS FINISHED ---------------------------------'))

    @classmethod
    def _runImportTs(cls,
                     filesPattern: str = DataSetRe4STATuto.tsPattern.value,
                     exclusionWords: str = DataSetRe4STATuto.exclusionWordsTs03ts54.value) -> SetOfTiltSeries:
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTs = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=filesPattern,
                                       exclusionWords=exclusionWords,
                                       anglesFrom=2,  # From tlt file
                                       voltage=DataSetRe4STATuto.voltage.value,
                                       magnification=DataSetRe4STATuto.magnification.value,
                                       sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                       amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                       samplingRate=cls.unbinnedSRate,
                                       doseInitial=DataSetRe4STATuto.initialDose.value,
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImgWithTltFile.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

        cls.launchProtocol(protImportTs)
        tsImported = getattr(protImportTs, protImportTs.OUTPUT_NAME, None)
        return tsImported

    @classmethod
    def _runImportCtf(cls, isTsSet: SetOfTiltSeries) -> SetOfCTFTomoSeries:
        print(magentaStr("\n==> Importing the CTFs:"))
        protImportCtf = cls.newProtocol(ProtImportTsCTF,
                                        filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                        filesPattern=DataSetRe4STATuto.ctfPattern.value,
                                        importFrom=ImportChoice.CTFFIND.value,
                                        inputSetOfTiltSeries=isTsSet)
        cls.launchProtocol(protImportCtf)
        importedCtfs = getattr(protImportCtf, protImportCtf._possibleOutputs.CTFs.name, None)
        return importedCtfs

    @classmethod
    def _runImportTrMatrix(cls,
                           inTsSet: SetOfTiltSeries,
                           binningTM: int = 1,
                           binningTS: int = 1,
                           filesPattern: str = DataSetRe4STATuto.transformPattern.value,
                           exclusionWords: str = None) -> SetOfTiltSeries:

        print(magentaStr("\n==> Importing the TS' transformation matrices with IMOD:"
                         f"\n\t- Files pattern = {filesPattern}"
                         f"\n\t- Excluded words = {exclusionWords}"
                         f"\n\t- Transformation matrix binning = {binningTM}"
                         f"\n\t- TS binning = {binningTS}"))
        protImportTrMatrix = cls.newProtocol(ProtImodImportTransformationMatrix,
                                             filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                             filesPattern=filesPattern,
                                             exclusionWords=exclusionWords,
                                             inputSetOfTiltSeries=inTsSet,
                                             binningTM=binningTM,
                                             binningTS=binningTS)
        protImportTrMatrix.setObjLabel(f'trMat_b{binningTM} ts_b{binningTS}')
        cls.launchProtocol(protImportTrMatrix)
        outTsSet = getattr(protImportTrMatrix, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runTsPreprocess(cls,
                         inTsSet: SetOfTiltSeries,
                         binning: int = 1,
                         densAdjustMode: int = 2,
                         excludedViews: bool = False, **kwargs) -> SetOfTiltSeries:

        print(magentaStr(f"\n==> Running the TS preprocessing:"
                         f"\n\t- Binning factor = {binning}"
                         f"\n\t- Excluded views = {excludedViews}"
                         f"\n\t- Adjust densities mode = {FLOAT_DENSITIES_CHOICES[densAdjustMode]}"))
        excludeViewsMsg = 'eV' if excludedViews else ''
        protTsNorm = cls.newProtocol(ProtImodTsNormalization,
                                     inputSetOfTiltSeries=inTsSet,
                                     binning=binning,
                                     floatDensities=densAdjustMode,
                                     **kwargs)

        protTsNorm.setObjLabel(f'Bin_{binning} Mode_{densAdjustMode} {excludeViewsMsg}')
        cls.launchProtocol(protTsNorm)
        tsPreprocessed = getattr(protTsNorm, OUTPUT_TILTSERIES_NAME, None)
        return tsPreprocessed

    @classmethod
    def _runPrevProts(cls, importCtf=True):
        importedCtfs = None
        if importCtf:
            importedCtfs = cls._runImportCtf(cls.importedTs)
        tsWithAlignment = cls._runImportTrMatrix(cls.importedTs)
        tsWithAliBin4 = cls._runTsPreprocess(tsWithAlignment, binning=4)
        return importedCtfs, tsWithAliBin4

    @classmethod
    def _gentestAcqObjDictReStacked(cls, ts03MinAngle=-999, ts03MaxAngle=999, ts54MinAngle=-999, ts54MaxAngle=999):
        # The angle min and angle max for the re-stacked TS, as these values may change if the
        # removed tilt-images are the first or the last, for example. FGor the CTF correction, the tilt
        # axis angle and the initial and accum dose will be always 0.
        testAcqObjDictReStacked = {}
        acq_TS_03 = cls.testAcqObjDict[TS_03].clone()
        acq_TS_03.setDoseInitial(0)
        acq_TS_03.setAccumDose(0)
        acq_TS_03.setTiltAxisAngle(0)
        acq_TS_03.setAngleMin(ts03MinAngle)
        acq_TS_03.setAngleMax(ts03MaxAngle)
        testAcqObjDictReStacked[TS_03] = acq_TS_03

        acq_TS_54 = cls.testAcqObjDict[TS_54].clone()
        acq_TS_54.setDoseInitial(0)
        acq_TS_54.setAccumDose(0)
        acq_TS_54.setTiltAxisAngle(0)
        acq_TS_54.setAngleMin(ts54MinAngle)
        acq_TS_54.setAngleMax(ts54MaxAngle)
        testAcqObjDictReStacked[TS_54] = acq_TS_54
        return testAcqObjDictReStacked

    @classmethod
    def _runCistemEstimateCtf(cls, inTsSet):
        print(magentaStr("\n==> Estimating the CTF with Cistem:"))
        protEstimateCtf = cls.newProtocol(CistemProtTsCtffind,
                                          inputTiltSeries=inTsSet,
                                          lowRes=50,
                                          highRes=5,
                                          minDefocus=5000,
                                          maxDefocus=50000)
        cls.launchProtocol(protEstimateCtf)
        ctfs = getattr(protEstimateCtf, CistemProtTsCtffind._possibleOutputs.CTFs.name, None)
        return ctfs

    @classmethod
    def _genReStackedCtf(cls):
        importedTs = cls._runImportTs()
        # Exclude some views from the TS at metadata level
        cls._excludeSetViews(importedTs, excludedViewsDict=cls.ctfExcludedViewsDict)
        # Re-stack that TS
        reStackedTsSet = cls._runExcludeViewsProt(importedTs)
        # Estimate the CTF using the re-stacked TS
        return cls._runCistemEstimateCtf(reStackedTsSet)

    @classmethod
    def _excludeSetViews(cls, inSet: Union[SetOfTiltSeries, SetOfCTFTomoSeries],
                         excludedViewsDict: Union[dict, None] = None) -> None:
        if not excludedViewsDict:
            excludedViewsDict = cls.excludedViewsDict
        objList = [obj.clone(ignoreAttrs=[]) for obj in inSet]
        for obj in objList:
            cls._excIntermediateSetViews(inSet, obj, excludedViewsDict[obj.getTsId()])

    @staticmethod
    def _excIntermediateSetViews(inSet, obj, excludedViewsList):
        tiList = [ti.clone() for ti in obj]
        for i, ti in enumerate(tiList):
            if i in excludedViewsList:
                ti._objEnabled = False
                obj.update(ti)
        obj.write()
        inSet.update(obj)
        inSet.write()
        inSet.close()

    @classmethod
    def _runExcludeViewsProt(cls, inTsSet, objLabel=None):
        print(magentaStr("\n==> Running the TS exclusion of views:"))
        protExcViews = cls.newProtocol(ProtImodExcludeViews, inputSetOfTiltSeries=inTsSet)
        if objLabel:
            protExcViews.setObjLabel(objLabel)
        cls.launchProtocol(protExcViews)
        outTsSet = getattr(protExcViews, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runNovaCtf(cls,
                    inTsSet: SetOfTiltSeries,
                    inCtfSet: SetOfCTFTomoSeries,
                    tsSetMsg: str,
                    ctfSetMsg: str, defocusTol=200, interpWidth=15):
        print(magentaStr(f"\n==> Running the CTF correction:"
                         f"\n\t- Tilt-series = {tsSetMsg}"
                         f"\n\t- CTFs: {ctfSetMsg}"))
        protNovaCtf = cls.newProtocol(ProtNovaCtf,
                                      inputSetOfTiltSeries=inTsSet,
                                      inputSetOfCtfTomoSeries=inCtfSet,

                                      )

    def testCtfCorrection01(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, importedCtfs,
                                              tsSetMsg=self.UNMODIFIED,
                                              ctfSetMsg=self.UNMODIFIED)
        self._checkInterpTs(tsSetCtfCorr,
                            testAcqObjDict=self.testInterpAcqObjDict,
                            anglesCountDict=self.anglesCountDict)

    def testCtfCorrection02(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        self._excludeSetViews(importedCtfs,
                              excludedViewsDict=self.ctfExcludedViewsDict)  # Excluded some views in the CTF at metadata level
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, importedCtfs,
                                              tsSetMsg=self.UNMODIFIED,
                                              ctfSetMsg=self.EXC_VIEWS)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-57,
                                                          ts54MaxAngle=54)
        self._checkInterpTs(tsSetCtfCorr,
                            testAcqObjDict=testAcqObjDict,
                            anglesCountDict=self.ctfAnglesCountDictExcluded)

    def testCtfCorrection03(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, importedCtfs,
                                              tsSetMsg=self.EXC_VIEWS,
                                              ctfSetMsg=self.UNMODIFIED)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-54,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTs(tsSetCtfCorr,
                            testAcqObjDict=testAcqObjDict,
                            anglesCountDict=self.anglesCountDictExcluded)

    def testCtfCorrection04(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        tsSetReStacked = self._runExcludeViewsProt(tsWithAliBin4)  # Re-stack the TS
        tsSetCtfCorr = self._runCtfCorrection(tsSetReStacked, importedCtfs,
                                              tsSetMsg=self.RE_STACKED,
                                              ctfSetMsg=self.UNMODIFIED)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-54,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTs(tsSetCtfCorr,
                            testAcqObjDict=testAcqObjDict,
                            anglesCountDict=self.anglesCountDictExcluded)

    def testCtfCorrection05(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        self._excludeSetViews(tsWithAliBin4,
                              excludedViewsDict=self.ctfExcludedViewsDict)  # Excluded some views in the CTF at metadata level
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, importedCtfs,
                                              tsSetMsg=self.EXC_VIEWS,
                                              ctfSetMsg=self.EXC_VIEWS)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTs(tsSetCtfCorr,
                            testAcqObjDict=testAcqObjDict,
                            anglesCountDict=self.intersectAnglesCountDictExcluded)

    def testCtfCorrection06(self):
        importedCtfs, tsWithAliBin4 = self._runPrevProts()
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        tsSetReStacked = self._runExcludeViewsProt(tsWithAliBin4)  # Re-stack the TS
        self._excludeSetViews(importedCtfs,
                              excludedViewsDict=self.ctfExcludedViewsDict)  # Excluded some views in the CTF at metadata level
        tsSetCtfCorr = self._runCtfCorrection(tsSetReStacked, importedCtfs,
                                              tsSetMsg=self.RE_STACKED,
                                              ctfSetMsg=self.EXC_VIEWS)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTs(tsSetCtfCorr,
                            testAcqObjDict=testAcqObjDict,
                            anglesCountDict=self.intersectAnglesCountDictExcluded)

    def testCtfCorrection07(self):
        _, tsWithAliBin4 = self._runPrevProts(importCtf=False)
        ctfSetReStacked = self._genReStackedCtf()  # Gen a CTF estimated on a re-stacked TS
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, ctfSetReStacked,
                                              tsSetMsg=self.UNMODIFIED,
                                              ctfSetMsg=self.RE_STACKED)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-57,
                                                          ts54MaxAngle=54)
        self._checkInterpTs(tsSetCtfCorr,
                            testAcqObjDict=testAcqObjDict,
                            anglesCountDict=self.ctfAnglesCountDictExcluded)

    def testCtfCorrection08(self):
        _, tsWithAliBin4 = self._runPrevProts(importCtf=False)
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        ctfSetReStacked = self._genReStackedCtf()  # Gen a CTF estimated on a re-stacked TS
        tsSetCtfCorr = self._runCtfCorrection(tsWithAliBin4, ctfSetReStacked,
                                              tsSetMsg=self.EXC_VIEWS,
                                              ctfSetMsg=self.RE_STACKED)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTs(tsSetCtfCorr,
                            testAcqObjDict=testAcqObjDict,
                            anglesCountDict=self.intersectAnglesCountDictExcluded)

    def testCtfCorrection09(self):
        _, tsWithAliBin4 = self._runPrevProts(importCtf=False)
        self._excludeSetViews(tsWithAliBin4)  # Excluded some views in the TS at metadata level
        tsSetReStacked = self._runExcludeViewsProt(tsWithAliBin4)  # Re-stack the TS
        ctfSetReStacked = self._genReStackedCtf()  # Gen a CTF estimated on a re-stacked TS
        tsSetCtfCorr = self._runCtfCorrection(tsSetReStacked, ctfSetReStacked,
                                              tsSetMsg=self.RE_STACKED,
                                              ctfSetMsg=self.RE_STACKED)
        # Check the results
        testAcqObjDict = self._gentestAcqObjDictReStacked(ts03MinAngle=-51,
                                                          ts03MaxAngle=54,
                                                          ts54MinAngle=-54,
                                                          ts54MaxAngle=51)
        self._checkInterpTs(tsSetCtfCorr,
                            testAcqObjDict=testAcqObjDict,
                            anglesCountDict=self.intersectAnglesCountDictExcluded)