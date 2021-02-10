# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import *
from pyworkflow.tests.test_utils import wait
import imod
from novactf.protocols import *
from pwem.emlib.image import ImageHandler
import tomo


class TestNovaCtfBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def _runImportTiltSeries(cls, filesPath, pattern, voltage, magnification, sphericalAberration, amplitudeContrast,
                             samplingRate, doseInitial, dosePerFrame, anglesFrom=0, minAngle=0.0, maxAngle=0.0,
                             stepAngle=1.0):
        cls.protImportTS = cls.newProtocol(tomo.protocols.ProtImportTs,
                                           filesPath=filesPath,
                                           filesPattern=pattern,
                                           voltage=voltage,
                                           anglesFrom=anglesFrom,
                                           magnification=magnification,
                                           sphericalAberration=sphericalAberration,
                                           amplitudeContrast=amplitudeContrast,
                                           samplingRate=samplingRate,
                                           doseInitial=doseInitial,
                                           dosePerFrame=dosePerFrame,
                                           minAngle=minAngle,
                                           maxAngle=maxAngle,
                                           stepAngle=stepAngle)
        cls.launchProtocol(cls.protImportTS)
        return cls.protImportTS

    @classmethod
    def _runCTFEstimation(cls, inputSoTS, defocusTol, expectedDefocusOrigin, expectedDefocusValue, expectedDefocusFile,
                          axisAngle, interactiveMode, leftDefTol, rightDefTol, tileSize, angleStep, angleRange,
                          startFreq, endFreq, extraZerosToFit, skipAstigmaticViews, searchAstigmatism,
                          findAstigPhaseCutonToggle, phaseShiftAstigmatism, cutOnFrequencyAstigmatism,
                          minimumViewsAstigmatism, minimumViewsPhaseShift, numberSectorsAstigmatism,
                          maximumAstigmatism):
        cls.protCTFEstimation = cls.newProtocol(imod.protocols.ProtImodCtfEstimation,
                                                inputSetOfTiltSeries=inputSoTS,
                                                defocusTol=defocusTol,
                                                expectedDefocusOrigin=expectedDefocusOrigin,
                                                expectedDefocusValue=expectedDefocusValue,
                                                expectedDefocusFile=expectedDefocusFile,
                                                axisAngle=axisAngle,
                                                interactiveMode=interactiveMode,
                                                leftDefTol=leftDefTol,
                                                rightDefTol=rightDefTol,
                                                tileSize=tileSize,
                                                angleStep=angleStep,
                                                angleRange=angleRange,
                                                startFreq=startFreq,
                                                endFreq=endFreq,
                                                extraZerosToFit=extraZerosToFit,
                                                skipAstigmaticViews=skipAstigmaticViews,
                                                searchAstigmatism=searchAstigmatism,
                                                findAstigPhaseCutonToggle=findAstigPhaseCutonToggle,
                                                phaseShiftAstigmatism=phaseShiftAstigmatism,
                                                cutOnFrequencyAstigmatism=cutOnFrequencyAstigmatism,
                                                minimumViewsAstigmatism=minimumViewsAstigmatism,
                                                minimumViewsPhaseShift=minimumViewsPhaseShift,
                                                numberSectorsAstigmatism=numberSectorsAstigmatism,
                                                maximumAstigmatism=maximumAstigmatism)
        cls.launchProtocol(cls.protCTFEstimation)
        return cls.protCTFEstimation

    @classmethod
    def _runCtfReconstruction(cls, inputSoTS, ctfEstimationType, protImodCtfEstimation, tomoThickness, tomoShift,
                              defocusStep, correctionType, radialFirstParameter, radialSecondParameter):
        cls.protCTFReconstruction = cls.newProtocol(ProtNovaCtfTomoDefocus,
                                                    inputSetOfTiltSeries=inputSoTS,
                                                    ctfEstimationType=ctfEstimationType,
                                                    protImodCtfEstimation=protImodCtfEstimation,
                                                    tomoThickness=tomoThickness,
                                                    tomoShift=tomoShift,
                                                    defocusStep=defocusStep,
                                                    correctionType=correctionType,
                                                    radialFirstParameter=radialFirstParameter,
                                                    radialSecondParameter=radialSecondParameter)
        cls.launchProtocol(cls.protCTFReconstruction)
        return cls.protCTFReconstruction


class TestNovaCtfReconstructionWorkflow(TestNovaCtfBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.inputDataSet = DataSet.getDataSet('novaCtfTestData')
        cls.inputSoTS = cls.inputDataSet.getFile('tsCtf')

        cls.thicknessTomo = 20.0

        cls.protImportTS = cls._runImportTiltSeries(filesPath=os.path.split(cls.inputSoTS)[0],
                                                    pattern="tomo1_bin4.mrc",
                                                    anglesFrom=0,
                                                    voltage=300,
                                                    magnification=50000,
                                                    sphericalAberration=2.7,
                                                    amplitudeContrast=0.07,
                                                    samplingRate=8.8,
                                                    doseInitial=0,
                                                    dosePerFrame=0.3,
                                                    minAngle=-60.0,
                                                    maxAngle=60.0,
                                                    stepAngle=3.0)

        cls.protCTFEstimation = cls._runCTFEstimation(inputSoTS=cls.protImportTS.outputTiltSeries,
                                                      defocusTol=200.0,
                                                      expectedDefocusOrigin=0,
                                                      expectedDefocusValue=6000,
                                                      expectedDefocusFile="",
                                                      axisAngle=0.0,
                                                      interactiveMode=1,
                                                      leftDefTol=2000.0,
                                                      rightDefTol=2000.0,
                                                      tileSize=256,
                                                      angleStep=2.0,
                                                      angleRange=20.0,
                                                      startFreq=0.0,
                                                      endFreq=0.0,
                                                      extraZerosToFit=0.0,
                                                      skipAstigmaticViews=1,
                                                      searchAstigmatism=1,
                                                      findAstigPhaseCutonToggle=1,
                                                      phaseShiftAstigmatism=0,
                                                      cutOnFrequencyAstigmatism=0,
                                                      minimumViewsAstigmatism=3,
                                                      minimumViewsPhaseShift=1,
                                                      numberSectorsAstigmatism=36,
                                                      maximumAstigmatism=1.2)

        cls.protCTFReconstruction = cls._runCtfReconstruction(inputSoTS=None,
                                                              ctfEstimationType=0,
                                                              protImodCtfEstimation=cls.protCTFEstimation,
                                                              tomoThickness=cls.thicknessTomo,
                                                              tomoShift=0,
                                                              defocusStep=15,
                                                              correctionType=0,
                                                              radialFirstParameter=0.3,
                                                              radialSecondParameter=0.05)

        wait(condition=lambda: not (cls.proj.getRuns() == 4 and cls.proj.getRuns()[3].isFinished()),
             timeout=600)

    def test_tomoReconstructionOutput(self):
        protTomoReconstruction = self.proj.getRuns()[3]
        self.assertIsNotNone(protTomoReconstruction.outputSetOfTomograms)

    def test_tomoReconstructionOutputSize(self):
        protTomoReconstruction = self.proj.getRuns()[3]
        self.assertTrue(protTomoReconstruction.outputSetOfTomograms.getSize() == 1)

    def test_tomoReconstructionOutputTomogramDimensions(self):
        ih = ImageHandler()
        protTomoReconstruction = self.proj.getRuns()[3]
        self.assertTrue(
            ih.getDimensions(protTomoReconstruction.outputSetOfTomograms.getFirstItem()) ==
            (960, 928, self.thicknessTomo, 1))
