# *****************************************************************************
#
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
# *****************************************************************************
import os

from pyworkflow.tests import setupTestProject, BaseTest, DataSet
from pwem.emlib.image import ImageHandler
import tomo.protocols
import imod.protocols

from ..protocols import *


class TestNovaCtfBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def _runImportTiltSeries(cls, filesPath, pattern, voltage, magnification, sphericalAberration, amplitudeContrast,
                             samplingRate, doseInitial, dosePerFrame, anglesFrom=0, minAngle=0.0, maxAngle=0.0,
                             stepAngle=1.0, tiltAxisAngle=87.2):
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
                                           stepAngle=stepAngle,
                                           tiltAxisAngle=tiltAxisAngle)

        cls.launchProtocol(cls.protImportTS)

        return cls.protImportTS

    @classmethod
    def _runCTFEstimation(cls, inputSoTS, defocusTol, expectedDefocusOrigin, expectedDefocusValue,
                          axisAngle, leftDefTol, rightDefTol, tileSize, angleStep, angleRange,
                          startFreq, endFreq, extraZerosToFit, skipAstigmaticViews, searchAstigmatism,
                          findAstigPhaseCutonToggle, phaseShiftAstigmatism, cutOnFrequencyAstigmatism,
                          minimumViewsAstigmatism, minimumViewsPhaseShift, numberSectorsAstigmatism,
                          maximumAstigmatism):
        cls.protCTFEstimation = cls.newProtocol(imod.protocols.ProtImodAutomaticCtfEstimation,
                                                inputSet=inputSoTS,
                                                defocusTol=defocusTol,
                                                expectedDefocusOrigin=expectedDefocusOrigin,
                                                expectedDefocusValue=expectedDefocusValue,
                                                axisAngle=axisAngle,
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
    def _runComputeCtfArray(cls, inputSetOfTiltSeries, inputSetOfCtfTomoSeries, tomoThickness, tomoShift,
                            defocusStep, correctionType, correctAstigmatism):
        cls.protCTFReconstruction = cls.newProtocol(ProtNovaCtfTomoDefocus,
                                                    inputSetOfTiltSeries=inputSetOfTiltSeries,
                                                    inputSetOfCtfTomoSeries=inputSetOfCtfTomoSeries,
                                                    tomoThickness=tomoThickness,
                                                    tomoShift=tomoShift,
                                                    defocusStep=defocusStep,
                                                    correctionType=correctionType,
                                                    correctAstigmatism=correctAstigmatism)

        cls.launchProtocol(cls.protCTFReconstruction)

        return cls.protCTFReconstruction


class TestNovaCtfReconstructionWorkflow(TestNovaCtfBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.inputDataSet = DataSet.getDataSet('novaCtfTestData')
        cls.inputSoTS = cls.inputDataSet.getFile('tsCtf')

        cls.protImportTS = cls._runImportTiltSeries(
            filesPath=os.path.dirname(cls.inputSoTS),
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

        cls.protCTFEstimation = cls._runCTFEstimation(
            inputSoTS=cls.protImportTS.outputTiltSeries,
            defocusTol=200.0,
            expectedDefocusOrigin=0,
            expectedDefocusValue=6000,
            axisAngle=0.0,
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

        cls.protCTFCompute = cls._runComputeCtfArray(
            inputSetOfTiltSeries=cls.protImportTS.outputTiltSeries,
            inputSetOfCtfTomoSeries=cls.protCTFEstimation.CTFTomoSeries,
            tomoThickness=20,
            tomoShift=0,
            defocusStep=15,
            correctionType=0,
            correctAstigmatism=1)

        cls.protReconstruct = cls.newProtocol(ProtNovaCtfTomoReconstruction,
                                              protTomoCtfDefocus=cls.protCTFCompute,
                                              applyAlignment=False)
        cls.launchProtocol(cls.protReconstruct)

        return cls.protReconstruct

    def test_tomoReconstructionOutput(self):
        self.assertIsNotNone(self.protReconstruct.outputSetOfTomograms)
        self.assertTrue(self.protReconstruct.outputSetOfTomograms.getSize() == 1)

        ih = ImageHandler()
        self.assertTrue(
            ih.getDimensions(self.protReconstruct.outputSetOfTomograms.getFirstItem()) ==
            (960, 928, 20, 1))
