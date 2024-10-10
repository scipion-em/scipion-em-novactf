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
from pyworkflow.utils import magentaStr
from pwem.emlib.image import ImageHandler
import tomo.protocols
import imod.protocols

from ..protocols import ProtNovaCtfReconstruction, ProtNovaCtfDefocus


class TestNovaCtfBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def _runImportTiltSeries(cls, filesPath, pattern, voltage, magnification,
                             samplingRate, dosePerFrame, anglesFrom=0, minAngle=0.0, maxAngle=0.0,
                             stepAngle=1.0, tiltAxisAngle=0.0):
        print(magentaStr("\n==> Importing data - tilt series:"))
        cls.protImportTS = cls.newProtocol(tomo.protocols.ProtImportTs,
                                           filesPath=filesPath,
                                           filesPattern=pattern,
                                           voltage=voltage,
                                           anglesFrom=anglesFrom,
                                           magnification=magnification,
                                           samplingRate=samplingRate,
                                           dosePerFrame=dosePerFrame,
                                           minAngle=minAngle,
                                           maxAngle=maxAngle,
                                           stepAngle=stepAngle,
                                           tiltAxisAngle=tiltAxisAngle)

        cls.launchProtocol(cls.protImportTS)

        return cls.protImportTS

    @classmethod
    def _runCTFEstimation(cls, inputSoTS, expectedDefocusOrigin, angleRange,
                          expectedDefocusValue, searchAstigmatism):
        print(magentaStr("\n==> Running CTF estimation with IMOD:"))
        cls.protCTFEstimation = cls.newProtocol(imod.protocols.ProtImodAutomaticCtfEstimation,
                                                inputSetOfTiltSeries=inputSoTS,
                                                expectedDefocusOrigin=expectedDefocusOrigin,
                                                expectedDefocusValue=expectedDefocusValue,
                                                angleRange=angleRange,
                                                searchAstigmatism=searchAstigmatism)

        cls.launchProtocol(cls.protCTFEstimation)

        return cls.protCTFEstimation

    @classmethod
    def _runComputeCtfArray(cls, inputSetOfTiltSeries, inputSetOfCtfTomoSeries, tomoThickness, tomoShift,
                            defocusStep, correctionType, correctAstigmatism):
        print(magentaStr("\n==> Testing novactf - compute defocus array:"))
        cls.protComputeDefocus = cls.newProtocol(ProtNovaCtfDefocus,
                                                 inputSetOfTiltSeries=inputSetOfTiltSeries,
                                                 inputSetOfCtfTomoSeries=inputSetOfCtfTomoSeries,
                                                 tomoThickness=tomoThickness,
                                                 tomoShift=tomoShift,
                                                 defocusStep=defocusStep,
                                                 correctionType=correctionType,
                                                 correctAstigmatism=correctAstigmatism)

        cls.launchProtocol(cls.protComputeDefocus)

        return cls.protComputeDefocus


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
            samplingRate=8.8,
            dosePerFrame=0.3,
            minAngle=-60.0,
            maxAngle=60.0,
            stepAngle=3.0,
            tiltAxisAngle=2.8)

        cls.protCTFEstimation = cls._runCTFEstimation(
            inputSoTS=cls.protImportTS.outputTiltSeries,
            expectedDefocusOrigin=0,
            expectedDefocusValue=6000,
            angleRange=20,
            searchAstigmatism=0)

        cls.protCTFCompute = cls._runComputeCtfArray(
            inputSetOfTiltSeries=cls.protImportTS.outputTiltSeries,
            inputSetOfCtfTomoSeries=cls.protCTFEstimation.CTFTomoSeries,
            tomoThickness=20,
            tomoShift=0,
            defocusStep=50,
            correctionType=0,
            correctAstigmatism=False)

        print(magentaStr("\n==> Testing novactf - 3D CTF correction and reconstruction:"))
        cls.protReconstruct = cls.newProtocol(ProtNovaCtfReconstruction,
                                              protNovaCtfDefocus=cls.protCTFCompute,
                                              applyAlignment=False)
        cls.launchProtocol(cls.protReconstruct)

        return cls.protReconstruct

    def test_tomoReconstructionOutput(self):
        outputName = ProtNovaCtfReconstruction._possibleOutputs.Tomograms.name
        output = getattr(self.protReconstruct, outputName)
        self.assertIsNotNone(output)
        self.assertTrue(output.getSize() == 1)

        ih = ImageHandler()
        self.assertTrue(
            ih.getDimensions(output.getFirstItem()) == (960, 928, 20, 1))
