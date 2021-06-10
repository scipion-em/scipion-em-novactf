# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
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

import os

from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.exceptions import PyworkflowException
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram, TomoAcquisition
from novactf import Plugin
from novactf.protocols import ProtNovaCtfTomoDefocus
from imod import Plugin as imodPlugin


class ProtNovaCtfTomoReconstruction(EMProtocol, ProtTomoBase):
    """
    Tomogram reconstruction with ctf correction procedure based on the novaCTF procedure.

    More info:
            https://github.com/turonova/novaCTF
    """

    _label = 'tomo ctf reconstruction'
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        ProtTomoBase.__init__(self)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('protTomoCtfDefocus',
                      params.PointerParam,
                      label="NovaCtf defocus estimation run",
                      pointerClass='ProtNovaCtfTomoDefocus',
                      help='Select the previous NovaCtf defocus estimation run.')

        ProtNovaCtfTomoDefocus.defineFilterParameters(form)

        form.addParallelSection(threads=4, mpi=1)

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):

        allCreateOutputId = []

        for index, ts in enumerate(self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get()):
            convertInputId = self._insertFunctionStep(self.convertInputStep,
                                                      ts.getObjId(),
                                                      prerequisites=[])

            intermediateStacksId = []

            for counter in range(0, self.protTomoCtfDefocus.get().numberOfIntermediateStacks[index].get()):
                ctfId = self._insertFunctionStep(self.processIntermediateStacksStep,
                                                 ts.getObjId(),
                                                 counter,
                                                 prerequisites=[convertInputId])
                intermediateStacksId.append(ctfId)

            reconstructionId = self._insertFunctionStep(self.computeReconstructionStep,
                                                        ts.getObjId(),
                                                        prerequisites=intermediateStacksId)

            createOutputId = self._insertFunctionStep(self.createOutputStep,
                                                      ts.getObjId(),
                                                      prerequisites=[reconstructionId])
            allCreateOutputId.append(createOutputId)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 prerequisites=allCreateOutputId)

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        with self._lock:
            ts = self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get()[tsObjId]
            ti = ts.getFirstItem()

        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        print(ti.getDim())

        """Apply the transformation form the input tilt-series"""
        outputTsFileName = os.path.join(tmpPrefix, ti.parseFileName())

        print(outputTsFileName)
        with self._lock:
            ts.applyTransform(outputTsFileName)
            """Generate angle file"""
            outputTltFileName = os.path.join(tmpPrefix, ti.parseFileName(extension=".tlt"))
            ts.generateTltFile(outputTltFileName)

    def processIntermediateStacksStep(self, tsObjId, counter):
        with self._lock:
            ts = self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get()[tsObjId]
            ti = ts.getFirstItem()

        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(tsId)
        extraPrefixPreviousProt = self.protTomoCtfDefocus.get()._getExtraPath(tsId)

        defocusFilePath = os.path.join(extraPrefixPreviousProt, ti.parseFileName(extension=".defocus_"))
        tltFilePath = os.path.join(tmpPrefix, ti.parseFileName(extension=".tlt"))
        outputFilePath = os.path.join(tmpPrefix, ti.parseFileName(extension=".st_"))

        # CTF correction step
        paramsCtfCorrection = {
            'Algorithm': "ctfCorrection",
            'InputProjections': os.path.join(tmpPrefix, ti.parseFileName()),
            'OutputFile': outputFilePath + str(counter),
            'DefocusFile': defocusFilePath + str(counter),
            'TiltFile': tltFilePath,
            'CorrectionType': self.getCorrectionType(),
            'DefocusFileFormat': "imod",
            'CorrectAstigmatism': self.protTomoCtfDefocus.get().correctAstigmatism.get(),
            'PixelSize': self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get().getSamplingRate() / 10,
            'AmplitudeContrast':
                self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get().getAcquisition().getAmplitudeContrast(),
            'SphericalAberration':
                self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get().getAcquisition().getSphericalAberration(),
            'Voltage': self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get().getAcquisition().getVoltage()
        }

        argsCtfCorrection = "-Algorithm %(Algorithm)s " \
                            "-InputProjections %(InputProjections)s " \
                            "-OutputFile %(OutputFile)s " \
                            "-DefocusFile %(DefocusFile)s " \
                            "-TILTFILE %(TiltFile)s " \
                            "-CorrectionType %(CorrectionType)s " \
                            "-DefocusFileFormat %(DefocusFileFormat)s " \
                            "-CorrectAstigmatism %(CorrectAstigmatism)s " \
                            "-PixelSize %(PixelSize)f " \
                            "-AmplitudeContrast %(AmplitudeContrast)f " \
                            "-Cs %(SphericalAberration)f " \
                            "-Volt %(Voltage)d"

        Plugin.runNovactf(self, 'novaCTF', argsCtfCorrection % paramsCtfCorrection)

        # Flipping step
        paramsClip = {
            'inputFilePath': os.path.join(tmpPrefix, ti.parseFileName(extension=".st_" + str(counter))),
            'outputFilePath': os.path.join(tmpPrefix, ti.parseFileName(suffix="_flip",
                                                                                      extension=".st_" + str(counter))),
        }

        argsClip = "flipyz " \
                   "%(inputFilePath)s " \
                   "%(outputFilePath)s"

        imodPlugin.runImod(self, 'clip', argsClip % paramsClip)

        # Filtering step
        flippedFilePath = os.path.join(tmpPrefix,
                                       ti.parseFileName(suffix="_flip", extension=".st_"))
        outputFilePath = os.path.join(tmpPrefix,
                                      ti.parseFileName(suffix="_flip_filter", extension=".st_"))
        tltFilePath = os.path.join(tmpPrefix, ti.parseFileName(extension=".tlt"))

        paramsFilterProjections = {
            'Algorithm': "filterProjections",
            'InputProjections': flippedFilePath + str(counter),
            'OutputFile': outputFilePath + str(counter),
            'TiltFile': tltFilePath,
            'StackOrientation': "xz",
            'Radial': str(self.radialFirstParameter.get()) + "," + str(self.radialSecondParameter.get())
        }

        argsFilterProjections = "-Algorithm %(Algorithm)s " \
                                "-InputProjections %(InputProjections)s " \
                                "-OutputFile %(OutputFile)s " \
                                "-TILTFILE %(TiltFile)s " \
                                "-StackOrientation %(StackOrientation)s " \
                                "-RADIAL %(Radial)s"

        Plugin.runNovactf(self, 'novaCTF', argsFilterProjections % paramsFilterProjections)

    def computeReconstructionStep(self, tsObjId):
        ts = self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = ts.getFirstItem()

        outputFileNameFlip = firstItem.parseFileName(suffix="_flip", extension=".mrc")
        outputFileName = firstItem.parseFileName(extension=".mrc")

        print("*** otuput file name")
        print(outputFileNameFlip)
        outputFilePathFlipped = os.path.join(tmpPrefix, outputFileNameFlip)
        print(outputFilePathFlipped)

        tltFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))

        params3dctf = {
            'Algorithm': "3dctf",
            'InputProjections': os.path.join(tmpPrefix, firstItem.parseFileName(suffix="_flip_filter",
                                                                                        extension=".st")),
            'OutputFile': outputFilePathFlipped,
            'TiltFile': tltFilePath,
            'Thickness': self.protTomoCtfDefocus.get().tomoThickness.get(),
            'FullImage': str(firstItem.getDim()[0]) + "," + str(firstItem.getDim()[1]),
            'Shift': "0.0," + str(self.protTomoCtfDefocus.get().tomoShift.get()),
            'PixelSize': self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get().getSamplingRate() / 10,
            'DefocusStep': self.protTomoCtfDefocus.get().defocusStep.get()
        }

        args3dctf = "-Algorithm %(Algorithm)s " \
                    "-InputProjections %(InputProjections)s " \
                    "-OutputFile %(OutputFile)s " \
                    "-TILTFILE %(TiltFile)s " \
                    "-THICKNESS %(Thickness)d " \
                    "-FULLIMAGE %(FullImage)s " \
                    "-SHIFT %(Shift)s " \
                    "-PixelSize %(PixelSize)f " \
                    "-DefocusStep %(DefocusStep)d"

        Plugin.runNovactf(self, 'novaCTF', args3dctf % params3dctf)

        paramsTrimvol = {
            'inputFilePath': outputFilePathFlipped,
            'outputFilePath': os.path.join(extraPrefix, outputFileName)
        }

        argsTrimvol = "%(inputFilePath)s " \
                      "%(outputFilePath)s " \
                      "-yz "

        imodPlugin.runImod(self, 'trimvol', argsTrimvol % paramsTrimvol)

    def createOutputStep(self, tsObjId):
        ts = self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get()[tsObjId]

        tsId = ts.getTsId()
        angleMax = ts[ts.getSize()].getTiltAngle()
        angleStepAverage = self.getAngleStepFromSeries(ts)

        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem().clone()

        """Remove intermediate files. Necessary for big sets of tilt-series"""
        # path.cleanPath(self._getTmpPath(tsId))

        """Generate output set"""
        outputSetOfTomograms = self.getOutputSetOfTomograms()

        newTomogram = Tomogram()
        newTomogram.setLocation(os.path.join(extraPrefix, firstItem.parseFileName(extension=".mrc")))

        if not os.path.exists(newTomogram.getFileName()):
            raise PyworkflowException("%s does not exist."
                                      "\ntsObj: %d"
                                      "\ntsId: %s"
                                      "\nobject tsId %s"
                                      % (newTomogram.getFileName(), tsObjId, tsId, firstItem.getTsId()))

        # Set tomogram origin
        origin = Transform()
        sr = self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get().getSamplingRate()
        origin.setShifts(firstItem.getXDim() / -2. * sr,
                         firstItem.getYDim() / -2. * sr,
                         self.protTomoCtfDefocus.get().tomoThickness.get() / -2 * sr)
        newTomogram.setOrigin(origin)

        # Set tomogram acquisition
        acquisition = TomoAcquisition()
        acquisition.setAngleMin(firstItem.getTiltAngle())
        acquisition.setAngleMax(angleMax)
        acquisition.setStep(angleStepAverage)
        newTomogram.setAcquisition(acquisition)

        outputSetOfTomograms.append(newTomogram)
        outputSetOfTomograms.write()
        self._store()

    def closeOutputSetsStep(self):
        self.getOutputSetOfTomograms().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getCorrectionType(self):
        if self.protTomoCtfDefocus.get().correctionType.get() == 0:
            correctionType = "phaseflip"
        elif self.protTomoCtfDefocus.get().correctionType.get() == 1:
            correctionType = "multiplication"

        return correctionType

    def getOutputSetOfTomograms(self):
        if hasattr(self, "outputSetOfTomograms"):
            self.outputSetOfTomograms.enableAppend()
        else:
            outputSetOfTomograms = self._createSetOfTomograms()
            outputSetOfTomograms.copyInfo(self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get())
            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTomograms=outputSetOfTomograms)
            self._defineSourceRelation(self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get(), outputSetOfTomograms)

        return self.outputSetOfTomograms

    @staticmethod
    def getAngleStepFromSeries(ts):
        """ This method return the average angles step from a series. """

        angleStepAverage = 0
        for i in range(1, ts.getSize()):
            angleStepAverage += abs(ts[i].getTiltAngle() - ts[i + 1].getTiltAngle())

        angleStepAverage /= ts.getSize() - 1

        return angleStepAverage

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfTomograms'):
            summary.append("Input Tilt-Series: %d.\nCTF corrected reconstructions calculated: %d.\n"
                           % (self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get().getSize(),
                              self.outputSetOfTomograms.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfTomograms'):
            methods.append("%d CTF corrected tomograms have been calculated using the NovaCtf software.\n"
                           % (self.outputSetOfTomograms.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
