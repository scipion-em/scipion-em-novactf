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

import os
import numpy as np
import math
import pyworkflow as pw
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem.protocols import EMProtocol
from pwem.emlib.image import ImageHandler
from tomo.protocols import ProtTomoBase
from tomo.convert import writeTiStack
from tomo.objects import Tomogram, TiltSeries
from novactf import Plugin
from imod import Plugin as imodPlugin


class ProtNovaCtfTomoReconstruction(EMProtocol, ProtTomoBase):
    """
    Tomogram reconstruction and ctf correction procedure based on the novaCTF procedure.

    More info:
            https://github.com/turonova/novaCTF
    """

    _label = 'tomo ctf reconstruction'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        ProtTomoBase.__init__(self)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('protTomoCtfDefocus',
                      params.PointerParam,
                      label="IMOD CTF estimation run",
                      pointerClass='ProtNovaCtfTomoDefocus',
                      help='Select the previous IMOD CTF estimation run.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):

        for ts in self.inputSetOfTiltSeries.get():
            allCtfId = []

            for counterCtf in range(0, numberOfIntermediateStacks + 1):
                ctfId = self._insertFunctionStep('computeCtfCorrectionStep',
                                                 ts.getObjId(),
                                                 counterCtf)
                allCtfId.append(ctfId)

            allFlipId = []

            for counterFlip in range(0, numberOfIntermediateStacks + 1):
                flipId = self._insertFunctionStep('computeFlipStep',
                                                  ts.getObjId(),
                                                  counterFlip,
                                                  prerequisites=allCtfId)
                allFlipId.append(flipId)

            allFilterId = []

            for counterFilter in range(0, numberOfIntermediateStacks + 1):
                filterId = self._insertFunctionStep('computeFilteringStep',
                                                    ts.getObjId(),
                                                    counterFilter,
                                                    prerequisites=allFlipId)
                allFilterId.append(filterId)

            reconstructionId = self._insertFunctionStep('computeReconstructionStep',
                                                        ts.getObjId(),
                                                        prerequisites=allFilterId)

            self._insertFunctionStep('createOutputStep',
                                     ts.getObjId(),
                                     prerequisites=[reconstructionId])

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)
        outputTsFileName = os.path.join(tmpPrefix, "%s.st" % tsId)

        """Apply the transformation form the input tilt-series"""
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, "%s.tlt" % tsId)
        ts.generateTltFile(angleFilePath)

    def generateImodDefocusFile(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        outputDefocusFile = self.getDefocusFile(ts)

        outputDefocusFilePrefix = self.protImodCtfEstimation.get()._getExtraPath(tsId)
        path.copyFile(os.path.join(outputDefocusFilePrefix, "%s.defocus" % tsId), outputDefocusFile)

    def generateCtffindDefocusFile(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        outputDefocusFile = self.getDefocusFile(ts)
        defocusInfo = []

        for ti in ts:
            ctfModel = ti.getCTF()
            defU, defV, defAngle = ctfModel.getDefocus()
            azimutAng = defAngle - 180  # Conversion to CTFFind angle
            phaseFlip = 0 if ctfModel.getPhaseShift() is None else ctfModel.getPhaseShift()
            fitQuality = ctfModel.getFitQuality()  # CTFFind cross correlation
            resolution = ctfModel.getResolution()  # CTFFind spacing up to which CTF rings were fit successfully
            tiDefocusInfoVector = [defU, defV, azimutAng, phaseFlip, fitQuality, resolution]
            defocusInfo.append(tiDefocusInfoVector)

        with open(outputDefocusFile, 'w') as f:
            f.writelines("# Columns: #1 - micrograph number; #2 - defocus 1 [Angstroms]; #3 - defocus 2; "
                         "#4 - azimuth of astigmatism; #5 - additional phase shift [radians]; #6 - cross correlation; "
                         "#7 - spacing (in Angstroms) up to which CTF rings were fit successfully\n")

            for counter, vector in enumerate(defocusInfo):
                line = "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n" % \
                       (counter + 1, vector[0], vector[1], vector[2], vector[3], vector[4], vector[5])
                f.writelines(line)

    def computeDefocusStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(ts.getTsId())

        paramsDefocus = {
            'Algorithm': "defocus",
            'InputProjections': os.path.join(tmpPrefix, "%s.st" % tsId),
            'FullImage': str(ts.getFirstItem().getDim()[0]) + "," + str(ts.getFirstItem().getDim()[1]),
            'Thickness': self.tomoThickness.get(),
            'TiltFile': os.path.join(tmpPrefix, "%s.tlt" % tsId),
            'Shift': "0.0," + str(self.tomoShift.get()),
            'CorrectionType': self.getCorrectionType(),
            'DefocusFileFormat': "ctffind4",
            'CorrectAstigmatism': 1,
            'DefocusFile': self.getDefocusFile(ts),
            'PixelSize': self.inputSetOfTiltSeries.get().getSamplingRate() / 10,
            'DefocusStep': self.defocusStep.get()
        }

        argsDefocus = "-Algorithm %(Algorithm)s " \
                      "-InputProjections %(InputProjections)s " \
                      "-FULLIMAGE %(FullImage)s " \
                      "-THICKNESS %(Thickness)d " \
                      "-TILTFILE %(TiltFile)s " \
                      "-SHIFT %(Shift)s " \
                      "-CorrectionType %(CorrectionType)s " \
                      "-DefocusFileFormat %(DefocusFileFormat)s " \
                      "-CorrectAstigmatism %(CorrectAstigmatism)d " \
                      "-DefocusFile %(DefocusFile)s " \
                      "-PixelSize %(PixelSize)s " \
                      "-DefocusStep %(DefocusStep)d"

        Plugin.runNovactf(self, 'novaCTF', argsDefocus % paramsDefocus)

    def computeCtfCorrectionStep(self, tsObjId, counter):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(ts.getTsId())
        tltFilePath = os.path.join(tmpPrefix, "%s.tlt" % tsId)
        outputFilePath = os.path.join(tmpPrefix, "%s.st_" % tsId)
        defocusFilePath = os.path.join(tmpPrefix, "%s.defocus_" % tsId)

        paramsCtfCorrection = {
            'Algorithm': "ctfCorrection",
            'InputProjections': os.path.join(tmpPrefix, "%s.st" % tsId),
            'OutputFile': outputFilePath + str(counter),
            'DefocusFile': defocusFilePath + str(counter),
            'TiltFile': tltFilePath,
            'CorrectionType': self.getCorrectionType(),
            'DefocusFileFormat': "ctffind4",
            'CorrectAstigmatism': 1,
            'PixelSize': self.inputSetOfTiltSeries.get().getSamplingRate() / 10,
            'AmplitudeContrast': self.inputSetOfTiltSeries.get().getAcquisition().getAmplitudeContrast(),
            'SphericalAberration': self.inputSetOfTiltSeries.get().getAcquisition().getSphericalAberration(),
            'Voltage': self.inputSetOfTiltSeries.get().getAcquisition().getVoltage()
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

    def computeFlipStep(self, tsObjId, counter):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(ts.getTsId())
        inputFilePath = os.path.join(tmpPrefix, "%s.st_" % tsId)
        outputFilePath = os.path.join(tmpPrefix, "%s_flip.st_" % tsId)

        argsFlip = "flipyz " + inputFilePath + str(counter) + " " + outputFilePath + str(counter)
        imodPlugin.runImod(self, 'clip', argsFlip)

    def computeFilteringStep(self, tsObjId, counter):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(ts.getTsId())
        flippedFilePath = os.path.join(tmpPrefix, "%s_flip.st_" % tsId)
        outputFilePath = os.path.join(tmpPrefix, "%s_flip_filter.st_" % tsId)
        tltFilePath = os.path.join(tmpPrefix, "%s.tlt" % tsId)

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
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(ts.getTsId())
        tmpPrefix = self._getTmpPath(ts.getTsId())
        outputFilePathFlipped = os.path.join(tmpPrefix, "%s.mrc" % tsId)
        tltFilePath = os.path.join(tmpPrefix, "%s.tlt" % tsId)

        params3dctf = {
            'Algorithm': "3dctf",
            'InputProjections': os.path.join(tmpPrefix, "%s_flip_filter.st" % tsId),
            'OutputFile': outputFilePathFlipped,
            'TiltFile': tltFilePath,
            'Thickness': self.tomoThickness.get(),
            'FullImage': str(ts.getFirstItem().getDim()[0]) + "," + str(ts.getFirstItem().getDim()[1]),
            'Shift': "0.0," + str(self.tomoShift.get()),
            'PixelSize': self.inputSetOfTiltSeries.get().getSamplingRate() / 10,
            'DefocusStep': self.defocusStep.get()
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

        outputFilePath = os.path.join(extraPrefix, "%s.mrc" % tsId)

        argsTrimvol = outputFilePathFlipped + " "
        argsTrimvol += outputFilePath + " "
        argsTrimvol += "-yz "

        imodPlugin.runImod(self, 'trimvol', argsTrimvol)

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        """Keep defocus file"""
        tmpDefocusFile = os.path.join(self._getTmpPath(tsId), tsId + ".defocus")
        extraDefocusFile = os.path.join(self._getExtraPath(tsId), tsId + ".defocus")
        path.moveFile(tmpDefocusFile, extraDefocusFile)

        """Remove intermediate files. Necessary for big sets of tilt-series"""
        path.cleanPath(self._getTmpPath(tsId))

        """Generate output set"""
        outputSetOfTomograms = self.getOutputSetOfTomograms()
        extraPrefix = self._getExtraPath(tsId)
        newTomogram = Tomogram()
        newTomogram.copyInfo(ts)
        newTomogram.setLocation(os.path.join(extraPrefix, "%s.mrc" % tsId))
        outputSetOfTomograms.append(newTomogram)
        outputSetOfTomograms.update(newTomogram)
        outputSetOfTomograms.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfTomograms(self):
        if not hasattr(self, "outputSetOfTomograms"):
            outputSetOfTomograms = self._createSetOfTomograms()
            outputSetOfTomograms.copyInfo(self.inputSetOfTiltSeries.get())
            self._defineOutputs(outputSetOfTomograms=outputSetOfTomograms)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTomograms)

        return self.outputSetOfTomograms

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        if self.ctfEstimationType.get() == 1 and \
                not self.inputSetOfTiltSeries.get().getFirstItem().getFirstItem().hasCTF():
            validateMsgs.append("You need to generate an estimation of the CTF associated to the tilt series to "
                                "calculate its corrected reconstruction")

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfTomograms'):
            summary.append("Input Tilt-Series: %d.\nCTF corrected reconstructions calculated: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
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
