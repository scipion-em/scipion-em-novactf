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
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram
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

        for index, ts in enumerate(self.protTomoCtfDefocus.get().getInputSetOfTiltSeries()):
            convertInputId = self._insertFunctionStep('convertInputStep',
                                                      ts.getObjId(),
                                                      prerequisites=[])

            allCtfId = []

            for counterCtf in range(0, self.protTomoCtfDefocus.get().numberOfIntermediateStacks[index].get()):
                ctfId = self._insertFunctionStep('computeCtfCorrectionStep',
                                                 ts.getObjId(),
                                                 counterCtf,
                                                 prerequisites=[convertInputId])
                allCtfId.append(ctfId)

            allFlipId = []

            for counterFlip in range(0, self.protTomoCtfDefocus.get().numberOfIntermediateStacks[index].get()):
                flipId = self._insertFunctionStep('computeFlipStep',
                                                  ts.getObjId(),
                                                  counterFlip,
                                                  prerequisites=allCtfId)
                allFlipId.append(flipId)

            allFilterId = []

            for counterFilter in range(0, self.protTomoCtfDefocus.get().numberOfIntermediateStacks[index].get()):
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
        ts = self.protTomoCtfDefocus.get().getInputSetOfTiltSeries()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        """Apply the transformation form the input tilt-series"""
        outputTsFileName = os.path.join(tmpPrefix, "%s.st" % tsId)
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        outputTltFileName = os.path.join(tmpPrefix, '%s.tlt' % tsId)
        ts.generateTltFile(outputTltFileName)

    def computeCtfCorrectionStep(self, tsObjId, counter):
        ts = self.protTomoCtfDefocus.get().getInputSetOfTiltSeries()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(tsId)
        extraPrefixPreviousProt = self.protTomoCtfDefocus.get()._getExtraPath(tsId)
        defocusFilePath = os.path.join(extraPrefixPreviousProt, "%s.defocus_" % tsId)
        tltFilePath = os.path.join(tmpPrefix, "%s.tlt" % tsId)
        outputFilePath = os.path.join(tmpPrefix, "%s.st_" % tsId)

        paramsCtfCorrection = {
            'Algorithm': "ctfCorrection",
            'InputProjections': os.path.join(tmpPrefix, "%s.st" % tsId),
            'OutputFile': outputFilePath + str(counter),
            'DefocusFile': defocusFilePath + str(counter),
            'TiltFile': tltFilePath,
            'CorrectionType': self.getCorrectionType(),
            'DefocusFileFormat': self.getDefocusFileFormat(),
            'CorrectAstigmatism': self.protTomoCtfDefocus.get().correctAstigmatism.get(),
            'PixelSize': self.protTomoCtfDefocus.get().getInputSetOfTiltSeries().getSamplingRate() / 10,
            'AmplitudeContrast':
                self.protTomoCtfDefocus.get().getInputSetOfTiltSeries().getAcquisition().getAmplitudeContrast(),
            'SphericalAberration':
                self.protTomoCtfDefocus.get().getInputSetOfTiltSeries().getAcquisition().getSphericalAberration(),
            'Voltage': self.protTomoCtfDefocus.get().getInputSetOfTiltSeries().getAcquisition().getVoltage()
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
        ts = self.protTomoCtfDefocus.get().getInputSetOfTiltSeries()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(ts.getTsId())
        inputFilePath = os.path.join(tmpPrefix, "%s.st_" % tsId)
        outputFilePath = os.path.join(tmpPrefix, "%s_flip.st_" % tsId)

        argsFlip = "flipyz " + inputFilePath + str(counter) + " " + outputFilePath + str(counter)
        imodPlugin.runImod(self, 'clip', argsFlip)

    def computeFilteringStep(self, tsObjId, counter):
        ts = self.protTomoCtfDefocus.get().getInputSetOfTiltSeries()[tsObjId]
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
        ts = self.protTomoCtfDefocus.get().getInputSetOfTiltSeries()[tsObjId]
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
            'Thickness': self.protTomoCtfDefocus.get().tomoThickness.get(),
            'FullImage': str(ts.getFirstItem().getDim()[0]) + "," + str(ts.getFirstItem().getDim()[1]),
            'Shift': "0.0," + str(self.protTomoCtfDefocus.get().tomoShift.get()),
            'PixelSize': self.protTomoCtfDefocus.get().getInputSetOfTiltSeries().getSamplingRate() / 10,
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

        outputFilePath = os.path.join(extraPrefix, "%s.mrc" % tsId)

        argsTrimvol = outputFilePathFlipped + " "
        argsTrimvol += outputFilePath + " "
        argsTrimvol += "-yz "

        imodPlugin.runImod(self, 'trimvol', argsTrimvol)

    def createOutputStep(self, tsObjId):
        ts = self.protTomoCtfDefocus.get().getInputSetOfTiltSeries()[tsObjId]
        tsId = ts.getTsId()

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
    def getCorrectionType(self):
        if self.protTomoCtfDefocus.get().correctionType.get() == 0:
            correctionType = "phaseflip"
        elif self.protTomoCtfDefocus.get().correctionType.get() == 1:
            correctionType = "multiplication"

        return correctionType

    def getDefocusFileFormat(self):
        if self.protTomoCtfDefocus.get().ctfEstimationType.get() == 0:
            outputDefocusFileFormat = "imod"
        if self.protTomoCtfDefocus.get().ctfEstimationType.get() == 1:
            outputDefocusFileFormat = "ctffind4"

        return outputDefocusFileFormat

    def getOutputSetOfTomograms(self):
        if hasattr(self, "outputSetOfTomograms"):
            self.outputSetOfTomograms.enableAppend()
        else:
            outputSetOfTomograms = self._createSetOfTomograms()
            outputSetOfTomograms.copyInfo(self.protTomoCtfDefocus.get().getInputSetOfTiltSeries())
            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTomograms=outputSetOfTomograms)
            self._defineSourceRelation(self.protTomoCtfDefocus.get().getInputSetOfTiltSeries(), outputSetOfTomograms)

        return self.outputSetOfTomograms

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfTomograms'):
            summary.append("Input Tilt-Series: %d.\nCTF corrected reconstructions calculated: %d.\n"
                           % (self.protTomoCtfDefocus.get().getInputSetOfTiltSeries().getSize(),
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
