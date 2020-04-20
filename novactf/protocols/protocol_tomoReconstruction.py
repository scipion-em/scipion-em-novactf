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
import pyworkflow as pw
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
from pwem.emlib.image import ImageHandler
from tomo.protocols import ProtTomoBase
from tomo.convert import writeTiStack
from tomo.objects import Tomogram, TiltSeries


class ProtTomoCtfReconstruction(EMProtocol, ProtTomoBase):
    """
    Tomogram reconstruction and ctf correction procedure based on the novaCTF procedure.

    More info:
            https://github.com/turonova/novaCTF
    """

    _label = 'tomo ctf reconstruction'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-Series')

        form.addParam('defocusFileFormat',
                      params.EnumParam,
                      choices=['Ctffind', 'IMOD'],
                      default=0,
                      label='Defocus file format',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Defocus file format from CTF estimation.')

        form.addParam('protCtffindCtfEstimation',
                      params.PointerParam,
                      label="Ctffind CTF estimation run",
                      condition='defocusFileFormat==0',
                      pointerClass='ProtTsCtffind',
                      help='Select the previous Ctffind CTF estimation run.')

        form.addParam('protImodCtfEstimation',
                      params.PointerParam,
                      label="IMOD CTF estimation run",
                      condition='defocusFileFormat==1',
                      pointerClass='ProtCtfEstimation',
                      help='Select the previous IMOD CTF estimation run.')

        form.addParam('tomoThickness', params.FloatParam,
                      default=100,
                      label='Tomogram thickness',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Size in pixels of the tomogram in the z axis (beam direction).')

        form.addParam('tomoShift',
                      params.FloatParam,
                      default=0,
                      label='Tomogram shift',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Shift in pixels of the tomogram in the z axis (beam direction).')

        form.addParam('defocusStep',
                      params.IntParam,
                      default=15,
                      label='Defocus step',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Minimum defocus difference used for reconstruction')

        form.addParam('correctionType',
                      params.EnumParam,
                      choices=['Phase filp', 'Multiplication'],
                      default=0,
                      label='Correction type',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Correction type to be applied for reconstruction')

        groupRadialFrequencies = form.addGroup('Radial filtering',
                                               help='This entry controls low-pass filtering with the radial weighting '
                                                    'function.  The radial weighting function is linear away from the '
                                                    'origin out to the distance in reciprocal space specified by the '
                                                    'first value, followed by a Gaussian fall-off determined by the '
                                                    'second value.',
                                               expertLevel=params.LEVEL_ADVANCED)

        groupRadialFrequencies.addParam('radialFirstParameter',
                                        params.FloatParam,
                                        default=0.3,
                                        label='First parameter',
                                        help='Linear region value')

        groupRadialFrequencies.addParam('radialSecondParameter',
                                        params.FloatParam,
                                        default=0.05,
                                        label='Second parameter',
                                        help='Gaussian fall-off parameter')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('computeDefocusStep', ts.getObjId())
            self._insertFunctionStep('computeCtfCorrectionStep', ts.getObjId())
            self._insertFunctionStep('computeFlipStep', ts.getObjId())
            self._insertFunctionStep('computeFilteringStep', ts.getObjId())
            self._insertFunctionStep('computeReconstructionStep', ts.getObjId())
            self._insertFunctionStep('createOutputStep', ts.getObjId())

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
        angleFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        ts.generateTltFile(angleFilePath)

    def computeDefocusStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(ts.getTsId())

        paramsDefocus = {
            'Algorithm': "defocus",
            'InputProjections': os.path.join(tmpPrefix, "%s.st" % tsId),
            'FullImage': str(ts.getFirstItem().getDim()[0]) + "," + str(ts.getFirstItem().getDim()[1]),
            'Thickness': self.tomoThickness.get(),
            'TiltFile': os.path.join(tmpPrefix, "%s.rawtlt" % tsId),
            'Shift': "0.0," + str(self.tomoShift.get()),
            'CorrectionType': self.getCorrectionType(),
            'DefocusFileFormat': self.getDefocusFileFormat(),
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

        self.runJob('/home/fede/novaCTF/novaCTF', argsDefocus % paramsDefocus)

    def computeCtfCorrectionStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(ts.getTsId())
        tltFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        outputFilePath = os.path.join(tmpPrefix, "%s.st_" % tsId)
        defocusFilePath = os.path.join(tmpPrefix, "%s.defocus_" % tsId)

        i = 0
        while os.path.exists(defocusFilePath + str(i)):
            paramsCtfCorrection = {
                'Algorithm': "ctfCorrection",
                'InputProjections': os.path.join(tmpPrefix, "%s.st" % tsId),
                'OutputFile': outputFilePath + str(i),
                'DefocusFile': defocusFilePath + str(i),
                'TiltFile': tltFilePath,
                'CorrectionType': self.getCorrectionType(),
                'DefocusFileFormat': self.getDefocusFileFormat(),
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

            self.runJob('/home/fede/novaCTF/novaCTF', argsCtfCorrection % paramsCtfCorrection)

            i += 1

    def computeFlipStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(ts.getTsId())
        inputFilePath = os.path.join(tmpPrefix, "%s.st_" % tsId)
        outputFilePath = os.path.join(tmpPrefix, "%s_flip.st_" % tsId)
        i = 0
        while os.path.exists(os.path.join(inputFilePath + str(i))):
            argsFlip = inputFilePath + str(i) + " " + outputFilePath + str(i)
            self.runJob('clip flipyz', argsFlip)
            i += 1

    def computeFilteringStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(ts.getTsId())
        outputFilePath = os.path.join(tmpPrefix, "%s_flip_filter.st_" % tsId)
        tltFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)

        i = 0
        while os.path.exists(os.path.join(tmpPrefix, "%s_flip.st_%d" % (tsId, i))):
            paramsFilterProjections = {
                'Algorithm': "filterProjections",
                'InputProjections': os.path.join(tmpPrefix, "%s_flip.st_%d" % (tsId, i)),
                'OutputFile': outputFilePath + str(i),
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
            self.runJob('/home/fede/novaCTF/novaCTF', argsFilterProjections % paramsFilterProjections)

            i += 1

    def computeReconstructionStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(ts.getTsId())
        tmpPrefix = self._getTmpPath(ts.getTsId())
        outputFilePath = os.path.join(extraPrefix, "%s.mrc" % tsId)
        tltFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)

        params3dctf = {
            'Algorithm': "3dctf",
            'InputProjections': os.path.join(tmpPrefix, "%s_flip_filter.st" % tsId),
            'OutputFile': outputFilePath,
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

        self.runJob('/home/fede/novaCTF/novaCTF', args3dctf % params3dctf)

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        """Keep defocus file"""
        tmpDefocusFile = os.path.join(self._getTmpPath(tsId), tsId + ".defocus")
        extraDefocusFile = os.path.join(self._getExtraPath(tsId), tsId + ".defocus")
        path.moveFile(tmpDefocusFile, extraDefocusFile)

        outputSetOfTomograms = self.getOutputSetOfTomograms()
        extraPrefix = self._getExtraPath(tsId)
        newTomogram = Tomogram()
        newTomogram.copyInfo(ts)
        outputSetOfTomograms.append(newTomogram)
        newTomogram.setLocation(os.path.join(extraPrefix, "%s.mrc" % tsId))
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

    def getCorrectionType(self):
        if self.correctionType.get() == 0:
            correctionType = "phaseflip"
        elif self.correctionType.get() == 1:
            correctionType = "multiplication"
        return correctionType

    def getDefocusFileFormat(self):
        if self.defocusFileFormat.get() == 0:
            defocusFileFormat = "ctffind4"
        elif self.defocusFileFormat.get() == 1:
            defocusFileFormat = "imod"
        return defocusFileFormat

    def getDefocusFile(self, ts):
        tsId = ts.getTsId()
        outputDefocusFile = os.path.join(self._getTmpPath(tsId), tsId + ".defocus")

        if self.defocusFileFormat.get() == 0:
            self.generateCtffindDefocusFile(ts, outputDefocusFile)
        elif self.defocusFileFormat.get() == 1:
            outputDefocusFilePrefix = self.protImodCtfEstimation.get()._getExtraPath(tsId)
            path.copyFile(os.path.join(outputDefocusFilePrefix, "%s.defocus" % tsId), outputDefocusFile)

        return outputDefocusFile

    def generateCtffindDefocusFile(self, ts, outputDefocusFile):
        tsId = ts.getTsId()
        defocusInfo = []

        for ti in ts:
            inputDefocusFileName = '%s_%03d_PSD.txt' % (tsId, ti.getObjId())
            inputDefocusFilePrefix = self.protCtffindCtfEstimation.get()._getExtraPath(tsId)
            inputDefocusFilePath = os.path.join(inputDefocusFilePrefix, inputDefocusFileName)
            with open(inputDefocusFilePath, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line[0] != '#':
                        defocusInfo.append(line.split(' '))

        with open(outputDefocusFile, 'w') as f:
            f.writelines("# Columns: #1 - micrograph number; #2 - defocus 1 [Angstroms]; #3 - defocus 2; "
                         "#4 - azimuth of astigmatism; #5 - additional phase shift [radians]; #6 - cross correlation; "
                         "#7 - spacing (in Angstroms) up to which CTF rings were fit successfully\n")

            for counter, vector in enumerate(defocusInfo):
                line = "%.6f %s %s %s %s %s %s" % \
                       (counter+1, vector[1], vector[2], vector[3], vector[4], vector[5], vector[6])
                f.writelines(line)
