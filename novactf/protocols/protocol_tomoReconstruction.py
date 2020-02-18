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
from tomo.objects import Tomogram, TomoAcquisition, TiltSeries


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

        form.addParam('inputSetOfTiltSeries', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-Series')

        form.addParam('tomoThickness', params.FloatParam,
                      default=100,
                      label='Tomogram thickness',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Size in pixels of the tomogram in the z axis (beam direction).')

        form.addParam('tomoShift', params.FloatParam,
                      default=0,
                      label='Tomogram shift',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Shift in pixels of the tomogram in the z axis (beam direction).')

        form.addParam('defocusStep', params.IntParam,
                      default=15,
                      label='Defocus step',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Minimum defocus difference used for reconstruction')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('computeDefocusStep', ts.getObjId())
            self._insertFunctionStep('computeCtfCorrectionStep', ts.getObjId())
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
        extraPrefix = self._getExtraPath(ts.getTsId())
        tltFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        for index, ti in enumerate(ts):
            defocusFilePath = os.path.join(extraPrefix, "%s.defocus" % tsId)
            path.copyFile("/home/fede/Downloads/paraFedeCTF.txt", defocusFilePath)
            paramsDefocus = {
                'Algorithm': "defocus",
                'InputProjections': str(ti.getLocation()[1]),
                'FullImage': str(ti.getDim()[0]) + "," + str(ti.getDim()[1]),
                'Thickness': self.tomoThickness.get(),
                'TiltFile': tltFilePath,
                'Shift': self.tomoShift.get(),
                'CorrectionType': "phaseflip",
                'DefocusFileFormat': "ctffind4",
                'CorrectAstigmatism': 1,
                'DefocusFile': defocusFilePath,
                'PixelSize': self.inputSetOfTiltSeries.get().getSamplingRate(),
                'DefocusStep': self.defocusStep.get()
            }
            argsDefocus = "-Algorithm %(Algorithm)s " \
                          "-InputProjections %(InputProjections)s " \
                          "-FULLIMAGE %(FullImage)s " \
                          "-THICKNESS %(Thickness)d " \
                          "-TILTFILE %(TiltFile)s " \
                          "-SHIFT %(Shift)d " \
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
        extraPrefix = self._getExtraPath(ts.getTsId())
        tmpPrefix = self._getTmpPath(ts.getTsId())
        tltFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        for index, ti in enumerate(ts):
            outputFilePath = os.path.join(extraPrefix, "%s_" % tsId + str(index) + ".ctfCorrection")
            defocusFilePath = os.path.join(extraPrefix, "%s.defocus" % tsId)
            paramsCtfCorrection = {
                'Algorithm': "ctfCorrection",
                'InputProjections': ti.getLocation()[1],
                'OutputFile': outputFilePath,
                'DefocusFile': defocusFilePath,
                'TiltFile': tltFilePath,
                'CorrectionType': "phaseflip",
                'DefocusFileFormat': "ctffind4",
                'CorrectAstigmatism': 1,
                'PixelSize': self.inputSetOfTiltSeries.get().getSamplingRate(),
                'AmplitudeContrast': self.inputSetOfTiltSeries.get().getAcquisition().getAmplitudeContrast(),
                'SphericalAberration': self.inputSetOfTiltSeries.get().getAcquisition().getSphericalAberration(),
                'Voltage': self.inputSetOfTiltSeries.get().getAcquisition().getVoltage()
            }
            print(paramsCtfCorrection)
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

    def computeFilteringStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(ts.getTsId())
        tmpPrefix = self._getTmpPath(ts.getTsId())
        outputFilePath = os.path.join(extraPrefix, "%s_out.filterProj" % tsId)
        tltFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        for ti in ts:
            paramsFilterProjections = {
                'Algorithm': "filterProjections",
                'InputProjections': ti.getLocation()[1],
                'OutputFile': outputFilePath,
                'TiltFile': tltFilePath,
                'StackOrientation': "xy",
                'Radial': "0.3,0.05"
            }
            argsFilterProjections = "-Algorithm %(Algorithm)s " \
                                    "-InputProjections %(InputProjections)s " \
                                    "-OutputFile %(OutputFile)s " \
                                    "-TILTFILE %(TiltFile)s " \
                                    "-StackOrientation %(StackOrientation)s " \
                                    "-RADIAL %(Radial)s"
            self.runJob('/home/fede/novaCTF/novaCTF', argsFilterProjections % paramsFilterProjections)

    def computeReconstructionStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(ts.getTsId())
        tmpPrefix = self._getTmpPath(ts.getTsId())
        outputFilePath = os.path.join(extraPrefix, "%s_out.3dctf" % tsId)
        tltFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        for ti in ts:
            params3dctf = {
                'Algorithm': "3dctf",
                'InputProjections': ti.getLocation()[1],
                'OutputFile': outputFilePath,
                'TiltFile': tltFilePath,
                'Thickness': self.tomoThickness.get(),
                'FullImage': str(ti.getDim()[0]) + "," + str(ti.getDim()[1]),
                'Shift': self.tomoShift.get(),
                'PixelSize': self.inputSetOfTiltSeries.get().getSamplingRate(),
                'DefocusStep': self.defocusStep.get()
            }
            args3dctf = "-Algorithm %(Algorithm)s " \
                        "-InputProjections %(InputProjections)s " \
                        "-OutputFile %(OutputFile)s " \
                        "-TILTFILE %(TiltFile)s " \
                        "-THICKNESS %(Thickness)d " \
                        "-FULLIMAGE %(FullImage)s " \
                        "-SHIFT %(Shift)d " \
                        "-PixelSize %(PixelSize)f " \
                        "-DefocusStep %(DefocusStep)d"
        self.runJob('/home/fede/novaCTF/novaCTF', args3dctf % params3dctf)

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        newTomogram = Tomogram()
        newTomogram.copyInfo(ts)
        newTomogram.setLocation(os.path.join(self._getExtraPath(tsId), '%s.mrc' % tsId))
        outputSetOfTomograms = self.getOutputSetOfTomograms()
        outputSetOfTomograms.append(newTomogram)

        """Debug code ***"""
        path.moveTree(self._getTmpPath(), self._getExtraPath())

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfTomograms(self):
        if not hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            outputSetOfTomograms = self._createSetOfTomograms()
            outputSetOfTomograms.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTomograms.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputSetOfTomograms=outputSetOfTomograms)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTomograms)
        return self.outputSetOfTomograms
