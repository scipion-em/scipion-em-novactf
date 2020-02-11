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
from tomo.objects import Tomogram, TomoAcquisition


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

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts)
            self._insertFunctionStep('computeDefocusStep', ts)
            self._insertFunctionStep('computeCtfCorrectionStep', ts)
            self._insertFunctionStep('computeFilteringStep', ts)
            self._insertFunctionStep('computeReconstructionStep', ts)
            self._insertFunctionStep('createOutputStep', ts)

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, ts):
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

    def computeDefocusStep(self, ts):
        tmpPrefix = self._getTmpPath(ts.gettsId())
        extraPrefix = self._getExtraPath(ts.gettsId())
        tltFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        defocusFilePath = os.path.join(extraPrefix, "%s.defocus" % tsId)
        for ti in ts:
            paramsDefocus = {
                'Algorithm': "defocus",
                'InputProjections': ti.getLocation(),
                'FullImage': ti.getDim(),
                'Thickness': self.tomoThickness.get(),
                'TiltFile': tltFilePath,
                'Shift': self.tomoShift.get(),
                'CorrectionType': "phaseflip",
                'DefocusFileFormat': "ctffind4",
                'CorrectAstigmatism': 1,
                'DefocusFile': defocusFilePath,
                'PixelSize': self.inputSetOfTiltSeries.get().getSamplingRate(),
                'DefocusStep': 15
            }
            argsDefocus = "-Algorithm %(Algorithm)s," \
                          "-InputProjections %(InputProjections)s," \
                          "-FULLIMAGE %(FullImage)s," \
                          "-THICKNESS %(Thickness)d," \
                          "-TILTFILE %(TiltFile)s," \
                          "-SHIFT %(Shift)d," \
                          "-CorrectionType %(CorrectionType)s," \
                          "-DefocusFileFormat %(DefocusFileFormat)," \
                          "-CorrectAstigmatism %(CorrectAstigmatism)s," \
                          "-DefocusFile %(DefocusFile)s," \
                          "-PixelSize %(PixelSize)s," \
                          "-DefocusStep %(DefocusStep)d"
            self.runJob('/usr/local/novaCTF/novaCTF', argsDefocus % paramsDefocus)

    def computeCtfCorrectionStep(self, ts):
        extraPrefix = self._getExtraPath(ts.gettsId())
        tmpPrefix = self._getTmpPath(ts.gettsId())
        outputFilePath = os.path.join(extraPrefix, "%s.ctfCorrection" % tsId)
        defocusFilePath = os.path.join(extraPrefix, "%s.defocus" % tsId)
        tltFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        for ti in ts:
            paramsCtfCorrection = {
                'Algorithm': "ctfCorrection",
                'InputProjections': ti.getLocation(),
                'OutputFile': outputFilePath,
                'DefocusFile': defocusFilePath,
                'TiltFile': tltFilePath,
                'CorrectionType': "phaseflip",
                'DefocusFileFormat': "ctffind4",
                'CorrectAstigmatism': 1,
                'PixelSize': self.inputSetOfTiltSeries.get().getSamplingRate(),
                'AmplitudeContrast': self.inputSetOfTiltSeries.get().getAmplitudeContrast(),
                'SphericalAberration': self.inputSetOfTiltSeries.get().getSphericalAberration(),
                'Voltage': self.inputSetOfTiltSeries.get().getVoltage()
            }
            argsCtfCorrection = "-Algorithm %(Algorithm)s," \
                                "-InputProjections %(InputProjections)s," \
                                "-OutputFile %(OutputFile)s," \
                                "-DefocusFile %(DefocusFile)s," \
                                "-TILTFILE %(TiltFile)s," \
                                "-CorrectionType %(CorrectionType)s," \
                                "-DefocusFileFormat %(DefocusFileFormat)," \
                                "-CorrectAstigmatism %(CorrectAstigmatism)s," \
                                "-PixelSize %(PixelSize)f," \
                                "-AmplitudeContrast %(AmplitudeContrast)f" \
                                "-Cs %(SphericalAberration)f" \
                                "-Volt %(Voltage)d"
            self.runJob('/usr/local/novaCTF/novaCTF', argsCtfCorrection % paramsCtfCorrection)

    def computeFilteringStep(self, ts):
        extraPrefix = self._getExtraPath(ts.gettsId())
        tmpPrefix = self._getTmpPath(ts.gettsId())
        outputFilePath = os.path.join(extraPrefix, "%s_out.filterProj" % tsId)
        tltFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        for ti in ts:
            paramsFilterProjections = {
                'Algorithm': "filterProjections",
                'InputProjections': ti.getLocation(),
                'OutputFile': outputFilePath,
                'TiltFile': tltFilePath,
                'StackOrientation': "xy",
                'Radial': "0.3,0.05"
            }
            argsFilterProjections = "-Algorithm %(Algorithm)s," \
                                    "-InputProjections %(InputProjections)s," \
                                    "-OutputFile %(OutputFile)s," \
                                    "-TILTFILE %(TiltFile)s," \
                                    "-StackOrientation %(StackOrientation)s," \
                                    "-RADIAL %(Radial)s"
            self.runJob('/usr/local/novaCTF/novaCTF', argsFilterProjections % paramsFilterProjections)

    def computeReconstructionStep(self, ts):
        for ts in self.inputSetOfTiltSeries.get():
            extraPrefix = self._getExtraPath(ts.gettsId())
            tmpPrefix = self._getTmpPath(ts.gettsId())
            outputFilePath = os.path.join(extraPrefix, "%s_out.3dctf" % tsId)
            tltFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
            for ti in ts:
                params3dctf = {
                    'Algorithm': "3dctf",
                    'InputProjections': ti.getLocation(),
                    'OutputFile': outputFilePath,
                    'TiltFile': tltFilePath,
                    'Thickness': self.tomoThickness.get(),
                    'FullImage': ti.getDim(),
                    'Shift': self.tomoShift.get(),
                    'PixelSize': self.inputSetOfTiltSeries.get().getSamplingRate(),
                    'DefocusStep': 15
                }
                args3dctf = "-Algorithm %(Algorithm)s," \
                                        "-InputProjections %(InputProjections)s," \
                                        "-OutputFile %(OutputFile)s," \
                                        "-TILTFILE %(TiltFile)s," \
                                        "-THICKNESS %(Thickness)d," \
                                        "-FULLIMAGE %(FullImage)s," \
                                        "-SHIFT %(Shift)d," \
                                        "-PixelSize %(PixelSize)f," \
                                        "-DefocusStep %(DefocusStep)d"
            self.runJob('/usr/local/novaCTF/novaCTF', args3dctf % params3dctf)

    def createOutputStep(self):
        self.outputSetOfTomograms = self._createSetOfTomograms()
        self.outputSetOfTomograms.setSamplingRate(self.inputSetOfTiltSeries.get().getSamplingRate())
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            newTomogram = Tomogram()
            newTomogram.setLocation(os.path.join(self._getExtraPath(tsId), '%s.mrc' % tsId))
            self.outputSetOfTomograms.append(newTomogram)
        self._defineOutputs(outputTomograms=self.outputSetOfTomograms)
        self._defineSourceRelation(self.inputSetOfTiltSeries, self.outputSetOfTomograms)

        """Debug code ***"""
        path.moveTree(self._getTmpPath(), self._getExtraPath())
