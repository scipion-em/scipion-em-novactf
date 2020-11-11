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
import pyworkflow.protocol.params as params
from pyworkflow.project import Manager
import pyworkflow.utils.path as path
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.object import Integer, List
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from novactf import Plugin


class ProtNovaCtfTomoDefocus(EMProtocol, ProtTomoBase):
    """
    Defocus estimation of each tilt-image procedure based on the novaCTF procedure.

    More info:
            https://github.com/turonova/novaCTF
    """

    _label = 'tomo ctf defocus'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        ProtTomoBase.__init__(self)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.numberOfIntermediateStacks = List([])

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('ctfEstimationType',
                      params.EnumParam,
                      choices=['IMOD', 'Other'],
                      default=1,
                      label='CTF estimation',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='CTF estimation origin.')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      condition='ctfEstimationType==1',
                      label='Input set of tilt-Series')

        form.addParam('protImodCtfEstimation',
                      params.PointerParam,
                      label="IMOD CTF estimation run",
                      condition='ctfEstimationType==0',
                      pointerClass='ProtImodCtfEstimation',
                      help='Select the previous IMOD CTF estimation run.')

        form.addParam('tomoThickness',
                      params.FloatParam,
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

        form.addParam('correctAstigmatism',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=0,
                      label='Correct astigmatism',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Correct for astigmatism in reconstruction')

        self.defineFilterParameters(form)

        form.addParallelSection(threads=8, mpi=1)

    @staticmethod
    def defineFilterParameters(form):
        groupRadialFrequencies = form.addGroup('Radial filtering',
                                               help='This entry controls low-pass filtering with the radial weighting '
                                                    'function.  The radial weighting function is linear away from the '
                                                    'origin out to the distance in reciprocal space specified by the '
                                                    'first value, followed by a Gaussian fall-off determined by the '
                                                    'second value.')
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

        for ts in self.getInputSetOfTiltSeries():
            self._insertFunctionStep('convertInputStep',
                                     ts.getObjId())

            if self.ctfEstimationType.get() == 0:
                self._insertFunctionStep('generateImodDefocusFileStep',
                                         ts.getObjId())

            elif self.ctfEstimationType.get() == 1:
                self._insertFunctionStep('generateCtffindDefocusFileStep',
                                         ts.getObjId())

            self._insertFunctionStep('computeDefocusStep',
                                     ts.getObjId())

            self._insertFunctionStep('getNumberOfIntermediateStacksStep',
                                     ts.getObjId())

        self._insertFunctionStep('triggerReconstructionProtocolStep')

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.getInputSetOfTiltSeries()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(tsId)
        extraPrefix = self._getExtraPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, "%s.tlt" % tsId)
        ts.generateTltFile(angleFilePath)

    def generateImodDefocusFileStep(self, tsObjId):
        ts = self.getInputSetOfTiltSeries()[tsObjId]
        tsId = ts.getTsId()
        outputDefocusFile = self.getDefocusFile(ts)

        outputDefocusFilePrefix = self.protImodCtfEstimation.get()._getExtraPath(tsId)
        path.copyFile(os.path.join(outputDefocusFilePrefix, "%s.defocus" % tsId), outputDefocusFile)

    def generateCtffindDefocusFileStep(self, tsObjId):
        ts = self.getInputSetOfTiltSeries()[tsObjId]
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
        ts = self.getInputSetOfTiltSeries()[tsObjId]
        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(ts.getTsId())

        paramsDefocus = {
            'Algorithm': "defocus",
            'InputProjections': ts.getFirstItem().getLocation()[1],
            'FullImage': str(ts.getFirstItem().getDim()[0]) + "," + str(ts.getFirstItem().getDim()[1]),
            'Thickness': self.tomoThickness.get(),
            'TiltFile': os.path.join(tmpPrefix, "%s.tlt" % tsId),
            'Shift': "0.0," + str(self.tomoShift.get()),
            'CorrectionType': self.getCorrectionType(),
            'DefocusFileFormat': self.getDefocusFileFormat(),
            'CorrectAstigmatism': self.correctAstigmatism.get(),
            'DefocusFile': self.getDefocusFile(ts),
            'PixelSize': self.getInputSetOfTiltSeries().getSamplingRate() / 10,
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

    def getNumberOfIntermediateStacksStep(self, tsObjId):
        ts = self.getInputSetOfTiltSeries()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(ts.getTsId())

        defocusFilePath = os.path.join(extraPrefix, "%s.defocus_" % tsId)
        numberOfIntermediateStacks = 0

        counter = 0
        while os.path.exists(defocusFilePath + str(counter)):
            numberOfIntermediateStacks += 1
            counter += 1

        self.numberOfIntermediateStacks.append(Integer(numberOfIntermediateStacks))

    def triggerReconstructionProtocolStep(self):
        # Local import to avoid looping
        from novactf.protocols import ProtNovaCtfTomoReconstruction

        manager = Manager()
        project = manager.loadProject(self.getProject().getName())

        protTomoReconstruction = ProtNovaCtfTomoReconstruction()
        protTomoReconstruction.setObjLabel('novactf - tomo ctf reconstruction')
        protTomoReconstruction.protTomoCtfDefocus.set(self)
        protTomoReconstruction.radialFirstParameter.set(self.radialFirstParameter.get())
        protTomoReconstruction.radialSecondParameter.set(self.radialSecondParameter.get())
        protTomoReconstruction.numberOfThreads.set(self.numberOfThreads.get())
        protTomoReconstruction.numberOfMpi.set(self.numberOfMpi.get())

        project.scheduleProtocol(protTomoReconstruction)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getInputSetOfTiltSeries(self):
        if self.ctfEstimationType.get() == 0:
            inputSetOfTiltSeries = self.protImodCtfEstimation.get().inputSetOfTiltSeries.get()
        elif self.ctfEstimationType.get() == 1:
            inputSetOfTiltSeries = self.inputSetOfTiltSeries.get()

        return inputSetOfTiltSeries

    def getCorrectionType(self):
        if self.correctionType.get() == 0:
            correctionType = "phaseflip"
        elif self.correctionType.get() == 1:
            correctionType = "multiplication"

        return correctionType

    def getDefocusFile(self, ts):
        tsId = ts.getTsId()
        outputDefocusFile = os.path.join(self._getExtraPath(tsId), tsId + ".defocus")

        return outputDefocusFile

    def getDefocusFileFormat(self):
        if self.ctfEstimationType.get() == 0:
            outputDefocusFileFormat = "imod"
        if self.ctfEstimationType.get() == 1:
            outputDefocusFileFormat = "ctffind4"

        return outputDefocusFileFormat

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

        counter = 0

        for ts in self.getInputSetOfTiltSeries():
            tsId = ts.getTsId()
            if os.path.exists(os.path.join(self._getExtraPath(tsId), tsId + ".defocus_0")):
                counter += 1

        if counter != 0:
            summary.append("Input Tilt-Series: %d.\n"
                           "Tilt-series defocus processed: %d.\n"
                           "Defocus files generated for each tilt-series: %d.\n"
                           % (self.getInputSetOfTiltSeries().getSize(),
                              counter,
                              self.numberOfIntermediateStacks[0]))
        else:
            summary.append("Output not ready yet.")
        return summary

    def _methods(self):
        methods = []

        counter = 0

        for ts in self.getInputSetOfTiltSeries():
            tsId = ts.getTsId()
            if os.path.exists(os.path.join(self._getExtraPath(tsId), tsId + ".defocus_0")):
                counter += 1

        if counter != 0:
            methods.append("%d defocus files have been generated for each of the %d tilt-series using the defocus "
                           "algorithm from novaCTF.\n"
                           % (self.numberOfIntermediateStacks[0],
                              counter))
        else:
            methods.append("Output not ready yet.")

        return methods
