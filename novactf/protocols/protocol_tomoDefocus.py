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
from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pyworkflow.project import Manager
import pyworkflow.utils.path as path
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.object import Integer, List
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from novactf import Plugin
from imod import utils as imodUtils


class ProtNovaCtfTomoDefocus(EMProtocol, ProtTomoBase):
    """
    Defocus estimation of each tilt-image procedure based on the novaCTF procedure.

    More info:
            https://github.com/turonova/novaCTF
    """

    _label = 'tomo ctf defocus'
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        ProtTomoBase.__init__(self)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.numberOfIntermediateStacks = List([])

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input set of tilt-Series')

        form.addParam('inputSetOfCtfTomoSeries',
                      params.PointerParam,
                      label="input tilt-series CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation from the set of tilt-series.')

        form.addParam('tomoThickness',
                      params.FloatParam,
                      default=100,
                      label='Tomogram thickness (voxels)',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Size in voxels of the tomogram in the z axis (beam direction).')

        form.addParam('tomoShift',
                      params.FloatParam,
                      default=0,
                      label='Tomogram shift (voxels)',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Shift in voxels of the tomogram in the z axis (beam direction).')

        form.addParam('defocusStep',
                      params.IntParam,
                      default=15,
                      label='Defocus step (nm)',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Minimum defocus difference used for reconstruction in nanometers.')

        form.addParam('correctionType',
                      params.EnumParam,
                      choices=['Phase flip', 'Multiplication'],
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

        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.convertInputStep,
                                     ts.getObjId())

            self._insertFunctionStep(self.computeDefocusStep,
                                     ts.getObjId())

            self._insertFunctionStep(self.getNumberOfIntermediateStacksStep,
                                     ts.getObjId())

        self._insertFunctionStep(self.triggerReconstructionProtocolStep)

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        ctfTomoSeries = self.inputSetOfCtfTomoSeries.get()[tsObjId]

        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(tsId)
        extraPrefix = self._getExtraPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        firstItem = ts.getFirstItem()

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))
        ts.generateTltFile(angleFilePath)

        """Generate defocus file"""
        defocusFilePath = os.path.join(extraPrefix, firstItem.parseFileName(extension=".defocus"))
        imodUtils.generateDefocusIMODFileFromObject(ctfTomoSeries, defocusFilePath)

    def computeDefocusStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        tmpPrefix = self._getTmpPath(ts.getTsId())
        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem()

        defocusFilePath = os.path.join(extraPrefix, firstItem.parseFileName(extension=".defocus"))

        paramsDefocus = {
            'Algorithm': "defocus",
            'InputProjections': firstItem.getLocation()[1],
            'FullImage': str(firstItem.getDim()[0]) + "," + str(firstItem.getDim()[1]),
            'Thickness': self.tomoThickness.get(),
            'TiltFile': os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt")),
            'Shift': "0.0," + str(self.tomoShift.get()),
            'CorrectionType': self.getCorrectionType(),
            'DefocusFileFormat': "imod",
            'CorrectAstigmatism': self.correctAstigmatism.get(),
            'DefocusFile': defocusFilePath,
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

    def getNumberOfIntermediateStacksStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]

        extraPrefix = self._getExtraPath(ts.getTsId())

        defocusFilePath = os.path.join(extraPrefix,
                                       ts.getFirstItem().parseFileName(extension=".defocus_"))
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
    def getCorrectionType(self):
        if self.correctionType.get() == 0:
            correctionType = "phaseflip"
        elif self.correctionType.get() == 1:
            correctionType = "multiplication"

        return correctionType

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        if self.inputSetOfTiltSeries.get().getSize() != self.inputSetOfCtfTomoSeries.get().getSize():
            validateMsgs.append("Input set of tilt-series and input set of CTF tomo estimations must contain the "
                                "same number of elements.")

        return validateMsgs

    def _summary(self):
        summary = []

        counter = 0

        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            if os.path.exists(os.path.join(self._getExtraPath(tsId), tsId + ".defocus_0")):
                counter += 1

        if counter != 0:
            summary.append("Input Tilt-Series: %d.\n"
                           "Tilt-series defocus processed: %d.\n"
                           "Defocus files generated for each tilt-series: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              counter,
                              self.numberOfIntermediateStacks[0]))
        else:
            summary.append("Output not ready yet.")
        return summary

    def _methods(self):
        methods = []

        counter = 0

        for ts in self.inputSetOfTiltSeries.get():
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
