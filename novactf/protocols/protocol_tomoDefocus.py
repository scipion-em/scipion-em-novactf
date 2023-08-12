# *****************************************************************************
# *
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

from glob import glob

from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pyworkflow.project import Manager
import pyworkflow.utils.path as path
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.object import Integer, List
from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase

from imod import utils as imodUtils
from .. import Plugin


class ProtNovaCtfTomoDefocus(EMProtocol, ProtTomoBase):
    """
    Compute defocus array for each tilt-image with novaCTF.

    More info:
            https://github.com/turonova/novaCTF
    """

    _label = 'compute defocus array'
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        ProtTomoBase.__init__(self)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input set of tilt-series')

        form.addParam('inputSetOfCtfTomoSeries',
                      params.PointerParam,
                      label="Input tilt-series CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation for the input tilt-series.')

        form.addParam('tomoThickness',
                      params.FloatParam,
                      default=400,
                      label='Tomogram thickness (voxels)',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Size of the tomogram in voxels in the Z direction.')

        form.addParam('tomoShift',
                      params.FloatParam,
                      default=0,
                      label='Tomogram shift (voxels)',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Shift of the tomogram in voxels in the Z direction. '
                           'The shift should be set to zero even if for '
                           'reconstruction we want to shift the tomogram in z! '
                           'We assume the defocus to be estimated at the center '
                           'of mass which should correspond to the shifted '
                           'tomogram and thus here the shift should be zero.')

        form.addParam('defocusStep',
                      params.IntParam,
                      default=15,
                      label='Defocus step (nm)',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='The space between min and max in Z is sliced '
                           'by defocus step. 15 nm is default step size. '
                           'See Fig. 2 in Turonova et al., 2017 for optimal number.')

        form.addParam('correctionType',
                      params.EnumParam,
                      choices=['Phase flip', 'Multiplication'],
                      default=0,
                      label='Correction type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='CTF correction type to be applied for the tilt-series.')

        form.addParam('correctAstigmatism',
                      params.EnumParam,
                      choices=['No', 'Yes'],
                      default=1,
                      label='Correct astigmatism',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Correct for astigmatism in reconstruction.')

        form.addParallelSection(threads=2, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self.numberOfIntermediateStacks = List([])
        for ts in self.getInputTs():
            self._insertFunctionStep(self.convertInputStep, ts.getObjId())
            self._insertFunctionStep(self.computeDefocusStep, ts.getObjId())
        self._insertFunctionStep(self.triggerReconstructionProtocolStep)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.getInputTs()[tsObjId]
        tsId = ts.getTsId()
        self.info("Generating tlt file and defocus file for %s (id %s)" % (tsId, tsObjId))

        ctfTomoSeries = self.getCtfTomoSeriesFromTsId(tsId)

        tmpPrefix = self._getTmpPath(tsId)
        extraPrefix = self._getExtraPath(tsId)

        # Create the folders for the tilt series
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        # Generate angle file
        angleFilePath = self.getTltFileName(tsId)

        self.info("Generating %s" % angleFilePath)
        ts.generateTltFile(angleFilePath)

        # Generate defocus file
        defocusFilePath = self.getDefocusFileName(tsId)

        imodUtils.generateDefocusIMODFileFromObject(ctfTomoSeries, defocusFilePath)

    def getTltFileName(self, tsId):
        return self._getTmpPath(tsId, tsId + ".tlt")

    def getDefocusFileName(self,tsId):
        return self._getExtraPath(tsId, tsId + ".defocus")

    def computeDefocusStep(self, tsObjId):
        ts = self.getInputTs()[tsObjId]
        tsId = ts.getTsId()

        firstItem = ts.getFirstItem()

        ih = ImageHandler()
        xDim, yDim, _, _ = ih.getDimensions(firstItem.getFileName() + ":mrc")

        defocusFilePath = self.getDefocusFileName(tsId)

        if self.tomoShift.get() > 0.01:
            defocusShift = self.tomoThickness.get() / 2 + self.tomoShift.get()
            with open(defocusFilePath.replace(".defocus", ".def_shift"), "w") as fn:
                fn.write(f"{defocusShift}")

        paramsDefocus = {
            'Algorithm': "defocus",
            'InputProjections': firstItem.getFileName(),
            'FullImage': str(xDim) + "," + str(yDim),
            'Thickness': self.tomoThickness.get(),
            'TiltFile': self.getTltFileName(tsId),
            'CorrectionType': "phaseflip" if self.correctionType.get() == 0 else "multiplication",
            'DefocusFileFormat': "imod",
            'CorrectAstigmatism': self.correctAstigmatism.get(),
            'DefocusFile': defocusFilePath,
            'PixelSize': self.getInputTs().getSamplingRate() / 10,
            'DefocusStep': self.defocusStep.get(),
            'DefocusShiftFile': defocusFilePath.replace(".defocus", ".def_shift")
        }

        argsDefocus = "-Algorithm %(Algorithm)s " \
                      "-InputProjections %(InputProjections)s " \
                      "-FULLIMAGE %(FullImage)s " \
                      "-THICKNESS %(Thickness)d " \
                      "-TILTFILE %(TiltFile)s " \
                      "-SHIFT 0.0,0.0 " \
                      "-CorrectionType %(CorrectionType)s " \
                      "-DefocusFileFormat %(DefocusFileFormat)s " \
                      "-CorrectAstigmatism %(CorrectAstigmatism)d " \
                      "-DefocusFile %(DefocusFile)s " \
                      "-PixelSize %(PixelSize)s " \
                      "-DefocusStep %(DefocusStep)d "

        if self.tomoShift.get() > 0.01:
            argsDefocus += "-DefocusShiftFile %(DefocusShiftFile)s "

        Plugin.runNovactf(self, 'novaCTF', argsDefocus % paramsDefocus)

        files = glob(defocusFilePath.replace(".defocus", ".defocus_*"))
        self.numberOfIntermediateStacks.append(Integer(len(files)))

        self._store()

    def triggerReconstructionProtocolStep(self):
        # Local import to avoid looping
        from . import ProtNovaCtfTomoReconstruction

        manager = Manager()
        project = manager.loadProject(self.getProject().getName())

        applyAlignment = self.getInputTs().getFirstItem().getFirstItem().hasTransform()

        label = self.getObjLabel() + " - REC"
        protTomoReconstruction = ProtNovaCtfTomoReconstruction(applyAlignment=applyAlignment,
                                                               objLabel=label)
        protTomoReconstruction.protTomoCtfDefocus.set(self)
        protTomoReconstruction.numberOfThreads.set(self.numberOfThreads.get())
        protTomoReconstruction.numberOfMpi.set(self.numberOfMpi.get())

        project.scheduleProtocol(protTomoReconstruction)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        if self.getInputTs().getSize() != self.inputSetOfCtfTomoSeries.get().getSize():
            validateMsgs.append("Input set of tilt-series and input set of CTFs "
                                " must contain the same number of items.")

        ctf = self.inputSetOfCtfTomoSeries.get().getFirstItem()
        if hasattr(ctf, "_IMODDefocusFileFlag"):
            defFlag = ctf.getIMODDefocusFileFlag()
            if defFlag in [0, 4] and self.correctAstigmatism.get() == 1:
                validateMsgs.append("CTF estimation does not have astigmatism values.")

        return validateMsgs

    def _summary(self):
        summary = []

        if len(self.numberOfIntermediateStacks):
            summary.append(f"Input tilt-series: {self.getInputTs().getSize()}\n"
                           f"Defocus files generated for each tilt-series: "
                           f"{self.numberOfIntermediateStacks[0]}\n\n"
                           "Use this protocol as input for novaCTF - 3D CTF "
                           "correction and reconstruction.")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    # --------------------------- UTILS functions -----------------------------
    def getInputTs(self):
        return self.inputSetOfTiltSeries.get()

    def getCtfTomoSeriesFromTsId(self, tsId):
        for ctfTomoSeries in self.inputSetOfCtfTomoSeries.get():
            if tsId == ctfTomoSeries.getTsId():
                return ctfTomoSeries
