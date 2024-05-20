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

from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.object import Integer, List
from pwem.protocols import EMProtocol

from imod import utils as imodUtils
from novactf import Plugin


class ProtNovaCtfDefocus(EMProtocol):
    """
    Compute defocus for each tilt-image with novaCTF. This is a metaprotocol
    and automatically will trigger novaCTf - 3D CTF correction and reconstruction

    More info:
            https://github.com/turonova/novaCTF

    NovaCTF is a tomogram reconstruction algorithm that allows a local CTF correction using
    the weighted back projection (WBP) algorithm. This local CTF correction enhances the
    quality of the tomogram leading to a higher resolution in subtomogram averaging.\n

    NovaCTF is an efficient implementation of G.J. Jensen, R.D. Kornberg, "Defocus-gradient
    corrected backprojection", Ultramicroscopy, 84, 57-64, (2000). This algorithm uses multiple CTF
    corrections to reconstruct the tomogram via WBP to have a local CTF corrected tomogram. To achieve
    this, the algorithm takes into account the gradient of defocus in the sample. This means that
    the top and the botton of the sample present different defocus values. A set of planes or
    heights are defined along the gradient of defocus to model the defocus gradient. The number of
    planes is determined by the parameter defocus step. By tilting the sample the defocus of a given
    point in the sample will change with the tilt angle. This is due to the change of position
    (height) of such point. Therefore, the defocus of the same voxel will be different in different
    tilt images. The algorithm of novaCTF carries out a multiple CTF correction per tilt image
    (as many as defocus steps will be defined). Then, A WBP will be carried out to reconstruct
    the tomogram, however, according to the position of the voxel to the reconstruction the
    corresponding CTF corrected image will be back projected. This ensures a local CTF correction
    in the reconstruction.
    """

    _label = 'compute defocus'
    _devStatus = PROD

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.numberOfIntermediateStacks = List()

    def _initialize(self):
        self._createFilenameTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
            'defocusFn': self._getExtraPath("%(tsId)s/%(tsId)s.defocus"),
            'stackDefocusFn': self._getExtraPath("%(tsId)s/%(tsId)s.defocus_%(counter)d"),
            'defocusShiftFn': self._getExtraPath("%(tsId)s/%(tsId)s.def_shift"),
            'tltFn': self._getTmpPath("%(tsId)s/%(tsId)s.tlt")
        }

        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series')

        form.addParam('inputSetOfCtfTomoSeries', params.PointerParam,
                      label="Associated CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      important=True,
                      help='CTF of the tilt-series.')

        form.addParam('tomoThickness', params.IntParam,
                      default=400,
                      label='Tomogram thickness (voxels)',
                      help='Size of the tomogram in voxels in the Z axis (beam direction).')

        form.addParam('tomoShift', params.IntParam,
                      default=0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Tomogram shift (voxels)',
                      help='Shift of the tomogram in voxels in the Z direction. '
                           'The shift should be set to zero even if for '
                           'reconstruction we want to shift the tomogram in z! '
                           'We assume the defocus to be estimated at the center '
                           'of mass which should correspond to the shifted '
                           'tomogram and thus here the shift should be zero.')

        form.addParam('defocusStep', params.IntParam,
                      default=15,
                      label='Defocus step (nm)',
                      help='The space between min and max in Z is sliced '
                           'by defocus step. 15 nm is default step size. '
                           'See Fig. 2 in Turonova et al., 2017 for optimal number.')

        form.addParam('correctionType', params.EnumParam,
                      choices=['Phase flip', 'Multiplication'],
                      default=0,
                      label='Correction type',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='CTF correction type to be applied for the tilt-series. \n'
                           ''
                           '_Phase Flipping Method_: The phase information of the Fourier transform '
                           'of the image is flipped by 180 degrees for those frequencies affected by '
                           'the CTF. Then, the inverse Fourier transform is applied to obtaing the'
                           'CTF corrected image.\n'
                           '_Multiplication Method_: The Fourier transform of the image is multiplied'
                           'by the CTF to attenuate the amplitudes of the frequencies affected by the '
                           'CTF. Finally, the inverse Fourier transform is applied to reconstruct the '
                           'corrected image.\n'
                           'Difference:\n'
                           '* Phase flipping directly modifies the phase information of the Fourier transform, '
                           'while multiplication modifies the amplitudes.\n'
                           '* Multiplication involves modeling and multiplication with a sinusoidal function, '
                           'which can be more complex computationally.')

        form.addParam('correctAstigmatism', params.BooleanParam,
                      default=True,
                      label='Correct astigmatism?',
                      help='Correct for astigmatism in reconstruction.')

        form.addParallelSection(threads=4)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()

        for ts in self.getInputTs().iterItems():
            tsId = ts.getTsId()
            objId = ts.getObjId()
            self._insertFunctionStep(self.convertInputStep, objId, tsId)
            self._insertFunctionStep(self.computeDefocusStep, objId, tsId)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, tsObjId, tsId):
        self.info(f"Generating tlt and defocus files for {tsId}")

        # Create the folders for the tilt series
        path.makePath(self._getTmpPath(tsId))
        path.makePath(self._getExtraPath(tsId))

        # Generate angle file
        with self._lock:
            ts = self.getInputTs()[tsObjId]
            ts.generateTltFile(self._getFileName("tltFn", tsId=tsId))

        # Generate defocus file
        ctfTomoSeries = self.getCtfTomoSeriesFromTsId(tsId)
        defocusFile = self._getFileName("defocusFn", tsId=tsId)
        imodUtils.generateDefocusIMODFileFromObject(ctfTomoSeries, defocusFile)

    def computeDefocusStep(self, tsObjId, tsId):
        with self._lock:
            ts = self.getInputTs()[tsObjId]
            firstItem = ts.getFirstItem()
            tsFn = firstItem.getFileName()
            xDim, yDim, _ = firstItem.getDim()

        paramsDefocus = {
            '-Algorithm': "defocus",
            '-InputProjections': tsFn,
            '-FULLIMAGE': f"{xDim},{yDim}",
            '-THICKNESS': self.tomoThickness.get(),
            '-TILTFILE': self._getFileName("tltFn", tsId=tsId),
            '-SHIFT': "0.0,0.0",
            '-CorrectionType': self.getCorrectionType(),
            '-DefocusFileFormat': "imod",
            '-CorrectAstigmatism': 1 if self.correctAstigmatism else 0,
            '-DefocusFile': self._getFileName("defocusFn", tsId=tsId),
            '-PixelSize': self.getInputSamplingRate() / 10,
            '-DefocusStep': self.defocusStep,
        }

        if self.tomoShift.get() > 0.01:
            defocusShiftFile = self._getFileName("defocusShiftFn", tsId=tsId)
            paramsDefocus["-DefocusShiftFile"] = defocusShiftFile
            defocusShift = self.tomoThickness.get() / 2 + self.tomoShift.get()
            with open(defocusShiftFile, "w") as fn:
                fn.write(f"{defocusShift}")

        Plugin.runNovactf(self, **paramsDefocus)

        nstacks = self.getNumberOfStacks(tsId)
        self.numberOfIntermediateStacks.append(Integer(nstacks))

        self._store()

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        if self.getInputTs().getSize() != self.getInputCtf().getSize():
            validateMsgs.append("Input set of tilt-series and input set of CTFs "
                                " must contain the same number of items.")

        ctf = self.getInputCtf().getFirstItem()
        if hasattr(ctf, "_IMODDefocusFileFlag"):
            defFlag = ctf.getIMODDefocusFileFlag()
            if defFlag in [0, 4] and self.correctAstigmatism:
                validateMsgs.append("CTF estimation does not have astigmatism values.")

        return validateMsgs

    def _summary(self):
        summary = []

        if len(self.numberOfIntermediateStacks):
            summary.append(f"Input tilt-series: {self.getInputTs().getSize()}\n"
                           "Defocus files generated for each tilt-series.\n"
                           "Use this protocol as input for novaCTF - 3D CTF "
                           "correction and reconstruction.")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    # --------------------------- UTILS functions -----------------------------
    def getInputTs(self):
        return self.inputSetOfTiltSeries.get()

    def getInputCtf(self):
        return self.inputSetOfCtfTomoSeries.get()

    def getCtfTomoSeriesFromTsId(self, tsId):
        for item in self.getInputCtf().iterItems(where=f"_tsId=='{tsId}'"):
            return item

    def getNumberOfStacks(self, tsId):
        defocusFilePath = self._getFileName("defocusFn", tsId=tsId)
        files = glob(defocusFilePath.replace(".defocus", ".defocus_*"))
        return len(files)

    def getInputSamplingRate(self):
        return self.getInputTs().getSamplingRate()

    def getCorrectionType(self):
        return "phaseflip" if self.correctionType.get() == 0 else "multiplication"
