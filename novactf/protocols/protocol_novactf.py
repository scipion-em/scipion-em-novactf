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
import logging
import traceback
from enum import Enum
from glob import glob
from os.path import join, exists
from imod import Plugin as ImodPlugin
from imod.convert import genXfFile
from pwem.emlib.image.image_readers import MRCImageReader
from pyworkflow.constants import PROD, SCIPION_DEBUG_NOCLEAN
import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.object import List, String, Set
from pwem.protocols import EMProtocol
from imod import utils as imodUtils
from novactf import Plugin
from pyworkflow.utils import Message, cyanStr, makePath, redStr, envVarOn, cleanPath
from tomo.objects import Tomogram, SetOfTomograms
from tomo.utils import getCommonTsAndCtfElements

logger = logging.getLogger(__name__)


class outputs(Enum):
    Tomograms = SetOfTomograms


class ProtNovaCtf(EMProtocol):
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

    _label = 'compute defocus and reconstruct tomogram'
    _devStatus = PROD
    _possibleOutputs = outputs
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **args):
        super().__init__(**args)
        self.numberOfIntermediateStacks = List()
        self.tsDict = None
        self.ctfDict = None
        self.failedItems = []
        self.nonMatchingTsIdsMsg = String()

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam('inputSetOfTiltSeries', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series')
        form.addParam('inputSetOfCtfTomoSeries', params.PointerParam,
                      label="Associated CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      important=True,
                      help='CTF of the tilt-series.')

        group = form.addGroup('CTF')
        group.addParam('defocusStep', params.IntParam,
                       default=15,
                       label='Defocus step (nm)',
                       help='The space between min and max in Z is sliced '
                            'by defocus step. 15 nm is default step size. '
                            'See Fig. 2 in Turonova et al., 2017 for optimal number.')
        group.addParam('correctionType', params.EnumParam,
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
        group.addParam('correctAstigmatism', params.BooleanParam,
                       default=True,
                       label='Correct astigmatism?',
                       help='Correct for astigmatism in reconstruction.')

        group = form.addGroup('Tomogram')
        group.addParam('tomoThickness', params.IntParam,
                       default=400,
                       label='Tomogram thickness (voxels)',
                       help='Size of the tomogram in voxels in the Z axis (beam direction).')
        group.addParam('tomoShift', params.IntParam,
                       default=0,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Tomogram shift (voxels)',
                       help='Shift of the tomogram in voxels in the Z direction. '
                            'The shift should be set to zero even if for '
                            'reconstruction we want to shift the tomogram in z! '
                            'We assume the defocus to be estimated at the center '
                            'of mass which should correspond to the shifted '
                            'tomogram and thus here the shift should be zero.')

        group = form.addGroup('Radial filtering',
                              help='This entry controls low-pass filtering with the radial weighting '
                                   'function. The radial weighting function is linear away from the '
                                   'origin out to the distance in reciprocal space specified by the '
                                   'first value, followed by a Gaussian fall-off determined by the '
                                   'second value. Expressed in digital units 0-0.5, where 0.5 means '
                                   'Nyquist frequency.')
        group.addParam('radialFirstParameter', params.FloatParam,
                       default=0.3, label='Linear region')
        group.addParam('radialSecondParameter', params.FloatParam,
                       default=0.05, label='Gaussian fall-off')

        form.addParallelSection(threads=4)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        closeSetDeps = []
        self._initialize()
        for tsId in self.tsDict.keys():
            cInId = self._insertFunctionStep(self.convertInputStep, tsId,
                                             prerequisites=[],
                                             needsGPU=False)
            cDefId = self._insertFunctionStep(self.computeDefocusStep, tsId,
                                              prerequisites=cInId,
                                              needsGPU=False)
            recId = self._insertFunctionStep(self.computeReconstructionStep, tsId,
                                             prerequisites=cDefId,
                                             needsGPU=False)
            cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                              prerequisites=recId,
                                              needsGPU=False)
            closeSetDeps.append(cOutId)
        self._insertFunctionStep(self.closeOutputSetStep,
                                 prerequisites=closeSetDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inTsSet = self.getInputTs()
        inCtfs = self.getInputCtf()
        tsIds = set(inTsSet.getTSIds())
        ctfTsIds = set(inCtfs.getTSIds())
        # Check the common elements
        matchingTsIds = tsIds & ctfTsIds
        nonMatchingTsIds = tsIds ^ ctfTsIds
        if not matchingTsIds:
            raise Exception('No matching tsIds were found among the given sets of tilt-series and CTFs.')
        if nonMatchingTsIds:
            msg = f'Some non-matching tsIds were found: {nonMatchingTsIds}'
            self.nonMatchingTsIdsMsg.set(msg)
            logger.info(cyanStr(msg))
            self._store(self.nonMatchingTsIdsMsg)
        self.tsDict = {ts.getTsId(): ts.clone() for ts in inTsSet if ts.getTsId() in matchingTsIds}
        self.ctfDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in inCtfs if ctf.getTsId() in matchingTsIds}

    def convertInputStep(self, tsId: str):
        logger.info(cyanStr(f"tsId = {tsId} Generating tlt and defocus files...'"))
        try:
            ts = self.tsDict[tsId]
            ctf = self.ctfDict[tsId]
            presentAcqOrders = getCommonTsAndCtfElements(ts, ctf)
            if len(presentAcqOrders) == 0:
                raise Exception(f'tsId = {tsId} -> No common acquisition orders found between the '
                                f'tilt-series and the CTF.')

            logger.info(cyanStr(f"tsId = {tsId} -> present acquisition orders in both "
                                f"the tilt-series and the CTF are {presentAcqOrders}.'"))
            with self._lock:
                firstItem = ts.getFirstItem()

            # Create the folders for the tilt series
            makePath(*[self._getTsIdResultsPath(tsId),
                       self._getTsIdTmpPath(tsId)])

            # Re-stack if there are excluded views
            ts.reStack(self._getTsTmpFileName(tsId), presentAcqOrders)

            # Generate the alignment file
            if firstItem.hasTransform():
                xfFn = self._getXfFn(tsId)
                genXfFile(ts, xfFn, presentAcqOrders=presentAcqOrders)

            # Generate angle file
            tltFn = self._getTltFn(tsId)
            ts.generateTltFile(tltFn, presentAcqOrders=presentAcqOrders)

            # Generate defocus file
            defocusFn = self._getDefocusFn(tsId)
            imodUtils.generateDefocusIMODFileFromObject(ctf, defocusFn, presentAcqOrders=presentAcqOrders)
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> input conversion failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def computeDefocusStep(self, tsId: str):
        if tsId in self.failedItems:
            return
        try:
            logger.info(cyanStr(f'tsId ={tsId} -> Computing the defocus...'))
            ts = self.tsDict[tsId]
            tsFn = self._getTsTmpFileName(tsId)
            xDim, yDim, _, _ = MRCImageReader.getDimensions(tsFn)
            tltFn = self._getTltFn(tsId)
            paramsDefocus = {
                '-Algorithm': "defocus",
                '-InputProjections': tsFn,
                '-FULLIMAGE': f"{xDim},{yDim}",
                '-THICKNESS': self.tomoThickness.get(),
                '-TILTFILE': tltFn,
                '-SHIFT': "0.0,0.0",
                '-CorrectionType': self.getCorrectionType(),
                '-DefocusFileFormat': "imod",
                '-CorrectAstigmatism': 1 if self.correctAstigmatism else 0,
                '-DefocusFile': self._getDefocusFn(tsId),
                '-PixelSize': self.getInputSamplingRate() / 10,
                '-DefocusStep': self.defocusStep,
            }

            if self.tomoShift.get() > 0.01:
                defocusShiftFile = self._getDefocusShiftFn(tsId)
                paramsDefocus["-DefocusShiftFile"] = defocusShiftFile
                defocusShift = self.tomoThickness.get() / 2 + self.tomoShift.get()
                with open(defocusShiftFile, "w") as fn:
                    fn.write(f"{defocusShift}")

            Plugin.runNovactf(self, **paramsDefocus)

            # ---------- CTF correction step --------------------------------------
            logger.info(cyanStr(f"tsId = {tsId} -> Processing the intermediate stacks..."))
            acq = ts.getAcquisition()
            rotationAngle = acq.getTiltAxisAngle()
            nstacks = self.getNumberOfStacks(tsId)
            paramsCtfCorrection = {
                '-Algorithm': "ctfCorrection",
                '-InputProjections': tsFn,
                '-TILTFILE': tltFn,
                '-CorrectionType': self.getCorrectionType(),
                '-DefocusFileFormat': "imod",
                '-CorrectAstigmatism': 1 if self.correctAstigmatism.get() else 0,
                '-PixelSize': self.getInputSamplingRate() / 10,
                '-AmplitudeContrast': acq.getAmplitudeContrast(),
                '-Cs': acq.getSphericalAberration(),
                '-Volt': acq.getVoltage()
            }

            for i in range(nstacks):
                outFn = self._getStackTsFnByIndex(tsId, i)
                paramsCtfCorrection['-OutputFile'] = outFn
                paramsCtfCorrection['-DefocusFile'] = self._getDefocusFnByIndex(tsId, i)
                Plugin.runNovactf(self, **paramsCtfCorrection)

                # --------- Alignment step --------------------------------------------
                if ts.hasAlignment():
                    outAliFn = self._getAliStackTsFnByIndex(tsId, i)
                    paramsAlignment = {
                        "-input": outFn,
                        "-output": outAliFn,
                        "-xform": self._getXfFn(tsId),
                        "-AdjustOrigin": "",
                        "-NearestNeighbor": "",
                        "-taper": "1,1"
                    }
                    outFn = outAliFn

                    # Check if rotation angle is greater than 45º.
                    # If so, swap x and y dimensions to adapt output image
                    # sizes to the final sample disposition.
                    if 45 < abs(rotationAngle) < 135:
                        paramsAlignment['-size'] = f"{yDim},{xDim}"

                    args = ' '.join([f"{k} {v}" for k, v in paramsAlignment.items()])
                    ImodPlugin.runImod(self, 'newstack', args)

                # ----------- Flipping step (XYZ to XZY) ------------------------------
                outFlipFn = self._getFlipStackFnByIndex(tsId, i)
                flipArgs = ["flipyz", outFn, outFlipFn]
                ImodPlugin.runImod(self, 'clip', " ".join(flipArgs))

                # ------------- Filtering step ----------------------------------------
                paramsFilter = {
                    "-Algorithm": "filterProjections",
                    "-InputProjections": outFlipFn,
                    "-OutputFile": self._getFilterStackFnByIndex(tsId, i),
                    "-TILTFILE": tltFn,
                    "-StackOrientation": "xz",
                    "-RADIAL": f"{self.radialFirstParameter.get()},{self.radialSecondParameter.get()}"
                }

                Plugin.runNovactf(self, **paramsFilter)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> NovaCTF execution failed computing the defocus or processing the '
                                f'intermediate stacks with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def computeReconstructionStep(self, tsId: str):
        if tsId in self.failedItems:
            return
        try:
            logger.info(f"tsId = {tsId} -> Reconstructing the tomogram...")
            ts = self.tsDict[tsId]
            with self._lock:
                firstItem = ts.getFirstItem()
            tsFn = firstItem.getFileName()
            xDim, yDim, _, _ = MRCImageReader.getDimensions(tsFn)
            acq = ts.getAcquisition()
            rotationAngle = acq.getTiltAxisAngle()

            # ----------- 3D CTF step ---------------------------------------------
            recTomoFn = self._getRecTomoTmpFn(tsId)
            params3dctf = {
                "-Algorithm": "3dctf",
                "-InputProjections": self._getFilterFn(tsId),
                "-OutputFile": recTomoFn,
                "-FULLIMAGE": f"{xDim},{yDim}",
                "-TILTFILE": self._getTltFn(tsId),
                "-THICKNESS": self.tomoThickness.get(),
                "-SHIFT": f"0.0,{self.tomoShift.get()}",
                "-PixelSize": self.getInputSamplingRate() / 10,
                "-NumberOfInputStacks": self.getNumberOfStacks(tsId),
                "-Use3DCTF": 1
            }

            # Check if rotation angle is greater than 45º.
            # If so, swap x and y dimensions to adapt output image
            # sizes to the final sample disposition.
            if 45 < abs(rotationAngle) < 135:
                params3dctf['-FULLIMAGE'] = f"{yDim},{xDim}"

            Plugin.runNovactf(self, **params3dctf)

            # ---------- Trim vol - rotate around X -------------------------------
            finalTomoFn = self._getFinalRecTomoFn(tsId)
            ImodPlugin.runImod(self, 'trimvol', " ".join(["-rx", recTomoFn, finalTomoFn]))

            # Remove intermediate files. Necessary for big sets of tilt-series
            if exists(finalTomoFn) and not envVarOn(SCIPION_DEBUG_NOCLEAN):
                cleanPath(self._getTsIdTmpPath(tsId))
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> NovaCTF execution failed reconstructing the tomogram '
                                f'with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def createOutputStep(self, tsId: str):
        if tsId is self.failedItems:
            return
        try:
            finalTomoFn = self._getFinalRecTomoFn(tsId)
            if exists(finalTomoFn):
                with self._lock:
                    inSRate = self.getInputSamplingRate()
                    ts = self.tsDict[tsId]
                    acq = ts.getAcquisition()
                    outputTomos = self._getOutputTomoSet()
                    newTomogram = Tomogram()
                    newTomogram.setFileName(finalTomoFn)
                    newTomogram.setTsId(tsId)
                    newTomogram.setSamplingRate(inSRate)
                    newTomogram.fixMRCVolume(inSRate)

                    # Set default tomogram origin
                    newTomogram.setOrigin(newOrigin=None)
                    newTomogram.setAcquisition(acq)

                    outputTomos.append(newTomogram)
                    outputTomos.write()
                    self._store()
            else:
                logger.error(redStr(f'tsId = {tsId} -> Output file {finalTomoFn} was not generated. Skipping... '))
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    def closeOutputSetStep(self):
        self._closeOutputSet()
        outTomos = getattr(self, self._possibleOutputs.Tomograms.name, None)
        if not outTomos or (outTomos and len(outTomos) == 0):
            raise Exception(f'No output/s {outTomos} were generated. Please check the '
                            f'Output Log > run.stdout and run.stderr')

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []
        ctf = self.getInputCtf().getFirstItem()
        if hasattr(ctf, "_IMODDefocusFileFlag"):
            defFlag = ctf.getIMODDefocusFileFlag()
            if defFlag in [0, 4] and self.correctAstigmatism:
                validateMsgs.append("CTF estimation does not have astigmatism values.")

        return validateMsgs

    def _warnings(self):
        warnMsg = []
        ts = self.getInputTs().getFirstItem()
        if not ts.getFirstItem().hasTransform():
            warnMsg.append("The introduced tilt-series do not have alignment information.")
        return warnMsg

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

    def getNumberOfStacks(self, tsId):
        defocusFilePath = self._getDefocusFn(tsId)
        files = glob(defocusFilePath.replace(".defocus", ".defocus_*"))
        return len(files)

    def getInputSamplingRate(self):
        return self.getInputTs().getSamplingRate()

    def getCorrectionType(self):
        return "phaseflip" if self.correctionType.get() == 0 else "multiplication"

    def _getTsIdTmpPath(self, tsId: str) -> str:
        return self._getTmpPath(tsId)

    def _getTsIdResultsPath(self, tsId: str) -> str:
        return self._getExtraPath(tsId)

    def _getTsTmpFileName(self, tsId: str) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}.mrc')

    def _getTltFn(self, tsId: str) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}.tlt')

    def _getDefocusFn(self, tsId: str) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}.defocus')

    def _getDefocusShiftFn(self, tsId: str) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}.def_shift')

    def _getXfFn(self, tsId: str) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}.xf')

    def _getDefocusFnByIndex(self, tsId: str, index: int) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}.defocus_{index}')

    def _getStackTsFnByIndex(self, tsId: str, index: int) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}.mrc_{index}')

    def _getAliStackTsFnByIndex(self, tsId: str, index: int) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}_ali.mrc_{index}')

    def _getFlipStackFnByIndex(self, tsId: str, index: int) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}_flip.mrc_{index}')

    def _getFilterFn(self, tsId: str) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}_filter.mrc')

    def _getFilterStackFnByIndex(self, tsId: str, index: int) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}_filter.mrc_{index}')

    def _getRecTomoTmpFn(self, tsId: str) -> str:
        return join(self._getTsIdTmpPath(tsId), f'{tsId}_rec.mrc')

    def _getFinalRecTomoFn(self, tsId: str) -> str:
        return join(self._getTsIdResultsPath(tsId), f'{tsId}.mrc')

    def _getOutputTomoSet(self) -> SetOfTomograms:
        outTomosAttribName = self._possibleOutputs.Tomograms.name
        outTomos = getattr(self, outTomosAttribName, None)
        if outTomos:
            outTomos.enableAppend()
        else:
            outTomos = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            outTomos.copyInfo(self.getInputTs())
            outTomos.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outTomosAttribName: outTomos})
            self._defineSourceRelation(self.inputSetOfTiltSeries, outTomos)

        return outTomos