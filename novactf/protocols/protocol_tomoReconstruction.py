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
import os
from enum import Enum

from pyworkflow.constants import PROD
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfTomograms, Tomogram

from imod import Plugin as imodPlugin
from imod import utils as imodUtils
from novactf import Plugin


class outputs(Enum):
    Tomograms = SetOfTomograms


class ProtNovaCtfReconstruction(EMProtocol, ProtTomoBase):
    """
    Tomogram reconstruction with 3D CTF correction by novaCTF.

    More info:
            https://github.com/turonova/novaCTF
    """

    _label = '3D CTF correction and reconstruction'
    _devStatus = PROD
    _possibleOutputs = outputs

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _initialize(self):
        self._createFilenameTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        tmpPath = lambda p: self._getTmpPath("%(tsId)s", "%(tsId)s" + p)
        myDict = {
            'tltFn': tmpPath(".tlt"),
            'xfFn': tmpPath(".xf"),
            'inputTsFn': tmpPath(".mrc"),
            'stackTsFn': tmpPath(".mrc_%(counter)d"),
            'stackAliFn': tmpPath("_ali.mrc_%(counter)d"),
            'stackEraseFn': tmpPath("_erase.mrc_%(counter)d"),
            'eraseFidFn': tmpPath("_erase.fid"),
            'stackFlipFn': tmpPath("_flip.mrc_%(counter)d"),
            'stackFilterFn': tmpPath("_filter.mrc_%(counter)d"),
            'filterFn': tmpPath("_filter.mrc"),
            'recFn': tmpPath("_rec.mrc"),
            'outputTsFn': self._getExtraPath("%(tsId)s", "%(tsId)s.mrc"),
        }

        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('protNovaCtfDefocus',
                      params.PointerParam,
                      label="NovaCTF compute defocus run",
                      pointerClass='ProtNovaCtfDefocus')

        form.addParam('applyAlignment', params.BooleanParam,
                      default=True,
                      label="Apply tilt-series alignment?")

        form.addSection("Erase gold beads")
        form.addParam('doEraseGold', params.BooleanParam,
                      default=False,
                      label='Erase gold beads',
                      help='Remove the gold beads from the tilt-series.')
        form.addParam('inputSetOfLandmarkModels',
                      params.PointerParam,
                      allowsNull=True,
                      condition='doEraseGold',
                      pointerClass='SetOfLandmarkModels',
                      label='Input set of fiducial models',
                      help='Set of fid. models with no gaps after alignment')
        form.addParam('goldDiam', params.IntParam,
                      condition='doEraseGold',
                      default=18,
                      label='Bead diameter (px)',
                      help="For circle objects, this entry "
                           "specifies a radius to use for points "
                           "without an individual point size "
                           "instead of the object's default sphere "
                           "radius. This entry is floating point "
                           "and can be used to overcome the "
                           "limitations of having an integer "
                           "default sphere radius. If there are "
                           "multiple circle objects, enter one "
                           "value to apply to all objects or a "
                           "value for each object.")

        form.addSection("Filtering")
        group = form.addGroup('Radial filtering',
                              help='This entry controls low-pass filtering with the radial weighting '
                                   'function. The radial weighting function is linear away from the '
                                   'origin out to the distance in reciprocal space specified by the '
                                   'first value, followed by a Gaussian fall-off determined by the '
                                   'second value.')
        group.addParam('radialFirstParameter', params.FloatParam,
                       default=0.3, label='Linear region')
        group.addParam('radialSecondParameter', params.FloatParam,
                       default=0.05, label='Gaussian fall-off')

        form.addParallelSection(threads=8)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        allCreateOutputId = []
        nstacks = self.getInputProt().numberOfIntermediateStacks

        for index, ts in enumerate(self.getInputTs().iterItems()):
            tsId = ts.getTsId()
            objId = ts.getObjId()
            convertInputId = self._insertFunctionStep(self.convertInputStep,
                                                      objId, tsId)

            intermediateStacksId = []
            for counter in range(nstacks[index].get()):
                ctfId = self._insertFunctionStep(self.processIntermediateStacksStep,
                                                 objId, tsId, counter,
                                                 prerequisites=[convertInputId])
                intermediateStacksId.append(ctfId)

            reconstructId = self._insertFunctionStep(self.computeReconstructionStep,
                                                     objId, tsId,
                                                     nstacks[index].get(),
                                                     prerequisites=intermediateStacksId)

            createOutputId = self._insertFunctionStep(self.createOutputStep,
                                                      objId, tsId,
                                                      prerequisites=[reconstructId])
            allCreateOutputId.append(createOutputId)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 prerequisites=allCreateOutputId)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, tsObjId, tsId):
        # Create the folders for the tilt series
        path.makePath(self._getTmpPath(tsId))
        path.makePath(self._getExtraPath(tsId))

        with self._lock:
            ts = self.getInputTs()[tsObjId]
            firstItem = ts.getFirstItem()
            tsFn = firstItem.getFileName()

            if firstItem.hasTransform():
                # Generate transformation matrices file
                outputTmFile = self._getFileName("xfFn", tsId=tsId)
                imodUtils.formatTransformFile(ts, outputTmFile)

            # Generate angle file
            ts.generateTltFile(self._getFileName("tltFn", tsId=tsId))

        # Link tilt series file
        path.createLink(tsFn,
                        self._getFileName("inputTsFn", tsId=tsId))

    def processIntermediateStacksStep(self, tsObjId, tsId, counter):
        self.info(f"Processing {tsId}, intermediate stack #{counter}")
        inputProt = self.getInputProt()
        inputProt._createFilenameTemplates()
        defocusFn = inputProt._getFileName("stackDefocusFn",
                                           tsId=tsId, counter=counter)

        with self._lock:
            ts = self.getInputTs()[tsObjId]
            rotationAngle = ts.getAcquisition().getTiltAxisAngle()
            firstItem = ts.getFirstItem()
            xDim, yDim, _ = firstItem.getDim()
            hasTransform = firstItem.hasTransform()

        # ---------- CTF correction step --------------------------------------
        paramsCtfCorrection = {
            '-Algorithm': "ctfCorrection",
            '-InputProjections': self._getFileName("inputTsFn", tsId=tsId),
            '-OutputFile': self._getFileName("stackTsFn",
                                             tsId=tsId, counter=counter),
            '-DefocusFile': defocusFn,
            '-TILTFILE': self._getFileName("tltFn", tsId=tsId),
            '-CorrectionType': self.getCorrectionType(),
            '-DefocusFileFormat': "imod",
            '-CorrectAstigmatism': 1 if inputProt.correctAstigmatism else 0,
            '-PixelSize': self.getInputSamplingRate() / 10,
            '-AmplitudeContrast':
                self.getInputAcquisition().getAmplitudeContrast(),
            '-Cs':
                self.getInputAcquisition().getSphericalAberration(),
            '-Volt': self.getInputAcquisition().getVoltage()
        }

        Plugin.runNovactf(self, **paramsCtfCorrection)

        currentFn = self._getFileName("stackTsFn",
                                      tsId=tsId, counter=counter)

        # --------- Alignment step --------------------------------------------
        if self.applyAlignment and hasTransform:
            paramsAlignment = {
                "-input": currentFn,
                "-output": self._getFileName("stackAliFn",
                                             tsId=tsId, counter=counter),
                "-xform": self._getFileName("xfFn", tsId=tsId),
                "-AdjustOrigin": "",
                "-NearestNeighbor": "",
                "-taper": "1,1"
            }

            # Check if rotation angle is greater than 45º.
            # If so, swap x and y dimensions to adapt output image
            # sizes to the final sample disposition.
            if 45 < abs(rotationAngle) < 135:
                paramsAlignment['-size'] = f"{yDim},{xDim}"

            args = ' '.join([f"{k} {v}" for k, v in paramsAlignment.items()])
            imodPlugin.runImod(self, 'newstack', args)

            currentFn = self._getFileName("stackAliFn",
                                          tsId=tsId, counter=counter)

        # ---------- Erase gold step ------------------------------------------
        # TODO: fixme
        if self.doEraseGold:
            lm = self.inputSetOfLandmarkModels.get().getLandmarkModelFromTsId(tsId=tsId)

            # apply alignment to fid. model
            paramsXfModel = {
                "-XformsToApply": outputTmFileName,
                lm.getModelName(): "",
                self._getFileName("eraseFidFn", tsId=tsId): ""
            }
            args = ' '.join([f"{k} {v}" for k, v in paramsXfModel.items()])
            imodPlugin.runImod(self, 'xfmodel', args)

            paramsCcderaser = {
                "-InputFile": currentFn,
                "-OutputFile": self._getFileName("stackEraseFn",
                                                 tsId=tsId, counter=counter),
                "-ModelFile": self._getFileName("eraseFidFn", tsId=tsId),
                "-BetterRadius": self.goldDiam.get() / 2,
                "-PolynomialOrder": 0,
                "-CircleObjects": "/",
                "-MergePatches": 1,
                "-ExcludeAdjacent": "",
                "-SkipTurnedOffPoints": 1,
                "-ExpandCircleIterations": 3
            }

            args = ' '.join([f"{k} {v}" for k, v in paramsCcderaser.items()])
            imodPlugin.runImod(self, 'ccderaser', args)

            currentFn = self._getFileName("stackEraseFn",
                                          tsId=tsId, counter=counter)

        # ----------- Flipping step (XYZ to XZY) ------------------------------
        flipArgs = [
            "flipyz",
            currentFn,
            self._getFileName("stackFlipFn", tsId=tsId, counter=counter)
        ]
        imodPlugin.runImod(self, 'clip', " ".join(flipArgs))

        currentFn = self._getFileName("stackFlipFn", tsId=tsId, counter=counter)

        # ------------- Filtering step ----------------------------------------
        paramsFilter = {
            "-Algorithm": "filterProjections",
            "-InputProjections": currentFn,
            "-OutputFile": self._getFileName("stackFilterFn",
                                             tsId=tsId, counter=counter),
            "-TILTFILE": self._getFileName("tltFn", tsId=tsId),
            "-StackOrientation": "xz",
            "-RADIAL":
                f"{self.radialFirstParameter.get()},{self.radialSecondParameter.get()}"
        }

        Plugin.runNovactf(self, **paramsFilter)

    def computeReconstructionStep(self, tsObjId, tsId, nstacks):
        with self._lock:
            ts = self.getInputTs()[tsObjId]
            firstItem = ts.getFirstItem()
            xDim, yDim, _ = firstItem.getDim()
            rotationAngle = ts.getAcquisition().getTiltAxisAngle()

        # ----------- 3D CTF step ---------------------------------------------
        params3dctf = {
            "-Algorithm": "3dctf",
            "-InputProjections": self._getFileName("filterFn", tsId=tsId),
            "-OutputFile": self._getFileName("recFn", tsId=tsId),
            "-FULLIMAGE": f"{xDim},{yDim}",
            "-TILTFILE": self._getFileName("tltFn", tsId=tsId),
            "-THICKNESS": self.getInputProt().tomoThickness,
            "-SHIFT": f"0.0,{self.getInputProt().tomoShift.get()}",
            "-PixelSize": self.getInputSamplingRate() / 10,
            "-NumberOfInputStacks": nstacks,
            "-Use3DCTF": 1
        }

        # Check if rotation angle is greater than 45º.
        # If so, swap x and y dimensions to adapt output image
        # sizes to the final sample disposition.
        if 45 < abs(rotationAngle) < 135:
            params3dctf['-FULLIMAGE'] = f"{yDim},{xDim}"

        Plugin.runNovactf(self, **params3dctf)

        # ---------- Trim vol - rotate around X -------------------------------
        inputFn = self._getFileName("recFn", tsId=tsId)
        outputFn = self._getFileName("outputTsFn", tsId=tsId)
        imodPlugin.runImod(self, 'trimvol',
                           " ".join(["-rx", inputFn, outputFn]))

        # Remove intermediate files. Necessary for big sets of tilt-series
        if os.path.exists(outputFn):
            path.cleanPath(self._getTmpPath(tsId))

    def createOutputStep(self, tsObjId, tsId):
        with self._lock:
            ts = self.getInputTs()[tsObjId]
            acq = ts.getAcquisition()

        outputTomos = self.getOutputSetOfTomograms()
        outputFn = self._getFileName("outputTsFn", tsId=tsId)

        if os.path.exists(outputFn):
            newTomogram = Tomogram()
            newTomogram.setLocation(outputFn)
            newTomogram.setTsId(tsId)
            newTomogram.setSamplingRate(self.getInputSamplingRate())

            # Set default tomogram origin
            newTomogram.setOrigin(newOrigin=None)
            newTomogram.setAcquisition(acq)

            outputTomos.append(newTomogram)
            outputTomos.write()
            self._store()

    def closeOutputSetsStep(self):
        self.getOutputSetOfTomograms().setStreamState(Set.STREAM_CLOSED)
        self._store()

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []
        ts = self.getInputTs()

        if self.applyAlignment and not ts.getFirstItem().getFirstItem().hasTransform():
            validateMsgs.append("Input tilt-series do not have alignment "
                                "information! You cannot apply alignment.")

        if self.doEraseGold and not self.inputSetOfLandmarkModels.hasValue():
            validateMsgs.append("You have to provide input set of landmarks to erase gold.")

        return validateMsgs

    def _summary(self):
        summary = []

        if hasattr(self, outputs.Tomograms.name):
            summary.append(f"Input tilt-series: {self.getInputTs().getSize()}\n"
                           "CTF corrected tomos calculated: "
                           f"{getattr(self, outputs.Tomograms.name).getSize()}")
        else:
            summary.append("Outputs are not ready yet.")

        return summary

    def _methods(self):
        methods = []

        if hasattr(self, outputs.Tomograms.name):
            methods.append(f"{getattr(self, outputs.Tomograms.name).getSize()} "
                           "3D CTF corrected tomograms have been calculated "
                           "with NovaCTF.")

        return methods

    # --------------------------- UTILS functions -----------------------------
    def getOutputSetOfTomograms(self):
        outputName = outputs.Tomograms.name
        if hasattr(self, outputName):
            getattr(self, outputName).enableAppend()
        else:
            outputSetOfTomograms = self._createSetOfTomograms()
            outputSetOfTomograms.copyInfo(self.getInputTs())
            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputName: outputSetOfTomograms})
            self._defineSourceRelation(self.getInputTs(pointer=True),
                                       outputSetOfTomograms)

        return getattr(self, outputName)

    def getInputTs(self, pointer=False):
        if pointer:
            return self.getInputProt().inputSetOfTiltSeries
        else:
            return self.getInputProt().inputSetOfTiltSeries.get()

    def getInputProt(self):
        return self.protNovaCtfDefocus.get()

    def getInputSamplingRate(self):
        return self.getInputTs().getSamplingRate()

    def getInputAcquisition(self):
        return self.getInputTs().getAcquisition()

    def getCorrectionType(self):
        if self.getInputProt().correctionType.get() == 0:
            correctionType = "phaseflip"
        else:
            correctionType = "multiplication"

        return correctionType
