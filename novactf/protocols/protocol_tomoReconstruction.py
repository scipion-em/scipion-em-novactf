# **************************************************************************
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
# **************************************************************************

import os

from pyworkflow import BETA
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram

from imod import Plugin as imodPlugin
import imod.utils as imodUtils

from .. import Plugin


class ProtNovaCtfTomoReconstruction(EMProtocol, ProtTomoBase):
    """
    Tomogram reconstruction with 3D CTF correction by novaCTF.

    More info:
            https://github.com/turonova/novaCTF
    """

    _label = '3D CTF correction and reconstruction'
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        ProtTomoBase.__init__(self)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('protTomoCtfDefocus',
                      params.PointerParam,
                      label="NovaCTF compute defocus run",
                      pointerClass='ProtNovaCtfTomoDefocus')

        form.addParam('applyAlignment', params.BooleanParam,
                      default=True,
                      label="Apply tilt-series alignment?")

        form.addSection("Erase gold beads")
        form.addParam('doEraseGold', params.BooleanParam,
                      default=False, label='Erase gold beads',
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
        group.addParam('radialFirstParameter',
                       params.FloatParam,
                       default=0.3,
                       label='Linear region')
        group.addParam('radialSecondParameter',
                       params.FloatParam,
                       default=0.05,
                       label='Gaussian fall-off')

        form.addParallelSection(threads=4, mpi=1)

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        allCreateOutputId = []
        nstacks = self.protTomoCtfDefocus.get().numberOfIntermediateStacks

        for index, ts in enumerate(self.getInputTs()):
            convertInputId = self._insertFunctionStep(self.convertInputStep,
                                                      ts.getObjId())

            intermediateStacksId = []

            for counter in range(nstacks[index].get()):
                ctfId = self._insertFunctionStep(self.processIntermediateStacksStep,
                                                 ts.getObjId(),
                                                 counter,
                                                 prerequisites=[convertInputId])
                intermediateStacksId.append(ctfId)

            reconstructionId = self._insertFunctionStep(self.computeReconstructionStep,
                                                        ts.getObjId(), nstacks[index].get(),
                                                        prerequisites=intermediateStacksId)

            createOutputId = self._insertFunctionStep(self.createOutputStep,
                                                      ts.getObjId(),
                                                      prerequisites=[reconstructionId])
            allCreateOutputId.append(createOutputId)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 prerequisites=allCreateOutputId)

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        with self._lock:
            ts = self.getInputTs()[tsObjId]
            firstItem = ts.getFirstItem()

        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        outputTsFileName = os.path.join(tmpPrefix,
                                        firstItem.parseFileName(extension=".mrc"))
        self.info("Tilt series %s linked." % tsId)
        path.createLink(firstItem.getFileName(), outputTsFileName)

        # Generate angle file
        outputTltFileName = os.path.join(tmpPrefix,
                                         firstItem.parseFileName(extension=".tlt"))
        ts.generateTltFile(outputTltFileName)

        if firstItem.hasTransform():
            # Generate transformation matrices file
            outputTmFileName = os.path.join(tmpPrefix,
                                            firstItem.parseFileName(extension=".xf"))
            imodUtils.formatTransformFile(ts, outputTmFileName)

    def processIntermediateStacksStep(self, tsObjId, counter):

        with self._lock:
            ts = self.getInputTs()[tsObjId]
            firstItem = ts.getFirstItem()

        tsId = ts.getTsId()
        tmpPrefix = self._getTmpPath(tsId)

        defocusFilePath = self.protTomoCtfDefocus.get().getDefocusFileName(tsId) + "_"
        tltFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))
        outputFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".mrc_"))

        # CTF correction step
        paramsCtfCorrection = {
            'Algorithm': "ctfCorrection",
            'InputProjections': os.path.join(tmpPrefix,
                                             firstItem.parseFileName(extension=".mrc")),
            'OutputFile': outputFilePath + str(counter),
            'DefocusFile': defocusFilePath + str(counter),
            'TiltFile': tltFilePath,
            'CorrectionType': self.getCorrectionType(),
            'DefocusFileFormat': "imod",
            'CorrectAstigmatism': self.protTomoCtfDefocus.get().correctAstigmatism.get(),
            'PixelSize': self.getInputTs().getSamplingRate() / 10,
            'AmplitudeContrast':
                self.getInputTs().getAcquisition().getAmplitudeContrast(),
            'SphericalAberration':
                self.getInputTs().getAcquisition().getSphericalAberration(),
            'Voltage': self.getInputTs().getAcquisition().getVoltage()
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
                            "-Volt %(Voltage)d "

        Plugin.runNovactf(self, 'novaCTF', argsCtfCorrection % paramsCtfCorrection)
        currentFn = outputFilePath + str(counter)

        # Alignment step
        if self.applyAlignment and firstItem.hasTransform():
            paramsAlignment = {
                'input': currentFn,
                'output': os.path.join(tmpPrefix,
                                       firstItem.parseFileName(suffix="_ali",
                                                               extension=".mrc_")) + str(counter),
                'xform': os.path.join(tmpPrefix,
                                      firstItem.parseFileName(extension=".xf"))
            }

            argsAlignment = "-input %(input)s " \
                            "-output %(output)s " \
                            "-xform %(xform)s -AdjustOrigin -NearestNeighbor -taper 1,1 "

            rotationAngle = ts.getAcquisition().getTiltAxisAngle()

            # Check if rotation angle is greater than 45º.
            # If so, swap x and y dimensions to adapt output image
            # sizes to the final sample disposition.
            if 45 < abs(rotationAngle) < 135:
                paramsAlignment.update({
                    'size': "%d,%d" % (firstItem.getYDim(), firstItem.getXDim())
                })

                argsAlignment += "-size %(size)s "

            imodPlugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)
            currentFn = os.path.join(tmpPrefix,
                                     firstItem.parseFileName(suffix="_ali",
                                                             extension=".mrc_")) + str(counter)

        # Erase gold step  # TODO: fixme
        if self.doEraseGold:
            lm = self.inputSetOfLandmarkModels.get().getLandmarkModelFromTsId(tsId=tsId)

            # apply alignment to fid. model
            imodPlugin.runImod(self, 'xfmodel',
                               f'-XformsToApply {outputTmFileName}'
                               f' {lm.getModelName()} '
                               f'{os.path.join(tmpPrefix, firstItem.parseFileName(suffix="_erase", extension=".fid"))}')

            paramsCcderaser = {
                'inputFile': currentFn,
                'outputFile': os.path.join(tmpPrefix,
                                           firstItem.parseFileName(suffix="_erase",
                                                                   extension=".mrc_" + str(counter))),
                'modelFile': os.path.join(tmpPrefix,
                                          firstItem.parseFileName(suffix="_erase", extension=".fid")),
                'betterRadius': self.goldDiam.get() / 2,
                'polynomialOrder': 0,
                'circleObjects': "/"
            }

            argsCcderaser = "-InputFile %(inputFile)s " \
                            "-OutputFile %(outputFile)s " \
                            "-ModelFile %(modelFile)s " \
                            "-BetterRadius %(betterRadius)f " \
                            "-PolynomialOrder %(polynomialOrder)d " \
                            "-CircleObjects %(circleObjects)s " \
                            "-MergePatches 1 " \
                            "-ExcludeAdjacent " \
                            "-SkipTurnedOffPoints 1 " \
                            "-ExpandCircleIterations 3 "

            imodPlugin.runImod(self, 'ccderaser', argsCcderaser % paramsCcderaser)
            currentFn = os.path.join(tmpPrefix,
                                     firstItem.parseFileName(suffix="_erase",
                                                             extension=".mrc_" + str(counter)))

        # Flipping step (XYZ to XZY)
        paramsClip = {
            'inputFilePath': currentFn,
            'outputFilePath': os.path.join(tmpPrefix,
                                           firstItem.parseFileName(suffix="_flip",
                                                                   extension=".mrc_" + str(counter))),
        }

        argsClip = "flipyz " \
                   "%(inputFilePath)s " \
                   "%(outputFilePath)s "

        imodPlugin.runImod(self, 'clip', argsClip % paramsClip)
        currentFn = os.path.join(tmpPrefix,
                                 firstItem.parseFileName(suffix="_flip",
                                                         extension=".mrc_" + str(counter)))

        # Filtering step
        outputFilePath = os.path.join(tmpPrefix,
                                      firstItem.parseFileName(suffix="_filter", extension=".mrc_"))

        paramsFilterProjections = {
            'Algorithm': "filterProjections",
            'InputProjections': currentFn,
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
                                "-RADIAL %(Radial)s "

        Plugin.runNovactf(self, 'novaCTF', argsFilterProjections % paramsFilterProjections)

    def computeReconstructionStep(self, tsObjId, nstacks):
        with self._lock:
            ts = self.getInputTs()[tsObjId]
            firstItem = ts.getFirstItem()

        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        outputFilePath = os.path.join(tmpPrefix,
                                      firstItem.parseFileName(suffix="_rec", extension=".mrc"))
        tltFilePath = os.path.join(tmpPrefix,
                                   firstItem.parseFileName(extension=".tlt"))

        params3dctf = {
            'Algorithm': "3dctf",
            'InputProjections': os.path.join(tmpPrefix, firstItem.parseFileName(suffix="_filter",
                                                                                extension=".mrc")),
            'OutputFile': outputFilePath,
            'FullImage': "%d,%d" % (firstItem.getXDim(), firstItem.getYDim()),
            'TiltFile': tltFilePath,
            'Thickness': self.protTomoCtfDefocus.get().tomoThickness.get(),
            'Shift': "0.0," + str(self.protTomoCtfDefocus.get().tomoShift.get()),
            'PixelSize': self.getInputTs().getSamplingRate() / 10,
            'NumberOfStacks': nstacks
        }

        args3dctf = "-Algorithm %(Algorithm)s " \
                    "-InputProjections %(InputProjections)s " \
                    "-OutputFile %(OutputFile)s " \
                    "-FULLIMAGE %(FullImage)s " \
                    "-TILTFILE %(TiltFile)s " \
                    "-THICKNESS %(Thickness)d " \
                    "-SHIFT %(Shift)s " \
                    "-PixelSize %(PixelSize)f " \
                    "-NumberOfInputStacks %(NumberOfStacks)d " \
                    "-Use3DCTF 1 "

        # Check if rotation angle is greater than 45º.
        # If so, swap x and y dimensions to adapt output image
        # sizes to the final sample disposition.

        rotationAngle = ts.getAcquisition().getTiltAxisAngle()

        if 45 < abs(rotationAngle) < 135:
            params3dctf['FullImage'] = "%d,%d" % (firstItem.getYDim(), firstItem.getXDim())

        Plugin.runNovactf(self, 'novaCTF', args3dctf % params3dctf)

        # Trim vol - rotate around X
        paramsTrimvol = {
            'inputFilePath': outputFilePath,
            'outputFilePath': os.path.join(extraPrefix,
                                           firstItem.parseFileName(extension=".mrc"))
        }

        argsTrimvol = "-rx %(inputFilePath)s " \
                      "%(outputFilePath)s "

        imodPlugin.runImod(self, 'trimvol', argsTrimvol % paramsTrimvol)

        # Remove intermediate files. Necessary for big sets of tilt-series
        path.cleanPath(self._getTmpPath(tsId))

    def createOutputStep(self, tsObjId):
        with self._lock:
            ts = self.getInputTs()[tsObjId]
            firstItem = ts.getFirstItem()
            acq = ts.getAcquisition()

        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        outputTomos = self.getOutputSetOfTomograms()
        outputFn = os.path.join(extraPrefix, firstItem.parseFileName(extension=".mrc"))

        if os.path.exists(outputFn):
            newTomogram = Tomogram()
            newTomogram.setLocation(outputFn)
            newTomogram.setTsId(tsId)
            newTomogram.setSamplingRate(ts.getSamplingRate())

            # Set default tomogram origin
            newTomogram.setOrigin(newOrigin=None)
            newTomogram.setAcquisition(acq)

            outputTomos.append(newTomogram)
            outputTomos.write()
            self._store()

    def closeOutputSetsStep(self):
        self.getOutputSetOfTomograms().setStreamState(Set.STREAM_CLOSED)
        self._store()

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsg = []

        ts = self.getInputTs()
        if self.applyAlignment and not ts.getFirstItem().getFirstItem().hasTransform():
            validateMsg.append("Input tilt-series do not have alignment "
                               "information! You cannot apply alignment.")

        if self.doEraseGold and not self.inputSetOfLandmarkModels.hasValue():
            validateMsg.append("You have to provide input set of landmarks to erase gold.")

        return validateMsg

    def _summary(self):
        summary = []

        if hasattr(self, 'outputSetOfTomograms'):
            summary.append(f"Input tilt-series: {self.getInputTs().getSize()}\n"
                           f"CTF corrected tomos calculated: "
                           f"{self.outputSetOfTomograms.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")

        return summary

    def _methods(self):
        methods = []

        if hasattr(self, 'outputSetOfTomograms'):
            methods.append(f"{self.outputSetOfTomograms.getSize()} "
                           "3D CTF corrected tomograms have been calculated "
                           "with NovaCTF.")

        return methods

    # --------------------------- UTILS functions ----------------------------
    def getInputTs(self, pointer=False):
        if pointer:
            return self.protTomoCtfDefocus.get().inputSetOfTiltSeries
        else:
            return self.protTomoCtfDefocus.get().inputSetOfTiltSeries.get()

    def getCorrectionType(self):
        if self.protTomoCtfDefocus.get().correctionType.get() == 0:
            correctionType = "phaseflip"
        else:
            correctionType = "multiplication"

        return correctionType

    def getOutputSetOfTomograms(self):
        if hasattr(self, "outputSetOfTomograms"):
            self.outputSetOfTomograms.enableAppend()
        else:
            outputSetOfTomograms = self._createSetOfTomograms()
            outputSetOfTomograms.copyInfo(self.getInputTs())
            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTomograms=outputSetOfTomograms)
            self._defineSourceRelation(self.getInputTs(pointer=True), outputSetOfTomograms)

        return self.outputSetOfTomograms
