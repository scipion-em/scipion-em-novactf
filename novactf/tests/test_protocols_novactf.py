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

from pyworkflow.tests import *
from imod.protocols import *
from novactf.protocols import protocol_tomoReconstruction
from pwem.emlib.image import ImageHandler
import tomo

class TestNovaCtfBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def _runImportTiltSeries(cls, filesPath, pattern, voltage, magnification, sphericalAberration, amplitudeContrast,
                             samplingRate, doseInitial, dosePerFrame, anglesFrom=0, minAngle=0.0, maxAngle=0.0,
                             stepAngle=1.0):
        cls.protImportTS = cls.newProtocol(tomo.protocols.ProtImportTs,
                                           filesPath=filesPath,
                                           filesPattern=pattern,
                                           voltage=voltage,
                                           anglesFrom=anglesFrom,
                                           magnification=magnification,
                                           sphericalAberration=sphericalAberration,
                                           amplitudeContrast=amplitudeContrast,
                                           samplingRate=samplingRate,
                                           doseInitial=doseInitial,
                                           dosePerFrame=dosePerFrame,
                                           minAngle=minAngle,
                                           maxAngle=maxAngle,
                                           stepAngle=stepAngle)
        cls.launchProtocol(cls.protImportTS)
        return cls.protImportTS

