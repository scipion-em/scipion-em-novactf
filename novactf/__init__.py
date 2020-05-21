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

import pwem

_logo = ""
_references = []


class Plugin(pwem.Plugin):
    @classmethod
    def _defineVariables(cls):
        pass

    @classmethod
    def getEnviron(cls):
        return None

    @classmethod
    def getDependencies(cls):
        neededPrograms = ['git', 'fftw3', 'fftw3f', 'lib64', 'gcc']

        return neededPrograms

    @classmethod
    def defineBinaries(cls, env):
        version = '1.0.0'
        NOVACTF_INSTALLED = 'novactf_%s_installed' % version
        warningMsg1 = "'WARNING: if program fails when executing reinstall in scipion3/software/em/novactf-1.0.0 " \
                      "setting the path to the include files from fftw3 as well as following fftw3 libraries: fftw3 " \
                      "fftw3f. Typically both libraries should be in one folder and thus one library path should be " \
                      "sufficient. Example command:'"
        warningMsg2 = "'make includepath = \"path_to_fftw_include_files\" libpath = \"path_to_fftw_libraries\"'"
        warningMsg3 = "'NovaCTF also requires standard libraries with the standard path being " \
                      "/usr/lib64 - if the libraries are elsewhere open makefile and change the path accordingly'"

        # Download git repo
        installationCmd = 'git clone https://github.com/turonova/novaCTF.git && '

        # Binaries compilation
        installationCmd += 'cd novaCTF && ' \
                           'make && ' \
                           'cd .. && '

        # Create installation finished flag file
        installationCmd += 'touch %s && ' % NOVACTF_INSTALLED

        # Display warning message
        installationCmd += 'echo %s && ' % warningMsg1
        installationCmd += 'echo %s && ' % warningMsg2
        installationCmd += 'echo %s' % warningMsg3

        env.addPackage('novactf',
                       version=version,
                       tar='void.tgz',
                       neededProgs=cls.getDependencies(),
                       commands=[(installationCmd, NOVACTF_INSTALLED)],
                       default=True)
