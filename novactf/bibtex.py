# coding: latin-1
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
"""

@article{Turonova2017,
title = "Efficient 3D-CTF correction for cryo-electron tomography using NovaCTF improves subtomogram averaging resolution to 3.4 Å",
journal = "Journal of Structural Biology",
volume = "199",
number = "1",
pages = "187 - 195",
year = "2017",
doi = "https://doi.org/10.1016/j.jsb.2017.07.007",
url = "https://www.sciencedirect.com/science/article/pii/S1047847717301272?via%3Dihub",
author = "Beata Turonova and Florian K.M.Schur and William Wan and John A.G. Briggsab",
abstract = "Cryo-electron tomography (cryo-ET) allows cellular ultrastructures and macromolecular complexes to be imaged in three-dimensions in their native environments. Cryo-electron tomograms are reconstructed from projection images taken at defined tilt-angles. In order to recover high-resolution information from cryo-electron tomograms, it is necessary to measure and correct for the contrast transfer function (CTF) of the microscope. Most commonly, this is performed using protocols that approximate the sample as a two-dimensional (2D) plane. This approximation accounts for differences in defocus and therefore CTF across the tilted sample. It does not account for differences in defocus of objects at different heights within the sample; instead, a 3D approach is required. Currently available approaches for 3D-CTF correction are computationally expensive and have not been widely implemented. Here we simulate the benefits of 3D-CTF correction for high-resolution subtomogram averaging, and present a user-friendly, computationally-efficient 3D-CTF correction tool, NovaCTF, that is compatible with standard tomogram reconstruction workflows in IMOD. We validate the approach on synthetic data and test it using subtomogram averaging of real data. Consistent with our simulations, we find that 3D-CTF correction allows high-resolution structures to be obtained with much smaller subtomogram averaging datasets than are required using 2D-CTF. We also show that using equivalent dataset sizes, 3D-CTF correction can be used to obtain higher-resolution structures. We present a 3.4 Å resolution structure determined by subtomogram averaging."
}



"""
