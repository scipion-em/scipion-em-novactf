==============
NovaCTF plugin
==============

This plugin provides wrappers for several programs of `NovaCTF <https://github.com/turonova/novaCTF>`_ software suite.

.. image:: https://img.shields.io/pypi/v/scipion-em-novactf.svg
        :target: https://pypi.python.org/pypi/scipion-em-novactf
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-novactf.svg
        :target: https://pypi.python.org/pypi/scipion-em-novactf
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-novactf.svg
        :target: https://pypi.python.org/pypi/scipion-em-novactf
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-novactf?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-novactf
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-novactf
        :target: https://pypi.python.org/pypi/scipion-em-novactf
        :alt: Downloads

Installation
------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

   .. code-block::

      scipion installp -p scipion-em-novactf

b) Developer's version

   * download repository

   .. code-block::

      git clone -b devel https://github.com/scipion-em/scipion-em-novactf.git

   * install

   .. code-block::

      scipion installp -p /path/to/scipion-em-novactf --devel

NovaCTF binaries will be downloaded and installed automatically with the plugin, but you can also link an existing installation. Default installation path assumed is ``software/em/novactf-master``, if you want to change it, set *NOVACTF_HOME* in ``scipion.conf`` file to the folder where the NovaCTF is installed.

To check the installation, simply run the test below:

``scipion3 tests novactf.tests.test_protocols_novactf.TestNovaCtfReconstructionWorkflow``

Supported versions
------------------

master (8 Nov 2018)

Protocols
---------

* compute defocus array
* 3D CTF correction and reconstruction

References
----------

1. Beata Turoňová and Florian K.M. Schur and William Wan and John A.G. Briggs. Efficient 3D-CTF correction for cryo-electron tomography using NovaCTF improves subtomogram averaging resolution to 3.4Å. Journal of Structural Biology 199 (3), p. 187-195, 2017. doi: https://doi.org/10.1016/j.jsb.2017.07.007
