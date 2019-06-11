molecbio
========

Molecular biology informatics functions

INSTALLATION
============

After cloning the repository, change directory into it and install the package with the command: 'python setup.py install'. Depending on your configuration, you may need to use: 'sudo python setup.py install'.

During installation, you may see the message "nwmodule.c was not correctly compiled, so the Python implementation will be used instead." It can safely be ignored. This is relating to the C implementation of the Needleman-Wunsch algorithm in the 'aligners' module, and the error is triggered by the import order during setup. If the C implementation was truly not compiled correctly, a related message will be printed when importing the 'aligners' module.
