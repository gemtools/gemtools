#!/usr/bin/env python
import sys
import os

# hack the python path
__gt_cwd__dir = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]
__gt_py_version = sys.version_info
sys.path.insert(0, "%s/lib64/python%d.%d/site-packages" % (__gt_cwd__dir, __gt_py_version[0], __gt_py_version[1]))

__requires__ = 'Gemtools'
import sys
from pkg_resources import load_entry_point

sys.exit(load_entry_point('Gemtools', 'console_scripts', 'gemtools')())
