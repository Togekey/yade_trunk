##########################################################################
#  2019        Janek Kozicki                                             #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################

"""
When yade uses high-precision number as ``Real`` type the usual (old) ``import minieigen``
has to be replaced with ``import yade.minieigenHP`` library which uses exacly the same number of decimal places
as yade is using everywhere else. Please note that such precision can be very arbitrary, because ``cpp_bin_float``
or ``mpfr`` take it as a compile-time argument. Hence such ``minieigenHP`` cannot be separately precompiled as a package.
Though it could be precompiled for some special types such as ``boost::multiprecision::float128``.
"""

# all C++ functions are accessible now:
import yade

if(yade.config.highPrecisionBits != 64):
	from yade._minieigenHP import *
else:
	from minieigen import *

