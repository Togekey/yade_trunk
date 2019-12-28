// 2019 Janek Kozicki
// This file is to maintain compatibility with external libraries that cannot take Real as an argument: VTK and Coinor.
// The macros perform conversion of arguments from Real to double.
// This file should be kept very small. If more than 50 lines must be added it means that there is some design issue.

#pragma once

#define ARRAY_2_DOUBLE(val1, val2)                                                                                                                             \
	{                                                                                                                                                      \
		(static_cast<double>(val1)), (static_cast<double>(val2))                                                                                       \
	}

#define ARRAY_3_DOUBLE(val1, val2, val3)                                                                                                                       \
	{                                                                                                                                                      \
		(static_cast<double>(val1)), (static_cast<double>(val2)), (static_cast<double>(val3))                                                          \
	}

#define ARRAY_4_DOUBLE(val1, val2, val3, val4)                                                                                                                 \
	{                                                                                                                                                      \
		(static_cast<double>(val1)), (static_cast<double>(val2)), (static_cast<double>(val3)), (static_cast<double>(val4))                             \
	}

#define ARRAY_9_DOUBLE(val1, val2, val3, val4, val5, val6, val7, val8, val9)                                                                                   \
	{                                                                                                                                                      \
		(static_cast<double>(val1)), (static_cast<double>(val2)), (static_cast<double>(val3)), (static_cast<double>(val4)),                            \
		        (static_cast<double>(val5)), (static_cast<double>(val6)), (static_cast<double>(val7)), (static_cast<double>(val8)),                    \
		        (static_cast<double>(val9)),                                                                                                           \
	}

