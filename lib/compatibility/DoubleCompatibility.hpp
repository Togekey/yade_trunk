// 2019 Janek Kozicki
// This file is to maintain compatibility with external libraries that cannot take Real as an argument: VTK and Coinor.
// The macros perform conversion of arguments from Real to double.
// This file should be kept very small. If more than 50 lines must be added it means that there is some design issue.

#pragma once

// These macros have such names on purpose: they create a comma separated list of static_casted arguments.
// When Real==double they evaluate to nothing and are removed during compilation.
#define ONE_DOUBLE(val) static_assert(false, "Better use static_cast<double>(val), because it should be very rare, and then it is obvious what's going on.");
#define TWO_DOUBLES(val1, val2) (static_cast<double>(val1)), (static_cast<double>(val2))
#define THREE_DOUBLES(val1, val2, val3) (static_cast<double>(val1)), (static_cast<double>(val2)), (static_cast<double>(val3))
#define FOUR_DOUBLES(val1, val2, val3, val4) (static_cast<double>(val1)), (static_cast<double>(val2)), (static_cast<double>(val3)), (static_cast<double>(val4))

// These macros create on ald C-array of doubles.
#define ARRAY_2_DOUBLE(val1, val2)                                                                                                                             \
	{                                                                                                                                                      \
		TWO_DOUBLES(val1, val2)                                                                                                                        \
	}

#define ARRAY_3_DOUBLE(val1, val2, val3)                                                                                                                       \
	{                                                                                                                                                      \
		THREE_DOUBLES(val1, val2, val3)                                                                                                                \
	}

#define ARRAY_4_DOUBLE(val1, val2, val3, val4)                                                                                                                 \
	{                                                                                                                                                      \
		FOUR_DOUBLES(val1, val2, val3, val4)                                                                                                           \
	}

#define ARRAY_9_DOUBLE(val1, val2, val3, val4, val5, val6, val7, val8, val9)                                                                                   \
	{                                                                                                                                                      \
		(static_cast<double>(val1)), (static_cast<double>(val2)), (static_cast<double>(val3)), (static_cast<double>(val4)),                            \
		        (static_cast<double>(val5)), (static_cast<double>(val6)), (static_cast<double>(val7)), (static_cast<double>(val8)),                    \
		        (static_cast<double>(val9)),                                                                                                           \
	}

#define VEC3_TO_ARRAY_DOUBLE(VectorName, ArrayName)                                                                                                            \
	static_assert(sizeof(VectorName) == sizeof(::yade::Vector3r));                                                                                         \
	const double ArrayName[3] = ARRAY_3_DOUBLE(VectorName[0], VectorName[1], VectorName[2])

#define VEC4_TO_ARRAY_DOUBLE(VectorName, ArrayName)                                                                                                            \
	static_assert(sizeof(VectorName) == sizeof(::yade::Vector4r));                                                                                         \
	const double ArrayName[4] = ARRAY_4_DOUBLE(VectorName[0], VectorName[1], VectorName[2], VectorName[3])

