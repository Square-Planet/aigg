#ifndef DEFINE_SCALAR_TYPES_H
#define DEFINE_SCALAR_TYPES_H

#include <string> // for FmtExportFloat

#ifdef SCALAR_IS_FLOAT_128
   #ifdef FLOAT_128_NATIVE_EXT
      // use __float128 as backend (compiler extension; e.g., in g++).
      #include <boost/multiprecision/float128.hpp>
      // bring in the math function overloads from boost::multiprecision into the
      // global namespace. We need that because depending on the environment, we
      // may have different namespaces in which the math functions reside (and I
      // have seen namespace aliases create all kinds of havoc on older compilers)
      using namespace boost::multiprecision;
      typedef float128
         FScalar;
      // ^- note: on gcc, this uses gcc's __float128 internally.
      // - This makes very large debug symbols, but apart from that the actual wrapper
      //   overhead is not large, and there is little point in using the __float128
      //   directly (I tried).
      // - this one has compiler-side named literals ending in "Q" (as in "f" for 32bit floats)
      // - unlike for the cpp_bin_float types below, the actual sizeof(T) of
      //   these wrapped float128 objects on g++ is indeed 16 bytes (i.e., without
      //   size overhead)
      #define FLOAT_LITERAL(z) z##Q
   #else
      #include <boost/multiprecision/cpp_bin_float.hpp>
      // bring in the math function overloads from boost::multiprecision into the
      // global namespace. We need that because depending on the environment, we
      // may have different namespaces in which the math functions reside (and I
      // have seen namespace aliases create all kinds of havoc on older compilers)
      using namespace boost::multiprecision;
      // Notes:
      // - For IEEE-float-like typedefs, see:
      //
      //   https://www.boost.org/doc/libs/1_66_0/libs/multiprecision/doc/html/boost_multiprecision/ref/cpp_bin_float_ref.html
      //
      //   for example, this has:
      //
      //   typedef number<backends::cpp_bin_float<113, backends::digit_base_2, void, boost::int16_t, -16382, 16383>, et_off>  cpp_bin_float_quad;
      //
      //   for the 128bit float type used here.
      //
      // - 113 mantiassa bits are good for approx np.log10(2.**113) = 34.02 decimal digits
      // - This coincides with the IEEE-754 definition described here: https://en.wikipedia.org/wiki/Quadruple-precision_floating-point_format:
      //   Sign bit: 1 bit; Exponent width: 15 bits; Significand precision: 113 bits (112 explicitly stored)
      //
      // - WARNING: due to internal storage format issues and alignment, the actual size
      //   (sizeof(T)) of a boost::multiprecision::cpp_bin_float in memory will generally be
      //   larger (often by up to factor of 2) than the formally expected size based on the number
      //   of bits. Also, https://www.boost.org/doc/libs/1_59_0/libs/multiprecision/doc/html/boost_multiprecision/tut/limits/constants.html indicates
      //   that the cpp_bin_float implementation does not actually use implicit mantissa bits (didn't check
      //   in detail. First scan showed conflicting info).
      // - Details on type sizes: see ~/dev/sphere_grids/test-scalars>
      typedef cpp_bin_float_quad
         FScalar;
      // this one has no named literal suffixes (I think?), but the objects can
      // be constructed from quoted strings.
      #define FLOAT_LITERAL(z) FScalar(#z)
   #endif
#elif defined(SCALAR_IS_FLOAT_512)
      #include <boost/multiprecision/cpp_bin_float.hpp>
      using namespace boost::multiprecision;
      // Notes:
      // - There does not appear to be an official IEEE binary512-float type definition
      //   (but, admittedly, I did not look very hard)
      // - As in the binary64 (=double) vs binary128 (=quad) definition, I do
      //   not extend the representable exponents beyond the binary256 float
      //   range here (don't need it)
      // - I get the number of mantissa bits as 237 + (512-256)
      //   (237 is what the still-defined binary256 type has, I also copied the exponent range from that)
      // - 493 mantiassa bits are good for approx np.log10(2.**493) = 148.41 decimal digits
      typedef number<backends::cpp_bin_float<493, backends::digit_base_2, void, boost::int32_t, -262142, 262143>, et_off>
         FScalar;
#ifdef INCLUDE_ABANDONED
//       // Notes:
//       // - one mantissa bit is implicit (which is why with 497 mantiassa bits we
//       //   request still have enough left for an int16_t)
//       // - As in the cpp_bin_float_quad typedef, I do not extend the representable
//       //   exponents beyond 80bit float range here (don't need it)
//       // - I get the number of mantissa bits as 113 + (512-128)
//       // - 497 mantiassa bits are good for approx np.log10(2.**497) = 149.61 decimal digits
//       typedef number<backends::cpp_bin_float<497, backends::digit_base_2, void, boost::int16_t, -16382, 16383>, et_off>
//          FScalar;
      // UPDATE:
      // - https://www.boost.org/doc/libs/1_59_0/libs/multiprecision/doc/html/boost_multiprecision/tut/limits/constants.html says:
      //   "The Boost.Multiprecision binary types do not use an implicit bit, so
      //    the digits member reflects exactly how many bits of precision were requested: ..."
      // - cpp_bin_float.hpp itself defines:
      //   typedef number<backends::cpp_bin_float<113, backends::digit_base_2, void, boost::int16_t, -16382, 16383>, et_off> cpp_bin_float_quad;
      //   typedef number<backends::cpp_bin_float<237, backends::digit_base_2, void, boost::int32_t, -262142, 262143>, et_off> cpp_bin_float_oct;
      //   However, I do not quite see how this works with 128bits of storage, for example. If that is really what it is. Should check.
#endif // INCLUDE_ABANDONED
      #define FLOAT_LITERAL(z) FScalar(#z)
#elif defined(SCALAR_IS_FLOAT_256)
      #include <boost/multiprecision/cpp_bin_float.hpp>
      using namespace boost::multiprecision;
      // https://en.wikipedia.org/wiki/Octuple-precision_floating-point_format says:
      // Sign bit: 1 bit; Exponent width: 19 bits; Significand precision: 237 bits (236 explicitly stored)
      // - 237 mantiassa bits are good for approx np.log10(2.**237) = 71.34 decimal digits
      // - This coincides with what cpp_bin_float.hpp itself defines in newer versions:
      //   typedef number<backends::cpp_bin_float<113, backends::digit_base_2, void, boost::int16_t, -16382, 16383>, et_off> cpp_bin_float_quad;
      //   typedef number<backends::cpp_bin_float<237, backends::digit_base_2, void, boost::int32_t, -262142, 262143>, et_off> cpp_bin_float_oct;
      // - For compatibility reasons we therefore switch to the 237 bit variant with larger exponent range.
      typedef number<backends::cpp_bin_float<237, backends::digit_base_2, void, boost::int32_t, -262142, 262143>, et_off>
         FScalar;
//       typedef number<backends::cpp_bin_float<241, backends::digit_base_2, void, boost::int16_t, -16382, 16383>, et_off>
//          FScalar;
      #define FLOAT_LITERAL(z) FScalar(#z)
#else
   #ifndef SCALAR_IS_FLOAT_64
      #define SCALAR_IS_FLOAT_64
   #endif
   typedef double
      FScalar;
   
   #define FLOAT_LITERAL(z) z
#endif

extern unsigned const
   // bit size of floating point data type behind FScalar
   // (e.g., 64 for the typical IEEE native double)
   g_ScalarFloatSizeBits,
   // number of bits used for the mantissa of the floating point numbers in
   // FScalar. This defines the representable relative precision of the scalars.
   g_ScalarFloatMantissaBits;
extern char const
   // textual description of the scalar data type implementation
   *g_pScalarImplDesc;

extern FScalar const
   // pi in sufficient precision for the given data type.
   X_PI,
   // Euler number (e = Exp[1]) in sufficient precision for the given data type.
   X_EULER;

extern FScalar const
   g_ThrVerySmall, // default ~1e-25 for 128bit floats
   g_ThrAlmostZero;  // default ~1e-20 for 128bit floats
bool IsAlmostZero(FScalar const &f, FScalar const &Thr=g_ThrAlmostZero);
inline bool IsAlmostEqual(FScalar const &f, FScalar const &g, FScalar const &Thr=g_ThrAlmostZero) { return IsAlmostZero(f - g, Thr); }


// format a scalar (double or float128) as string, for sake of final data
// export, with extended precision. Bypass around hacking float128 support into
// cppformat.
// Precision: max number of decimal digits to export. -1 means "auto".
std::string FmtExportFloat(FScalar f, int Precision=-1);

FScalar ParseFloat(std::string const &s);


#endif // DEFINE_SCALAR_TYPES_H
