#include "Aigg.h"
#include "ScalarTypes.h"
#include <sstream>

// next: setup description data and default thresholds for and test functions
// for configured floating point data type.

#ifdef SCALAR_IS_FLOAT_128
   unsigned const
      g_ScalarFloatSizeBits = 128,
      g_ScalarFloatMantissaBits = 113;
   #ifdef FLOAT_128_NATIVE_EXT
   char const
      *g_pScalarImplDesc = "boost::multiprecision on __float128";
   #else
   char const
      *g_pScalarImplDesc = "boost::multiprecision::cpp_bin_float";
   #endif

   FScalar const
      g_ThrVerySmall = FLOAT_LITERAL(1e-25),
      g_ThrAlmostZero = FLOAT_LITERAL(1e-20);

#elif defined(SCALAR_IS_FLOAT_256)
   unsigned const
      g_ScalarFloatSizeBits = 256,
//       g_ScalarFloatMantissaBits = 241;
      g_ScalarFloatMantissaBits = 237;
   char const
      *g_pScalarImplDesc = "boost::multiprecision::cpp_bin_float";

   FScalar const
      g_ThrVerySmall = FLOAT_LITERAL(1e-50),
      g_ThrAlmostZero = FLOAT_LITERAL(1e-40);

#elif defined(SCALAR_IS_FLOAT_512)
   unsigned const
      g_ScalarFloatSizeBits = 512,
//       g_ScalarFloatMantissaBits = 497;
      g_ScalarFloatMantissaBits = 493;
   char const
      *g_pScalarImplDesc = "boost::multiprecision::cpp_bin_float";

   FScalar const
      g_ThrVerySmall = FLOAT_LITERAL(1e-100),
      g_ThrAlmostZero = FLOAT_LITERAL(1e-80);

#elif defined(SCALAR_IS_FLOAT_64)
   unsigned const
      g_ScalarFloatSizeBits = 64,
      g_ScalarFloatMantissaBits = 53;
   char const
      *g_pScalarImplDesc = "native double";
   FScalar const
      g_ThrVerySmall = 1e-12,
      g_ThrAlmostZero = 1e-10;
#else
   #error "Configuration error: None of the SCALAR_IS_FLOAT_xx flags was set in DefineScalarTypes.cpp"
#endif // SCALAR_IS_FLOAT_xx


#ifdef SCALAR_IS_FLOAT_64
   bool IsAlmostZero(FScalar const &f, FScalar const &Thr) {
      return std::abs(f) < Thr;
   }
#else
   bool IsAlmostZero(FScalar const &f, FScalar const &Thr) {
//       return abs(f) < 1e-24Q;
      return abs(f) < Thr;
   }
#endif // SCALAR_IS_FLOAT_xx


FScalar const
//    X_PI(FLOAT_LITERAL(3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825));
   // ^- not enough for the 512bit floats... they have ~150 decimal digits.
   // This one comes from N[Pi,200]:
   X_PI(FLOAT_LITERAL(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303820)),
   // ...should probably use boost::math::constants<ScalarT> instead... especially since
   // we anyway assume boost to be around, don't we?
   // UPDATE: """All the constants are accurate to at least the 34 decimal digits required for 128-bit long doubles, and most are accurate to 100 digits or more..."""
   // so that won't be enough :(.
   // ...anyway...here's N[E,200]:
   X_EULER(FLOAT_LITERAL(2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952605956307381323286279434907632338298807531952510190));


std::string FmtExportFloat(FScalar f, int Precision)
{
   if (!(Precision == -1 || Precision > 0))
      throw std::runtime_error(fmt::format("FmtExportFloat({},{}): decimal precision specifier must be -1 or positive", double(f), Precision));
   unsigned nDigits10 = std::numeric_limits<FScalar>::max_digits10;
   // clamp precision to max number of decimal digits which are meaningful
   // for current scalar type
   if (Precision > 0)
      Precision = std::min(unsigned(Precision), nDigits10);

#ifdef SCALAR_IS_FLOAT_64
   if (Precision == -1) {
      return fmt::format("{:24.16e}", f);
   } else {
      std::string sFmt = fmt::format("{{:{}.{}e}}", 8+Precision, Precision);
      return fmt::format(sFmt, f);
   }
#else
   std::stringstream ss;
   // Notes:
   // - The boost::multiprecision types can be streamed from and to std streams says:
   //    https://www.boost.org/doc/libs/1_58_0/libs/multiprecision/doc/html/boost_multiprecision/tut/conversions.html
   // - We could use std::numeric_limits<FScalar>::max_digits10 to export with full precision.
   //   ...but I guess 24 decimal digits should be enough, shouldn't they?
   //   We only make double precision grids, after all.
   if (Precision == -1) {
      if (0) {
         ss.precision(24);
         ss.width(32);
      } else {
         nDigits10 = std::min(nDigits10, unsigned(32));
   //       nDigits10 = std::min(nDigits10, unsigned(64));
         ss.precision(nDigits10);
         ss.width(10+nDigits10);
      }
   } else {
      ss.precision(Precision);
      ss.width(10+Precision);
   }
   ss << std::scientific;
   ss << std::right;
   ss << f;
//    {:>32}
   // TODO: apparently there is a .str() method on boost::multiprecision
   // numbers. Would support things like `val.str(10,
   // std::ios_base::scientific)` or `val.str(10, std::ios_base::fixed)` or
   // similar.
   return ss.str();
#endif // SCALAR_IS_FLOAT_64
}

FScalar ParseFloat(std::string const &s)
{
   FScalar
      r(0);
   bool
      Failed = false;
#ifdef SCALAR_IS_FLOAT_64
   // c/p'd from CxParse1.cpp
   char const
      *pEnd = 0;
   r = strtod(&*s.begin(), const_cast<char**>(&pEnd));
   if (pEnd != &*(s.end()-1) + 1)
      Failed = true;
#else
   std::stringstream
      ss(s);
   ss >> r;
   if (!ss.eof() || ss.bad())
      Failed = true;
//    if (Failed) {
//       std::cout << fmt::format("ParseFloat(): failed to parse '{}' as a floating point number. r = {}  ss.eof = {}  ss.bad = {}", s, double(r), ss.eof(), ss.bad()) << std::endl;
//    }
#endif // SCALAR_IS_FLOAT_64
   if (Failed) {
      throw std::runtime_error(fmt::format("ParseFloat(): failed to parse '{}' as a floating point number", s));
   }
   return r;
}
