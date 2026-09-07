#include "Aigg.h"
#include <sstream>
#include <fstream>
#include <ctime>
#include <stdlib.h> // for atoi

#include "CxIo.h" // for GetExePath()
#include "CxOsInt.h" // for GetExePath()

#include "SearchSphereGrid.h"
#include "TargetSpace.h"


ct::FLogStdStream io(std::cout); // <-- main console output goes there.


#ifdef INCLUDE_ABANDONED
// char const
//    *p1ResultFmt = " {:<32}{:18.12f}\n",
//    *p1ResultFmtAnnotated = " {:<32}{:18.12f}  ({})\n",
//    *p1CountFmt = " {:<32}{:5}\n",
//    *p1CountFmtAnnotated =" {:<32}{:5}  ({})\n",
//    *p1InfoFmt = " {:<32}{:>18}\n",
//    *p1InfoFmtAnnotated = " {:<32}{:>18}  ({})\n",
//    *p1TimingFmt = "{:<50}{:10.2f} sec\n",
//    *p1TimingFmtAnnotated = "{:<50}{:10.2f} sec  (n={})\n",
//    *p1InfoExpfFmt = " {:<32}{:18.4e}\n",
//    *p1InfoExpfFmtAnnotated = " {:<32}{:18.4e}  ({})\n";
// // static char const
// //    *p1TimingFmt = "{:<50}{:10.2f} sec\n";
// 
// namespace io_detail {
//    std::ostream &g_OutputStream = std::cout;
//    fmt::MemoryWriter g_IoWriteProxy;
// }
// 
// void IoFlush() {
//    g_OutputStream << io_detail::g_IoWriteProxy.c_str();
//    io_detail::g_IoWriteProxy.clear();
//    g_OutputStream.flush();
// }
//    
// 
// template<class TValue>
// std::string WriteValueT(char const *pFmt, char const *pFmtAnnotated, fmt::BasicStringRef<char> Name, TValue const &Value, fmt::BasicStringRef<char> Annotation) {
//    if (Annotation == "")
//       Write(pFmt, Name, iValue);
//    else
//       Write(pFmtAnnotated, Name, iValue, Annotation);
// }
// 
// void WriteCount(std::string const &Name, ptrdiff_t Value, fmt::BasicStringRef<char> Annotation) {
//    WriteValueT(p1CountFmt, p1CountFmtAnnotated, Name, Value, Annotation);
// }
// 
// void WriteInfo(std::string const &Name, fmt::BasicStringRef<char> Value, fmt::BasicStringRef<char> Annotation) {
//    WriteValueT(p1InfoFmt, p1InfoFmtAnnotated, Name, Value, Annotation);
// }
// 
// void WriteInfoExpf(std::string const &Name, FScalar const &Value, fmt::BasicStringRef<char> Annotation) {
//    WriteValueT(p1InfoExpfFmt, p1InfoExpfFmtAnnotated, Name, Value, Annotation);
// }
// 
// void WriteResult(std::string const &Name, FScalar const &Value, fmt::BasicStringRef<char> Annotation) {
//    WriteValueT(FormatValueT(p1ResultFmt, p1ResultFmtAnnotated, Name, Value, Annotation);
// }
// 
// void WriteTiming(std::string const &Name, FScalar const &TimeInSec, size_t nCount) {
//    if (nCount >= 1)
//       Write(p1TimingFmt, Name, TimeInSec/nCount);
//    else
//       Write(pFmtAnnotated, Name, TimeInSec);
// }
#endif // INCLUDE_ABANDONED
   
#ifdef INCLUDE_ABANDONED
// template<class TValue>
// std::string FormatValueT(char const *pFmt, char const *pFmtAnnotated, fmt::BasicStringRef<char> Name, TValue const &Value, fmt::BasicStringRef<char> Annotation) {
//    if (Annotation == "")
//       return fmt::format(pFmt, Name, iValue);
//    else
//       return fmt::format(pFmtAnnotated, Name, iValue, Annotation);
// }
// 
// void WriteCount(std::string const &Name, ptrdiff_t Value, fmt::BasicStringRef<char> Annotation) {
//    std::cout << FormatValueT(p1CountFmt, p1CountFmtAnnotated, Name, Value, Annotation); std::cout.flush();
// }
// 
// void WriteInfo(std::string const &Name, fmt::BasicStringRef<char> Value, fmt::BasicStringRef<char> Annotation) {
//    std::cout << FormatValueT(p1InfoFmt, p1InfoFmtAnnotated, Name, Value, Annotation); std::cout.flush();
// }
// 
// void WriteInfoExpf(std::string const &Name, FScalar const &Value, fmt::BasicStringRef<char> Annotation) {
//    std::cout << FormatValueT(p1InfoExpfFmt, p1InfoExpfFmtAnnotated, Name, Value, Annotation); std::cout.flush();
// }
// 
// void WriteResult(std::string const &Name, FScalar const &Value, fmt::BasicStringRef<char> Annotation) {
//    std::cout << FormatValueT(p1ResultFmt, p1ResultFmtAnnotated, Name, Value, Annotation); std::cout.flush();
// }
// 
// void WriteTiming(std::string const &Name, FScalar const &TimeInSec, size_t nCount) {
//    if (nCount >= 1)
//       return fmt::format(p1TimingFmt, Name, TimeInSec/nCount);
//    else
//       return fmt::format(pFmtAnnotated, Name, TimeInSec);
// }
#endif // INCLUDE_ABANDONED

#ifdef INCLUDE_ABANDONED
// std::string FormatCount(fmt::BasicStringRef<char> Name, ptrdiff_t iValue, fmt::BasicStringRef<char> Annotation) {
//    if (Annotation == "")
//       return fmt::format(p1CountFmt, Name, iValue);
//    else
//       return fmt::format(p1CountFmtAnnotated, Name, iValue, Annotation);
// }
// 
// std::string FormatInfo(fmt::BasicStringRef<char> Name, fmt::BasicStringRef<char> Value, fmt::BasicStringRef<char> Annotation) {
//    if (Annotation == "")
//       return fmt::format(p1InfoFmt, Name, Value);
//    else
//       return fmt::format(p1InfoFmtAnnotated, Name, Value, Annotation);
// }
// 
// std::string FormatInfoExpf(fmt::BasicStringRef<char> Name, FScalar const &Value, fmt::BasicStringRef<char> Annotation) {
//    if (Annotation == "")
//       return fmt::format(p1InfoExpfFmt, Name, Value);
//    else
//       return fmt::format(p1InfoExpfFmtAnnotated, Name, Value, Annotation);
// }
// 
// void WriteCount(std::string const &Name, ptrdiff_t iValue, fmt::BasicStringRef<char> Annotation) {
//    std::cout << FormatCount(Name, iValue, Annotation); std::cout.flush();
// }
// 
// void WriteInfo(std::string const &Name, fmt::BasicStringRef<char> Value, fmt::BasicStringRef<char> Annotation) {
//    std::cout << FormatCount(Name, iValue, Annotation); std::cout.flush();
// }
// 
// void WriteInfoExpf(std::string const &Name, FScalar const &Value, fmt::BasicStringRef<char> Annotation) {
//    std::cout << FormatInfoExpf(Name, Value, Annotation); std::cout.flush();
// }
// 
#endif // INCLUDE_ABANDONED



// void CxAssertFail(char const *pExpr, char const *pFile, int iLine)
// {
//    throw std::runtime_error(fmt::format("assertion failed:\n{}:{}: {}", pFile, iLine, pExpr));
// }



int main_direct_invoke(int argc, char *argv[])
{
   if (!((argc == 2 && argv[1][0] == '{') || 
         (argc == 3 && (strcmp(argv[1], "--config") == 0 || strcmp(argv[1], "-c") == 0)) ||
         (argc == 3 && (strcmp(argv[1], "--preprocess") == 0 || strcmp(argv[1], "-p") == 0))
        )) {
      std::cout <<
         "\n :::: aigg [v@Pr0gVeR] -- (A)ngular (I)ntegration (G)rid (G)enerator ::::"
         "\n"
         "\n      ...The Knizia Group's program for computing optimized integration grids"
         "\n         for functions defined on the surface of the 3D unit sphere (S_2)."
         "\n "
         "\n                           Authors: Gerald Knizia, Mieke Peels, 2016--2021"
         "\n "
         "\n Invocation:"
         "\n "
         "\n   aigg '{<options>}'"
         "\n "
         "\n or"
         "\n "
         "\n   aigg --config <input-filename>"
         "\n   aigg -c <input-filename>"
         "\n "
         "\n or"
         "\n "
         "\n   aigg --preprocess '{point-group: ..., order: ...}'"
         "\n   aigg -p '{point-group: ..., order: ...}'"
         "\n "
         "\n For usage examples, and possible options, see example input files in"
         "\n the inputs/ directory."
         "\n "
         "\n Alternatively to input files, options can be given on the command line"
         "\n directly (in the same format, including the curly braces). For example,"
         "\n "
         "\n   aigg '{point-group: icosahedral; degree: [4,step:+1,-1]; initial-points: hex-grid{5;1;01c}}'"
         "\n "
         "\n is an equivalent of inputs/1-incremental-search-up.ini. Note that the"
         "\n outer curly braces { } are REQUIRED and might need to be enclosed in"
         "\n single quotes (as above) in order to be passed to the program and not"
         "\n be interpreted by the shell (e.g., bash)."
         "\n "
         "\n See input examples or manual for an explanation of the options and"
         "\n usage of this program."
         "\n "
         "\n The --preprocess options may be used to generate data files (stored in"
         "\n <aigg-exe-path>/data) which encode minimal representations of the target"
         "\n function space for a given point group and a given target integration order;"
         "\n pre-generating optimal function sets is not required, but can speed up and"
         "\n improve subsequent regular grid optimizations. It should only be done with"
         "\n aigg compiled for high floating point accuracy (at least 256bits)."
         << std::endl;
      return 0;
   }
   if (1) {
      std::cout << "-- Configuration\n\n";
//       char const *pInfoFmt = "  {:<25} {}\n";
//       char const *pInfoFmt = " {:<32}{:}\n";
      char const *pInfoFmt = " {:<25}{:}\n";
      std::cout << fmt::format(pInfoFmt, "Version:", MakeAiggVersionString());
      std::cout << fmt::format(pInfoFmt, "Scalar backend driver:", g_pScalarImplDesc);
      std::cout << fmt::format(pInfoFmt, "Scalar size in bits:", g_ScalarFloatSizeBits);
      double fDecDigits = double(log(FScalar(2))/log(FScalar(10))) * g_ScalarFloatMantissaBits;
      std::cout << fmt::format(pInfoFmt, "Scalar precision:", fmt::format("{} bits (~{:.2f} decimal digits)", g_ScalarFloatMantissaBits, fDecDigits));
      std::cout << std::endl;
   }
   std::cout << "-- Main process\n" << std::endl;
//    return 0;

   
   if (argc == 3 && (strcmp(argv[1], "--preprocess") == 0 || strcmp(argv[1], "-p") == 0)) {
      FTargetSpacePreprocessOptionsPtr
         pOptions;
      std::string
         sOptions = std::string(argv[2]);
      pOptions = FTargetSpacePreprocessOptionsPtr(new FTargetSpacePreprocessOptions(sOptions));
      
      // pass control to preprocess driver routine
      RunTargetSpacePreprocess(*pOptions);
   } else {
      // read input options and task description
      FGridSearchOptionsPtr
         pOptions;
      std::string
         sOptions;
      if (argc == 2) {
         // read options from command line directly. They still need to come as a formatted property list, however...
         sOptions = std::string(argv[1]);
      } else {
         assert(argc == 3 && (strcmp(argv[1], "--config") == 0 || strcmp(argv[1], "-c") == 0));
         // options come from a input file. attempt to read it... and convert it to a str.
         char const
            *pFileName = argv[2];
         std::cout << fmt::format("! Reading configuration from file '{}'", pFileName) << std::endl;
         std::ifstream
            inp(pFileName);
         std::string
            line;
         while (std::getline(inp, line)) {
            if (!line.empty() && line[0] == '#')
               continue;
            sOptions += line; // totally elegant...
         }
      }
      pOptions = new FGridSearchOptions(sOptions);

      // pass control to main driver routine
      RunGridSearchSchedule(*pOptions);
   }
   return 0;
}

std::string MakeAiggVersionString()
{
//    return fmt::format("kgaigg v@Pr0gVeR (compiled: " __DATE__ ", fprec: {})", g_ScalarFloatSizeBits);
   return fmt::format("aigg v@Pr0gVeR (compiled: " __DATE__ ", fprec: {})", g_ScalarFloatSizeBits);
}


std::string MakeCurrentTimeString()
{
   std::time_t t = time(0);
   struct std::tm *now = std::gmtime(&t);
   if (now) {
      return fmt::format("{:04}-{:02}-{:02} {:02}:{:02} (UTC)", now->tm_year+1900, now->tm_mon+1, now->tm_mday, now->tm_hour, now->tm_min);
   } else {
      return "unknown-time";
   }
}


// concatenate a string `pattern` together for `count` times (like Python `count * pattern` would do);
// if given, append `left` to the left and `right` to the right.
std::string Repeat(size_t count, std::string const &pattern, std::string const &left, std::string const right)
{
   std::string out;
   out.reserve(left.size() + pattern.size() * count + right.size());
   out.append(left);
   for (size_t ipat = 0; ipat < count; ++ ipat)
      out.append(pattern);
   out.append(right);
   return out;
}


int SymFnTest();


int main(int argc, char *argv[])
{
   // src/Aigg.cpp:303:23: warning: 'void Eigen::initParallel()' is deprecated: Initialization is no longer needed. [-Wdeprecated-declarations]
   // Eigen::initParallel();
   ct::FFileLocator::SetBasePath(ct::GetExePath());
   
//    return LeakTest();
   if (0 != SymFnTest())
      return -1;
   
   bool RunUnderExceptionHandler = true;
#ifdef _DEBUG
   RunUnderExceptionHandler = false;
#endif // _DEBUG
   if (RunUnderExceptionHandler) {
      try {
         return main_direct_invoke(argc, argv);
      } catch (std::exception &e) {
         std::cerr << fmt::format("\n\n** SPHERE GRID GENERATOR CRASHED. Uncaught exception:\n{}", e.what()) << std::endl;
      }
   } else {
      return main_direct_invoke(argc, argv);
   }
   return -1;
}
