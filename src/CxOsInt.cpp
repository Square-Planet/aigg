#ifndef _WIN32
    #include <unistd.h> // for getpid and pid_t. Anyone got a better idea to obtain the base path name?
#else
    #include <stdlib.h> // for malloc/free
    #define NOMINMAX
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#endif

#include "CxDefs.h"
#include "CxOsInt.h"
#include "format.h"
#include <stdexcept>
#include <cstdlib> // for wcstombs
#include <cstring> // for len()
#include <algorithm> // for std::min() and std::reverse
#include <sstream> // for path split/join stuff
#include <fstream> // for ifstream in file_readable_q
// #include <list> // for split_path

namespace ct {

#ifdef _WIN32

// #ifdef UNICODE
//    #error "Sorry, this stuff assumes 8byte chars. Please compile with ASCII API bindings."
// #endif


std::string FmtLastError()
{
   // Retrieve the system error message for the last-error code
   LPTSTR lpMsgBuf;
   DWORD dw = GetLastError();

   FormatMessage(
      FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
      NULL,
      dw,
      MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
      (LPTSTR) &lpMsgBuf,
      0, NULL );

   // convert to std::string.
   std::string
      s(lpMsgBuf, lpMsgBuf + lstrlen(lpMsgBuf));
   LocalFree(lpMsgBuf);
   return s;
}

std::string GetExePath()
{
   // _CxAssertRt(sizeof(TCHAR) == sizeof(char)); // this may become a serious problem on non-european systems...
   // return std::string("C:\\Linux\\dev\\microscf.20180116");
   DWORD
      nBuf = 512,
      nChars;
   TCHAR
      *pwBuf = (TCHAR*)::malloc(sizeof(TCHAR) * nBuf);
   char
      *pcBuf = (char*)::malloc(sizeof(char) * nBuf);
   nChars = GetModuleFileName(0, pwBuf, nBuf);

   if (nChars == 0)
      throw std::runtime_error("Failed to obtain current path name. GetLastError() says: " + FmtLastError() );

   if (sizeof(TCHAR) == sizeof(wchar_t)) {
      std::wcstombs(pcBuf, (wchar_t*)pwBuf, std::min(size_t(nBuf), size_t(1 + std::wcslen((wchar_t*)pwBuf))));
   } else if (sizeof(TCHAR) == sizeof(char)) {
      std::memcpy(pcBuf, pwBuf, sizeof(char) * std::min(size_t(nBuf), size_t(1 + std::strlen((char*)pwBuf))));
   } else {
      throw std::runtime_error("GetExePath(): broken character format.");
   }

   // convert to std::string.
   std::string
      s(pcBuf, pcBuf + nChars);

   ::free(pcBuf);
   ::free(pwBuf);

   return s;
}


#else

// return the absolute path of the currently executed binary on
// a linux machine. Required to access our compiled algorithm files,
// because working directory is usually changed to some place else.
// if anyone got a better idea on how to do it ... [what follows now
// appears to be not some kind of obscure hack, but the actual linux
// standard way of doing this]
std::string GetExePath()
{
    pid_t
        ProcessId = getpid();
    std::string
        ExeSymLink(fmt::format("/proc/{}/exe", ProcessId)),
        ExeName(513,'\0');
    readlink(ExeSymLink.c_str(), &ExeName[0], ExeName.size() - 1);
    int
        iEnd = ExeName.find('\0');
    while ( iEnd > 0 && ExeName[iEnd] != '/' && ExeName[iEnd] != '\\' && ExeName[iEnd] != ':' )
        -- iEnd;
    ++ iEnd;
    return ExeName.substr(0, iEnd);
}

#endif

} // namespace ct



namespace ct {
   namespace path {
      // returns 0 if c is not one of the dir separators, otherwise returns the dir separator
      char dir_seperator_q(char c) {
         for (size_t iDirSep = 0; iDirSep < CX_PATH_NUM_SEPARATORS; ++ iDirSep) {
            char cDirSepToTest = CX_PATH_ALL_SEPARATORS[iDirSep];
            _CxAssert(cDirSepToTest != 0); // 0 is then obviously not supported as directory separator...
            if (c == cDirSepToTest)
               return c;
         }
         return 0;
      }


      size_t _iBaseNameStartQ(FFileName const &FileName) {
         size_t n = FileName.size();
         size_t nBase = 0;
         for ( ; (nBase < n) && (dir_seperator_q(FileName[n - nBase - 1]) == 0); ++ nBase) {
         }
         size_t nPath = n - nBase;
//          fmt::print(" {:20} n = {:<3} np = {:<3} nb = {:<3}\n\n", "'"+FileName+"'", n, nPath, nBase);
         return nPath;
      }


      FFileName_PathAndBase split(FFileName const &FileName) {
         size_t nPath = _iBaseNameStartQ(FileName);
         return FFileName_PathAndBase(FileName.substr(0, nPath), FileName.substr(nPath, FileName.size() - nPath));
      }

      FFileName_StemAndExt split_ext(FFileName const &FileName) {
         // extract section of FileName corresponding to the base name of `FileName`
         size_t iBase = _iBaseNameStartQ(FileName);

         // ignore leading periods in the base name (just as Python's os.path.splitext).
         // This takes care of file names like ".", "..", or ".bashrc", which
         // should be returned unchanged (the latter indicates a "hidden" file
         // in unix-like/linux-like operating systems)
         size_t iExt0 = iBase;
         while (iExt0 < FileName.size() && FileName[iExt0] == '.')
            iExt0 += 1;

         // now scan for periods (from right to left, but not further than iExt0)
         for (size_t ii = size_t(FileName.size()) - size_t(1); ii < FileName.size() && ii >= iExt0; -- ii) {
            if (FileName[ii] == '.')
               // found the right-most period >= iExt0 at ii.
               return FFileName_StemAndExt(FileName.substr(0,ii), FileName.substr(ii));
            // we should not reach a path separator before an extension
            // separator, because for ii >= iExt0 there should not be any path
            // separators (it wouldn't be part of the base name!)
            _CxAssert(dir_seperator_q(FileName[ii]) == 0);
         }
         return FFileName_StemAndExt(FileName, std::string());
      }

#ifdef INCLUDE_ABANDONED
//       FFileName_StemAndExt split_ext(FFileName const &FileName) {
//          size_t n = FileName.size();
//          size_t nBase = 0;
//          // extract section corresponding to base name of `FileName`
//
//          // ignore leading periods (just as Python's os.path.splitext)
//          for ( ; nBase < n; ++ nBase) {
//             size_t ii = n - nBase - 1;
// //             fmt::print(" [{:4}]  n = {:<4}  nb = {:<4}  c = '{}'   '{}'\n", n - nBase - 1, n, nBase, FileName[n - nBase - 1], FileName);
//             if (FileName[ii] == '.')
//                // ^-- not sure what to do with "file" names like "." and "..". I.e., what the right result would even be.
//                return FFileName_StemAndExt(FileName.substr(0,ii), FileName.substr(ii));
//
//             if (dir_seperator_q(FileName[ii]) != 0)
//                // reached path separator before an extension separator -> no extension
//                break;
//          }
//          return FFileName_StemAndExt(FileName, std::string());
//       }
//          for ( ; nBase < n; ++ nBase) {
// //             fmt::print(" [{:4}]  n = {:<4}  nb = {:<4}  c = '{}'   '{}'\n", n - nBase - 1, n, nBase, FileName[n - nBase - 1], FileName);
//             if (dir_seperator_q(FileName[n - nBase - 1]) != 0)
//                break;
//          }
//       FFileName_PathAndBase split(FFileName const &FileName) {
//          size_t n = FileName.size();
//          size_t nBase = 0;
//          char cDirSep = 0; // set from 0 to one of CX_PATH_ALL_SEPARATORS if one is found.
//          for ( ; nBase < n && cDirSep == 0; ++ nBase) {
// //             fmt::print(" [{:4}]  n = {:<4}  nb = {:<4}  c = '{}'   '{}'\n", n - nBase - 1, n, nBase, FileName[n - nBase - 1], FileName);
//             for (size_t iDirSep = 0; iDirSep < CX_PATH_NUM_SEPARATORS; ++ iDirSep) {
//                char cDirSepToTest = CX_PATH_ALL_SEPARATORS[iDirSep];
//                if (FileName[n - nBase - 1] == cDirSepToTest) {
//                   cDirSep = cDirSepToTest;
//                   break;
//                }
//             }
//             if (cDirSep != 0) break;
//          }
//          size_t nPath = n - nBase;
// //          fmt::print(" {:20} n = {:<3} np = {:<3} nb = {:<3}\n\n", "'"+FileName+"'", n, nPath, nBase);
//          return FFileName_PathAndBase(FileName.substr(0, nPath), FileName.substr(nPath, nBase));
//       }
#endif // INCLUDE_ABANDONED

      FFileName join(FFileName const &a, FFileName const &b) {
         if (!a.empty() && 0 != dir_seperator_q(a.back()))
            // `a` already ends in a dir separator. Don't need to add another one in-between
            return a + b;
         else
            return fmt::format("{}" CX_PATH_SEPARATOR "{}", a, b);
//             return a + std::string(CX_PATH_SEPARATOR) + b;
      }

      FFileName join(FFileNames const &L) {
         std::stringstream ss;
         bool NeedSep = false;
         for (FFileNames::const_iterator it = L.begin(); it != L.end(); ++ it) {
            if (it->empty())
               continue;
            if (NeedSep)
               ss << CX_PATH_SEPARATOR;
            ss << *it;
            NeedSep = (0 == dir_seperator_q(it->back()));
            // ^-- if *it already ends in a dir separator, we don't need to add another one
            // before inserting the next argument
//             fmt::print(" frag = {:20}  NeedSep? {}\n", fmt::format("'{}'", *it), NeedSep);
         }
         return ss.str();
      }

      static bool _CanPathBeSplit(FFileName const &Path) {
         // check if Path is a special such as c:, \\computer (for network shares in Windows)
         // or /usr for a path specified relative to the root directory.
         if (Path.empty())
            return false;
         size_t
            nPathSepsAtStart = 0;
         for (size_t i = 0; i < Path.size() && dir_seperator_q(Path[i]); ++ i) {
            nPathSepsAtStart = i;
         }
         for (size_t i = nPathSepsAtStart; i < Path.size() && dir_seperator_q(Path[i]); ++ i) {
            // found another path separator *not* right at the start --> path can be split at that one.
            // (e.g., "\\storm\netshare" can be split (into "\\storm\" and "netshare"), but "\\storm"
            // cannot be split)
            return true;
         }
         return false;
      }

      FFileNames split_path(FFileName const &FileName, size_t MaxSplit)
      {  // WARNING: I think this code may be 100% untested!
         FFileNames Fragments;
         {
            std::string
               Rest = FileName;
            while (!Rest.empty() && Fragments.size() <= MaxSplit) { // <- n splits make n+1 fragments.
               // split path/base into 'path/' and 'base'
               size_t iBase = _iBaseNameStartQ(Rest);
               FFileName_PathAndBase rs{Rest.substr(0,iBase), Rest.substr(iBase, std::string::npos)};

               // remove trailing path separator from path, but only if there is anything left in it
               // which is not a path separator (otherwise things like network shares (\\computer\bla)
               // or root-relative files (/opt/bla) would cause problems).
               Fragments.push_back(without_trailing_path_separator(std::move(rs.second)));
               Rest = std::move(rs.first);
            }
            if (!Rest.empty())
               Fragments.push_back(std::move(Rest));
         }
         std::reverse(Fragments.begin(), Fragments.end());
         return Fragments;
      }
#ifdef INCLUDE_ABANDONED
//       FFileNames split_path(FFileName const &FileName) {
//          using FFileNameLinkedList = std::list<FFileName>;
//          FFileNameLinkedList
//             FragmentList;
//          {
//             std::string
//                Rest = FileName;
//             while (!Rest.empty()) {
//                // split path/base into 'path/' and 'base'
//                FFileName_PathAndBase rs = split(Rest);
//                // remove trailing path separator from path, but only if there is anything left in it
//                // which is not a path separator (otherwise things like network shares (\\computer\bla)
//                // or root-relative files (/opt/bla) would cause problems).
// #ifdef INCLUDE_ABANDONED
// //                if (!rs.first.empty() && dir_seperator_q(rs.first.back()) != 0 && _CanPathBeSplit(rs.first))
// //                   rs.first.pop_back();
// //                FragmentList.emplace_back(std::move(rs.second));
// #endif // INCLUDE_ABANDONED
//                FragmentList.emplace_back(without_trailing_path_separator(std::move(rs.second)));
//                Rest = std::move(rs.first);
//             }
//          }
//
//          FFileNames out;
//          out.reserve(FragmentList.size());
//          for (FFileNameLinkedList::reverse_iterator it = FragmentList.rbegin(); it != FragmentList.rend(); ++ it) {
//             out.emplace_back(std::move(*it));
//          }
//          return out;
//       }
#endif // INCLUDE_ABANDONED


      bool file_readable_q(char const *pFileName) {
         std::ifstream si(pFileName);
         return si.good(); // <-- slow but simple...
      }

      bool file_readable_q(std::string const &FileName) { return file_readable_q(FileName.c_str()); }


      FFileName with_trailing_path_separator(FFileName &&PathName) {
         if (PathName.empty())
            return PathName;
         if (!dir_seperator_q(PathName.back()))
            PathName.push_back(CX_PATH_ALL_SEPARATORS[0]);
         return PathName;
      }

      FFileName with_trailing_path_separator(FFileName const &PathName) { return with_trailing_path_separator(std::move(FFileName(PathName))); }


      FFileName without_trailing_path_separator(FFileName &&PathName) {
         // remove trailing path separator from path, but only if there is anything left in it
         // which is not a path separator (otherwise things like network shares (\\computer\bla)
         // or root-relative files (/opt/bla) would cause problems).
         if (!PathName.empty() && dir_seperator_q(PathName.back()) != 0 && _CanPathBeSplit(PathName))
            PathName.pop_back();
         return PathName;
      }

      FFileName without_trailing_path_separator(FFileName const &PathName) { return without_trailing_path_separator(std::move(FFileName(PathName))); }

   } // namespace path

} // namespace ct
