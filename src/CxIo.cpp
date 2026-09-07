#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <cctype> // for tolower()
#include <cstdlib> // for getenv

#include <list> // used in FFileLocator
#include <map> // used in FFileLocator
#include <set> // used in FFileLocator

#ifdef QT_CORE_LIB
   #include <QString>
   #include <QByteArray>
   #include <QFile>
   #include <QApplication>
   #include <QClipboard>
#endif

#if (defined(CX_UNICODE) && CX_UNICODE >= 1)
   #include "CxUnicode.h"
   namespace ct {
      ptrdiff_t TextWidthQ(char const *p, size_t n) { return static_cast<ptrdiff_t>(ct::unicode::ColumnWidthQ(p, n)); }
   }
#else
   namespace ct {
      ptrdiff_t TextWidthQ(char const *, size_t n) { return static_cast<ptrdiff_t>(n); }
   }
#endif

// #define USE_INDENT_STREAM
#ifdef USE_INDENT_STREAM
   #define INDENTSTREAM_IMPL
   #include "CxIndentStream.h"
   namespace ct {
      static fmt::FIndentStream1
         s_OutputStreamAdp(std::cout,true,1);
      std::ostream
         &xout = s_OutputStreamAdp.stream,
         &xerr = std::cerr;
   } // namespace ct
#else
   namespace ct {
      std::ostream
         &xout = std::cout,
         &xerr = std::cerr;
   } // namespace ct
#endif // USE_INDENT_STREAM

#include "CxParse1.h" // used in FFileLocator
#include "CxIo.h"


namespace ct {

bool g_FileLocator_DebugPrints = false;

// char const
//    *pResultFmt = " %-32s%18.12f\n",
//    *pResultFmtAnnoted = " %-32s%18.12f  (%s)\n",
//    *pResultFmtI = " %-32s%5i\n",
//    *pResultFmtIAnnoted = " %-32s%5i  (%s)\n",
//    *pTimingFmt = " Time for %s:%40t%10.2f sec\n",
//    *pTimingPerFmt = " Time per %s:%40t%10.2f sec  (n=%i)\n";
//    *pTimingFmt = " Time for %s:%40t%10.4f sec\n",
//    *pTimingPerFmt = " Time per %s:%40t%10.4f sec  (n=%i)\n";

// int
//    Verbosity = 0;
// ^-- still used anywhere?

void FatalError( std::string const &Message,
    char const *pFromWhere, int nLine )
{
   std::stringstream str;
   str << Message;
   if (pFromWhere) {
      str << " at " << pFromWhere << ":" << nLine;
   };
   throw std::runtime_error(str.str());
}



// makes a string to be output into the log as the first line of a major
// program component.
void MajorProgramIntro( std::ostream &out, const std::string &Name, const std::string &Version )
{
//    out << fmt::unind();
   out << "\n*** " << Name;
   if (Version != "")
      out << " [Ver. " << Version << "]";
   out << " ***\n" << std::endl;
//    out << fmt::eind();
}


bool LoadFileIntoMemory(TArray<char> &pFileContent,
      std::string const &FileName, unsigned *pFileLength)
{
//    std::cout << "in LoadFileIntoMemory. Name: '" << FileName << "'." << std::endl;
#ifdef QT_CORE_LIB
   // Use QT functions to load the file. This gives us access to files
   // stored in resources and the application clipboard.
   QByteArray
      Text;
   if (FileName == ":/!clipboard!") {
      QClipboard
         *clipboard = QApplication::clipboard();
      QString
         Subtype = "plain",
         sText = clipboard->text(Subtype);
      Text = sText.toUtf8();
   } else {
      QFile
         StyleFile(QString(FileName.c_str()));
      if (!StyleFile.open(QFile::ReadOnly))
         return false;
      Text = StyleFile.readAll();
   }
   pFileContent.resize(Text.size() + 2);
   // 0-terminate the string.
   pFileContent[Text.size()] = 0;
   pFileContent[Text.size()+1] = 0;
   memcpy(&pFileContent[0], Text.data(), Text.size());
   if (pFileLength)
       // if requested, store actual file length at target location.
      *pFileLength = Text.size();
   return true;
#else
   // read file data with C++ standard library routines
   std::ifstream
      File(FileName.c_str());
   std::size_t
      FileLength;
   if (false == File.good())
      return false;
   File.seekg(0, std::ios::end);
   FileLength = File.tellg();
   if (0 != pFileLength)
      // if requested, store actual file length at target location.
      *pFileLength = FileLength;
   pFileContent.resize(2 + FileLength);
   // ^- cgk 2020-02-27: hm... what's with these +2 things? I guess one might
   // be for allowing interpretation as 0-terminated string (but since we don't
   // scan for other 0s in the file content, that would not be safe, either).
   // But two?
   // Update: apparently the result of this is attached to std::stringstreams
   // in several places, which require 0-termination. Still doesn't explain the
   // two 0s.
   memset(&pFileContent[0], 0, 2 + FileLength);
   File.seekg(0, std::ios::beg);
   File.read(&pFileContent[0], FileLength);
   return true;
#endif
}




std::string FormatConvergenceError(std::string const &CalcType, ptrdiff_t iFinalIt, double fFinalResidual, double fFinalEnergyChange)
{
   fmt::MemoryWriter w;
   w.write("{} failed to converge", CalcType);
   if (iFinalIt != -1)
      w.write(" within {} iterations", iFinalIt);
   w << ".";
   if (fFinalResidual != -1. || fFinalEnergyChange != -1.) {
      w << " (";
      bool bFirst = false;
      if (fFinalResidual != -1.) {
         bFirst = false;
         w.write("Final residual: {:8.2e}", fFinalResidual);
      }
      if (fFinalEnergyChange != -1.) {
         if (!bFirst) w << "; ";
         bFirst = false;
         w.write("Final energy change: {:16.8f}", fFinalEnergyChange);
      }
      w << ")";
   }
   return w.str();
}

FConvergenceError::FConvergenceError(std::string const &CalcType, ptrdiff_t iFinalIt, double fFinalResidual, double fFinalEnergy)
   : FBase(FormatConvergenceError(CalcType, iFinalIt, fFinalResidual, fFinalEnergy))
{}



FLogError::FLogError(std::string const &Reason)
   : FBase(Reason)
{}


FLog::FLog()
   : m_IoFlags(IOFLAG_FlushAfterWrite)
{
}


FLog::~FLog()
{
}

#ifdef INCLUDE_ABANDONED
// char const
//    *p1ResultFmtForTabulation = "${:<32}| {:24.16e} |  ({})\n",
//    // ^- used to mark quantities intended for tabulation by a script analyzing the program output
//    *p1ResultFmt =            " {:<32}{:18.12f}\n",
//    *p1ResultFmtAnnotated =   " {:<32}{:18.12f}  ({})\n",
//    *p1CountFmt =             " {:<32}{:5}\n",
//    *p1CountFmtAnnotated =    " {:<32}{:5}  ({})\n",
//    *p1InfoFmt =              " {:<32}{:>18}\n",
//    *p1InfoFmtAnnotated =     " {:<32}{:>18}  ({})\n",
//    *p1InfoExpfFmt =          " {:<32}{:18.4e}\n",
//    *p1InfoExpfFmtAnnotated = " {:<32}{:18.4e}  ({})\n",
//    *p1TimingFmt =            "{:<40}{:10.2f} sec\n",
//    *p1TimingFmtAnnotated =   "{:<40}{:10.2f} sec  ({}{})\n";
//
// ^-- before 2021-06-27 update.
//     Note: Apart from CtEmb code, MicroScf and wmme apparently only use Writecount in two instances:
//
//         CtMain.cpp:   m_pLog->WriteCount("Allocated workspace (MB)", m_WorkSpaceMb);
//         CtMain.cpp:   m_pLog->WriteCount("Number of CPU-threads", omp_get_max_threads());
//         CtMain.cpp:   Log.WriteCount("Allocated Workspace (MB)", WorkSpaceMb);
//
//     So at least this particular problem can probably be fixed freely.
//     The timing fmt without annotation appears to be off by one char, just as
//     it looks it should be from the numbers (normally it is meant to align the right
//     border of the times with the right border of result outputs, I think,
//     but it presently doesn't:
//
//             10     -266.96663500    -0.00000011     0.00000000   1.10e-05      0.73  9  9
//             11     -266.96663502    -0.00000001     0.00000000   3.78e-06      0.78  0 10
//
//          Time per iteration:                          0.05 sec  (n=11)
//
//          Nuclear repulsion energy          163.607581115976
//          1-electron energy                -687.675690850437
//          Band energy                      -148.707968220842
//
//          Lost electrons                         -2.7699e-05
//          .../sqrt(nAtoms)                       -9.2329e-06
//
//     So I'll move this by one character. Once done, the timing annotation (for n = ...) needs to
//     stay where it is, however. Moving it would impact some things.
#endif // INCLUDE_ABANDONED

char const
   *p1ResultFmtForTabulation = "${:<32}| {:24.16e} |  ({})\n",
   // ^- used to mark quantities intended for tabulation by a script analyzing the program output
   *p1ResultFmt =            " {:<32}{:18.12f}\n",
   *p1ResultFmtAnnotated =   " {:<32}{:18.12f}  ({})\n",
   *p1CountFmt =             " {:<38}{:12}\n",
   *p1CountFmtAnnotated =    " {:<38}{:12}  ({})\n",
   *p1InfoFmt =              " {:<32}{:>18}\n",
   *p1InfoFmtAnnotated =     " {:<32}{:>18}  ({})\n",
   *p1InfoExpfFmt =          " {:<32}{:18.4e}\n",
   *p1InfoExpfFmtAnnotated = " {:<32}{:18.4e}  ({})\n",
   *p1TimingFmt =            "{:<41}{:10.2f} sec\n",        // <-- note: gets an extra leading " " before "Time for..."
   *p1TimingFmtAnnotated =   "{:<41}{:10.2f} sec  ({}{})\n";


std::string FmtInfo(fmt::BasicStringRef<char> Name, fmt::BasicStringRef<char> Value, fmt::BasicStringRef<char> Annotation) {
   return Annotation.empty()?
        fmt::format(p1InfoFmt, Name, Value)
      : fmt::format(p1InfoFmtAnnotated, Name, Value, Annotation);
}

std::string FmtExpfInfo(fmt::BasicStringRef<char> Name, double Value, fmt::BasicStringRef<char> Annotation) {
   return Annotation.empty()?
        fmt::format(p1InfoExpfFmt, Name, Value)
      : fmt::format(p1InfoExpfFmtAnnotated, Name, Value, Annotation);
}

std::string FmtResult(fmt::BasicStringRef<char> Name, double Value, fmt::BasicStringRef<char> Annotation) {
   return Annotation.empty()?
        fmt::format(p1ResultFmt, Name, Value)
      : fmt::format(p1ResultFmtAnnotated, Name, Value, Annotation);
}

std::string FmtCount(fmt::BasicStringRef<char> Name, ptrdiff_t Value, fmt::BasicStringRef<char> Annotation) {
   return Annotation.empty()?
        fmt::format(p1CountFmt, Name, Value)
      : fmt::format(p1CountFmtAnnotated, Name, Value, Annotation);
}

std::string FmtTiming(fmt::BasicStringRef<char> Name, double fTimeInSeconds, size_t nTasks) {
   if (nTasks == size_t(-1))
      return fmt::format(p1TimingFmt, fmt::format(" Time for {}:", Name), fTimeInSeconds);
   else
      return fmt::format(p1TimingFmtAnnotated, fmt::format(" Time per {}:", Name), fTimeInSeconds/double(nTasks), "n=", nTasks);
}



void FLog::WriteResult(fmt::BasicStringRef<Char> Name, double fValue) {
   w.write(p1ResultFmt, Name, fValue);
   Flush();
}


void FLog::WriteResult(fmt::BasicStringRef<Char> Name, double fValue, fmt::BasicStringRef<Char> Annotation) {
   w.write(p1ResultFmtAnnotated, Name, fValue, Annotation);
   Flush();
}


void FLog::WriteResult(fmt::BasicStringRef<Char> Name, double fValue, unsigned Flags) {
   if ((Flags & OUTPUT_Tabulate) != 0) {
      WriteResult(Name, fValue); // write it normally once
      w.write(p1ResultFmtForTabulation, Name, fValue, "->!Tabulate-Average");
      Flush();
   } else {
      WriteResult(Name, fValue);
   }
}

#if defined(USE_GENERAL_SCALARS) && !defined(SCALAR_IS_FLOAT_64)
void FLog::WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue) {
   w.write(p1ResultFmt, Name, fValue); Flush();
}

void FLog::WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue, fmt::BasicStringRef<Char> Annotation) {
   w.write(p1ResultFmtAnnotated, Name, fValue, Annotation); Flush();
}


void FLog::WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue, unsigned Flags) {
   if ((Flags & OUTPUT_Tabulate) != 0) {
      WriteResult(Name, fValue); // write it normally once
      w.write(p1ResultFmtForTabulation, Name, fValue, "->!Tabulate-Average");
      Flush();
   } else {
      WriteResult(Name, fValue);
   }
}
#endif // USE_GENERAL_SCALARS



void FLog::WriteCount(fmt::BasicStringRef<Char> Name, ptrdiff_t iValue) {
   w.write(p1CountFmt, Name, iValue); Flush();
}

void FLog::WriteCount(fmt::BasicStringRef<Char> Name, ptrdiff_t iValue, fmt::BasicStringRef<Char> Annotation) {
   w.write(p1CountFmtAnnotated, Name, iValue, Annotation); Flush();
}

void FLog::WriteInfo(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Value) {
   w.write(p1InfoFmt, Name, Value); Flush();
}

void FLog::WriteInfo(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Value, fmt::BasicStringRef<Char> Annotation) {
   w.write(p1InfoFmtAnnotated, Name, Value, Annotation); Flush();
}

void FLog::WriteInfoExpf(fmt::BasicStringRef<Char> Name, double Value) {
   w.write(p1InfoExpfFmt, Name, Value); Flush();
}


void FLog::WriteInfoExpf(fmt::BasicStringRef<Char> Name, double Value, fmt::BasicStringRef<Char> Annotation) {
   if (Annotation.empty())
      return WriteInfoExpf(Name, Value);
   w.write(p1InfoExpfFmtAnnotated, Name, Value, Annotation);
   Flush();
}


void FLog::WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds) {
   w.write(p1TimingFmt, fmt::format(" Time for {}:", Name), fTimeInSeconds); Flush();
}


void FLog::WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds, fmt::BasicStringRef<Char> Annotation) {
   if (!Annotation.empty()) {
      w.write(p1TimingFmtAnnotated, fmt::format(" Time for {}:", Name), fTimeInSeconds, Annotation, "");
      Flush();
   } else {
      return this->WriteTiming(Name, fTimeInSeconds);
   }
}


void FLog::WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds, size_t nTasks) {
   w.write(p1TimingFmtAnnotated, fmt::format(" Time per {}:", Name), fTimeInSeconds/double(nTasks), "n=", nTasks);
   Flush();
}


void FLog::WriteProgramIntro(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Version) {
   w << "\n*** " << Name;
   if (Version.size() != 0u)
      w << " [Ver. " << Version << "]";
   w << " ***\n\n";
   Flush();
}


FLog::FRunStatus FLog::GetStatus() {
   // do nothing. derived classes might do more involved things.
   return STATUS_Okay;
}


void FLog::CheckStatus() {
   // do nothing. derived classes might, however.
   FRunStatus
      Status = GetStatus();
   if (Status == STATUS_Okay)
      return;
   if (Status == STATUS_AbortSignalled)
      throw ct::FLogError("Abort signalled.");
   throw ct::FLogError("FLog signalled unrecognized error.");
}


FLogStdStream::FLogStdStream(std::ostream &TargetStream)
   : m_TargetStream(TargetStream)
{
}


void FLog::EmitWarning(fmt::BasicStringRef<Char> Name) {
   w << " WARNING: " << Name << "\n";
   Flush();
}


void FLog::EmitError(fmt::BasicStringRef<Char> Name) {
   w << " ERROR: " << Name << "\n";
   Flush();
}


void FLogStdStream::Flush() {
   m_TargetStream << w.c_str();
   w.clear();
   m_TargetStream.flush();
}




FMemoryLog::FMemoryLog()
   : FLogStdStream(m_TextStream)
{}


FMemoryLog::~FMemoryLog()
{}


std::string FMemoryLog::GetText(unsigned Flags) {
   Flush();
   std::string Text = m_TextStream.str();
   if (bool(Flags & GETTEXT_ClearLog))
      Clear();
   return Text;
}


void FMemoryLog::Clear() {
   Flush();
   m_TextStream.str(std::string());
}



} // namespace ct


namespace ct {
   char const *FIoError::_FormatOpType(FOpType OpType) {
      switch (OpType) {
         case IOOP_Read:    return "read from";
         case IOOP_Write:   return "write to";
         case IOOP_Open:    return "open";
         case IOOP_Sync:    return "sync/flush";
         case IOOP_Create:  return "create";
         case IOOP_Delete:  return "delete";
         default:           { _CxAssert(OpType == IOOP_NotSpecified); return 0; }
      }
   }


   char const *FIoError::_FormatFileType(FFileType OpType) {
      switch (OpType) {
         case IOFILE_RegularFile:  return "file";
         case IOFILE_Directory:    return "directory";
         default:                  { _CxAssert(OpType == IOFILE_NotSpecified); return "file/path"; }
      }
   }


   std::string FIoError::_FormatMessage(std::string const &Reason, char const *pAffectedFileOrPathName, FOpType OpType, FFileType FileType) {
      fmt::MemoryWriter w;
      w << "IO error";
      char const *pOpFmt = _FormatOpType(OpType);
      if (pOpFmt != 0 || pAffectedFileOrPathName != 0) {
         w << " (@";
         size_t iArg = 0;
         if (pOpFmt != 0) {
            if (iArg != 0) w << " ";
            w << pOpFmt;
            iArg += 1;
         }
         if (FileType != IOFILE_NotSpecified) {
            if (iArg != 0) w << " ";
            w << _FormatFileType(FileType);
            iArg += 1;
         }
         if (pAffectedFileOrPathName != 0) {
            if (iArg != 0) w << " ";
            w.write("'{}'", pAffectedFileOrPathName);
            iArg += 1;
         }
         w << ")";
      }
      if (!Reason.empty())
         w.write(": {}", Reason);
      return w.str();
   }


   FIoError::FIoError(std::string const &Reason, char const *pAffectedFileOrPathName, FOpType OpType, FFileType FileType)
      : FBase(_FormatMessage(Reason, pAffectedFileOrPathName, OpType, FileType)) {
   }
}


namespace ct {

static bool
   // to allow for empty base paths...
   g_FileLocator_BasePathIsSet = false;
static std::string
   // FFileLocator objects will attempt to resolve file names relative to this.
   // By default, that would be what is returned by GetExePath()
   g_FileLocator_BasePath;


void FFileLocator::SetBasePath(std::string const &BasePath)
{
   g_FileLocator_BasePath = BasePath;
   g_FileLocator_BasePathIsSet = true;
}



enum FStrSubstFlags {
   STRSUBST_MatchAtFrontOnly = 0x0001
};

struct FFileLocatorState : public FIntrusivePtrDest
{
   typedef std::list<std::string>
      FDirectoryList;
   FDirectoryList
      m_SearchPaths;
   typedef std::set<std::string>
      FDirectorySet;
   FDirectorySet
      m_SearchPaths_Unordeded;

   struct FSubstEntry {
      std::string Target;
      unsigned Flags;
      FSubstEntry() : Flags(0) {};
      FSubstEntry(std::string Target_, unsigned Flags_) : Target(Target_), Flags(Flags_) {}
   };
   typedef std::map<std::string, FSubstEntry>
      FSubstMap;
   FSubstMap
      m_Substitutions;
};
typedef TIntrusivePtr<FFileLocatorState>
   FFileLocatorStatePtr;


struct FFileLocatorImpl
{
   void AddSearchPath(std::string const &PathName);

   FFileLocatorImpl();

   std::string FindFile(std::string const &FileName, unsigned Flags) const;

   void PushState();
   void PopState();

   void AddHomePath();
   void ImportEnvironmentVar(std::string VarName);
   void AddSubstitution(std::string const &From, std::string const &To, unsigned Flags);
protected:
   typedef FFileLocatorState::FDirectoryList
      FDirectoryList;
   typedef FFileLocatorState::FDirectorySet
      FDirectorySet;
   typedef std::list<FFileLocatorStatePtr>
      FStateStack;
   typedef std::list<bool>
      FStateStackFlags;
   typedef FFileLocatorState::FSubstEntry
      FSubstEntry;
   typedef FFileLocatorState::FSubstMap
      FSubstMap;

   FStateStack
      m_StateStack;
   FStateStackFlags
      // we only copy the actual state stack contents when they are
      // modified. Before that, the current level is aliased to the
      // higher one. This list keeps track of whether or not the current
      // state is (still) aliased to the previous one.
      m_StateAliased;
   size_t
      m_StateStackDepth;

   void UnaliasState(); // if current stack level is aliased to the last one, make a copy and un-alias it.
   FFileLocatorState &State() { _CxAssert(!m_StateStack.empty()); return *m_StateStack.front(); }
   FFileLocatorState const &State() const { _CxAssert(!m_StateStack.empty()); return *m_StateStack.front(); }
   FDirectoryList &m_SearchPaths() { return State().m_SearchPaths; }
   FDirectoryList const &m_SearchPaths() const { return State().m_SearchPaths; }
   FDirectorySet &m_SearchPaths_Unordeded() { return State().m_SearchPaths_Unordeded; }
   FDirectorySet const &m_SearchPaths_Unordeded() const { return State().m_SearchPaths_Unordeded; }
   FSubstMap &m_Substitutions() { return State().m_Substitutions; }
   FSubstMap const &m_Substitutions() const { return State().m_Substitutions; }
   bool &_StateAliased() { _CxAssert(!m_StateAliased.empty()); return m_StateAliased.front(); }

   // apply substitutions and remove trailing '/' / '\' characters.
   std::string ProcessFileName(std::string const &FileName) const;
private:
   void operator = (FFileLocatorImpl const &); // not implemented
   FFileLocatorImpl(FFileLocatorImpl const &); // not implemented
};


FFileLocator::FFileLocator()
   : p(new FFileLocatorImpl)
{
}

FFileLocator::~FFileLocator() {
   delete p;
   p = 0;
}


FFileLocatorImpl::FFileLocatorImpl() {
   m_StateStack.push_front(FFileLocatorStatePtr(new FFileLocatorState()));
   m_StateStackDepth = 0; // no "extra" layers of state (only the current top layer, but we don't count this here).
   m_StateAliased.push_front(false); // original state is never aliased.
   if (!g_FileLocator_BasePathIsSet)
      throw std::runtime_error("FFileLocator constructor: Use SetBasePath() before instanciating FFileLocator objects.");
   std::string
      // note: overrides in the base path are used in the Debug/Release configurations of VC
      // on windows. Otherwise the relative directory stuff does not work.
      sBasePath = g_FileLocator_BasePath;
//    if (sBasePath.size() == 0)
//       sBasePath = GetExePath();
   if (sBasePath.empty())
      AddSubstitution("$BASEPATH", ".", 0);
   else
      AddSubstitution("$BASEPATH", ProcessFileName(sBasePath), 0);
}


void FFileLocatorImpl::AddSubstitution(std::string const &From, std::string const &To, unsigned Flags) {
   if (g_FileLocator_DebugPrints)
      std::cout << fmt::format(": register file path substitution: '{}' --> '{}' (flags: {})\n", From, To, Flags);
   FSubstMap::iterator
      it = m_Substitutions().find(From);
   if (it == m_Substitutions().end() || it->second.Flags != Flags || it->second.Target != To) {
      UnaliasState();
      m_Substitutions()[From] = FSubstEntry(To, Flags);
   }
}


static std::string GetEnvironmentVar(char const *pName) {
   _CxAssert(pName != 0);
   if (pName[0] == '$')
      // skip over leading "$" if in the variable name.
      pName += 1;
   return std::string(std::getenv(pName));
}


void FFileLocatorImpl::AddHomePath() {
   std::string
      sHomePath;
#ifdef _WIN32
// "%HOMEDRIVE%%HOMEPATH%"
   sHomePath = FFileLocator::JoinPath(GetEnvironmentVar("HOMEDRIVE"), GetEnvironmentVar("HOMEPATH"));
#else
   sHomePath = GetEnvironmentVar("HOME");
#endif
   AddSubstitution("$HOME", sHomePath, 0);
   AddSubstitution("~/", sHomePath + "/", STRSUBST_MatchAtFrontOnly);
   // ^- the extra '/' is necessary because we remove the '/' of '~/'.
}


void FFileLocator::AddHomePath() {
   return p->AddHomePath();
}


void FFileLocatorImpl::ImportEnvironmentVar(std::string VarName) {
   if (VarName.empty())
      return;
   if (VarName[0] == '$')
      // remove leading "$" if supplied (e.g., turn "$HOME" into "HOME")
      VarName.erase(0, 1);
   AddSubstitution("$" + VarName, GetEnvironmentVar(VarName.c_str()), 0);
}


void FFileLocator::ImportEnvironmentVar(std::string const &VarName) {
   return p->ImportEnvironmentVar(VarName);
}


void FFileLocatorImpl::PushState() {
   m_StateStack.push_front(m_StateStack.front());
   // mark that current stack level is aliased to the last one
   // (the line above copies the *pointer* to the lower level state).
   m_StateAliased.push_front(true);
   m_StateStackDepth += 1;
}


void FFileLocatorImpl::PopState() {
   if (m_StateStackDepth == 0)
      throw std::runtime_error("FFileLocatorImpl::PopState(): attempted to restore a non-existent state");
   m_StateStackDepth -= 1;
   m_StateStack.pop_front();
   m_StateAliased.pop_front();
}


void FFileLocatorImpl::UnaliasState()
{
   _CxAssert(!m_StateAliased.empty());
   if (_StateAliased()) {
      // make a copy of the (aliased) state of the current level, and replace current
      // level state pointer by a pointer to this copy.
      FFileLocatorStatePtr
         pCopyOfState(new FFileLocatorState(State()));

      m_StateStack.front() = pCopyOfState;
      _StateAliased() = false;
   }
}


void FFileLocator::PushState() {
   p->PushState();
}


void FFileLocator::PopState() {
   p->PopState();
}


std::string FFileLocator::JoinPath(std::string const &Path, std::string const &FileName)
{
   string_slice
      slPath(Path);
   if (slPath.empty() || slPath == "." || slPath == "./")
      return FileName;
   if (slPath.endswith("/")
#ifdef _WIN32
         || slPath.endswith("\\")
#endif
      )
      return fmt::format("{}{}", Path, FileName);
   else {
#ifdef _WIN32
      return fmt::format("{}\\{}", Path, FileName);
#else
      return fmt::format("{}/{}", Path, FileName);
#endif
   }
}


static void ReplaceSubStr(std::string &inout, std::string const &key, std::string const &value, unsigned Flags)
{
   size_t ipos = std::string::npos;
   if (bool(Flags & STRSUBST_MatchAtFrontOnly)) {
      if (StartsWith(inout, string_slice(key)))
         ipos = 0;
   } else {
      ipos = inout.find(key);
   }
   if (ipos != std::string::npos) {
      inout.replace(ipos, key.size(), value);
   }
}


std::string FFileLocatorImpl::ProcessFileName(std::string const &FileName_) const
{
   std::string
      FileName = FileName_;
   FSubstMap::const_iterator
      itSubst;
   // apply directory substitutions
   for (itSubst = m_Substitutions().begin(); itSubst != m_Substitutions().end(); ++ itSubst) {
      ReplaceSubStr(FileName, itSubst->first, itSubst->second.Target, itSubst->second.Flags);
   }

   // get rid of spaces to either side and trailing back- and front-slashes.
   string_slice
      slName(FileName);
   slName.trim();
   while (slName.first != slName.last && (slName.endswith("/") || slName.endswith("\\")))
      slName.last -= 1;
   return slName.to_str();
}


void FFileLocatorImpl::AddSearchPath(std::string const &PathName)
{
   if (PathName.empty())
      return;
   std::string
      s = ProcessFileName(PathName);
   // add the string to our search path list... unless we already have it there.
   if (m_SearchPaths_Unordeded().find(s) == m_SearchPaths_Unordeded().end()) {
      // we'll modify the state. Unlink it from higher level states now, unless already done.
      UnaliasState();
      m_SearchPaths().push_front(ProcessFileName(PathName));
      m_SearchPaths_Unordeded().insert(s);
   }
}


static bool IsFileReadable(std::string const &FileName) {
   // well... slow, but simple.
   std::ifstream
      File(FileName.c_str());
   return File.good();
}


std::string FFileLocatorImpl::FindFile(std::string const &FileName_, unsigned Flags) const
{
   if (g_FileLocator_DebugPrints) {
      std::cout << fmt::format(": searching for {} (#paths: {})\n", FileName_, m_SearchPaths().size());
   }
   std::string
      FileName;
   if (bool(Flags & FFileLocator::FINDFILE_NoSubstitutions)) {
      FileName = FileName_;
   } else {
      FileName = ProcessFileName(FileName_);
   }

   FDirectoryList::const_iterator
      itPath;
   for (itPath = m_SearchPaths().begin(); itPath != m_SearchPaths().end(); ++ itPath) {
      std::string
         s = FFileLocator::JoinPath(*itPath, FileName);
      if (IsFileReadable(s)) {
         if (g_FileLocator_DebugPrints) {
            std::cout << fmt::format(": -> found '{}'\n", s);
         }
         return s;
      }
   }
   // couldn't find it.
   if (bool(Flags & FFileLocator::FINDFILE_RaiseErrorIfNotFound)) {
      fmt::MemoryWriter w;
      w << "! Search paths in FindFile():\n";
      for (itPath = m_SearchPaths().begin(); itPath != m_SearchPaths().end(); ++ itPath) {
         w << fmt::format("!: '{}'\n", *itPath);
      }
      w << fmt::format("FindFile: Failed to locate file '{}'", FileName);
//       throw std::runtime_error(w.str());
      throw FIoError(w.str(), FileName_.c_str(), FIoError::IOOP_Find, FIoError::IOFILE_RegularFile);
   } else {
      return FileName;
   }
}


std::string FFileLocator::FindFile(std::string const &FileName, unsigned Flags) const {
   return p->FindFile(FileName, Flags);
}


void FFileLocator::AddSearchPath(std::string const &PathName) {
   return p->AddSearchPath(PathName);
}


static size_t FindPathFileSeparator(std::string const &PathAndFileName)
{
   size_t
      iSlash = PathAndFileName.rfind('/');
#ifdef _WIN32
   size_t
      iBackSlash = PathAndFileName.rfind('\\');
   if (iBackSlash != std::string::npos && (iSlash == std::string::npos || iBackSlash > iSlash))
      iSlash = iBackSlash;
#endif
   _CxAssert(iSlash == std::string::npos || iSlash < PathAndFileName.size());
   return iSlash;
}


std::string FFileLocator::BaseName(std::string const &PathAndFileName)
{
   size_t
      iSlash = FindPathFileSeparator(PathAndFileName);
   if (iSlash != std::string::npos) {
      if (iSlash + 1 >= PathAndFileName.size())
         return "";
      else
         return PathAndFileName.substr(iSlash+1);
   } else {
      // no directory separator. Return original file name.
      return PathAndFileName;
   }
}


std::string FFileLocator::PathName(std::string const &PathAndFileName)
{
   size_t
      iSlash = FindPathFileSeparator(PathAndFileName);
   if (iSlash != std::string::npos) {
      // str contains a directory separator. Return input str up to
      // this separator (excluding the separator itsef).
      // (note: That will prevent adding just "/" as a search path,
      // but I do not think this would be a great idea in any case...)
      return PathAndFileName.substr(0, iSlash);
   } else {
      // no directory separator. Return empty string.
      return std::string();
   }
}


namespace io {
#ifdef INCLUDE_ABANDONED
//    void wspace(std::ostream &out, ptrdiff_t count) {
//       if (count > 0) {
//          // special case for spaces: use ostream's own filling algo to fill it up.
//          size_t w0 = out.width(_len); out << ""; out.width(w0);
//       }
//    }
#endif // INCLUDE_ABANDONED
   //                              1234567890123456
   static char const *_wspace16 = "                ";
   void wspace(std::ostream &out, ptrdiff_t count) {
      while (count > 16) { out.write(_wspace16, 16); count -= 16; }
      if (count > 0) out.write(_wspace16, size_t(count));
   }

}

namespace obj_io {
   node_t::node_t(config_t cfg, osc_t osc)
      : _parent(0), _config(cfg), _osc{osc}, _idepth(0), _istate(0)
   { _objOpen(); };


   node_t::node_t(node_t *outer, osc_t osc)
      : _parent(outer), _config(*_parent->config()), _osc(osc), _idepth(_parent->_idepth + 1), _istate(0)
   { _CxAssert(!_parent->_zombieQ()); _objOpen(); }

   bool node_t::_zombieQ() const { return _idepth == -1; }
   void node_t::_zombify() { _idepth = -1; _config = config_t{0,0}; }

   node_t::~node_t() { if (!_zombieQ()) _objClose(); }
   node_t::node_t(node_t &&src) { *this = std::move(src); }

   auto node_t::operator = (node_t &&src) -> node_t& {
      this->_parent = src._parent;
      this->_config = src._config;
      this->_osc = src._osc;
      this->_idepth = src._idepth;
      this->_istate = src._istate;
      src._zombify();
      return *this;
   }

#ifdef INCLUDE_ABANDONED
//    void node_t::_wspace(ptrdiff_t w) { xout() << ct::io::wspace(w); }
//    void node_t::_text(char const *p) { if (p) { xout() << p; } }
#endif // INCLUDE_ABANDONED

   void node_t::_text(char const *p, size_t n) {
      xout().write(p, n);
   }
   void node_t::_wspace(ptrdiff_t w) {
      ct::io::wspace(xout(), w);
   }
   void node_t::_wspace() {
      this->_text(" ", 1);
   }

   void node_t::_text(int quote, char const *p, size_t n) {
      if (quote == 0) return this->_text(p,n);
      // TODO: escape the text! especially quotes and newlines!
      if (quote == 2) quote = int('"');
      if (quote == 1) quote = int('\'');
      auto &_out = xout();
      _out.put(static_cast<char>(quote));
      _out.write(p,n);
      _out.put(static_cast<char>(quote));
   }

   void node_t::_text(int quote, char const *p) {
      this->_text(quote, p, std::char_traits<char>::length(p));
   }
   void node_t::_text(int quote, std::string const &s) {
      this->_text(quote, s.data(), s.size());
   }
   void node_t::_text(char const *p) {
      if (p) { this->_text(p, std::char_traits<char>::length(p)); }
   }

   void node_t::_value1(char const *v) { _text(2, v); }
   void node_t::_value1(std::string const &v) { _text(2, v); }

   void node_t::_objOpen() {
      if (_parent)
         _parent->_elemSep();
      if(!_parent && _osc.linesQ()) lineInd();
      // ⬑ hm… not sure if that really is quite right. indents etc. Should cross-check another time when it comes up.
      _text(_osc.open);
      if(_osc.linesQ()) { _idepth += 1; newLine(); }
      _istate = 0;
   }

   void node_t::_elemSep() {
      _CxAssert(!_zombieQ());
      if (_istate > 0) {
         _text(_osc.sep);
         if(_osc.linesQ()) newLine();
         if(_osc.normalQ()) _wspace();
      }
      _istate += 1;
   }

   void node_t::_objClose() {
      _CxAssert(!_zombieQ());
      if(_osc.linesQ()) { _idepth -= 1; newLine(); }
      _text(_osc.close);
      if (_parent) _parent->_istate += 1;
   }


   void node_t::lineInd() { _text(ind0()); _wspace(2*_idepth); }
   void node_t::newLine() { _text("\n"); lineInd(); }

   auto node_t::list(int spacing) -> node_t { return {this, osc_t{"[", ",", "]", ":", spacing}}; }
   auto node_t::obj(int spacing) -> node_t { return {this, osc_t{"{", ",", "}", ":", spacing}}; }


   auto node_t::field(char const *name) -> node_t& {
      _CxAssert(!_zombieQ());
      _elemSep();
      _istate = 0; // ⟵ prevent emission of element-separator at next `value()` call
      _text(int(config()->quoteFieldsQ()? int('"') : 0), name);
      _text(_osc.assoc);
      if (_osc.normalQ()) {
         _wspace();
      } else if (_osc.linesQ()) {
         ptrdiff_t w = 16;
         _wspace(w - ptrdiff_t(TextWidthQ(name)));
      }
      return *this;
   }

   void _xDataOpen(std::ostream &xout, char const *name, config_t const &config) {
      xout <<   "$ ##!beg:" << name << (config.quoteFieldsQ()? "/json/" : "/js/") <<"\n";
   }

   void _xDataClose(std::ostream &xout, char const *name, config_t const &config) {
      xout << "\n$ ##!end:" << name << (config.quoteFieldsQ()? "/json/" : "/js/") <<"\n";
   }
} // namespace obj_io




} // namespace ct

