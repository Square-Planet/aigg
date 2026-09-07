#ifndef CT_IO_H
#define CT_IO_H

// TODO: port ~/dev/geocow/src/io.h 's IndentedLog to FIndentLog here.

#include <string>
#include <stdexcept>
#include <sstream>

#include "CxDefs.h" // for assert
#if defined(USE_GENERAL_SCALARS) && !defined(SCALAR_IS_FLOAT_64)
   #include "CxScalarTypes.h"
#endif

#include "CxPodArray.h"
#include "format.h" // brings in namespace "fmt", including its format function (which makes std::strings)

namespace ct {
#ifdef INCLUDE_ABANDONED
// extern char const
//    *pResultFmt,
//    *pResultFmtAnnotated,
//    *pResultFmtI,
//    *pResultFmtIAnnotated,
//    *pTimingFmt,
//    *pTimingPerFmt;
#endif // INCLUDE_ABANDONED
   std::string FmtInfo(fmt::BasicStringRef<char> Name, fmt::BasicStringRef<char> Value, fmt::BasicStringRef<char> Annotation = fmt::StringRef("",0));
   std::string FmtExpfInfo(fmt::BasicStringRef<char> Name, double Value, fmt::BasicStringRef<char> Annotation = fmt::StringRef("",0));
   std::string FmtResult(fmt::BasicStringRef<char> Name, double Value, fmt::BasicStringRef<char> Annotation = fmt::StringRef("",0));
   std::string FmtCount(fmt::BasicStringRef<char> Name, ptrdiff_t Value, fmt::BasicStringRef<char> Annotation = fmt::StringRef("",0));
   std::string FmtTiming(fmt::BasicStringRef<char> Name, double fTimeInSeconds, size_t nTasks = size_t(-1));
} // namespace ct


namespace ct {
   extern std::ostream
      &xerr, &xout;
   extern int
      Verbosity;
   // ^-- TODO: I think the global `Verbosity` is no longer used anywhere, and the xerr/xout probably also
   // are only retained for some legacy code. Consider removing.
}

namespace ct {

// If compiled in unicode mode, measures the physical output width (in terms of monospace text
// columns) of an UTF-8-encoded string `p` (i.e., does the same as CxUnicode.h's `ColumnWidthQ()`).
// If not compiled in unicode mode, just returns `n`.
ptrdiff_t TextWidthQ(char const *p, size_t n);
inline ptrdiff_t TextWidthQ(char const *p) { return TextWidthQ(p, std::char_traits<char>::length(p)); }
inline ptrdiff_t TextWidthQ(std::string const &s) { return TextWidthQ(s.data(), s.size()); }



// makes a string to be output into the log as the first line of a major
// program component.
void MajorProgramIntro(std::ostream &out, const std::string &Name, const std::string &Version="");

// will load an entire file into a TArray<char>. Adds a 0-terminator at end.
// Returns false if failed.
// if 0 != pFileLength, *pFileLength will receive the number of loaded chars.
bool LoadFileIntoMemory(TArray<char> &pFileContent,
        std::string const &FileName, unsigned *pFileLength = 0);


// // error in some sort of input -- (e.g., syntax, brken file formats).
// class FInputError : public std::runtime_error
// {
// public:
//    typedef std::runtime_error
//       FBase;
//    explicit FInputError(std::string const &Reason);
// };
//
// ^-- moved to CxParse1.h

// // some unsupported calculation or combination of options was asked for (e.g., a tau-type meta-GGA with DFXC=1)
// class FUnsupportedError : public std::runtime_error
// {
// public:
//    typedef std::runtime_error
//       FBase;
//    explicit FUnsupportedError(std::string const &Reason);
// };
//
// ^-- moved to CxParse1.h


/// Error in an input/output operation (e.g., attempt to load a non-existing
/// file, failed file write due to inaccessible path, full disk, exceeded quota,
/// insufficient permissions, etc)
class FIoError : public std::runtime_error {
public:
   enum FOpType {
      IOOP_NotSpecified, ///< no specific operation type specified
      IOOP_Read,   ///< operation: attempted reading of a file
      IOOP_Write,  ///< operation: attempted writing to a file
      IOOP_Open,   ///< operation: attempted opening of a file or directory
      IOOP_Find,   ///< operation: find the physical location of a file by name
      IOOP_Sync,   ///< operation: require for explicit flush/sync of queued operations
      IOOP_Create, ///< operation: attempted creation of a file or directory
      IOOP_Delete  ///< operation: attempted deletion of a file or directory
   };
   enum FFileType {
      IOFILE_NotSpecified,
      IOFILE_RegularFile,
      IOFILE_Directory
   };

   typedef std::runtime_error FBase;
   /// Instanciate FIoError to indicate a failed input/output operation
   ///
   /// @param pAffectedFileOrPathName if != 0, indicates the name of the file/path being processed when the error occurred
   /// @param OpType If != IOOP_NotSpecified, indicates the type of IO operation attempted when the error occurred
   explicit FIoError(std::string const &Reason, char const *pAffectedFileOrPathName = 0, FOpType OpType = IOOP_NotSpecified, FFileType FileType = IOFILE_NotSpecified);

protected:
   static char const *_FormatOpType(FOpType OpType);
   static char const *_FormatFileType(FFileType OpType);
   static std::string _FormatMessage(std::string const &Reason, char const *pAffectedFileOrPathName, FOpType OpType, FFileType FileType);
};


class FLogError : public std::runtime_error {
public:
   typedef std::runtime_error
      FBase;
   explicit FLogError(std::string const &Reason);
};


// error in converging some sort of calculation for which there is no workaround declared
class FConvergenceError : public std::runtime_error {
public:
   typedef std::runtime_error
      FBase;
   explicit FConvergenceError(std::string const &CalcType, ptrdiff_t iFinalIt = -1, double fFinalResidual = -1., double fFinalEnergy = -1);
};
std::string FormatConvergenceError(std::string const &CalcType, ptrdiff_t iFinalIt = -1, double fFinalResidual = -1., double fFinalEnergyChange = -1);



// a class encapsulating output streams and associated messaging.
//
// Main point about not simply using std::cout or some other global thing
// is that computations might be done on distributed machines or in some
// other non-local/non-sequential way.
class FLog {
protected:
public:
   fmt::MemoryWriter
      w;
   typedef char
      Char;

   FLog();
   virtual ~FLog();

   enum FWriteResultFlags{
      OUTPUT_Tabulate = 0x02 // used with WriteResult(...). Marks an output meant for tabulation by an external program
   };

   void Write(fmt::BasicStringRef<Char> format) {
      w << format << '\n';
      if (GetFlag(IOFLAG_FlushAfterWrite))
         Flush();
   }
   void WriteNoNl(fmt::BasicStringRef<Char> format) {
      w << format;
      if (GetFlag(IOFLAG_FlushAfterWrite))
         Flush();
   }
   void Write(fmt::BasicStringRef<Char> format, fmt::ArgList args) {
      fmt::BasicFormatter<Char>(w).format(format, args);
      w << '\n';
      if (GetFlag(IOFLAG_FlushAfterWrite))
         Flush();
   }
   void WriteNoNl(fmt::BasicStringRef<Char> format, fmt::ArgList args) {
      fmt::BasicFormatter<Char>(w).format(format, args);
      if (GetFlag(IOFLAG_FlushAfterWrite))
         Flush();
   }
   FMT_VARIADIC_VOID(Write, fmt::BasicStringRef<Char>)
   FMT_VARIADIC_VOID(WriteNoNl, fmt::BasicStringRef<Char>)

   virtual void Flush() = 0;

   virtual void EmitWarning(fmt::BasicStringRef<Char> Name);
   virtual void EmitError(fmt::BasicStringRef<Char> Name);
   void WriteResult(fmt::BasicStringRef<Char> Name, double fValue);
   void WriteResult(fmt::BasicStringRef<Char> Name, double fValue, fmt::BasicStringRef<Char> Annotation);
   void WriteResult(fmt::BasicStringRef<Char> Name, double fValue, unsigned Flags);
#if defined(USE_GENERAL_SCALARS) && !defined(SCALAR_IS_FLOAT_64)
   void WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue);
   void WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue, fmt::BasicStringRef<Char> Annotation);
   void WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue, unsigned Flags);
#endif // USE_GENERAL_SCALARS
   void WriteInfo(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Value);
   void WriteInfo(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Value, fmt::BasicStringRef<Char> Annotation);
   void WriteInfoExpf(fmt::BasicStringRef<Char> Name, double fValue);
   void WriteInfoExpf(fmt::BasicStringRef<Char> Name, double fValue, fmt::BasicStringRef<Char> Annotation);
   void WriteCount(fmt::BasicStringRef<Char> Name, ptrdiff_t iValue);
   void WriteCount(fmt::BasicStringRef<Char> Name, ptrdiff_t iValue, fmt::BasicStringRef<Char> Annotation);
   void WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds);
   void WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds, fmt::BasicStringRef<Char> Annotation);
   // print relative timing per task.
   void WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds, size_t nTasks);

   void WriteLine() { w << "\n"; };
   void WriteProgramIntro(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Version="");

   // this is our mechanism for communicating with the outside world (the other side of the log).
   // GetStatus may return something else than STATUS_Okay, for example, if an abort request was
   // placed on the other side by the user.
   // Default implementation does nothing (returns STATUS_Okay).
   enum FRunStatus {
      STATUS_Okay = 0,
      STATUS_AbortSignalled = 1
   };
   virtual FRunStatus GetStatus();
   // call GetStatus; throw FLogError if status says something differnt from STATUS_Okay.
   void CheckStatus();
   bool StatusOkay() { return GetStatus() == STATUS_Okay; };

public:
   enum FIoFlag {
      IOFLAG_FlushAfterWrite = 0x0001,
      IOFLAG_Max = IOFLAG_FlushAfterWrite
   };

   inline bool GetFlag(FIoFlag Flag) {
      _CxAssert(Flag <= IOFLAG_Max);
      return bool(m_IoFlags & Flag);
   }
   void SetFlag(FIoFlag Flag, bool NewState) {
      _CxAssert(Flag <= IOFLAG_Max);
      if (NewState)
         m_IoFlags |= Flag;
      else
         m_IoFlags &= ~Flag;
   }
   // meant for flags with more than two states; e.g., could be used for
   // local indent level.
   typedef unsigned
      FIoFlagValue;
   inline FIoFlagValue GetFlagValue(FIoFlag Flag) {
      _CxAssert(Flag <= IOFLAG_Max);
      return m_IoFlags & FIoFlagValue(Flag);
   }
   void SetFlagValue(FIoFlag Flag, FIoFlagValue NewValue) {
      _CxAssert(Flag <= IOFLAG_Max);
      m_IoFlags = (m_IoFlags & (~Flag)) | NewValue;
   }
protected:
   FIoFlagValue
      m_IoFlags;
private:
   void operator = (FLog const &); // not implemented
   FLog(FLog const &); // not implemented.
};


class FLogIoFlagOverride {
public:
   FLogIoFlagOverride(FLog *pLog, FLog::FIoFlag Flag, bool NewState) {
      m_pLog = pLog;
      m_Flag = Flag;
      m_NewState = NewState;
      if (m_pLog) {
         m_OriginalState = m_pLog->GetFlag(m_Flag);
         m_pLog->SetFlag(m_Flag, m_NewState);
      } else {
         m_OriginalState = false;
      }
   }
   ~FLogIoFlagOverride() {
      if (m_pLog) {
         m_pLog->SetFlag(m_Flag, m_OriginalState);
         if (m_Flag == FLog::IOFLAG_FlushAfterWrite && m_NewState == false && m_OriginalState == true)
            m_pLog->Flush();
      }
   }
protected:
   FLog
      *m_pLog;
   FLog::FIoFlag
      m_Flag;
   bool
      m_OriginalState, m_NewState;
};


/// log of which the output goes to a generic std::ostream.
class FLogStdStream : public FLog {
public:
   explicit FLogStdStream(std::ostream &TargetStream);
   void Flush(); // override
protected:
   std::ostream
      &m_TargetStream;
};



/// log which keeps the text in a memory object, without actually emitting it
/// anywhere actively.
class FMemoryLog : public FLogStdStream {
public:
   FMemoryLog();
   ~FMemoryLog();

   enum FGetTextFlags {
      GETTEXT_RetainLog = 0x0000,
      // if set, the log buffer will be cleared out after GetText() is called.
      GETTEXT_ClearLog = 0x0001
   };

   /// assemble a single std::string containing the entire log content, and return it.
   /// Note: this does *not* clear the current text.
   std::string GetText(unsigned Flags = GETTEXT_ClearLog);
   /// flush and clera out the entire controlled log text
   void Clear();
protected:
   std::stringstream
      m_TextStream;
};


struct FFileLocatorImpl;



/// A support structure holding functionality for locating physical files (e.g.,
/// data files from the application, such as basis set libraries) relative to
/// the application directory, rather than the current working directory.
///
/// Used to allow executing programs like MicroScf to locate its auxiliary data
/// if executed via command line from other directories.
///
/// Note:
/// - The BasePath of the file search environment has to be globally set once
///   via SetBasePath().
/// - Typically, one would use the location of the executable:
///
///    FFileLocator::SetBasePath(GetExePath()); // with GetExePath() from CxOsInt.h
///
/// - Afterwards, files and additional search directories can be specified
///   relative to this base path. E.g., in AddSearchPath, '$BASEPATH' is replaced
///   by the base path of the application
/// - The implementation in principle supports additional path substitutions
///   (e.g., for HOME directories or CONFIG/PROFILE directories). But there
///   is no interface to those yet, because it is not yet used in my programs.
struct FFileLocator {

   FFileLocator();
   ~FFileLocator();
   /// Store application base path in a global variable subsequently accessed
   /// in newly constructed FFileLocator objects.
   ///
   /// Notes:
   /// - Global to allow multiple instanciations of the object in lower level
   ///   structures of the code.
   /// - Typically this would be the application exe directory on linux, or
   ///   the application base directory on Windows.
   /// - Must be called before first invocation of FFileLocator ("." is okay
   ///   if no additional info is available).
   static void SetBasePath(std::string const &BasePath);

   static std::string JoinPath(std::string const &Path, std::string const &FileName);
   static std::string BaseName(std::string const &PathAndFileName);
   static std::string PathName(std::string const &PathAndFileName);

   enum FFindFileFlags {
      FINDFILE_RaiseErrorIfNotFound = 0x0001,
      // if set, the substitutions in m_Substitutions will not be applied
      // to the file name. These are otherwise used to register substitutions
      // like $BASEPATH or $HOME.
      FINDFILE_NoSubstitutions = 0x0002,
      FINDFILE_Default = FINDFILE_RaiseErrorIfNotFound
   };

   std::string FindFile(std::string const &FileName, unsigned Flags=FINDFILE_Default) const;
   /// note: '$BASEPATH' is replaced by the base path of the application.
   void AddSearchPath(std::string const &PathName);

   /// adds a substitution of $HOME and ~/ (at start) to the current user's
   /// home directory, if possible.
   void AddHomePath();

   /// imports environment variable sVarName as a substitution into *this.
   void ImportEnvironmentVar(std::string const &VarName);


   void PushState();
   void PopState();
protected:
   FFileLocatorImpl
      *p;
   void operator = (FFileLocator const &); // not implemented
   FFileLocator(FFileLocator const &); // not implemented
};


/// On construction, saves the current state of a FFileLocator object;
/// on destruction, restores the original state. This can be used to
/// savely add local search paths (e.g., when processing input files
/// in different directories) without compromising the structure once
/// done.
struct FFileLocatorStateGuard {
   FFileLocatorStateGuard(FFileLocator &FileLocator_)
      : m_FileLocator(FileLocator_)
   {
      m_FileLocator.PushState();
   }
   ~FFileLocatorStateGuard()
   {
      m_FileLocator.PopState();
   }
protected:
   FFileLocator
      &m_FileLocator;
private:
   void operator = (FFileLocatorStateGuard const&); // not implemented.
   FFileLocatorStateGuard(FFileLocatorStateGuard const&); // not implemented.
};


namespace io {

   template<class string_t>
   struct FRepeatProxy {
      ptrdiff_t _count;
      string_t _unit;
      void _emitUnit(std::ostream &out) const { out << this->_unit; }
   };

   template<>
   struct FRepeatProxy<char const*> {
      ptrdiff_t _count;
      char const *_unit;
      size_t _len;
      void _emitUnit(std::ostream &out) const { if (_len != 0) out.write(_unit, _len); }
   };
#ifdef INCLUDE_ABANDONED
//       explicit FRepeatProxy(size_t count, char const *unit, size_t len = size_t(-1))
//          : _count(count), _unit(unit),
//            _len((len != size_t(-1))? len : ((unit == 0)? 1 : std::char_traits<char>::length(unit)))
//       {}
//       void _emitUnit(std::ostream &out) const {
//          if (_len == 0) return;
//          // actual string set? ⟶ emit whatever it is.
//          if (_unit) return (void)out.write(_unit, _len);
//          // special case: (unit == 0 && len != 0) ⟶ emit spaces.
//          size_t w0 = out.width(_len); out << ""; out.width(w0);
//       }
#endif // INCLUDE_ABANDONED

   template<class string_t>
   std::ostream &operator << (std::ostream &out, FRepeatProxy<string_t> const &rp) {
      for (ptrdiff_t i = 0; i < rp._count; ++ i)
         rp._emitUnit(out);
      return out;
   }

   template<class string_t>
   auto repeat(ptrdiff_t count, string_t &&unit) -> FRepeatProxy<string_t> {
      return {count, std::forward<string_t>(unit)};
   }

   inline auto repeat(ptrdiff_t count, char const *unit) -> FRepeatProxy<char const *> {
      return {count, unit, bool(unit)? std::char_traits<char>::length(unit) : 0u};
   }

   struct FWhiteSpaceProxy { ptrdiff_t _count; };
   void wspace(std::ostream &out, ptrdiff_t count);
   inline auto wspace(ptrdiff_t count) -> FWhiteSpaceProxy { return { count }; }

   inline std::ostream &operator << (std::ostream &out, FWhiteSpaceProxy rp) { wspace(out, rp._count); return out; }


}


namespace obj_io {
   // “osc” ⟶ “open, separator, close” — controls translation of data nodes
   struct osc_t {
      // strings to emit before/after the node (e.g., `[`, `]`), between its elements (e.g., `,`),
      // and between field-names and field-values (e.g., `:`)
      char const *open, *sep, *close, *assoc;

      // set to 1 to request fields/values being placed on separate lines;
      // set to 0 for default spacing;
      // set to -1 or -2 to request that space compression (emit fewer/no unnecessary spaces)
      int spacing;
      constexpr bool linesQ() const { return spacing > 0; }
      constexpr bool normalQ() const { return spacing == 0; }
      constexpr bool compressQ() const { return spacing < 0; }
   };

   struct config_t {
      std::ostream
         // target output stream. Null if `*this` has been moved from, otherwise must be valid stream.
         *out;
      // base indent string supplied from outside — supposed to be emitted left of any new line
      char const
         *ind0;
      bool quoteFieldsQ() const { return true; }
   };

   struct node_t {
      auto config() const -> config_t const *{ return &_config; }
      auto xout() -> std::ostream& { return *this->config()->out; }
      auto ind0() const -> char const *{ return this->config()->ind0; }

      node_t(config_t config, osc_t osc);
      node_t(node_t *outer, osc_t osc);
      ~node_t();

      void newLine();
      void lineInd();

      // adds a field to an object-valued node (should be followed by a value to assign it to)
      auto field(char const *name) -> node_t&; // ⟵ returns `*this`

      // push a new scalar value to `*this`. This is specialized for
      // strings (to quote them, and later hopefully properly escape them)
      template<class T>
      auto value(T const &t) -> node_t& { _elemSep(); _value1(t); return *this; }

      // shortcut: `node(t)` translates to `node.value(t)`
      template<class T>
      auto operator()(T const &t) -> node_t& { return this->value(t); }

      // push a new list-valued item to `*this`; create and return a child `node_t` object
      // representing its content
      auto list(int spacing=0) -> node_t;
      // push a new object-valued item to `*this` (i.e., a sequence of (field, value)-associations);
      // create and return a child `node_t` object representing its content
      auto obj(int spacing=1) -> node_t;

      // “reflective map”: calls `fn(*this, e)` for each element of `seq`
      template<class Seq, class Fn>
      void rmap(Seq &&seq, Fn &&fn) { for (auto &&e : seq) fn(*this, e); }

      // “value map”: calls `this->value(fn(e))` for each element of `seq`
      template<class Seq, class Fn>
      void vmap(Seq &&seq, Fn &&fn) { for (auto &&e : seq) this->value(fn(e)); }

      // “value push”: calls `this->value(e)` for each element of `seq`
      template<class Seq, class Fn>
      void values(Seq &&seq) { for (auto &&e : seq) this->value(e); }

      // just a convenience function which combines `this->field(name)` with `this->value(value)`.
      template<class T>
      auto field(char const *name, T const &value) -> node_t& { this->field(name); this->value(value); return *this; }

      // combines `this->field(name)` with `this->list().vmap(seq, fn)`
      template<class Seq, class Fn>
      auto field(char const *name, Seq const &seq, Fn &&fn) -> node_t& { this->field(name); this->list().vmap(seq, fn); return *this; }

      node_t(node_t &&);
      auto operator = (node_t &&) -> node_t&;
   protected:
      // node this object is a child of, if any
      node_t *_parent; // ⟵ may be 0
      config_t _config;
      // open/separator/close tags: to strings to emit before the node content (e.g., '[' for a
      // list), to separate node-fields/values (e.g., ', '), and after the node content (e.g., ']')
      osc_t _osc;
      // dynamic indent depth; increases if this node is a child object of another
      int _idepth;
      // counts items emitted for this node. mostly used to track whether or not
      // a given item is the first one, and therefore requires an object separator.
      int _istate;

      // returns whether `*this` is an “undead” object — the leftover of a move c'tor.
      inline bool _zombieQ() const;
      // marks this node object as a zombie, to prevent its destructor from
      // taking actions (e.g., ending lines or closing the object)
      void _zombify();

      void _text(int quote, char const *p, size_t n);
      inline void _text(char const *p, size_t n);
      inline void _text(char const *p);
      inline void _text(int quote, char const *p);
      inline void _text(int quote, std::string const &s);
      inline void _wspace(ptrdiff_t w);
      inline void _wspace();

      template<class T>
      void _value1(T const &v) { xout() << v; }
      void _value1(char const *v);
      void _value1(std::string const &v);

      void _objOpen();
      void _elemSep();
      void _objClose();

   private:
      node_t(node_t const &) = delete;
      void operator = (node_t const &) = delete;
   };

   void _xDataOpen(std::ostream &xout, char const *name, config_t const &config);
   void _xDataClose(std::ostream &xout, char const *name, config_t const &config);

   // export a machine-readable object/list definition to the xout stream (at top-level)
   template<class Fn>
   auto xdata(std::ostream &xout, char const *name, osc_t osc, Fn fn) -> std::ostream& {
      config_t cfg{&xout, "  "};
      _xDataOpen(xout, name, cfg);
      { node_t root{cfg, osc}; fn(root); }
      _xDataClose(xout, name, cfg);
      return xout;
   }

   template<class Fn> inline auto xobj(std::ostream &xout, int spacing, char const *name, Fn fn) -> std::ostream& {
      return xdata(xout, name, osc_t{"{",",","}",":",spacing}, fn);
   }

   template<class Fn> inline auto xlist(std::ostream &xout, int spacing, char const *name, Fn fn) -> std::ostream& {
      return xdata(xout, name, osc_t{"[",",","]",":",spacing}, fn);
   }

} // namespace obj_io




} // namespace ct.

#endif // CT_IO_H
