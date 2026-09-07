from socket import gethostname
from os import path, listdir
import os

# Set up floating point environments. Comments:
# - The floating point types control the precision with which aigg internally
#   represents intermediate results.
#
# - In addition to using native "double precision" (on typical machines these
#   are standard IEEE 64bit floats), it allows using various different high-
#   precision floating point types in which the floating point arithmetic is
#   emulated.
#
# - Larger floating point types allow constructing grids of higher order, but
#   the arithmetic does become slower and the memory demand becomes higher.
#
#   + Native 64bit arithmetic is *BY FAR* the fastest option. However,
#     we recommend using this ONLY for pre-optimization of grids, because
#     when using only 64bits to represent intermediate results, it is clearly
#     not possible to generate full 64bit precision (53 mantissa bits) in the
#     final output data.
#
#   + Apart from this, high precision floats may also be needed because aigg
#     currently uses straight non-orthogonal Cartesian monomials to span the
#     target space, and the overlap matrix of these can acquire a very large
#     condition number.
#
#   + Emulated 128bit floats should be good to handle up to ~L=70 grids reliably
#     (for full 64bit output precision), and up to ~L=100 with some possible
#     hickups (e.g., slow or erratic convergence).
#
# - This script allows compiling multiple aigg executables with different
#   floating point types in parallel. This will result in different executables
#   distinguished by a file name suffix.
class FFloatEnv(object):
   def __init__(self, Suffix, CcFlags, LinkFlags, Libs=[]):
      """Defines parameters of the build environment associated
      with a given floating point type used to represent intermediate
      numbers inside aigg.
      
      Arguments:
      - Suffix: will be appended to exe name 
      - CcFlags: list of command line arguments to pass to the compiler (e.g., g++)
      - LinkFlags: list of command line arguments to pass to the linker
      - Libs: list of additional library files to include
        (same as passing -l arguments to LinkFlags)
      """
      self.Suffix = Suffix
      self.CcFlags = CcFlags
      self.LinkFlags = LinkFlags
      self.Libs = Libs

#BaseCxxStd = "-std=gnu++11"
def BaseCxxStd(s):
   #BaseVersion = "11"
   BaseVersion = "17"
   return "--std={}{}".format(s, BaseVersion)


# this one uses C++ "double" types for scalars
FloatEnv_NativeDouble = FFloatEnv("_64", CcFlags=["-DSCALAR_IS_FLOAT_64", BaseCxxStd("c++")], LinkFlags=[], Libs=[])

# this one uses the non-standard C++-extension "__float128" types for scalars
# (now wrapped via boost::multiprecision), which some compilers provide. In
# particular, this applies to GNU g++, which is what the present code assumes to
# be used.
FloatEnv_GccFloat128 = FFloatEnv("_128", CcFlags=["-DSCALAR_IS_FLOAT_128", "-DFLOAT_128_NATIVE_EXT", BaseCxxStd("gnu++"), "-fext-numeric-literals"], LinkFlags=[], Libs=["quadmath"])

# this one is also using 128bit floats, but employs
# boost::multiprecision::cpp_bin_float as computational backend. This is more
# portable (does not rely on a gcc-like __float128 compiler extension), but may
# be slower (on my system by ~factor two for larger optimizations).
FloatEnv_CppFloat128 = FFloatEnv("_128", CcFlags=["-DSCALAR_IS_FLOAT_128", "-faligned-new", "-std=c++17"], LinkFlags=[])
# 256bit floats on cpp_bin_float. In my inits tests, a run with this
# required about two times as long as with 128bit cpp_bin_float
FloatEnv_CppFloat256 = FFloatEnv("_256", CcFlags=["-DSCALAR_IS_FLOAT_256", "-faligned-new", "-std=c++17"], LinkFlags=[])
# 512bit floats on cpp_bin_float. In my inits tests, a run with this
# required about two times as long as with 256bit cpp_bin_float and
# four times as long as with 128bit cpp_bin_float
FloatEnv_CppFloat512 = FFloatEnv("_512", CcFlags=["-DSCALAR_IS_FLOAT_512", "-faligned-new", "-std=c++17"], LinkFlags=[])


# configure list of different floating point types to build the
# release (=production) version of aigg for.
FloatTypesToBuild_Release = []
FloatTypesToBuild_Release += [FloatEnv_NativeDouble]
if 1:
   # use this if you have gcc or something else which can deal with
   # __float128 and the g++ arguments.
   FloatTypesToBuild_Release += [FloatEnv_GccFloat128]
else:
   # use this for compatible c++-only version.
   FloatTypesToBuild_Release += [FloatEnv_CppFloat128]
# note: if you want to do both 128bit versions (e.g., for performance comparison)
# you can enable both, provided one of them gets a different suffix from the other one.
FloatTypesToBuild_Release += [FloatEnv_CppFloat256]
FloatTypesToBuild_Release += [FloatEnv_CppFloat512]

# by default, build only one debug version of aigg (whatever is
# first in the release variants).
FloatTypesToBuild_Debug = FloatTypesToBuild_Release[:1]


# make a list of compiler and linker flags which will be common
# to all build variants (debug/release and different float types)
CcFlagsCommon = []
LinkFlagsCommon = []
LibsCommon = []

# use commands like this to overwrite boost or eigen installation directories.
# Note that neither of them needs to be compiled for this program---the header-
# only versions of the boost libraries are sufficient (and eigen is normally not
# separately compiled anyway)
if 0:
   CcFlagsCommon += ["-I/opt/prog/boost_1_69_0"]

if gethostname() in ["mirage", "charge", "storm", "tempest"]:
   # ^- these are my (cgk)'s computers. I like keeping a git clone
   #    of the master branch of eigen.
   CcFlagsCommon += ["-I/home/cgk/Programs/eigen"]


if 1:
   # this turns on OpenMP parallelization. However, a large part of the
   # computational time is spent inside eigen's linear algebra subroutines, and
   # almost none of those are actually parallelized. So... in general you should
   # not assume this to do much.
   CcFlagsCommon += ["-fopenmp"]
   LinkFlagsCommon += ["-fopenmp"]


CcFlagsDebug = ["-g", "-O0", "-D_DEBUG"]
CcFlagsRelease = []
CcFlagsRelease += ["-DNDEBUG"] # switch off assertions
if 1:
   CcFlagsRelease += ["-O2"]
else:
   # that's -O3 with -ffast-math. On my system this doesn't make a
   # noticeable difference to -O2 with -ffast-math.
   CcFlagsRelease += ["-Ofast"]

if 1:
   # set of details of compiler configuration. This one is for g++ and compilers
   # which pretend to be g++ (e.g., clang or icc under linux)
   CcFlagsCommon += "-fmax-errors=2 -Wall -Wno-unknown-pragmas -Wno-unused-function -Wno-ignored-attributes".split()
   CcFlagsDebug = ["-g", "-O0", "-D_DEBUG"]
   # - setting -ffast-math allows the compiler to do some transformations
   #   to floating point expressions during optimization passes which may
   #   not necessarily preserve the exact value of the original expressions.
   #   For example, it allows the compiler treat floating point arithmetic
   #   expressions AS IF they were commutative (strictly speaking, they are not),
   #   and then do re-orderings of summations like
   #
   #       (a + b) + (c + d) -> (a + c) + (b + d)
   #
   #   or simple factorizations of expressions like
   #
   #       a * b + c * b = (a+c)*b.
   #
   # - This will not affect the emulated high-precision floats, but generally does
   #   help quite a bit when using native floating point arithmetic
   CcFlagsRelease += ["-ffast-math"]
   # - march=native tells the compiler to optimize for the CPU/architecture the
   #   program is compiled on
   CcFlagsRelease += ["-march=native"]

if 0:
   # if this set, include debug symbol names also in the (optimized) release
   # builds. note: this can make VERY large executables due to all the floats
   # being wrapping with boost multiprecision classes (which results in very
   # long symbol names)
   CcFlagsRelease += ["-g"]

# collect list of target files to compile (all cpp files in source directory)
SourceDir = "src"
#FileList = [path.join(SourceDir,s) for s in listdir(SourceDir) if ".cpp" in path.splitext(s)[1]]
# FileList = [s for s in listdir(SourceDir) if ".cpp" in path.splitext(s)[1]]
# Use SCons's built-in Glob to find files inside the src directory safely
FileList = Glob(f'{SourceDir}/*.cpp')


# set up actual build environments
def MakeTarget(NameOut, ObjDir, Env ):
    Env.Append(CPPPATH=['/usr/include/eigen3'])
    Env.VariantDir(ObjDir, SourceDir, duplicate=0)
    obj_sources = [os.path.join(ObjDir, os.path.basename(str(f))) for f in FileList]    
    return Env.Program(NameOut, source=obj_sources)
 
def MakeTargetForScalarType(ExeName, BuildDir, CcFlags, LinkFlags, Libs):
   # This program requires at least C++11 standard (or its gnu-extensions
   # variant) check if a C++ standard parameter was provided in CcFlags.
   # If not, add one.
   CxxStandardGiven = False
   for CcFlag in CcFlags:
      if "-std=" in CcFlag:
         CxxStandardGiven = True
   if not CxxStandardGiven:
      CcFlags = [o for o in CcFlags] + ["-std=c++11"]      

   target = MakeTarget(ExeName, BuildDir,
      Environment(CCFLAGS = " ".join(CcFlags), LIBS = Libs, LINKFLAGS = " ".join(LinkFlags)))
   return target

def MakeDebugTargetForScalarType(FloatEnv):
   ExeName = "aigg{}_d".format(FloatEnv.Suffix)
   BuildDir = path.join(*["build", "debug{}".format(FloatEnv.Suffix)])
   CcFlags = CcFlagsCommon + CcFlagsDebug + FloatEnv.CcFlags
   LinkFlags = LinkFlagsCommon + FloatEnv.LinkFlags
   Libs = LibsCommon + FloatEnv.Libs
   return MakeTargetForScalarType(ExeName, BuildDir, CcFlags, LinkFlags, Libs)

def MakeReleaseTargetForScalarType(FloatEnv):
   ExeName = "aigg{}".format(FloatEnv.Suffix)
   BuildDir = path.join(*["build", "release{}".format(FloatEnv.Suffix)])
   CcFlags = CcFlagsCommon + CcFlagsRelease + FloatEnv.CcFlags
   LinkFlags = LinkFlagsCommon + FloatEnv.LinkFlags
   Libs = LibsCommon + FloatEnv.Libs
   return MakeTargetForScalarType(ExeName, BuildDir, CcFlags, LinkFlags, Libs)

#def MakeReleaseTargetForScalarType(NameOut, ObjDir):
   #opt = MakeTarget( "aigg", "build/release/", Environment(CCFLAGS = "-g -Ofast -ffast-math -march=native -DNDEBUG " + CcFlagsCommon, LIBS = Libs, LINKFLAGS = " ".join(LinkFlags)))

#FloatTypesToBuild_Release = []

for FloatTypeToBuild in FloatTypesToBuild_Debug:
   tgt = MakeDebugTargetForScalarType(FloatTypeToBuild)
   Default(tgt)

for FloatTypeToBuild in FloatTypesToBuild_Release:
   tgt = MakeReleaseTargetForScalarType(FloatTypeToBuild)
   Default(tgt)


#Default(dbg)
#Default(opt)
