all:
	scons -Q --jobs=16
opt:
	scons -Q opt
dbg:
	scons -Q dbg
clean:
	scons --clean
