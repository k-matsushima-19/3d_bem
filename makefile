TARGETDIR = ./bin
TARGETNAME = a.out

# インストール先
DEST = ../../bin

# 再帰的にmakeを実行するディレクトリ (引数は継承されない)
SUBDIRS = 

# Library paths (すでにpathが通っていないもの)
LIB_PATHS = 

# Library names (右は左を必要としない順で並べる)
LIB_NAMES = gsl gslcblas lapack blas miniball

# OPENMPを使うか
OPENMP = 1

# 最適化をon
#OPT = 1

# Fortran
FC = gfortran
FFLAGS +=

# C
CC = gcc
CFLAGS += 
INCLUDES = ./

# Linker
LINK = $(FC)
LDFLAGS +=


#
# Compile options
#

# Fortran
ifeq ($(FC),gfortran)
  FFLAGS += -I$(OBJDIR) -J$(OBJDIR) -cpp $(FDEBUGS) -DGFORTRAN
  ifdef OPENMP
    FFLAGS += -fopenmp
    LDFLAGS += -fopenmp
  endif
  FOPTS = -O3
  FDEBUGS = -O0 -g -fbounds-check -Wuninitialized -Wconversion -ffpe-trap=invalid,zero,overflow -fbacktrace -DDEBUG
endif

ifeq ($(FC),ifort)
  FFLAGS += -I$(OBJDIR) -module $(OBJDIR) -fpp $(FDEBUGS) -DIFORT
  ifdef OPENMP
    FFLAGS += -qopenmp
    LDFLAGS += -qopenmp
  endif
  FOPTS = -fast
  FDEBUGS = -O0 -g -traceback -DDEBUG -CB -check uninit -fpe0
  LIB_NAMES := $(filter-out lapack blas,$(LIB_NAMES)) # 通常のlapack, blasをリンクしないようにする
  LDFLAGS += -lgfortran -mkl # 代わりにmklを使う
endif

ifeq ($(FC),frtpx)
  FFLAGS += -I$(OBJDIR) -M$(OBJDIR) -Cpp $(FDEBUGS)
  ifdef OPENMP
    FFLAGS += -Kopenmp
  endif
  FOPTS = -Kfast
  FDEBUGS = -O0 -g -DDEBUG
  LDFLAGS += -Kopenmp
endif

# C
ifeq "$(CC)" "gcc"
  CFLAGS += $(CDEBUGS)
  COPTS = -O3 -funroll-loops -funswitch-loops -march=native
  CDEBUGS = -O0 -g -fbounds-check -Wuninitialized -Wunused -DDEBUG
endif

# Directory
SRCDIR = ./src
OBJDIR = ./obj


#
# 以下たぶん変更の必要なし
#
LDFLAGS += $(addprefix -L,$(LIB_PATHS)) $(addprefix -l,$(LIB_NAMES))

TARGET = $(TARGETDIR)/$(TARGETNAME)
OBJS  = $(addprefix $(OBJDIR)/, $(patsubst %.f90,%.o,$(notdir $(wildcard $(SRCDIR)/*.f90))))
OBJS += $(addprefix $(OBJDIR)/, $(patsubst %.f,%.o,$(notdir $(wildcard $(SRCDIR)/*.f))))
OBJS += $(addprefix $(OBJDIR)/, $(patsubst %.c,%.o,$(notdir $(wildcard $(SRCDIR)/*.c))))

DEPS  = $(addprefix $(OBJDIR)/, $(patsubst %.f90,%.d,$(notdir $(wildcard $(SRCDIR)/*.f90))))
DEPS += $(addprefix $(OBJDIR)/, $(patsubst %.c,%.d,$(notdir $(wildcard $(SRCDIR)/*.c))))

ifdef OPT
  FFLAGS := $(filter-out $(FDEBUGS),$(FFLAGS))
  FFLAGS += $(FOPTS)
  CFLAGS := $(filter-out $(CDEBUGS),$(CFLAGS))
  CFLAGS += $(COPTS)
endif

.PHONY: $(SUBDIRS) all clean c prof install doxygen
all: $(SUBDIRS) $(TARGET) doxygen

$(SUBDIRS):
	@$(MAKE) -C $@

$(TARGET): $(OBJS)
	@-mkdir -p $(TARGETDIR)
	$(LINK) $(OBJS) $(LDFLAGS) -o $@

# Fortran
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@-mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) -o $@ -c $<

$(OBJDIR)/%.o: $(SRCDIR)/%.f
	@-mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) -o $@ -c $<

$(OBJDIR)/%.mod: $(SRCDIR)/%.f90 $(OBJDIR)/%.o
	@:
# C
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(addprefix -I,$(INCLUDES)) $(CFLAGS) -o $@ -c $<

# SUBDIRSを再帰的にclean
clean: c
	@for SUBDIR in $(SUBDIRS); do $(MAKE) clean -C $$SUBDIR ; done

# このディレクトリだけclean
c:
	@-rm -rf $(OBJS) $(DEPS) $(TARGET) $(OBJDIR)/*.mod

# doxygen
doxygen:
	@-doxygen > /dev/null

prof: FFLAGS += -pg
prof: LDFLAGS += -pg
prof: all

install:
	@mkdir -p $(DEST)
	@cp $(TARGET) $(DEST)

-include $(DEPS)

# dependency fileのrule
$(OBJDIR)/%.d: $(SRCDIR)/%.f90
	@-mkdir -p $(OBJDIR)
	@awk -v obj=$* -v objdir=$(OBJDIR) -v srcdir=$(SRCDIR) '/^[ \t \r\f]*[uU][sS][eE][ \t \r\f]+([^ \t\n\r\f]+)/{sub(/^[ \t \r\f]*[uU][sS][eE][ \t \r\f]+/,"");sub(/[ \t \r\f]*,.*/,"");if(toupper($$0)!~/^(ISO_FORTRAN_ENV|ISO_C_BINDING|IEEE_EXCEPTIONS|IEEE_ARITHMETIC|IEEE_FEATURES|OMP_LIB|OPENACC)$$/){use[$$0];}}/^[ \t \r\f]*[mM][oO][dD][uU][lL][eE][ \t \r\f]+([^ \t\n\r\f]+)/&&!/^[ \t \r\f]*[mM][oO][dD][uU][lL][eE][ \t \r\f]+[pP][rR][oO][cC][eE][dD][uU][rR][eE][ \t \r\f]+/{sub(/^[ \t \r\f]*[mM][oO][dD][uU][lL][eE][ \t \r\f]+/,"");sub(/[ \t \r\f]+$$/,"");mod[$$0];} END{printf objdir "/" obj ".o";printf ": " srcdir "/" obj ".f90";for(key in use){printf " " objdir "/" key ".mod";}printf "\n";for(key in mod){print objdir "/" key ".mod: " srcdir "/" obj ".f90 " objdir "/" obj ".o\n";}}' $< > $(OBJDIR)/$*.d

$(OBJDIR)/%.d: $(SRCDIR)/%.c
	@-mkdir -p $(OBJDIR)
	@gcc $(addprefix -I,$(INCLUDES)) -MM -MP $< > $(OBJDIR)/$*.d
