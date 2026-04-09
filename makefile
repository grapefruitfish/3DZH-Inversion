# ---------------- 配置 ----------------
CMD      ?= 3DZHTomo
FC       ?= gfortran
OPENMP   ?= -fopenmp

ROOTDIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

SRCDIR  := $(ROOTDIR)/src
OBJDIR  := $(ROOTDIR)/obj
MODDIR  := $(OBJDIR)/mod
BINDIR  := $(ROOTDIR)/bin

# 通用编译选项
FFLAGS   ?= -O -ffloat-store -g -fbacktrace -Wall
#FFLAGS   ?= -O0 -g -fcheck=all -fbacktrace -Wall -Wextra \
#         -finit-real=snan -finit-integer=-999999
# Fortran77 fixed-form
FFLAGS_F := $(FFLAGS) -std=legacy -ffixed-form -ffixed-line-length-none

# module 输出/搜索路径
MODFLAGS := -J$(MODDIR) -I$(MODDIR)

# ---------------- 源文件 ----------------
FSRCS  = zhFwd.f sdisp96_nodat.f surfdisp96.f

F90SRC = gaussian.f90 Def_zh.f90 Def_model.f90 \
         aprod.f90 lsmrDataModule.f90 lsmrblasInterface.f90 \
         lsmrblas.f90 lsmrModule.f90 Read_Data.f90 \
         Def_variable.f90 utils.f90 CalZHsG.f90 main.f90

# 对象文件（统一放入 OBJDIR）
OBJS = $(addprefix $(OBJDIR)/, \
        $(FSRCS:%.f=%.o) \
        $(F90SRC:%.f90=%.o))

.PHONY: all clean prepare

# 默认目标
all: prepare $(BINDIR)/$(CMD)

# 目录准备
prepare:
	@mkdir -p $(BINDIR) $(OBJDIR) $(MODDIR)

# ---------------- 链接 ----------------
$(BINDIR)/$(CMD): $(OBJS)
	$(FC) $(OPENMP) $^ -o $@

# ---------------- 编译规则 ----------------

# f90 -> o （注意 SRCDIR）
$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | prepare
	$(FC) $(OPENMP) $(FFLAGS) $(MODFLAGS) -c $< -o $@

# f -> o （fixed form）
$(OBJDIR)/%.o: $(SRCDIR)/%.f | prepare
	$(FC) $(OPENMP) $(FFLAGS_F) $(MODFLAGS) -c $< -o $@

# ---------------- 清理 ----------------
clean:
	rm -rf $(OBJDIR) $(BINDIR)
