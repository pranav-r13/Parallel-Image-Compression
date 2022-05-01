SEQ_MPI_BIN=seq-mpi-bin
OMP_BIN=omp-bin
SEQ_MPI_OBJDIR=objs-seq-mpi
OMP_OBJDIR=objs-omp
IMGDIR=images
COMPDIR=compressed
PNGDIR=lodepng
SEQ_MPI_CXX=mpic++ -m64
OMP_CXX=g++ -m64 -fopenmp
CXXFLAGS=-O3 -std=c++11

SEQ_MPI_OBJS=$(SEQ_MPI_OBJDIR)/seq-mpi.o $(SEQ_MPI_OBJDIR)/$(PNGDIR)/lodepng.o $(SEQ_MPI_OBJDIR)/dct.o $(SEQ_MPI_OBJDIR)/image.o $(SEQ_MPI_OBJDIR)/quantize.o $(SEQ_MPI_OBJDIR)/rle.o $(SEQ_MPI_OBJDIR)/dpcm.o
OMP_OBJS=$(OMP_OBJDIR)/omp.o $(OMP_OBJDIR)/$(PNGDIR)/lodepng.o $(OMP_OBJDIR)/dct.o $(OMP_OBJDIR)/image.o $(OMP_OBJDIR)/quantize.o $(OMP_OBJDIR)/rle.o $(OMP_OBJDIR)/dpcm.o


.PHONY: default dirs clean

default: dirs $(SEQ_MPI_BIN) $(OMP_BIN)

dirs:
		mkdir -p $(SEQ_MPI_OBJDIR) $(SEQ_MPI_OBJDIR)/$(PNGDIR) $(OMP_OBJDIR) $(OMP_OBJDIR)/$(PNGDIR) $(IMGDIR) $(COMPDIR)

clean:
		rm -rf $(SEQ_MPI_OBJDIR) $(OMP_OBJDIR) $(IMGDIR) $(COMPDIR) *~ $(SEQ_MPI_BIN) $(OMP_BIN) $(LOGS)

$(SEQ_MPI_BIN): $(SEQ_MPI_OBJS)
		$(SEQ_MPI_CXX) $(CXXFLAGS) -o $@ $(SEQ_MPI_OBJS) $(OBJS) $(LDFLAGS)

$(OMP_BIN): $(OMP_OBJS)
		$(OMP_CXX) $(CXXFLAGS) -o $@ $(OMP_OBJS) $(OBJS) $(LDFLAGS)

$(SEQ_MPI_OBJDIR)/%.o: src/%.cpp
		$(SEQ_MPI_CXX) $< $(CXXFLAGS) -c -o $@

$(SEQ_MPI_OBJDIR)/$(PNGDIR)/%.o: src/$(PNGDIR)/%.cpp
		$(SEQ_MPI_CXX) $< $(CXXFLAGS) -c -o $@

$(OMP_OBJDIR)/%.o: src/%.cpp
		$(OMP_CXX) $< $(CXXFLAGS) -c -o $@

$(OMP_OBJDIR)/$(PNGDIR)/%.o: src/$(PNGDIR)/%.cpp
		$(OMP_CXX) $< $(CXXFLAGS) -c -o $@

