CXX=mpiicpc
CXXFLAGS+=-I./include
CXXFLAGS+=-std=c++11
LIBPATH+=
vpath %.cpp src
vpath %.h include
vpath %.o obj
INJECTION: calculatecorrelation.o interface.o
	mkdir -p obj bin
	$(CXX) -o CORRELATION $(LIBPATH) calculatecorrelation.o interface.o
	mv *o bin/
%.o:%.c $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $^
clean:
	rm -rf bin obj *.o
