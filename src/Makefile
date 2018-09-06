# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
# -o3 for most optimized
# -g3.. for debugging information, -fbounds-check.. for check all array index are proper
LDLIBS= -L/disk/dawn-1/gragh/project_parameter_estimation/ANRAG.S.WRAP1/ -lANRAGfortlib
FCFLAGS = -O3 #-fbounds-check -mcmodel=large -ffpe-trap=invalid

#Param
#Put the variable modules

Param= 

# object files
OBJS = project.o

run: AN.exe


AN.exe: $(Param) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $@ $(LDLIBS)

$(OBJS): $(Param)

%.o: %.f90
	$(FC) $(FCFLAGS)  -c $< $(LDLIBS)

clean:
	rm -f *.o *~ *.mod exe *.nfs0*

