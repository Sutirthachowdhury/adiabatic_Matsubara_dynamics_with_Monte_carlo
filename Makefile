FC=gfortran
FLAGS = -lfftw3 -llapack -lblas -fbounds-check #-DDEBUG -D"ASSERT(x)=call assert(x, 'x', __LINE__, __FILE__)"           

SRC= global.f90 gran.f90 parameter_input.f90 initialize_mat_pos.f90 nm_t.f90 momentum_sample.f90 init_pot.f90 forces.f90 forces_t0.f90 evolve.f90 transform_matrix.f90 back_transform_matrix.f90 monte_carlo.f90 

OBJS  =  ${SRC:.f90=.o}

all:  $(OBJS)
	$(FC) $(OBJS) mat_1d_quartic.f90 -o mat_1d.exe $(FLAGS)


%.o : %.f90
	$(FC) -c $(FLAGS) $< -o $@

clean:
	 rm -rf *.o
