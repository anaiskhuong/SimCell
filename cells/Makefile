### LINUX
simcell  : driver.o domain.o main/data.hpp main/geom.hpp main/Vec.hpp main/Mat.hpp
	c++ -ggdb main/main.cpp -o simcell driver.o domain.o main/potter.o main/params.o main/models.o main/viewer.o main/shaker.o main/keeper.o main/linker.o main/bumper.o main/sorter.o main/fixers.o main/cell.o main/data.o main/geom.o main/util.o -lGL -lGLU -lglut -lpthread -lm 

### MAC
### simcell  : driver.o domain.o main/data.hpp main/geom.hpp main/Vec.hpp main/Mat.hpp
###	c++ -ggdb main/main.cpp -o simcell driver.o domain.o main/potter.o main/params.o main/models.o main/viewer.o main/shaker.o main/keeper.o main/linker.o main/bumper.o main/sorter.o main/fixers.o main/cell.o main/data.o main/geom.o main/util.o -lm -framework OpenGL -framework GLUT -L/opt/local/lib  -lIL -lILU -lILUT



driver.o : driver.cpp
	c++ -c -ggdb driver.cpp

domain.o : domain.cpp
	c++ -c -ggdb domain.cpp
