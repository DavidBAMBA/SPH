# Definir el compilador
CC=g++ 

# Definir el nombre del archivo ejecutable
EXECUTABLE=a.out

# Definir el script de Python
PYTHON_SCRIPT=/home/yo/Documents/Tesis/codes/SPH/shock-Tube/c++/plot_mesh-2d.py

# La regla 'all' se ejecutará por defecto si solo se escribe 'make' en la terminal.
# Agregamos un argumento 'n' que se puede pasar al llamar a make.
all: $(EXECUTABLE)
	./$(EXECUTABLE) $(n)
	/bin/python3 $(PYTHON_SCRIPT) $(n)

# Regla para construir el ejecutable
$(EXECUTABLE): shock-tube-2d.cpp
	$(CC) shock-tube-2d.cpp -o $(EXECUTABLE)

# Regla para limpiar los archivos compilados (opcional)
clean:
	rm -f $(EXECUTABLE)
