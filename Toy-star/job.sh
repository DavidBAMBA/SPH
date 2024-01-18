#!/bin/bash

for n in {1..300} # Esto iterará desde 1 hasta 100
do
    echo "Ejecutando simulación con n=$n"
    make n=$n
    # Aquí puedes agregar comandos adicionales si es necesario
done
