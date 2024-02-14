#!/bin/bash
for n in $(seq 1 1 300) # Esto iterará desde 3 hasta 300 en pasos de 3
do
    echo "Ejecutando simulación con n=$n"
    make n=$n
    # Aquí puedes agregar comandos adicionales si es necesario
done
