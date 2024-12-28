Fiber-Mesh Intersection
======================

Ejecutar el archivo "interx.py" para obtener la intersección entre los fascículos con el mallado cortical. Dentro del archivo se encuentra la ruta de los fascículos y del mallado en esta carpeta, los cuales pueden ser modificados para probar con otro sujeto. El "interx.py" llama al Makefile para realizar el cálculo de intersección en C++ de forma paralela.

Luego se genera la carpeta "/intersection/" que contiene los datos de intersección de cada fascículo con cada hemisferio almacenados en archivos binarios, los cuales son:

- Índice del triángulo que intersecta con el extremo inicial de la fibra.
- Índice del triángulo que intersecta con el extremo final de la fibra.
- Punto exacto de intersección con el extremo inicial de la fibra.
- Punto exacto de intersección con el extremo final de la fibra.
- Índice de la fibra que intersecta con los triángulos.

## Referencias
<a id="1">[1]</a>
F. Silva, M. Guevara, C. Poupon, J.-F. Mangin, C. Hernandez, and P. Guevara, “Cortical Surface Parcellation Based on Graph Representation of Short Fiber Bundle Connections,” in 2019 IEEE 16th International Symposium on Biomedical Imaging (ISBI 2019), 2019.
