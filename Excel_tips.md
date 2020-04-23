# Algunos tips que he recopilado para Microsoft Excel #

**1. Comando para buscar metadata específica de un ID (C1) asumiendo que la información de todos los IDs posibles está en el mismo archivo. Una columna tiene los IDs ordenados (Columna A) y la columna contigua tiene la información asociada a cada ID (Columna B)**

    =BUSCAR(Celda_con_ID_a_buscar;Columna_con_IDs;Columna_con_info)
	=BUSCAR(C1;A:A;B:B)

**2. Unir elementos de celda A1 y B1 separadas por un espacio**

    =A1&" "&B1

## Recomendaciones usando el paquete Real Statistics [http://www.real-statistics.com/](http://www.real-statistics.com/ "Real Statistics Using Excel") ##

**1. Mann Whitney Test:**
[http://www.real-statistics.com/non-parametric-tests/mann-whitney-test/](http://www.real-statistics.com/non-parametric-tests/mann-whitney-test/)

Ejemplo:

 Conocer el **pvalue** de una comparación usando MWTEST con 2 colas (tails), correción de vínculos (ties) y corrección de continuidad (cont) [parámetros default].

    =MWTEST(A1:C1;D1:F1;2;1;1) # con ambas correcciones
    =MWTEST(A1:C1;D1:F1;2;0;0) # sin correcciones

También se puede usar el MW\_EXACT(A1:C1;D1:F1;colas) y el MW\_TEST(resultado full) pero lo recomendado es MWTEST

Explicación rápida sobre uso de 1 o 2 colas
https://help.xlstat.com/s/article/cual-es-la-diferencia-entre-una-prueba-de-dos-colas-bilateral-y-de-una-cola-unilateral?language=es

...o en este video
https://es.coursera.org/lecture/datos-tecnicas/pruebas-de-hipotesis-de-dos-colas-y-de-una-cola-5bWWz


<p align="right">by Jonathan Maldonado<br>
https://github.com/jomaldon
</p>
