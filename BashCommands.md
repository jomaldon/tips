# Lista de Comandos útiles para la terminal de Linux/BASH
#### En esta lista encontrarán una recopilación de comandos que he guardado en mi experiencia con la terminal de Linux/Bash. La mayoría funcionarán en cualquier distribución de Linux (incluso Unix), pero unos cuantos pueden ser particulares de Ubuntu.

>En los programas de la terminal, el caracter **-** se usa para agregar "modificadores" para los programas o también conocidos como "opciones".


## Ejemplos para listar contenido
Lista el contenido de un directorio

    ls

Lista el contenido mostrando fecha y hora

    ls -l

Lista el contenido mostrando el tamaño en "unidades human-readable" 

    ls -lh

Ordena los archivos según fecha, pone el mas reciente arriba

    ls -lht

Ordena según tamaño de los archivos; mayor a menor

    ls -lhS

Ordena según tamaño de los archivos; menor a mayor

    ls -lhSr

Ordena el contenido por fecha; la mas reciente queda abajo [mi favorito]

    ls -lhtr

Ordena el contenido por tamaño de los archivos; el más grande queda abajo

    ls -lhSr

Ayuda del comando. El modificador *--help* es un comodín que sirve para lo mismo en la mayoría de los comandos

    ls --help

Muestra la info de ayuda pero en secciones. Al presionar *espacio* se muestra la sección siguiente. Al presionar la tecla *q* se sale.

    ls --help | more

Otra forma de obtener ayuda es invocar el comando *man* que invoca al manual del programa que uno desee.

    man ls

>En la propiedades de un archivo se ve algo similar lo siguiente: **drwxrwxrwx**  
>Corresponde a los **permisos** que tienen los usuarios sobre los archivos.  
Se dividen en 4 bloques: d rwx rwx rwx  
El primer bloque indica si el archivo es un directorio  
El segundo bloque son los permisos del dueño del archivo  
El tercer bloque son los permisos para el grupo  
El cuarto bloque son los permisos para todos los usuarios  
Las letras significan lo siguiente:
>>d: directorio, r: lectura, w: escritura, x: ejecutar. 

## Otros comandos de uso general

Mostrar en qué carpeta estoy trabajando

    pwd 

**Crear** directorio

    mkdir nombre_directorio

**Entrar** en directorio (cambiar de directorio

    cd nombre_directorio

**Salir** del directorio (cambia al directorio anterior)

    cd ..

**Borrar** directorio vacío

    rmdir nombre_directorio

**Borrar** archivo

    rm nombre_archivo

**Borrar** directorio con archivos en su interior **[usar con precaución!!]**

	rm -r nombre_directorio

**Mover** archivo desde el directorio1 al directorio2

    mv nombre directorio1/nombre_archivo nombre_directorio2/nombre_archivo

**Renombrar** archivo desde archivo1 a archivo2

    mv nombre directorio1/nombre_archivo nombre_directorio2/nombre_archivo

**Copiar** archivo a un directorio nuevo

    cp nombre_archivo nombre_directorio_destino
Mostrar todo el contenido de un archivo (sólo hacerlo con archivos de texto)

    cat nombre_archivo

Mostrar partes del contenido de un archivo. Se avanza presionando *espacio* y se sale presionando tecla *q*

    more nombre archivo

Mostrar el **comienzo** del archivo. Por default muestra las primeras 10 líneas pero si se usa el modificador *-n* se puede especificar el número de líneas

    head nombre_archivo

Mostrar el **final** del archivo. Por default muestra las últimas 10 líneas pero si se usa el modificador *-n* se puede especificar el número de líneas

    tail nombre_archivo
  
Mostrar **todo** el contenido del archivo. Sólo recomendable cuando se usan archivos con poco contenido.

    cat nombre_archivo

Encontrar algún archivo. Al usar el modificador -size se puede especificar un tamaño

    find nombre_archivo

Limpiar la terminal

    clear

Buscar un texto dentro de un archivo de texto

    grep texto nombre_archivo



Mostrar variable de entorno

    echo $nombre_de_variable
    echo $PATH

Buscar explicación *simple* de algún comando

    apropos comando  

## Editores de texto desde la línea de comando

En bioinformática generalmente debemos revisar archivos de texto grandes. Esto es un desafío para un editor de texto gráfico pues se requieren muchos recursos para desplegarlo. Sin embargo, es posible revisar e incluso editar estos archivos desde la línea de comandos usando mucho menos recursos.
Si sólo se desea revisar el contenido lo recomendable es usar comandos como **head**, **tail** o **more** explicados anteriormente.
Por otro lado, si deseamos editar texto simple usando la línea de comandos lo típico es usar **vi** o su versión algo mas amigable **vim**. El problema es que no son amigables y se requiere práctica para dominar su uso.  
Una altenativa mucho mas amigable es **nano** (un derivado del editor **pico**) el cual recomiendo para quienes estén iniciando en el uso de bash.  

>La gracia de **vi** es que siempre viene instalado en cualquier distribución de linux (y unix) mientras que nano generalmente se debe instalar. Por ello recomiendo eventualmente aprender a usar **vi**

En Ubuntu estos son los comandos para instalar nano
 
    sudo apt-get update  
    sudo apt-get install nano  

Para usarlo simplemente se ejecuta

    nano nombre_archivo.txt

Ejemplo de pantalla del programa *nano*
![Nano example](images/Nano.PNG "Nano Example")

Los controles aparecen abajo y se ejecutan con la tecla ***control*** y luego la letra del comando deseado. Por ejemplo, salir es **Control+X**.

## Editores de texto simple usando interfáz gráfica

Si la distribución de linux utilizada posee interfaz gráfica y/o el sistema está configurado para ejecutar programas la interfaz entonces una opción simple es **mousepad** (liviano) o **gedit**.


## Ejemplo de comparación de listas
Si A.txt y B.txt son archivos de texto con información que deseo comparar  

ejemplo de archivo A.txt

    cat A.txt  
>sec1  
sec2  
sec3  
sec4  

ejemplo de archivo B.txt  

    cat B.txt
>sec2  
sec3  
sec5  
sec6  

Los siguientes comandos permiten distintos tipos de comparaciones

    cat A B | sort | uniq -c |  grep  ' 2 ' > DUPLICADOS  

    cat A B | sort | uniq -c | awk '{if($1>1) print $2}' > INT  
    cat A B | sort | uniq -c | awk '{if($1==1) print $2}' > UNICOS  

    cat A INT | sort | uniq -c | awk '{if($1==1) print $2}' > UNICOS_A  
    cat A UNICOS | sort | uniq -c | awk '{if($1>1) print $2}' > UNICOS_A  

Donde:  
DUPLICADOS = archivo que contendrá los datos duplicados entre A y B  
INT = archivo que contendrá intercepción de los archivos A y B  
UNICOS =  archivo que contendrá los datos que son únicos entre A y B  
UNICOS_A = archivo que contendrá los datos que “sólo” pertenecen al archivo A  



## Modificadores
**\*** = Cualquier texto  
**?** = un carácter cualquiera  
**\.** = también puede funcionar como 1 caracter cualquiera en expr reg  
**\[\<caracteres\>\]** = cualquiera de los caracteres ingresados  
**comando1 | comando2**	= tubería envía el resultado de un comando1 a otro  
**;** = Separa dos comandos, como si fuesen independientes  



## Comandos útiles de sistema

Mostrar **discos** y sus propiedades 

    df -h  
 
Mostrar las propiedades del disco al que pertenece la carpeta de trabajo actual

    df -h .  

Mostrar **info** del de sistema instalado  

    uname -a  

Mostrar **versión** del sistema operativo (sólo en Ubuntu)

    lsb_release -a  

Memoria **RAM**  

    grep MemTotal /proc/meminfo  

Memoria **swap**

    grep SwapTotal /proc/meminfo  
  
Mostrar **particiones** de discos (detalles) 

    sudo fdisk -l
 
Datos sobre memoria ram y swap usada y libre  

    free -h  

Mostrar info de la CPU

    grep "model name" /proc/cpuinfo  

Mostrar info del último reboot

    last reboot  
	  
Mostrar el tiempo de uso

    uptime  

Apagar o reiniciar el sistema. 0: apaga el equipo, 6: reinicia el sistema  

    init <numero>  


## Administración de procesos
Mostrar los procesos activos y su codigo PID

    ps -A  

Mostrar los procesos que posean una palabra clave <texto>
    ps -A | grep <texto>  
    ps -A | grep cron

Mostrar los procesos activos, ordenado según uso de RAM de menor a mayor. Lo último que salga en la lista lo que usa mayor RAM

    ps -e -orss=,args= | sort -b -k1,1n | pr -TW$COLUMNS  
    ps -e -orss=,args= | sort -b -k1,1n | pr -TW$COLUMNS | tail  

Mostrar la dependencia de los procesos, como un árbol

    pstree  

Cerrar el proceso de codigo <PID>

     kill -9 <PID>  

Información a tiempo real del uso de los procesos en el sistema

    top  

Crear un alias como nombrealias del comando en archivo .bashrc o .bash_aliases

    alias <nombrealias>='<comando>'  
    alias micomando = ls -ltSr

Crear un acceso a directorio simbólico

    ls -s <directorio> <newpath>  

Elimnar acceso directo

    unlink <directorio> <newpath>  

Ejecutar comando en background (independiente de la consola)

    <comando> &  
  
Abortar proceso un proceso que está en ejecución

    CRTL + c  

Pausa el proceso actual  

    CRTL + z  

Muestra todos los procesos parados

    jobs  
	  
Continúa en primer plano el proceso detenido de número **N** según información obtenida con comando **jobs** muestra  

    fg N  

Continúa en segundo plano el proceso detenido de número **N**

    bg N  

## SCRIPTING

Entregar atributos al archivo para ser ejecutable

    chmod +x script.sh

Ejecutar el script en el archivo <file.sh

    bash <file.sh>

Mostrar valor de variable en consola

    echo $<Variable>

Ejecutar repetidamente el mismo comando cada 2s (default)

    watch <comando>

Leer la siguiente línea que se escriba y la almacena en una variable de entorno

    read mi_variable  
        esto es un ejemplo (presionar enter)
    echo $mi_variable

Generar una lista de números por ejemplo de 1 a 100

    seq 1 100

Almacenar el resultado de comando en una variable de entorno

    mi_variable=$(comando)
    mi_variable=$(ls)
    echo $mi_variable

Expresion matemática + - * / sin decimales % decimales ( ) anidar

    $((expresión))
    $(( (4-1)*3 ))

## REDES

**Conectar** con servidor mediante ssh

    ssh nombre_usuario@IP_servidor
    ssh nombre_usuario@dominio_servidor

**Copiar** archivo o carpeta a servidor mediante ssh

    scp archivo_local nombre_usuario@IP_servidor:carpeta_destino
    scp -r carpeta_local nombre_usuario@IP_servidor:carpeta_destino

**Traer** archivo o carpeta desde servidor mediante ssh

    scp nombre_usuario@IP_servidor:archivo_remoto .
    scp -r nombre_usuario@IP_servidor:carpeta_remota .

**Copiar** archivos o carpetas usando interfaz gráfica. Recomiendo el uso de la aplicación Filezilla.

    sudo apt install filezilla
    filezilla

**Establecer puente** hacia servidor interno pasando por servidor externo usando el puerto XX (por ejemplo, 1000). Suponemos que el servidor externo se conecta con el *usuario\_ex* y el servidor interno se conecta con el *usuario\_in*. Ojo que el firewall correspondiente debe permitir el paso de datos por el puerto escogido.

    ssh -L XX:IP_servidor_interno:22 usuario_ex@IP_servidor_externo -N -f

**Usar puente** para conectarse a servidor interno o para copiar archivo (ojo que la letra p es diferente). La premisa es la misma que para una conexión ssh normal pero debemos usar el puerto adecuado (XX) y la IP interna (127.0.0.1) que dirige el puente al lugar correcto.

    ssh -p XX usuario_in@127.0.0.1 # p minúscula
    scp -P XX usuario_in@127.0.0.1:carpeta_destino archivo # P mayúscula

**Cerrar puente** obteniendo el PID y luego usando kill

    ps aux | grep ssh # para ver el proceso y rescatar el PID
    kill -2 PID

**Cerrar puente** terminando todos los procesos ssh

    killall ssh

## Búsqueda de archivos

Buscar archivo segun tag

    find -name <tag>

Buscar archivos de mayor tamaño, con -<tamaño> de menor

    find -size +<tamaño>
    find -size +1M

Ejecutar un comando en los archivos encontrados que poseen el tag. El comando se ejecuta tantas veces como archivos encontrados.

    find / -name <tag>  -exec <comando> {} \;  

Ejecutar un comando en los archivos encontrados que poseen el tag. El comando se ejecuta sólo una vez considerando todos los archivos encontrados.

    find / -name <tag>  -exec <comando> {} \+

## Archivos y Directorios

Tamaño de directorio y sub directorio

    du -s -h  <nombredirectorio>
    du -sh  <nombredirectorio> # es lo mismo

## Compresión
Descomprimir el contenido de un archivo tipo tgz en la carpeta actual

    tar -xfzv <archivo.tgz> 

Comprimir un directorio en el archivo salida.tgz

    tar -cfzv <salida.tgz> <directorio> 

Mostrar el contenido del archivo

    tar -tzvf <archivo.tgz>

Dividir el archivo en varios archivos de tamaño “bytes” (no borra el original)

    split -b <bytes> <file>

Unir los archivos que responden al tag, antes separados con split

    cat <tag>* > <file>


## EDICION DE TEXTO

Cambiar las mayusculas a minusculas del texto dentro file

    tr 'A-Z' 'a-z' < <file>

Agregar el texto en la última línea de file

    echo "<text>" >> <file> 

Cambiar las letras de un texto por otras por ejemplo “a” por “A”

    echo <text> | sed 's/a/A/g'

## Filtro de información

Mostrar las diferencias de ambos archivos

    diff fichero1 fichero2 

Mostrar las coincidencias

    join fichero1 fichero2 

Mostrar las líneas que contengan la solicitud tag dentro del archivo file o de los archivos en el directorio

    grep <tag> <file/directorio> 

Idéntico pero en todos los archivos o directorios

    grep <tag> <file/directorio2> <file/directorio3> <file/directorio4>

>grep *opción* \<tag\> \<file/directorio\>  
	-i 	Ignora mayusculas o minusculas al realizar la búsqueda  
	-r	Búsqueda en todos los subdirectorios  
	-c	Cuenta cuantas veces se encontró el tag  
	--color	Destaca de color la palabra encontrada  
	-n 	Indica la linea donde esta la palabra  
	-v	entrega lo que NO hizo match  
	-f <file>	usa como patrón de búsqueda la info contenida en el archivo file  
	-w	búsqueda de la palabra completa (ie, contig1 es diferente de contig10)  

Busqueda de ambas palabras

    egrep -w 'word1|word2' <file>

Buscar tag pero excluye otros_tags

    grep <tag> <file> | grep -v <otros_tags>

Mostrar las líneas alrededor del tag

    grep –context=6 <tag> <file>

Mostrar las líneas que comienzan por tag

    grep '^<tag>' <file>

Ordenar archivo de texto

    sort <file>


<p style="text-align: right;">by Jonathan Maldonado<br>
[https://github.com/jomaldon](https://github.com/jomaldon "GitHub jomaldon")

<p align="right">by Jonathan Maldonado<br>
[https://github.com/jomaldon](https://github.com/jomaldon "GitHub jomaldon")

<p style="text-align: right;"> by Jonathan Maldonado</p>

<p align="right">by Jonathan Maldonado</p>