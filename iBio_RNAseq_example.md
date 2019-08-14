# Sesión 5 | Ejemplo básico de cuantificación para análisis de RNAseq

#### La siguiente guía presenta los pasos básicos para obtener un archivo de cuantificación compatible con el paquete DESeq2 de R.
Todos los pasos pueden ser ejecutados en su computador personal si cuenta con *macOS*, *Linux* (preferencia ubuntu 18.04) o *Windows con Subsistema de Linux* **(WSL)**.

## 1. Descarga e instalación de programas

Recomiendo generar una carpeta de trabajo particular para todo lo que corresponde a esta guía

    mkdir test1
    cd test1

Recomiendo que las descargas de programas queden en una carpeta llamada Downloads dentro de test1.

    mkdir Downloads
    cd Downloads

### 1.1 Descarga y habilitación de Trimmomatic (útil para trimming de secuencias por calidad y adaptadores)

    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
    unzip Trimmomatic/Trimmomatic-0.39.zip

### 1.2 Instalación de Anaconda
El programa Anaconda es un sistema que permite instalar casi sin problemas software que posee dependencias de otros programas y sin necesidad de permisos de administrador. Una vez instalado se debe "crear un ambiente" en el cual se el usuario instala los programas de interés. Es posible generar diversos ambientes lo cual permite aislar programas cuyas dependencias son incompatibles.

>Si Ud. ya posee instalado Anaconda o Miniconda, puede proceder con la instalación del **ambiente** de trabajo que usaremos en este ejemplo (Sección 1.3).

Los links de descarga puede encontrarlos directamente en [https://www.anaconda.com/distribution/](https://www.anaconda.com/distribution/) donde debe seleccionar *macOS* o *Linux* según corresponda. **No use la versión para Windows pues recuerde que en Windows usaremos el WSL**.

* **Usaremos la versión Python 3.7**

Con el botón derecho del ratón sobre el "botón Download" del sitio web seleccione "copiar ubicación del enlace".

Ahora, en su terminal posicionada en la carpeta Downloads ejecute el siguiente comando para descargar el archivo:


    wget <pegar dirección>

Ejemplo para Linux

    wget https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh

Para instalar el programa se debe ejecutar el siguiente comando y seguir las instrucciones en pantalla  

    bash Anaconda3-2019.07-Linux-x86_64.sh

Una vez terminado, se debe **salir del sistema y volver a entrar** para que los cambios en las variables de entorno actualizadas. Esto significa cerrar la terminal y volver a abrirla o, deslogearse del sistema y volver a entrar según el caso. [comando *exit*]

Una vez adentro, el siguiente comando previene que el sistema Anaconda se inicie de forma automática 

    conda config --set auto_activate_base false

Salir del sistema y volver a entrar

### 1.3 Instalación de los programas necesarios usando Anaconda

Ahora vamos a crear un ambiente Anaconda (conda) y lo llamaremos "ambiente1". En este ambiente instalaremos algunos paquetes base que nos permitirán trabajar con los mapeadores. 

    conda create -n ambiente1 -c conda-forge jemalloc icu zlib  tbb libstdcxx-ng libboost bzip2 libgcc-ng

Ahora vamos a agregar a nuestro "ambiente1" los programas que usaremos en esta guía:  

- Instalamos **fastqc** y **multiqc** que nos permiten revisar la calidad de los reads.
- Distintos **mapeadores** y **cuantificadores** de uso común y otros mas nuevos.
- También instalaremos la herramienta **BUSCO** que es útil para evaluar la *calidad* de los transcriptomas de referencia.
#
    conda install -n ambiente1 -c bioconda fastqc multiqc bowtie2 minimap last star hisat2 salmon kallisto busco 

>Para información sobre las diferencias entre los mapeadores recomiendo **comenzar** por la discusión de este foro:
>[https://bioinformatics.stackexchange.com/questions/4507/better-aligner-than-bowtie2](https://bioinformatics.stackexchange.com/questions/4507/better-aligner-than-bowtie2)

Finalmente **activamos** el ambiente con el siguiente comando

    conda activate ambiente1

Una vez terminado todo el trabajo podemos salir del ambiente cerrando la terminal o ejecutando el siguiente comando

    conda deactivate

## 1.4 Descarga de reads necesarios para este ejemplo

Utilizaremos los datos publicados por el siguiente artículo [https://www.ncbi.nlm.nih.gov/pubmed/26022254](https://www.ncbi.nlm.nih.gov/pubmed/26022254). Les pongo a disposición los archivos *originales* y también una versión *reducida* que les permitirá procesar los datos más rápido.

Para descargar los archivos *reducidos* se puede ejecutar los siguientes comandos

    mkdir raw && cd raw
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811524_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811525_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811526_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811527_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811528_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811529_1M.fastq.gz
    cd ..

Para descargar los archivos *originales* se puede ejecutar los siguientes comandos

    mkdir raw && cd raw
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR181/004/SRR1811524/SRR1811524.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR181/005/SRR1811525/SRR1811525.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR181/006/SRR1811526/SRR1811526.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR181/007/SRR1811527/SRR1811527.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR181/008/SRR1811528/SRR1811528.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR181/009/SRR1811529/SRR1811529.fastq.gz
    cd ..



Además, descargaremos desde TAIR la información del genoma de *Arabidopsis thaliana*

    mkdir databases && cd databases
    wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_blastsets/TAIR10_cds_20101214_updated
    wget https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt.gz
    wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
    cd ..

## 2. Control de calidad
### 2.1 Revisión de Calidad de los reads

Siempre se debe revisar la calidad de los reads, aún cuando los autores digan que la información entregada ya fue revisada pues, en mi experiencia, la información de los autores no siempre es confiable. Para ello utilizaremos la herramienta **fastqc**. Los resultados serán comparados en un solo reporte gracias a la herramienta **multiqc**. El modificador **-t** permite indicar el número de procesadores que se utilizarán.

    mkdir raw_fastqc
    fastqc -t 12 raw/*.gz -o raw_fastqc
    multiqc raw_fastqc -o raw_fastqc

Dentro de la carpeta **raw\_fastqc** encontrarán el archivo **multiqc\_report.html** el cual pueden abrir en cualquier explorador de internet.

Dentro de la carpeta **raw\_fastqc/multiqc\_data** encontrarán los resultados en diferentes archivos de texto.


### 2.2 Trimming de adaptadores y de segmentos de mala calidad

Asumiendo que es necesario hacer un trimming de los reads, vamos a utilizar el programa Trimmomatic.  
Primero entramos en la carpeta donde están los reads y copiamos en ella el archivo que contiene la secuencia de los adaptadores

    cd raw
    cp ~/Downloads/Trimmomatic-0.39/adapters/TruSeq3-PE.fa .

El proceso de trimming se debe realizar a todos los archivos, uno a uno, usando el siguiente comando, *modificando los nombres* cada vez que corresponda donde **SRR1811524\_1M\_trimmed.fastq.gz** es el nombre que he asignado al archivo de salida (el cual lógicamente también hay que ir modificando según el caso)

    java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 SRR1811524_1M.fastq.gz SRR1811524_1M_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:30 MINLEN:50 AVGQUAL:25

Una alternativa al paso anterior es crear un script en bash para hacer esta tarea sobre todos los archivos de forma **automática**. Para ello utilizaremos el string **1M.fastq.gz** como "palabra clave" para que bash encuentre todos los archivos deseados. Esta palabra clave la deben modificar si utilizan otras secuencias.

Pueden usar el editor **nano** para crear el archivo

    nano doTrimming.sh

Pegar en el archivo el siguiente texto, **guardar y salir**

    for f in *1M.fastq.gz; do
    	bash java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 $f $(echo $f | sed s/1M\.fastq\.gz/1M_trimmed\.fastq\.gz/) ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:30 MINLEN:50 AVGQUAL:25;
    done

Ejecutar el script de la siguiente manera

    bash doTrimming.sh

Esto realizará el trimming en todos los archivos cuyo nombre termina en **1M.fastq.gz** y generará nuevos archivos cuyo nombre termina en **1M_trimmed.fastq.gz**

### 2.3 Revisión de resultados del trimming
Una vez terminado el proceso es recomendable revisar los resultados lo cual se puede realizar de dos formas
#### 2.3.1 Módulo de multiqc para Trimmomatic
 El programa **multiqc** tiene un módulo especial para leer los archivos log producidos por Trimmomatic.

    cd ..
    mkdir trimmed_multiqc
    multiqc -m trimmomatic raw/ -o trimmed_multiqc

Dentro de la carpeta **trimmed\_fastqc** encontrarán el archivo **multiqc\_report.html** el cual pueden abrir en cualquier explorador de internet.

Dentro de la carpeta **trimmed\_fastqc/multiqc\_data** encontrarán los resultados en diferentes archivos de texto.

#### 2.3.2 fastqc+multiqc
De forma alternativa, se puede realizar el **procedimiento descrito en 2.1** para ver en detalle la calidad de los rads post-trimming

    cd ..
    mkdir trimmed_fastqc
    fastqc -t 12 raw/*trimmed.fastq.gz -o trimmed_fastqc
    multiqc trimmed_fastqc -o trimmed_fastqc

Dentro de la carpeta **trimmed\_fastqc** encontrarán el archivo **multiqc\_report.html** el cual pueden abrir en cualquier explorador de internet.

Dentro de la carpeta **trimmed\_fastqc/multiqc\_data** encontrarán los resultados en diferentes archivos de texto.

## 3. Mapeo y cuantificación de reads

Existen distintos métodos para mapear reads y cuantificar la expresión de genes según lo explicado en la clase teórica. A modo de ejemplo, utilizaremos el programa **salmon** que es similar en su etrategia a la utilizada por **kallisto**.

### 3.1. Salmon: index

Primero es necesario **indexar** el transcriptoma de referencia en un formato compatible con el programa. Note que el comando **-p** permite configurar el número de núcleos de CPU que utilizará el proceso.

    salmon index -t Downloads/TAIR10_cds_20101214_updated -i At_index -p 12

### 3.2 Salmon: mapeo y cuantificación

En el caso de salmon, el mapeo y la cuantificación se realizan en un único paso. Note que **At_index** es el directorio donde quedó el *índice* construido en el paso anterior. Este comando se debe repetir por cada librería modificando los nombres *correspondientes*

    salmon quant -i At_index -l A -r raw/SRR1811524_1M_trimmed.fastq.gz -o SRR1811524_1M_quants --validateMappings -p 12

### 3.3 Salmon: mapeo y cuantificación múltiple

El paso anterior se puede automatizar creando un script bash como el siguiente:

Pueden usar el editor **nano** para crear el archivo

    nano doSalmonQuant.sh

Pegar en el archivo el siguiente texto, **guardar y salir**

    for fn in raw/*_1M_trimmed.fastq.gz;
    do
        samp=`basename ${fn}`
        echo "Processing sample ${samp} from ${samp%_1M_trimmed.fastq.gz}"
        echo ""
        #echo "full path ${fn}"
        #echo "${samp%_1M_trimmed.fastq.gz}"
        salmon quant -i At_index -l A \
        -r ${fn}\
        -p 12 --validateMappings -o quants/${samp%_1M_trimmed.fastq.gz}_quant
    done

Ejecutar el script de la siguiente manera

    bash doSalmonQuant.sh

Una vez terminado el proceso de mapeo y cuantificación encontrará una nueva carpeta llamada **quants**. 

    cd quants && ls -lt

Dentro de ella encontrará una carpeta por cada archivo de entrada. Dentro de cada **subcarpeta** habrá un archivo **quant.sf** el cual contiene toda la cuantificación del archivo de entrada respectivo. Es un archivo de texto tabular.

Puede ver el encabezado de alguno de ellos con el comando *head*

    head SRR1811524_quant/quant.sf

### 3.3 Juntar los resultados en una sola tabla

Más adelante necesitaremos cargar los datos de cuantificación en R y así utilizar el paquete *DESeq2*. 
Para ello primero debemos generar una **tabla consolidada** de la cuantificación de las muestras.

El primer paso es renombrar cada archivo quant.sf a **sample**_quant.sf (donde **sample** es el nombre de la librería), y luego poner todos los archivos en un mismo directorio lo cual se puede hacer uno a uno o con el script de bash siguiente

    cd quants

    nano doRenameAndMove.sh

Pegar en el archivo el siguiente texto, **guardar y salir**

    for fn in */*.sf;
    do
        echo "Input=${fn} | Output=${fn%_quant/quant.sf}.quant.sf"
        cp ${fn} ${fn%\_quant/quant.sf}.quant.sf
    done

Ejecutar el script de la siguiente manera

    bash doRenameAndMove.sh

    ls -ltr

Si todo salió bien, tendremos en el directorio de trabajo quants una lista de archivos .sf cuyo nombre indica la librería de origen

Finalmente, vamos a unir todos los archivos utilizando un script en Python que he adaptado para los archivos de cuantificación tipo *salmon* y que pueden descargar desde mi github

Los comandos son los siguientes:
 
    curl -LJO https://raw.githubusercontent.com/jomaldon/scripts_bioinfo/master/merge_quants.py
    python merge_quants.py *.sf > At.separated.quants.txt

    head At.separated.quants.txt

Si todo sale bien, tendremos un archivo tabular (At.separated.quants.txt), cuya primera columna es el nombre del gen y las columnas consecutivas serán la cuantificación de dichos genes en las diferentes librerás.

#por ahora... FIN


<p align="right">by Jonathan Maldonado<br>
https://github.com/jomaldon
</p>