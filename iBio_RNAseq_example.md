# Sesión 5 | Ejemplo básico de cuantificación para análisis de RNAseq

#### La siguiente guía presenta los pasos básicos para obtener un archivo de cuantificación compatible con el paquete DESeq2 de R usando la herramienta Salmon.

A modo de referencia, en los siguientes enlaces pueden encontrar información y ejemplos de uso de Salmon y DESeq2

- [https://combine-lab.github.io/salmon/](https://combine-lab.github.io/salmon/)
- [https://combine-lab.github.io/salmon/getting_started/](https://combine-lab.github.io/salmon/getting_started)  
- [https://salmon.readthedocs.io/en/latest/salmon.html](https://salmon.readthedocs.io/en/latest/salmon.html)
- [http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix)
- [https://bioconductor.github.io/BiocWorkshops/rna-seq-data-analysis-with-deseq2.html](https://bioconductor.github.io/BiocWorkshops/rna-seq-data-analysis-with-deseq2.html)
- [ftp://ftp.jax.org/dgatti/MouseGen2016/DifferentialExpression.html](ftp://ftp.jax.org/dgatti/MouseGen2016/DifferentialExpression.html)

Todos los pasos que describo a continuación pueden ser ejecutados en su computador personal si cuenta con *macOS*, *Linux* (preferencia ubuntu 18.04) o *Windows con Subsistema de Linux* **(WSL)**.

## 1. Descarga e instalación de programas

Recomiendo generar una carpeta de trabajo particular para todo lo que corresponde a esta guía

    mkdir test1
    cd test1

Recomiendo que las descargas de programas queden en una carpeta llamada Downloads dentro de test1.

    mkdir Downloads
    cd Downloads

### 1.1 Descarga y habilitación de Trimmomatic (útil para trimming de secuencias por calidad y adaptadores)

    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
    unzip Trimmomatic-0.39.zip

### 1.2 Instalación de Anaconda
El programa Anaconda es un sistema que permite instalar casi sin problemas software que posee dependencias de otros programas y sin necesidad de permisos de administrador. Una vez instalado se debe "crear un ambiente" en el cual se el usuario instala los programas de interés. Es posible generar diversos ambientes lo cual permite aislar programas cuyas dependencias son incompatibles.  

>Si Ud. ya posee instalado Anaconda o Miniconda, puede proceder con la instalación del **ambiente** de trabajo que usaremos en este ejemplo (Sección 1.3).
#
>**Si Ud. utilizará el servidor proporcionado para el curso puede saltar a la Sección 1.4** pero nuestra recomendación es utilizar su propio computador o acceso personalizado a servidor y seguir toda la guía.

Los links de descarga puede encontrarlos directamente en [https://www.anaconda.com/distribution/](https://www.anaconda.com/distribution/) donde debe seleccionar *macOS* o *Linux* según corresponda. **No use la versión para Windows pues recuerde que en Windows usaremos el WSL**.

* **Usaremos la versión Python 3.7**

Con el botón derecho del ratón sobre el "botón Download" del sitio web seleccione "copiar ubicación del enlace".

Ahora, en su terminal posicionada en la carpeta Downloads ejecute el siguiente comando para descargar el archivo:


    wget <pegar dirección>

Ejemplo para Linux

    wget https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh

Ejemplo para MacOS

    wget https://repo.anaconda.com/archive/Anaconda3-2019.07-MacOSX-x86_64.sh

Para instalar el programa se debe ejecutar el siguiente comando y seguir las instrucciones en pantalla  

    bash Anaconda3-2019.07-Linux-x86_64.sh

Una vez terminada la instalación se debe responder "SI" (yes) al último mensaje con lo cual el sistema los dejará dentro del "ambiente1" de conda.

El siguiente comando previene que el sistema Anaconda se inicie de forma automática 

    conda config --set auto_activate_base false

Finalmente, **salir del sistema y volver a entrar** para que los cambios en las variables de entorno sean actualizadas. Esto significa cerrar la terminal y volver a abrirla o deslogearse del sistema y volver a entrar según el caso. [comando *exit*]


### 1.3 Instalación de los programas necesarios usando Anaconda

Ahora vamos a crear un ambiente Anaconda (conda) y lo llamaremos "ambiente1". En este ambiente instalaremos algunos paquetes base que nos permitirán trabajar con los mapeadores. La instalación de *java* es opcional pues la mayoría de los sistemas ya incluyen este programa.

    conda create -n ambiente1 -c conda-forge jemalloc icu zlib tbb libboost bzip2 libcxx
    conda install -n ambiente1 -c cyclus java-jre ## sólo requerido si su sistema no posee java

Ahora vamos a agregar a nuestro "ambiente1" los programas que usaremos en esta guía:  

- Instalamos **fastqc** y **multiqc** que nos permiten revisar la calidad de los reads.
- Distintos **mapeadores** y **cuantificadores** de uso común y otros mas nuevos.
- También instalaremos la herramienta **BUSCO** que es útil para evaluar la *calidad* de los transcriptomas de referencia.
#
    conda install -n ambiente1 -c bioconda fastqc multiqc salmon

**Como alternativa**, si desea una instalación que le permita probar y comparar diversos "mapeadores" puede ejecutar el siguiente comando (tomará más tiempo y ocupará más espacio de disco).

    conda install -n ambiente2 -c bioconda fastqc multiqc bowtie2 minimap last star hisat2 salmon kallisto busco

>Para información sobre las diferencias entre los mapeadores recomiendo **comenzar** por la discusión de este foro:
>[https://bioinformatics.stackexchange.com/questions/4507/better-aligner-than-bowtie2](https://bioinformatics.stackexchange.com/questions/4507/better-aligner-than-bowtie2)

Finalmente **activamos** el ambiente Anaconda con el siguiente comando

    conda activate ambiente1

Una vez terminado todo el trabajo podemos salir del ambiente cerrando la terminal o ejecutando el siguiente comando

    conda deactivate

## 1.4 Descarga de reads necesarios para este ejemplo

Utilizaremos los datos publicados por el siguiente artículo [https://www.ncbi.nlm.nih.gov/pubmed/26022254](https://www.ncbi.nlm.nih.gov/pubmed/26022254). Les pongo a disposición los archivos *originales* y también una versión *reducida* que les permitirá procesar los datos más rápido.

Para descargar los archivos *reducidos* se puede ejecutar los siguientes comandos

    cd test1
    mkdir raw && cd raw
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811524_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811525_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811526_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811527_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811528_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR1811529_1M.fastq.gz
    cd ..

Para descargar los archivos *originales* se puede ejecutar los siguientes comandos

    cd test1
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
    cp ~/test1/Downloads/Trimmomatic-0.39/adapters/TruSeq3-SE.fa .

El proceso de trimming se debe realizar a todos los archivos, uno a uno, usando el siguiente comando, *modificando los nombres* cada vez que corresponda donde **SRR1811524\_1M\_trimmed.fastq.gz** es el nombre que he asignado al archivo de salida (el cual lógicamente también hay que ir modificando según el caso)

    java -jar ~/test1/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 SRR1811524_1M.fastq.gz SRR1811524_1M_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:30 MINLEN:50 AVGQUAL:25

Una alternativa al paso anterior es crear un script en bash para hacer esta tarea sobre todos los archivos de forma **automática**. Para ello utilizaremos el string **1M.fastq.gz** como "palabra clave" para que bash encuentre todos los archivos deseados. Esta palabra clave la deben modificar si utilizan otras secuencias.

Pueden usar el editor **nano** para crear el archivo

    nano doTrimming.sh

Pegar en el archivo el siguiente texto, **guardar y salir**

    [ -f TruSeq3-SE.fa ] || cp ~/test1/Downloads/Trimmomatic-0.39/adapters/TruSeq3-SE.fa .
    for f in *1M.fastq.gz; do
    	java -jar ~/test1/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 $f $(echo $f | sed s/1M\.fastq\.gz/1M_trimmed\.fastq\.gz/) ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:30 MINLEN:50 AVGQUAL:25;
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

    salmon index -t databases/TAIR10_cds_20101214_updated -i At_index -p 12

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
    	cp ${fn} ${fn%_quant/quant.sf}.quant.sf
    done

> TIP: Estos comandos de bash permiten buscar y reemplazar un texto dentro de una variable:  
> %patrón/reemplazo  
> permite buscar "desde atrás" la palabra "patrón"   y la cambia por "reemplazo"
> \#patrón/reemplazo   
> permite buscar "desde el principio" la palabra "patrón" y la cambia por "reemplazo"

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



## 4. Segundo set de datos basado en *Saccharomyces cerevisiae*

Pruebe lo que ha aprendido realizand los pasos anteriores con este otro set de datos de *S. cerevisiae* cuya información podrá encontrar en el siguiente link [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA386475](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA386475)

donde:

SRR5534511 = R3 SA\_YPD  
SRR5534512 = R2 SA\_YPD  
SRR5534513 = R1 SA\_YPD  
SRR5534514 = R2 SA\_MPA  
SRR5534515 = R1 SA\_MPA  
SRR5534519 = R3 SA\_MPA  

SA = genotipo,
YPD = control,
MPA = tratamiento

    mkdir test2

Para descargar los archivos *reducidos* se puede ejecutar los siguientes comandos

    cd test2
    mkdir raw && cd raw
    wget http://genius.bio.puc.cl/genius/workshop/SRR5534511_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR5534512_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR5534513_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR5534514_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR5534515_1M.fastq.gz
    wget http://genius.bio.puc.cl/genius/workshop/SRR5534519_1M.fastq.gz

    cd ..

Para descargar los archivos *originales* se puede ejecutar los siguientes comandos

    cd test2
    mkdir raw && cd raw
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/001/SRR5534511/SRR5534511.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/002/SRR5534512/SRR5534512.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/003/SRR5534513/SRR5534513.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/004/SRR5534514/SRR5534514.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/005/SRR5534515/SRR5534515.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/009/SRR5534519/SRR5534519.fastq.gz
    cd ..

Además, descargaremos desde YeastGenome la información del genoma de *Saccharomyces cerevisiae*

    mkdir databases && cd databases
    wget https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz
    wget https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-2-1_20150113.tgz
    wget https://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff
    cd ..

> NOTA1: Estos reads tienen un tamaño de 50pb, por lo cual al hacer el trimming el largo mínimo se debe cambiar a un valor compatible, por ejemplo **MINLEN:20**
#
> NOTA2: Al crear el índice no olvide cambiar el nombre **At\_index** por uno adecuado, por ejemplo **Cs\_index** y no olvide cambiar el nombre del índice cuando realice el segundo paso de *Salmon* (mapeo). 


#por ahora... FIN

<p align="right">by Jonathan Maldonado<br>
https://github.com/jomaldon
</p>