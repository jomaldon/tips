# Comandos para administrar repositorio Git desde ubuntu 16.04 en Windows10 #
## Configurar credenciales github ##
    git config --global user.name "username"
    git config --global user.email "NN@users.noreply.github.com"

## Entrar en la carpeta local ##
    cd ~/Documents/git

## Iniciar un repositorio nuevo localmente ##
    mkdir tips
    git init tips

## Hacer cambios en archivos ##
    cd tips
    vi README.md

## Añadir o eliminar del índice ##
    git add README.md
    git rm README.md

## borrar del índice dejando la versión local ##

    git rm --cached README.md

## Publicar cambios añadiendo un comentario ##
    git commit -m "this is a test"

## Vincular repositorio remoto con el local (nombre exacto) ##
    git remote add origin https://github.com/jomaldon/tips.git
    git remote add origin https://github.com/jomaldon/scripts_bioinfo

## Bajar web a local ##
    git pull origin master

## Subir local a web ##
    git push origin master

## Turn off showUntrackedFiles ##
    git config --global status.showUntrackedFiles no

### test ###


<p align="right">by Jonathan Maldonado<br>
https://github.com/jomaldon
</p>