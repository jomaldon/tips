# Descargar archivo desde github #

1. Click the file name in a GitHub repo.

2. Click **Raw** to display the file contents.

3. Copy the URL in your browser.

4. In the command line, run either:

  *      `wget --no-check-certificate --content-disposition https://*URL-from-step3*/`
  *      `curl -LJO **https://URL-from-step3/**`

