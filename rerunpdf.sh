#!/bin/bash
# Delete everything except the .tex file, then runs pdflatex and bibtex the amount of
# times necessary to produce a proper pdf file.
usage="
$(basename "$0") [-h] [-d] -- program to compile your pdf file
out of a LaTeX program. The syntax we use as default is the following:

>  pdflatex -synctex=1 -interaction=nonstopmode filename

Note that there is no '.tex' at the end of <filename>.
Please don't use it.
"
if [ -z $1 ]; then
    echo "$usage"
    exit 1
fi
debugmode=0
while getopts ':dh' option; do
  case "$option" in
    h) echo "$usage"
    exit
    ;;
    d) debugmode=1
    echo "DEBUG MODE ON"
    ;;
    \?) echo "$usage" >&2
    exit 1
    ;;
  esac
done
shift $((OPTIND-1))
mv $1.tex ../
rm $1.*
mv ../$1.tex ./

if [ "$debugmode" == "0" ]; then
    pdflatex -synctex=1 -interaction=nonstopmode $1
    # Uncomment both lines if you have a .bib file available in this directory.
    bibtex $1
    pdflatex -synctex=1 -interaction=nonstopmode $1
    pdflatex -synctex=1 -interaction=nonstopmode $1
else
    pdflatex -synctex=1 $1
    # Uncomment the two following lines if you have a .bib file
    #available in this directory.
    bibtex $1
    pdflatex -synctex=1 $1
    pdflatex -synctex=1 $1
fi

clear

printf  "\nCompiling done!"
printf "\nIf you had a lot of errors, you can try the -d (debug) option.\n"
