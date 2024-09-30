# use XeTeX as compiler for UTF-8 support
$pdf_mode = 5;
$latex = 'pdflatex %O %S';
$pdflatex = 'pdflatex %O %S';

# for bibtex, uncomment
$bibtex = 'bibtex %O %B';
$bibtex_use = 2;

$dvi_mode = 0; # disable .dvi generation
$postscript_mode = 0;	# no postscript files

@default_files = ('gradus.tex')
