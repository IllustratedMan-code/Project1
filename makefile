
Project1.pdf: Project1.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='pdf_document')"

index.html: Project1.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='html_document', output_file='index.html')"

Project1.docx: Project1.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='word_document')"

all: Project1.pdf index.html Project1.docx