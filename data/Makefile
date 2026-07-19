# data/Makefile: download and decompress astronomical catalogs

DATA_DIR = .
CURL = curl -L --create-dirs

.PHONY: all clean hipparcos yale hyg tycho-2

all: hipparcos yale hyg tycho-2

clean:
	rm -rf hipparcos yale-bright-star hyg tycho-2

hipparcos: $(DATA_DIR)/hipparcos/ReadMe $(DATA_DIR)/hipparcos/hip_main.dat

$(DATA_DIR)/hipparcos/ReadMe:
	$(CURL) -o $@ https://cdsarc.cds.unistra.fr/ftp/I/239/ReadMe

$(DATA_DIR)/hipparcos/hip_main.dat:
	$(CURL) -o $@ https://cdsarc.cds.unistra.fr/ftp/I/239/hip_main.dat

yale: $(DATA_DIR)/yale-bright-star/ReadMe $(DATA_DIR)/yale-bright-star/catalog.dat $(DATA_DIR)/yale-bright-star/notes.dat

$(DATA_DIR)/yale-bright-star/ReadMe:
	$(CURL) -o $@ https://cdsarc.cds.unistra.fr/ftp/V/50/ReadMe

$(DATA_DIR)/yale-bright-star/catalog.dat:
	$(CURL) -o $(DATA_DIR)/yale-bright-star/catalog.gz https://cdsarc.cds.unistra.fr/ftp/V/50/catalog.gz
	gunzip -f $(DATA_DIR)/yale-bright-star/catalog.gz

$(DATA_DIR)/yale-bright-star/notes.dat:
	$(CURL) -o $(DATA_DIR)/yale-bright-star/notes.gz https://cdsarc.cds.unistra.fr/ftp/V/50/notes.gz
	gunzip -f $(DATA_DIR)/yale-bright-star/notes.gz

hyg: $(DATA_DIR)/hyg/hyg_v44.csv

$(DATA_DIR)/hyg/hyg_v44.csv:
	$(CURL) -o $(DATA_DIR)/hyg/hyg_v44.csv.gz https://codeberg.org/astronexus/hyg/media/branch/main/data/hyg/CURRENT/hyg_v44.csv.gz
	gunzip -f $(DATA_DIR)/hyg/hyg_v44.csv.gz

tycho-2: $(DATA_DIR)/tycho-2/ReadMe $(DATA_DIR)/tycho-2/tyc2.dat

$(DATA_DIR)/tycho-2/ReadMe:
	$(CURL) -o $@ https://cdsarc.cds.unistra.fr/ftp/I/259/ReadMe

$(DATA_DIR)/tycho-2/tyc2.dat:
	@mkdir -p $(DATA_DIR)/tycho-2
	@rm -f $@.tmp
	@for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19; do \
		echo "Downloading tyc2.dat.$$i.gz..."; \
		$(CURL) -o $(DATA_DIR)/tycho-2/tyc2.dat.$$i.gz https://cdsarc.cds.unistra.fr/ftp/I/259/tyc2.dat.$$i.gz; \
		echo "Appending tyc2.dat.$$i.gz to tyc2.dat..."; \
		gunzip -c $(DATA_DIR)/tycho-2/tyc2.dat.$$i.gz >> $@.tmp; \
		rm -f $(DATA_DIR)/tycho-2/tyc2.dat.$$i.gz; \
	done
	mv $@.tmp $@
