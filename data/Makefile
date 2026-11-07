# data/Makefile: download and decompress astronomical catalog

DATA_DIR = .
CURL = curl -L --create-dirs

.PHONY: all clean hyg

all: hyg

clean:
	rm -rf hyg

hyg: $(DATA_DIR)/hyg/hyg_v44.csv

$(DATA_DIR)/hyg/hyg_v44.csv:
	$(CURL) -o $(DATA_DIR)/hyg/hyg_v44.csv.gz https://codeberg.org/astronexus/hyg/media/branch/main/data/hyg/CURRENT/hyg_v44.csv.gz
	gunzip -f $(DATA_DIR)/hyg/hyg_v44.csv.gz
