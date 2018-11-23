
DATA = ../data/2018-10-20-orm_measurements/M8-E108-F1-I9-G1_

model: model_base model_params model_ealign
plots:            plots_params plots_ealign

model_base:
	python model_orbits.py -i init.txt -e base.list

model_params:
	python model_orbits.py -i init.txt -e params/params.pos.list
	python model_orbits.py -i init.txt -e params/params.neg.list
	python model_orbits.py -i init.txt -e params/params.inv.list

model_attrs:
	python model_orbits.py -i init.txt -e params/attrs.pos.list
	python model_orbits.py -i init.txt -e params/attrs.neg.list
	python model_orbits.py -i init.txt -e params/attrs.inv.list
	python model_orbits.py -i init.txt -e params/attrs.off.list

model_ealign:
	python model_orbits.py -i init.txt -e params/ealign.pos.list
	python model_orbits.py -i init.txt -e params/ealign.neg.list

plots_params:
	python plot_folders.py -b $(DATA)/base/model_0 -m $(DATA)/ $(DATA)/params.pos/*(/)
	python plot_folders.py -b $(DATA)/base/model_0 -m $(DATA)/ $(DATA)/params.neg/*(/)
	python plot_folders.py -b $(DATA)/base/model_0 -m $(DATA)/ $(DATA)/params.inv/*(/)

plots_ealign:
	python plot_folders.py -b $(DATA)/base/model_0 -m $(DATA)/ $(DATA)/ealign.pos/*(/)
	python plot_folders.py -b $(DATA)/base/model_0 -m $(DATA)/ $(DATA)/ealign.neg/*(/)
