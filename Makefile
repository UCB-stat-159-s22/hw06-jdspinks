#Makefile

html:
	jupyter-book build .
    
conf.py : _config.yml _toc.yml
	jupyter-book config sphinx .

html-hub: conf.py
	sphinx-build  . _build/html -D html_baseurl=${JUPYTERHUB_SERVICE_PREFIX}/proxy/absolute/8000
	@echo "Start the Python http server and visit:
	@echo "https://stat159.datahub.berkeley.edu/user-redirect/proxy/8000/index.html"

.PHONY: env
env: 
	mamba env create -f environment.yml 
	base -ic 'conda activate ligo;python -m ipykernel install --user --name ligo --display-name "IPython - ligo"'
	
.PHONY: clean
clean:
	rm -f figures/*.png
	rm -f audio/*.wav
	rm -rf _build/*

   