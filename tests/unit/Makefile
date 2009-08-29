TOP=../../..
RUNIT_DIR=$(shell pwd)
INST_PATH=$(shell dirname $(RUNIT_DIR))
PKG_PATH=$(shell dirname $(INST_PATH))

all: install test

full: build install test

install:
	cd ${TOP}; \
	R CMD INSTALL --clean ${PKG_PATH} && \
	cd ${RUNIT_DIR}

build:
	cd ${TOP}; \
	R CMD build ${PKG_PATH} && \
	cd ${RUNIT_DIR}

test:
	cd ..; \
	R --slave < runalltests.R

%.R: FORCE
	cd ..; \
	./runfile.sh $@

FORCE:
