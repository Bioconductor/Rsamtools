#!/bin/bash

TEST_FILE=$1

R --slave <<-EOF
	library('RUnit')
	library('Rsamtools')
	res <- runTestFile('unit/${TEST_FILE}',
		rngKind='default', rngNormalKind='default')
	printTextProtocol(res, showDetails=FALSE)
EOF
