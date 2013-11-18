#!/bin/bash

PREFIX=test.SAM.00

# run map 2 sam
../../bin/gt.map2sam -p -i ../../datasets/${PREFIX}.map > $TEST_DIR/result.sam || exit 1

# remove the time stamp
grep -v "^@CO" ../../datasets/${PREFIX}.sam > $TEST_DIR/expected
grep -v "^@CO" $TEST_DIR/result.sam > $TEST_DIR/result

diff $TEST_DIR/result $TEST_DIR/expected
