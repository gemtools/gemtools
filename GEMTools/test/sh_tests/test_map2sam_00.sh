#!/bin/bash

PREFIX=../../datasets/test.SAM.00.map

# run map 2 sam
../../bin/gt.map2sam -i ../../datasets/${PREFIX}.map > $TEST_DIR/result.sam

# remove the time stamp
grep -v "^@CO" ../../datasets/${PREFIX}.sam > expected
grep -v "^@CO" result.sam > result

diff result expected

