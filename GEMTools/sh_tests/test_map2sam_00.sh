#!/bin/bash

# test file prefix
PREFIX=test.SAM.00

../bin/gt.map2sam -i ../datasets/$PREFIX.map > $TEST_DIR/result.sam

# remove the time stamp from the header
grep -v '^@CO' $TEST_DIR/result.sam > $TEST_DIR/created.sam
grep -v '^@CO' ../datasets/$PREFIX.sam > $TEST_DIR/existing.sam

# compare
diff $TEST_DIR/created.sam $TEST_DIR/existing.sam

