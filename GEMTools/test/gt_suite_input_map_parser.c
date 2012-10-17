/*
 * PROJECT: GEM-Tools library
 * FILE: gt_suite_input_map_parser.c
 * DATE: 03/10/2012
 * DESCRIPTION: // TODO
 */

#include "gt_test.h"

gt_map* map;

void gt_input_map_parser_setup(void) {
  map = gt_map_new();
}

void gt_input_map_parser_teardown(void) {
  gt_map_delete(map);
}

START_TEST(gt_test_imp_string_map)
{
  // Parse old map
  fail_unless(gt_input_map_parse_map("chr7:F127708134G27<+5>30A34<-5>40T88",map)==0,"Failed parsing old map");
  fail_unless(gt_strcmp(gt_map_get_seq_name(map),"chr7"),"Failed parsing old map. Seq Name");
  fail_unless(gt_map_get_strand(map)==FORWARD,"Failed parsing old map. Strand");
  // TODO

  // Parse old split-map
  fail_unless(gt_input_map_parse_map("[26]=chr7:R1203797~chr7:R1203108",map)==0,"Failed parsing old split-map");

  // Parse new map
  fail_unless(gt_input_map_parse_map("chr12:+:9570521:6C2>1+1>1-3T>1-2T12>4-8T23T8",map)==0,"Failed parsing old split-map");
}
END_TEST

Suite *gt_input_map_parser_suite(void) {
  Suite *s = suite_create("gt_input_map_parser");

  /* String parsers test case */
  TCase *tc_map_string_parser = tcase_create("MAP parser. String parsers");
  tcase_add_checked_fixture(tc_map_string_parser,gt_input_map_parser_setup,gt_input_map_parser_teardown);
  tcase_add_test(tc_map_string_parser,gt_test_imp_string_map);
  suite_add_tcase(s,tc_map_string_parser);

  return s;
}
