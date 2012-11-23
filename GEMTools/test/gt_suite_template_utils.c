/*
 * PROJECT: GEM-Tools library
 * FILE: gt_suite_input_map_parser.c
 * DATE: 03/10/2012
 * DESCRIPTION: // TODO
 */

#include "gt_test.h"

gt_template* source;
gt_template* target;

void gt_template_utils_setup(void) {
  source = gt_template_new();
  target = gt_template_new();
}

void gt_template_utils_teardown(void) {
  gt_template_delete(source, true, true);
  gt_template_delete(target, true, true);
}

START_TEST(gt_test_template_merge)
{

  fail_unless(gt_input_map_parse_template(
      "ID\tACGT\t####\t1\tchr1:-:20:4",source)==0);
  fail_unless(gt_input_map_parse_template(
      "ID\tACGT\t####\t1:1\tchr1:-:20:4,chr9:+:50:2C1",target)==0);
  // merge into source
  gt_template_merge_template_mmaps(source, target);


}
END_TEST

START_TEST(gt_test_template_to_string)
{

  fail_unless(gt_input_map_parse_template(
      "ID\tACGT\t####\t1\tchr1:-:20:4",source)==0);
  gt_string* string = gt_string_new(1024);
  gt_output_map_sprint_template(string, source, GT_ALL, true);
  // convert to string
  char * line = gt_string_get_string(string);
  printf("\n\n111'%s'222\n\n", line);
}
END_TEST

Suite *gt_template_utils_suite(void) {
  Suite *s = suite_create("gt_template_utils");

  /* String parsers test case */
  TCase *test_case = tcase_create("Template Utils");
  tcase_add_checked_fixture(test_case,gt_template_utils_setup,gt_template_utils_teardown);
  tcase_add_test(test_case,gt_test_template_merge);
  tcase_add_test(test_case,gt_test_template_to_string);
  suite_add_tcase(s,test_case);

  return s;
}
