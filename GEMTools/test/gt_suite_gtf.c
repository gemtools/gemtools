/*
 * PROJECT: GEM-Tools library
 * FILE: gt_suite_alignment.c
 * DATE: 03/10/2012
 * DESCRIPTION: // TODO
 */

#include "gt_test.h"
#include "gt_gtf.h"


void gt_gtf_setup(void) {
}

void gt_gtf_teardown(void) {
}

START_TEST(gt_test_gtf_entry_create)
{
	gt_string* type = gt_string_new(10);
	gt_string_set_string(type, "exon");
	gt_gtf_entry* e = gt_gtf_entry_new(1, 100, FORWARD, type);

    fail_unless(e->start==1,"Start not set properly");
    fail_unless(e->end==100,"End not set properly");
    fail_unless(e->strand==FORWARD,"Strand not set properly");
	fail_unless(e->type==type, "Type not assigned");

	gt_gtf_entry_delete(e);
	gt_string_delete(type);

}
END_TEST

START_TEST(gt_test_gtf_read_line)
{
	gt_string* l = gt_string_new(1024);
	gt_string_set_string(l, "chr1\thg19_refGene\texon\t11874\t12227\t0.000000\t+\t.\tgene_id \"NR_046018\"; transcript_id \"NR_046018\";");
	gt_gtf* gtf = gt_gtf_new();
	gt_gtf_read_line(gt_string_get_string(l), gtf);
	fail_unless(gt_shash_get_num_elements(gtf->types)==1, "Type not iserted");

	gt_gtf_ref* ref = gt_gtf_get_ref(gtf, "chr1");
	fail_unless(gt_vector_get_used(ref->entries)==1, "Reference entry no added");
	
	gt_string_set_string(l, "chr1\thg19_refGene\texon\t11880\t12235\t0.000000\t+\t.\tgene_id \"NR_046018\"; transcript_id \"NR_046018\";");
	gt_gtf_read_line(gt_string_get_string(l), gtf);
	fail_unless(gt_shash_get_num_elements(gtf->types)==1, "Type not iserted");
	fail_unless(gt_vector_get_used(ref->entries)==2, "Reference entry no added");
	
	gt_string_set_string(l, "chr2\thg19_refGene\tintron\t11880\t12235\t0.000000\t+\t.\tgene_id \"NR_046018\"; transcript_id \"NR_046018\";");
	gt_gtf_read_line(gt_string_get_string(l), gtf);
	fail_unless(gt_shash_get_num_elements(gtf->refs)==2, "Reference not inserted");
	fail_unless(gt_shash_get_num_elements(gtf->types)==2, "Second type not iserted");
	fail_unless(gt_vector_get_used(ref->entries)==2, "Reference entry no added");
	

	gt_gtf_delete(gtf);
	gt_string_delete(l);

}
END_TEST

START_TEST(gt_test_gtf_read)
{
    FILE* fp = fopen("testdata/chr1.gtf", "r");
    gt_gtf* gtf =  gt_gtf_read_from_stream(fp, 1);
	fclose(fp);
	fail_unless(gt_shash_get_num_elements(gtf->types)==4, "Not all types inserted");
	
//	gt_gtf_ref* ref_1 = gt_gtf_get_ref(gtf, "chr1");
//	gt_gtf_ref* ref_2 = gt_gtf_get_ref(gtf, "chr2");
//	fail_unless(gt_vector_get_used(ref_1->entries)==5, "Not all chr1 entries created");
//	fail_unless(gt_vector_get_used(ref_2->entries)==7, "Not all chr2 entries created");

	gt_gtf_delete(gtf);
}
END_TEST


START_TEST(gt_test_gtf_search)
{
  FILE* fp = fopen("testdata/chr1.gtf", "r");
  gt_gtf* gtf =  gt_gtf_read_from_stream(fp, 1);
  fclose(fp);
//	gt_vector* entries = gt_vector_new(10, sizeof(gt_gtf_entry*));
//
//	uint64_t i = 10;
//	// start
//	i = gt_gtf_search(gtf, entries, 'chr1', 0, 0);
//	fail_unless(i == 0 ,"Search failed for chr1 and 0");
//  i = gt_gtf_search(gtf, entries, 'chr2', 69091, 69091);
//  fail_unless(i == 3 ,"Search failed for chr1 and 69091");
//
	// end
//	i = gt_gtf_bin_search(entries, 69091);
//	fail_unless(i == 3 ,"Search failed for chr1 and 69091");
//	i = gt_gtf_bin_search(entries, 69092);
//	fail_unless(i == 3 ,"Search failed for chr1 and 69092");
//	// other values
//	i = gt_gtf_bin_search(entries, 11874);
//	fail_unless(i == 0 ,"Search failed for chr1 and 11874");
//	i = gt_gtf_bin_search(entries, 11875);
//	fail_unless(i == 0 ,"Search failed for chr1 and 11875");
//
//	i = gt_gtf_bin_search(entries, 12613);
//	fail_unless(i == 1 ,"Search failed for chr1 and 12613");
//	i = gt_gtf_bin_search(entries, 12614);
//	fail_unless(i == 1 ,"Search failed for chr1 and 12614");
//	i = gt_gtf_bin_search(entries, 69091);
//	fail_unless(i == 3 ,"Search failed for chr1 and 69091");
//	i = gt_gtf_bin_search(entries, 700000203);
//	fail_unless(i == 4 ,"Search failed for chr1 and 700000203");


	
	//GT_VECTOR_ITERATE(ref_1->entries,element,counter,gt_gtf_entry*) 
	//{
	//	fail_unless(s <= (*element)->start, "Not sorted");
	//	s = (*element)->start;
	//}

	gt_gtf_delete(gtf);
}
END_TEST

START_TEST(gt_test_gtf_find_matches)
{
    FILE* fp = fopen("testdata/chr1.gtf", "r");
    gt_gtf* gtf =  gt_gtf_read_from_stream(fp, 1);
	fclose(fp);
	gt_vector* target = gt_vector_new(5, sizeof(gt_gtf_entry*));	
	gt_gtf_entry* e;
	gt_gtf_search(gtf, target, "chr1", 1, 100);
	fail_unless(gt_vector_get_used(target) == 0, "Found something :(");

	// find exact
	gt_gtf_search(gtf, target, "chr1", 11874, 12227);
	fail_unless(gt_vector_get_used(target) == 1, "Found nothing");
	e = *gt_vector_get_elm(target, 0, gt_gtf_entry*);
	fail_unless( e->start == 11874, "Wrong Entry found");
	
	// find start < s
	gt_gtf_search(gtf, target, "chr1", 1000, 12227);
	fail_unless(gt_vector_get_used(target) == 1, "Found nothing");
	e = *gt_vector_get_elm(target, 0, gt_gtf_entry*);
	fail_unless( e->start == 11874, "Wrong Entry found");
	
	// find start > s
	gt_gtf_search(gtf, target, "chr1", 11900, 12227);
	fail_unless(gt_vector_get_used(target) == 1, "Found nothing");
	e = *gt_vector_get_elm(target, 0, gt_gtf_entry*);
	fail_unless( e->start == 11874, "Wrong Entry found");
	
	// find start > s end < e
	gt_gtf_search(gtf, target, "chr1", 11900, 12200);
	fail_unless(gt_vector_get_used(target) == 1, "Found nothing");
	e = *gt_vector_get_elm(target, 0, gt_gtf_entry*);
	fail_unless( e->start == 11874, "Wrong Entry found");
	
	// find start > s end > e
	gt_gtf_search(gtf, target, "chr1", 11900, 12230);
	fail_unless(gt_vector_get_used(target) == 2, "Found nothing");
	e = *gt_vector_get_elm(target, 1, gt_gtf_entry*);
	fail_unless( e->start == 11874, "Wrong Entry found");
	
	
	// find all
	gt_gtf_search(gtf, target, "chr1", 1, 1000000000);
	fail_unless(gt_vector_get_used(target) == 7, "Found nothing");
	
	// find tails
	gt_gtf_search(gtf, target, "chr1", 14409, 69092);
	fail_unless(gt_vector_get_used(target) == 2, "Found nothing");
	e = *gt_vector_get_elm(target, 1, gt_gtf_entry*);
	fail_unless( e->end == 69093, "Wrong Entry found");
	e = *gt_vector_get_elm(target, 0, gt_gtf_entry*);
	fail_unless( e->end == 70005, "Wrong Entry found");

	gt_gtf_delete(gtf);
}
END_TEST

Suite *gt_gtf_suite(void) {
  Suite *s = suite_create("gt_gtf");

  /* Core test case */
  TCase *tc_core = tcase_create("gt gtf");
  tcase_add_checked_fixture(tc_core,gt_gtf_setup,gt_gtf_teardown);
  tcase_add_test(tc_core,gt_test_gtf_entry_create);
  tcase_add_test(tc_core,gt_test_gtf_read_line);
  tcase_add_test(tc_core,gt_test_gtf_read);
  tcase_add_test(tc_core,gt_test_gtf_search);
  tcase_add_test(tc_core,gt_test_gtf_find_matches);
  suite_add_tcase(s,tc_core);

  return s;
}
