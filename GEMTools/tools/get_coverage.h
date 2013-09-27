/*
 * get_coverage.h
 *
 *  Created on: 25 Sep 2013
 *      Author: heath
 */

#ifndef GET_COVERAGE_H_
#define GET_COVERAGE_H_

#define RANGE_BLK_SIZE 1024
#define NO_COUNT 0xffff
#define MAX_COUNT 0xfffe
#define OUT_WIDTH 80
#define ELAND_FMT 1
#define GEM_FMT 0

typedef u_int16_t count;

struct range_blk {
	struct range_blk *next;
	int idx;
	u_int32_t x1[RANGE_BLK_SIZE];
	u_int32_t x2[RANGE_BLK_SIZE];
	char *name[RANGE_BLK_SIZE];
};

struct contig {
	char *name;
	u_int32_t size;
	u_int32_t tsize;
	u_int32_t tot_bases;
	u_int32_t tot_tbases;
	u_int64_t tot_count;
	u_int64_t tot_tcount;
	count *counts;
	count *tcounts;
	struct range_blk *ranges;
	struct range_blk *tranges;
	pthread_mutex_t mut;
	UT_hash_handle hh;
};

struct match {
	char *ctg;
	char *cigar;
	u_int32_t pos;
	int orientation;
	int n_miss;
	int miss_score;
};

#endif /* GET_COVERAGE_H_ */
